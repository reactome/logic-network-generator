"""Virtual reactions from variant keys (deltasignal specs/046, LNG_VARIANT_NODES).

Replaces, under the flag, the two structures the emitter is built from:

- ``reaction_id_map``: one row per reaction copy (``uid``, ``reactome_id``);
- ``vr_entities``: ``uid -> (input_ids, output_ids, input_stoich, output_stoich)``,
  where the ids are now variant keys (src/variant_keys.py);

and the catalyst / regulator maps, whose ``entity_id`` per copy is the
participant's variant key under that copy's choice. Everything after that
(Phase 1 sharing per (key, reaction, role), Phase 2 merging across connected
reaction pairs, Phase 3 edges, boundary layer) is unchanged, which is how
decision D1 (identity across connected pairs only) is kept.

Copies (decision D2): one choice per set per reaction, over inputs, outputs,
catalysts and regulators (D4: negative regulators too). An output slot no
input-side participant holds is BOUND to an input-side slot whose members
correspond one to one by isoform-level reference identity; a tie is broken by
sorted rank and counted, never by list position (replaces finding F9). An
output slot that binds to nothing fans out.

Past the cap (more than LNG_MAX_VARIANTS copies) decision D6 applies, in order:
1. a participant none of whose slots reaches an output is read as ONE pool of
   its variants (``<stId>::pool``, wired by the generator, combined by OR:
   decision D5), and the reaction is enumerated over the rest;
2. if still over, every input-side participant with a slot outside the output
   variants is pooled, giving one copy per output variant;
3. if still over, ONE copy reads every slot-bearing input-side participant as
   a pool and writes its outputs by plain stable id (counted).
The first catalog arm used only a plain-id single copy; nothing produced those
plain ids, so a capped reaction was cut from everything upstream (IFN alpha/beta
"Expression of IFN-induced genes": 8 perturbations x 25 readouts).
"""
import uuid
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd

from src import variant_keys as vk

# `participant` is the curated catalyst/regulator the row came from (a bare
# set's copies name its MEMBERS in entity_id), read by depleter pooling.
_CAT_REG_COLUMNS = ["reaction_id", "entity_id", "edge_type", "uuid", "reaction_uuid", "participant"]

STATS: Counter = Counter()

_PARTICIPANT_CYPHER = """
UNWIND $rids AS rid
MATCH (r:ReactionLikeEvent {stId: rid})
OPTIONAL MATCH (r)-[i:input]->(ie:PhysicalEntity)
WITH r, collect(DISTINCT [ie.stId, coalesce(i.stoichiometry, 1)]) AS ins
OPTIONAL MATCH (r)-[o:output]->(oe:PhysicalEntity)
WITH r, ins, collect(DISTINCT [oe.stId, coalesce(o.stoichiometry, 1)]) AS outs
OPTIONAL MATCH (r)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(c:PhysicalEntity)
WITH r, ins, outs, collect(DISTINCT c.stId) AS cats
OPTIONAL MATCH (r)-[:regulatedBy]->(:PositiveRegulation)-[:regulator]->(p:PhysicalEntity)
WITH r, ins, outs, cats, collect(DISTINCT p.stId) AS pos
OPTIONAL MATCH (r)-[:regulatedBy]->(:NegativeRegulation)-[:regulator]->(n:PhysicalEntity)
RETURN r.stId AS rid, ins, outs, cats, pos, collect(DISTINCT n.stId) AS neg
"""


def fetch_participants(graph, reaction_ids: Sequence[str]) -> Dict[str, dict]:
    """reaction stId -> {"in": {stId: stoich}, "out": {...}, "cat": [...],
    "pos": [...], "neg": [...]}, everything sorted for determinism."""
    out: Dict[str, dict] = {}
    rids = sorted(set(str(r) for r in reaction_ids))
    for k in range(0, len(rids), 500):
        for row in graph.run(_PARTICIPANT_CYPHER, rids=rids[k:k + 500]).data():
            def stoich(pairs):
                d: Dict[str, int] = {}
                for s, n in pairs:
                    if s:
                        d[str(s)] = d.get(str(s), 0) + int(n or 1)
                return dict(sorted(d.items()))
            out[row["rid"]] = {
                "in": stoich(row["ins"]), "out": stoich(row["outs"]),
                "cat": sorted(x for x in row["cats"] if x),
                "pos": sorted(x for x in row["pos"] if x),
                "neg": sorted(x for x in row["neg"] if x),
            }
    return out


def _slots_of(participants: Sequence[str]) -> Set[str]:
    out: Set[str] = set()
    for p in participants:
        if vk.is_set(p):
            out.add(p)
        elif vk.is_complex(p) and not vk.is_capped(p):
            out |= set(vk.own_slots(p))
    return out


_sig_cache: Dict[str, Tuple[str, ...]] = {}


def _reference_signature(member: str) -> Tuple[str, ...]:
    """Isoform-level identity of a set member: the sorted multiset of its
    leaves' reference entities (a member complex keeps its leaf multiset)."""
    if member in _sig_cache:
        return _sig_cache[member]
    from src.neo4j_connector import get_reference_entity_id
    from src.reaction_generator import get_terminal_components
    leaves = sorted(get_terminal_components(member)) if vk.is_complex(member) else [member]
    sig = tuple(sorted(str(get_reference_entity_id(x) or x) for x in leaves))
    _sig_cache[member] = sig
    return sig


def bind_output_slots(input_side: Set[str], output_slots: Set[str]) -> Dict[str, Tuple[str, Dict[str, str]]]:
    """Output-only slot -> (input-side slot, {input member: output member}).

    A binding exists when every member of the output slot has a member of the
    input slot with the same reference signature, one to one. A signature
    shared by several members is paired by sorted rank (counted). The first
    input-side slot, in sorted order, that binds is used."""
    bound: Dict[str, Tuple[str, Dict[str, str]]] = {}
    free_in = sorted(s for s in input_side if s not in output_slots)
    for w in sorted(output_slots - input_side):
        wm = vk.flat_members(w)
        for u in free_in:
            um = vk.flat_members(u)
            if len(um) != len(wm):
                continue
            by_sig_u, by_sig_w = defaultdict(list), defaultdict(list)
            for m in um:
                by_sig_u[_reference_signature(m)].append(m)
            for m in wm:
                by_sig_w[_reference_signature(m)].append(m)
            if sorted((k, len(v)) for k, v in by_sig_u.items()) != \
               sorted((k, len(v)) for k, v in by_sig_w.items()):
                continue
            mapping: Dict[str, str] = {}
            for k, ms in by_sig_u.items():
                if len(ms) > 1:
                    STATS["slot_binding_ties"] += 1
                for a, b in zip(sorted(ms), sorted(by_sig_w[k])):
                    mapping[a] = b
            bound[w] = (u, mapping)
            STATS["output_slots_bound"] += 1
            break
        else:
            STATS["output_slots_free"] += 1
    return bound


def _choices(parts: dict, limit: int) -> Tuple[List[Dict[str, str]], bool]:
    """All copies of one reaction (choices), and whether the cap was hit."""
    input_side = list(parts["in"]) + parts["cat"] + parts["pos"] + parts["neg"]
    outputs = list(parts["out"])
    in_slots, out_slots = _slots_of(input_side), _slots_of(outputs)
    bound = bind_output_slots(in_slots, out_slots)
    # outputs whose slots are all bound contribute no free variable
    free_outputs = [o for o in outputs if not (_slots_of([o]) and _slots_of([o]) <= set(bound))]
    base = list(vk.reaction_choices(input_side + free_outputs, limit=limit + 1))
    result: List[Dict[str, str]] = []
    for sig in base:
        sig = dict(sig)
        for w, (u, mapping) in bound.items():
            if u in sig and sig[u] in mapping:
                sig[w] = mapping[sig[u]]
        # a bound member complex may open slots of its own: fan them out
        pending = [sig]
        while pending:
            cur = pending.pop()
            extra = _open_slots(outputs, cur)
            if not extra:
                result.append(cur)
                continue
            for more in vk.reaction_choices(extra, limit=limit + 1):
                pending.append({**cur, **more})
            if len(result) + len(pending) > limit:
                return result + pending, True
        if len(result) > limit:
            return result, True
    return result, len(result) > limit


def _missing_slots(e: str, sig: Dict[str, str], depth: int = 0) -> Set[str]:
    """Slots `e` needs under `sig` that `sig` does not assign, followed
    through chosen members recursively (a bound set's chosen member complex
    can open slots of its own)."""
    if depth > 12:
        return set()
    if vk.is_set(e):
        if e not in sig:
            return {e}
        return _missing_slots(sig[e], sig, depth + 1)
    if vk.is_complex(e) and not vk.is_capped(e):
        return {x for x in vk._complex_slots(e, sig) if x not in sig}
    return set()


def _open_slots(entities: Sequence[str], sig: Dict[str, str]) -> List[str]:
    """Sets still unassigned in `entities` under `sig` (to be fanned out)."""
    out: Set[str] = set()
    for e in entities:
        out |= _missing_slots(e, sig)
    return sorted(out)


POOL_SUFFIX = vk.POOL_SUFFIX


def deep_slots(p: str) -> Set[str]:
    """Every slot a participant can open, through chosen member complexes."""
    out: Set[str] = set()
    todo = list(_slots_of([p]))
    while todo:
        x = todo.pop()
        if x in out:
            continue
        out.add(x)
        for m in vk.flat_members(x):
            if vk.is_complex(m) and not vk.is_capped(m):
                todo += vk.own_slots(m)
    return out


def _without(parts: dict, pooled: Set[str]) -> dict:
    return {"in": {e: n for e, n in parts["in"].items() if e not in pooled},
            "out": parts["out"],
            "cat": [e for e in parts["cat"] if e not in pooled],
            "pos": [e for e in parts["pos"] if e not in pooled],
            "neg": [e for e in parts["neg"] if e not in pooled]}


def capped_fallback(parts: dict, limit: int) -> Tuple[List[Optional[Dict[str, str]]], Set[str], int]:
    """(choices, pooled participants, D6 step) for a reaction over the cap."""
    inside = list(parts["in"]) + parts["cat"] + parts["pos"] + parts["neg"]
    outputs = list(parts["out"])
    out_slots: Set[str] = set()
    for o in outputs:
        out_slots |= deep_slots(o)
    bound = bind_output_slots(_slots_of(inside), _slots_of(outputs))
    reach = out_slots | {u for (u, _) in bound.values()}
    slots = {x: deep_slots(x) for x in inside}
    # Step 1 pools only participants that share no slot with a kept one:
    # pooling A while a kept B fixes A's slot S would let each copy read every
    # S variant of A (decision D2, one choice per set per reaction; review of
    # vn6). Shrink to a fixed point.
    step1 = {x for x in inside if slots[x] and not (slots[x] & reach)}
    while True:
        kept_slots = set().union(*[slots[x] for x in inside if x not in step1]) if inside else set()
        shrink = {x for x in step1 if slots[x] & kept_slots}
        if not shrink:
            break
        step1 -= shrink
    step2 = {x for x in inside if slots[x] and not (slots[x] <= reach)}
    for step, pooled in ((1, step1), (2, step2)):
        choices, over = _choices(_without(parts, pooled), limit)
        if not over and choices:
            return choices, pooled, step
    return [None], {x for x in inside if slots[x]}, 3


def build_variant_reactions(graph, reaction_ids: Sequence[str]
                            ) -> Tuple[pd.DataFrame, Dict[str, tuple], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(reaction_id_map, vr_entities, catalyst_map, negative_regulator_map,
    positive_regulator_map) for LNG_VARIANT_NODES=1."""
    STATS.clear()
    cap = vk._max_variants()
    parts_by_rx = fetch_participants(graph, reaction_ids)
    rows: List[dict] = []
    vr_entities: Dict[str, tuple] = {}
    cat_rows: List[dict] = []; neg_rows: List[dict] = []; pos_rows: List[dict] = []
    for rx in sorted(parts_by_rx):
        parts = parts_by_rx[rx]
        choices, over = _choices(parts, cap if cap > 0 else 10 ** 9)
        if not choices and not over:
            # Never drop a reaction silently: zero copies means a set with no
            # members reached the enumeration, i.e. the structure is wrong.
            raise ValueError(f"{rx}: the variant enumeration produced no copy")
        pooled: Set[str] = set()
        if over:
            STATS["reactions_over_cap"] += 1
            choices, pooled, step = capped_fallback(parts, cap)
            STATS[f"capped_step{step}"] += 1
            STATS["pooled_participants"] += len(pooled)
        STATS["copies"] += len(choices)
        for sig in choices:
            uid = str(uuid.uuid4())

            def name(e, s=sig, pooled=pooled, output=False):
                # pooling applies to the input side only; an output is written
                # by its key (or plain stId in step 3)
                if e in pooled and not output:
                    return e + POOL_SUFFIX
                return e if s is None else vk.vkey(e, s)
            ins, outs = Counter(), Counter()
            for e, n in parts["in"].items():
                ins[name(e)] += n
            for e, n in parts["out"].items():
                outs[name(e, output=True)] += n
            rows.append({"uid": uid, "reactome_id": rx, "input_hash": None, "output_hash": None})
            vr_entities[uid] = (sorted(ins), sorted(outs), dict(ins), dict(outs))
            for lst, target in ((parts["cat"], cat_rows), (parts["neg"], neg_rows), (parts["pos"], pos_rows)):
                for e in lst:
                    target.append({"reaction_id": rx, "entity_id": name(e),
                                   "edge_type": "catalyst" if target is cat_rows else "regulator",
                                   "uuid": str(uuid.uuid4()), "reaction_uuid": uid,
                                   "participant": e})
    rid_map = pd.DataFrame(rows, columns=["uid", "reactome_id", "input_hash", "output_hash"])
    mk = lambda r: pd.DataFrame(r, columns=_CAT_REG_COLUMNS)  # noqa: E731
    return rid_map, vr_entities, mk(cat_rows), mk(neg_rows), mk(pos_rows)
