"""Specialization bridges (deltasignal specs/067, ``LNG_SPECIALIZATION_BRIDGE``).

Reactome often curates one species under two ids: a GENERIC complex whose
components are sets (PDGF:Phospho-PDGF receptor dimer = Phospho-PDGF receptor
dimer set : Active PDGF dimers set, R-HSA-186811) and a SPECIFIC complex
whose components are members of those sets (Phospho PDGF beta receptor:PDGF
chain B homodimer, R-HSA-389074). One reaction produces one form and the next
consumes the other; the two share no entity, and only a ``precedingEvent``
links the reactions. In Signaling by PDGF this cuts every receptor
perturbation off from every signalling readout.

The generator's variant nodes (deltasignal specs/046) already give the generic
complex one node per member combination, so a copy of the specific complex
and the copy of the generic one with the same member leaves are the same
species. This module finds the pairs; ``apply_specialization_bridges`` joins
the copies.

Modes: ``add`` (default since deltasignal specs/067 was adopted, 2026-10-17:
curator held-out +59, experimental 0; emit the bridge), ``replace`` (emit it
and drop the boundary assembly/composition in-edges of each bridged target
copy; measured, failed the gate: DNA Damage Bypass -14), ``off``.
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

_MODES = ("off", "add", "replace")
_CLOSED_MODES = ("off", "twin")
EDGE_TYPE = "specialization"
_MAX_SLOTS = 12


def specialization_bridge_mode() -> str:
    """LNG_SPECIALIZATION_BRIDGE: ``add`` (default, deltasignal specs/067) |
    ``replace`` | ``off``. Any other value is an error at startup."""
    value = os.environ.get("LNG_SPECIALIZATION_BRIDGE", "add")
    if value not in _MODES:
        raise ValueError(f"LNG_SPECIALIZATION_BRIDGE must be one of {_MODES}, got {value!r}")
    return value


def specialization_closed_mode() -> str:
    """LNG_SPECIALIZATION_CLOSED (deltasignal specs/068): ``off`` (default) |
    ``twin``. Under ``twin`` the carrier-loop closedness test ignores
    ``specialization`` in-edges: the bridge joins two copies of one species,
    as ``closed_twin``'s root-copy supply does, so it does not open the loop
    (PLCG1's release of the generic PDGF receptor, R-HSA-186811). The edges
    stay in the network; only ``carrier_loops.csv`` changes."""
    value = os.environ.get("LNG_SPECIALIZATION_CLOSED", "off")
    if value not in _CLOSED_MODES:
        raise ValueError(f"LNG_SPECIALIZATION_CLOSED must be one of {_CLOSED_MODES}, got {value!r}")
    return value


class Structure:
    """Complex components (a multiset, by stoichiometry) and set members (one
    level), as read from Neo4j or given by a test."""

    def __init__(self, components: Mapping[str, Sequence[str]], members: Mapping[str, Iterable[str]]):
        self.components = {k: list(v) for k, v in components.items()}
        self.members = {k: set(v) for k, v in members.items()}
        self._star: Dict[str, Set[str]] = {}
        self._spec: Dict[Tuple[str, str], bool] = {}

    def members_star(self, x: str) -> Set[str]:
        """Every member of set `x`, transitively (a set of sets included)."""
        if x not in self._star:
            self._star[x] = set()      # guards a cycle
            out: Set[str] = set()
            for m in self.members.get(x, ()):
                out.add(m)
                out |= self.members_star(m)
            self._star[x] = out
        return self._star[x]

    def fits(self, s: str, g: str) -> bool:
        """Slot `s` of the specific complex can fill slot `g` of the generic."""
        if s == g:
            return True
        mg = self.members_star(g)
        if mg:
            if s in mg:
                return True
            ms = self.members_star(s)
            if ms and ms <= mg:
                return True
            if any(self.specializes(s, m) for m in mg if m in self.components):
                return True
        return s in self.components and g in self.components and self.specializes(s, g)

    def specializes(self, s: str, g: str) -> bool:
        """Complex `s` is one member combination (or a narrower set choice) of
        complex `g`: components match one to one, with stoichiometry, each
        `s` component fitting its `g` component."""
        if s == g or s not in self.components or g not in self.components:
            return False
        key = (s, g)
        if key in self._spec:
            return self._spec[key]
        self._spec[key] = False        # guards recursion through nesting
        a, b = self.components[s], self.components[g]
        if len(a) != len(b) or not a or len(a) > _MAX_SLOTS:
            return False
        adj = [[j for j, y in enumerate(b) if self.fits(x, y)] for x in a]
        match = [-1] * len(b)

        def augment(i: int, seen: Set[int]) -> bool:
            for j in adj[i]:
                if j in seen:
                    continue
                seen.add(j)
                if match[j] < 0 or augment(match[j], seen):
                    match[j] = i
                    return True
            return False

        ok = all(augment(i, set()) for i in range(len(a)))
        self._spec[key] = ok
        return ok


def linked_pairs(reactions: Sequence[Tuple[str, Set[str], Set[str], Set[str], Set[str]]],
                 structure: Structure) -> List[Tuple[str, str]]:
    """(produced, consumed) complex pairs to bridge.

    `reactions` rows are (reaction, inputs, outputs, catalysts and regulators,
    preceding reactions). A pair qualifies when one complex specializes the
    other and some reaction consuming the second lists some reaction producing
    the first as ``precedingEvent``: the curators' own statement of order, so
    a bridge never runs against it (and a mutual pair is not welded)."""
    produced: Dict[str, Set[str]] = defaultdict(set)
    consumed: Dict[str, Set[str]] = defaultdict(set)
    preceding: Dict[str, Set[str]] = {}
    for rx, ins, outs, mods, pre in reactions:
        preceding[rx] = set(pre)
        for x in outs:
            produced[x].add(rx)
        for x in set(ins) | set(mods):
            consumed[x].add(rx)
    out = []
    for a in sorted(x for x in produced if x in structure.components):
        for b in sorted(x for x in consumed if x in structure.components):
            if a == b:
                continue
            if not (structure.specializes(a, b) or structure.specializes(b, a)):
                continue
            if any(p in preceding[c] for c in consumed[b] for p in produced[a]):
                out.append((a, b))
    return out


def copy_bridges(pairs: Sequence[Tuple[str, str]],
                 node_rows: Sequence[Mapping[str, str]]) -> List[Tuple[str, str]]:
    """(source uuid, target uuid) for each pair: every copy of the produced
    entity to every copy of the consumed one with identical member leaves.
    Pool nodes (the entity as a whole) are not copies of a species."""
    by_entity: Dict[str, List[Tuple[str, frozenset]]] = defaultdict(list)
    for r in node_rows:
        if r.get("node_kind") in ("variant_pool", "set_pool", "reaction"):
            continue
        ent = r.get("diagram_entity_id")
        if not ent:
            continue
        leaves = frozenset(x for x in (r.get("member_leaves") or "").split("|") if x)
        by_entity[ent].append((str(r["uuid"]), leaves))
    out = []
    for a, b in pairs:
        for u, la in by_entity.get(a, ()):
            for v, lb in by_entity.get(b, ()):
                if u != v and la == lb:
                    out.append((u, v))
    return sorted(set(out))


def apply_specialization_bridges(network, boundary_mask: List[bool],
                                 bridges: Sequence[Tuple[str, str]], mode: str):
    """Return (network, boundary_mask) with the bridges added (curated, not
    boundary) and, under ``replace``, the boundary assembly/composition
    in-edges of each bridged target removed."""
    import pandas as pd

    if mode == "off" or not bridges:
        return network, boundary_mask
    mask = list(boundary_mask)
    if mode == "replace":
        targets = {v for _, v in bridges}
        keep = [not (m and str(t) in targets and et in ("assembly", "composition"))
                for t, et, m in zip(network["target_id"], network["edge_type"], mask)]
        network = network[keep].reset_index(drop=True)
        mask = [m for m, k in zip(mask, keep) if k]
    existing = set(zip(network["source_id"].astype(str), network["target_id"].astype(str),
                       network["edge_type"].astype(str)))
    rows = [{"source_id": u, "target_id": v, "pos_neg": "pos", "and_or": "or",
             "edge_type": EDGE_TYPE, "stoichiometry": 1, "edge_reaction_id": None}
            for u, v in bridges if (u, v, EDGE_TYPE) not in existing]
    if rows:
        add = pd.DataFrame(rows, columns=list(network.columns))
        if "stoichiometry" in network.columns:
            add["stoichiometry"] = add["stoichiometry"].astype(network["stoichiometry"].dtype)
        network = pd.concat([network, add], ignore_index=True)
        mask = mask + [False] * len(rows)
    return network, mask


def fetch_reactions(pathway_id: str):
    """Every reaction-like event of the pathway: inputs, outputs, catalysts
    and regulators, and its preceding events."""
    from src.neo4j_connector import get_graph
    rows = get_graph().run(
        """
        MATCH (p:Pathway {stId: $pid})-[:hasEvent*0..]->(rx:ReactionLikeEvent)
        WITH DISTINCT rx
        OPTIONAL MATCH (rx)-[:input]->(i)
        OPTIONAL MATCH (rx)-[:output]->(o)
        OPTIONAL MATCH (rx)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(c)
        OPTIONAL MATCH (rx)-[:regulatedBy]->(:Regulation)-[:regulator]->(g)
        OPTIONAL MATCH (rx)-[:precedingEvent]->(pe)
        RETURN rx.stId AS rx, collect(DISTINCT i.stId) AS ins, collect(DISTINCT o.stId) AS outs,
               collect(DISTINCT c.stId) + collect(DISTINCT g.stId) AS mods,
               collect(DISTINCT pe.stId) AS pre
        """,
        pid=pathway_id,
    ).data()
    return [(r["rx"], set(r["ins"]), set(r["outs"]), set(r["mods"]), set(r["pre"]))
            for r in sorted(rows, key=lambda r: r["rx"])]


def fetch_structure(entity_ids: Iterable[str]) -> Structure:
    """Components (with stoichiometry) and set members of `entity_ids` and
    everything they contain."""
    from src.neo4j_connector import get_graph
    ids = sorted({e for e in entity_ids if e})
    if not ids:
        return Structure({}, {})
    rows = get_graph().run(
        """
        UNWIND $ids AS s
        MATCH (e:PhysicalEntity {stId: s})-[:hasComponent|hasMember|hasCandidate*0..12]->(x:PhysicalEntity)
        WITH DISTINCT x
        OPTIONAL MATCH (x)-[r:hasComponent]->(c)
        WITH x, collect([c.stId, coalesce(r.stoichiometry, 1)]) AS comps
        OPTIONAL MATCH (x)-[:hasMember|hasCandidate]->(m)
        RETURN x.stId AS x, 'Complex' IN labels(x) AS is_complex, comps, collect(DISTINCT m.stId) AS mems
        """,
        ids=ids,
    ).data()
    components: Dict[str, List[str]] = {}
    members: Dict[str, Set[str]] = {}
    for r in rows:
        if r["is_complex"]:
            components[r["x"]] = sorted(c for c, st in r["comps"] if c
                                        for _ in range(int(st or 1)))
        if r["mems"]:
            members[r["x"]] = set(r["mems"])
    return Structure(components, members)


def specialization_bridges(pathway_id: str, node_rows: Sequence[Mapping[str, str]]):
    """(pairs, copy bridges) for a pathway, read from Neo4j."""
    reactions = fetch_reactions(pathway_id)
    entities = set()
    for _rx, ins, outs, mods, _pre in reactions:
        entities |= ins | outs | mods
    structure = fetch_structure(entities)
    pairs = linked_pairs(reactions, structure)
    return pairs, copy_bridges(pairs, node_rows)
