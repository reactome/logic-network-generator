"""Variant keys: one canonical name per variant of an entity (deltasignal specs/046).

A *choice* sigma maps slots to chosen members. A slot is a set occurrence: the
OUTERMOST set on a hasComponent path from a participant, named by its stId.
Its chosen member is always a non-set entity (nested sets are flattened, as
break_apart_entity does). A chosen member that is itself a complex containing
sets contributes its own slots, so choices are made recursively.

    VK(e, sigma)
      leaf, modifier-isoform set, set-free complex, capped complex -> e
      bare set S                                    -> VK(sigma[S], sigma)
      complex C otherwise -> C + "::variant::" + "_".join(sorted(
                                 f"{slot}={VK(sigma[slot], sigma)}" ...))

Member complexes keep their identity (their stId, not their leaves), which is
what lets a producer that outputs an INNER set and a consumer that inputs an
OUTER set name the same node (the F10 fix). Keys are built from sorted
strings only: no uuid, dict or hash-seed order reaches them.

Design: specs/046 spec.md "Merged design" (Adam's decisions D1-D6) and
derivation-fable.md section 1; decision D2 is "one choice per set per
reaction", which `reaction_choices` implements by sharing slots across all of
a reaction's participants.
"""
from itertools import product as _product  # noqa: F401  (kept for callers)
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

SET_LABELS = ("EntitySet", "DefinedSet", "CandidateSet", "OpenSet")
VARIANT_SEP = "::variant::"

# Injected lookups, so the module is testable without a database. The
# defaults are the generator's cached Neo4j accessors.
_lookups: Dict[str, object] = {}
# Per-process memo of structure lookups (one Neo4j round trip per entity).
_memo: Dict[Tuple[str, str], List[str]] = {}


def reset_caches() -> None:
    """Forget memoised structure (tests, or a new database)."""
    _memo.clear()
    _variants_cache.clear()


def _labels(s: str) -> List[str]:
    if ("labels", s) in _memo:
        return _memo[("labels", s)]
    _memo[("labels", s)] = r = _labels_uncached(s)
    return r


def _labels_uncached(s: str) -> List[str]:
    f = _lookups.get("labels")
    if f is None:
        from src.neo4j_connector import get_labels as f
    try:
        return list(f(s) or [])
    except IndexError:
        return []


def _components(s: str) -> List[str]:
    if ("components", s) in _memo:
        return _memo[("components", s)]
    _memo[("components", s)] = r = _components_uncached(s)
    return r


def _components_uncached(s: str) -> List[str]:
    f = _lookups.get("components")
    if f is None:
        from src.neo4j_connector import get_complex_components as f
    return sorted(str(x) for x in (f(s) or {}))


def _members(s: str) -> List[str]:
    if ("members", s) in _memo:
        return _memo[("members", s)]
    _memo[("members", s)] = r = _members_uncached(s)
    return r


def _members_uncached(s: str) -> List[str]:
    f = _lookups.get("members")
    if f is None:
        from src.neo4j_connector import get_set_members as f
    return sorted(str(x) for x in (f(s) or {}))


def _atomic_sets() -> Set[str]:
    f = _lookups.get("atomic_sets")
    if f is None:
        from src.reaction_generator import modifier_isoform_set_ids as f
    return set(f())


def _max_variants() -> int:
    v = _lookups.get("max_variants")
    if v is None:
        from src.reaction_generator import MAX_VARIANTS as v
    return int(v() if callable(v) else v)


def is_set(s: str) -> bool:
    return any(t in _labels(s) for t in SET_LABELS) and s not in _atomic_sets()


def is_complex(s: str) -> bool:
    return "Complex" in _labels(s)


def flat_members(s: str, _depth: int = 0) -> List[str]:
    """A set's alternatives with nested sets flattened: non-set entities only."""
    out: Set[str] = set()
    if _depth > 12:
        return []
    for m in _members(s):
        if is_set(m):
            out |= set(flat_members(m, _depth + 1))
        else:
            out.add(m)
    return sorted(out)


def own_slots(e: str, _depth: int = 0) -> List[str]:
    """The set occurrences directly reachable from `e` (a complex or a bare
    set) without passing through another set: the slots `e` itself opens.
    A bare set is its own (single) slot. Sorted."""
    if _depth > 12:
        return []
    if is_set(e):
        return [e]
    if not is_complex(e):
        return []
    out: Set[str] = set()
    for c in _components(e):
        out |= set(own_slots(c, _depth + 1))
    return sorted(out)


_variants_cache: Dict[str, int] = {}


def variant_count(e: str, _depth: int = 0) -> int:
    """How many variants `e` has (product over its slots of the summed
    variants of each slot's members). Used to decide whether `e` is capped,
    identically wherever it appears."""
    if e in _variants_cache:
        return _variants_cache[e]
    n = 1
    if _depth <= 12:
        if is_set(e):
            n = sum(variant_count(m, _depth + 1) for m in flat_members(e)) or 1
        elif is_complex(e):
            for s in own_slots(e):
                n *= variant_count(s, _depth + 1)
    _variants_cache[e] = n
    return n


def is_capped(e: str) -> bool:
    mv = _max_variants()
    return mv > 0 and is_complex(e) and variant_count(e) > mv


def _complex_slots(e: str, sigma: Dict[str, str], _depth: int = 0) -> List[str]:
    """Every slot of complex `e` under choice `sigma`, the chosen member
    complexes' own slots hoisted in. Sorted, unique."""
    out: Set[str] = set()
    if _depth > 12:
        return []
    for s in own_slots(e):
        out.add(s)
        m = sigma.get(s)
        if m is not None and is_complex(m) and not is_capped(m):
            out |= set(_complex_slots(m, sigma, _depth + 1))
    return sorted(out)


def vkey(e: str, sigma: Dict[str, str]) -> str:
    """The variant key of `e` under choice `sigma` (every slot `e` opens,
    recursively, must be assigned)."""
    e = str(e)
    if e in _atomic_sets():
        return e
    if is_set(e):
        if e not in sigma:
            raise KeyError(f"slot {e} has no chosen member")
        return vkey(sigma[e], sigma)
    if not is_complex(e) or is_capped(e):
        return e
    slots = _complex_slots(e, sigma)
    if not slots:
        return e
    missing = [s for s in slots if s not in sigma]
    if missing:
        raise KeyError(f"{e}: slots {missing} have no chosen member")
    tokens = []
    for s in slots:
        m = sigma[s]
        # The member complex keeps its identity; its own slots are already
        # tokens of this key (hoisted), so name it by its plain stId here.
        tokens.append(f"{s}={m}")
    return e + VARIANT_SEP + "_".join(sorted(tokens))


def parse_variant_key(key: str) -> Tuple[str, Dict[str, str]]:
    """(parent stId, {slot: chosen member}) — the one parser every consumer
    of variant ids should use."""
    if VARIANT_SEP not in key:
        return key, {}
    parent, tail = key.split(VARIANT_SEP, 1)
    choice: Dict[str, str] = {}
    for tok in tail.split("_"):
        if not tok:
            continue
        if "=" not in tok:
            raise ValueError(f"malformed variant token {tok!r} in {key!r}")
        s, m = tok.split("=", 1)
        choice[s] = m
    return parent, choice


def variant_leaves(key: str, _depth: int = 0) -> Set[str]:
    """The leaf entities a variant is made of: the parent's structure with
    every slot replaced by its chosen member."""
    parent, sigma = parse_variant_key(key)
    return _leaves(parent, sigma, _depth)


def _leaves(e: str, sigma: Dict[str, str], _depth: int) -> Set[str]:
    if _depth > 12:
        return {e}
    if e in _atomic_sets():
        return {e}
    if is_set(e):
        m = sigma.get(e)
        return _leaves(m, sigma, _depth + 1) if m else {e}
    if is_complex(e) and not is_capped(e):
        out: Set[str] = set()
        for c in _components(e):
            out |= _leaves(c, sigma, _depth + 1)
        return out or {e}
    return {e}


def reaction_choices(participants: Sequence[str],
                     limit: Optional[int] = None) -> Iterator[Dict[str, str]]:
    """Every choice for a reaction whose (non-negative-regulator) participants
    are `participants`, with ONE variable per slot shared across all of them
    (decision D2). A chosen member complex opens its own slots. Yields
    choices in a deterministic order; stops after `limit` (the caller treats
    reaching it as "over the cap"). Capped complexes open no slots."""
    roots: Set[str] = set()
    for p in participants:
        if is_set(p):
            roots.add(p)
        elif is_complex(p) and not is_capped(p):
            roots |= set(own_slots(p))
    count = [0]

    def rec(sigma: Dict[str, str], open_slots: List[str]) -> Iterator[Dict[str, str]]:
        pending = sorted(s for s in open_slots if s not in sigma)
        if not pending:
            count[0] += 1
            yield dict(sigma)
            return
        s = pending[0]
        for m in flat_members(s):
            if limit is not None and count[0] >= limit:
                return
            sigma[s] = m
            new = set(open_slots)
            if is_complex(m) and not is_capped(m):
                new |= set(own_slots(m))
            yield from rec(sigma, sorted(new))
            del sigma[s]

    yield from rec({}, sorted(roots))


def variant_parts(node_id: str) -> Tuple[str, List[str]]:
    """(parent stId, chosen members) of a node id in EITHER format: the
    specs/046 key (``parent::variant::slot=member_...``) or the legacy
    ``parent::variant::m1_m2`` id. A plain stId gives (stId, []). Every
    consumer that needs a variant's members goes through this, so the two
    formats cannot be misread as each other."""
    if VARIANT_SEP not in node_id:
        return node_id, []
    parent, tail = node_id.split(VARIANT_SEP, 1)
    toks = [t for t in tail.split("_") if t]
    if any("=" in t for t in toks):
        return parent, sorted({t.split("=", 1)[1] for t in toks})
    return parent, toks
