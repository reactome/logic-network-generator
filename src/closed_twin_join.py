"""Join a closed recycling copy to its one fed identity twin (deltasignal specs/053).

The generator's uuids are position-aware: an entity appearing in two places
gets one uuid only where curated connectivity links the places. A cycle such as
bind → act → release then has its own copy of the released protein, fed by
nothing but the release step, even when the pathway also makes that protein
elsewhere. specs/053 traced the RAS GAP: SPRED dimer:NF1 made from NF1
(R-HSA-5658438) and SPRED dimer:NF1 released by GAP hydrolysis (R-HSA-5658231)
are two uuids with an identical nodes.csv identity, so an NF1 knockout never
reaches the GAP.

A node X is CLOSED when it has at least one in-edge, every in-edge is an
`output` edge, and each producing reaction has an `input` from a node whose
member leaves strictly contain X's. A FED TWIN is a node with X's identity (every
nodes.csv column but the uuid) that is not closed and has an in-edge. X is merged
into W when W is its only fed twin. With several fed twins X is left alone:
picking one would depend on random uuids, and merging all would bridge fed
copies to each other (the silo bridge, measured at 1:3.9 against).
"""
from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import pandas as pd

IDENTITY = ("node_kind", "diagram_entity_id", "compartment", "member_leaves",
            "source_sets", "chosen_members")


def _leaves(row: Mapping) -> frozenset:
    members = [m for m in str(row.get("member_leaves") or "").split("|") if m]
    if members:
        return frozenset(members)
    d = row.get("diagram_entity_id")
    return frozenset([str(d)]) if d else frozenset()


def _identity(row: Mapping) -> tuple:
    return tuple("" if row.get(k) is None or (isinstance(row.get(k), float) and pd.isna(row.get(k)))
                 else str(row.get(k)) for k in IDENTITY)


def find_joins(edges: Iterable[Tuple[str, str, str]], rows: Sequence[Mapping]
               ) -> Tuple[Dict[str, str], int]:
    """({closed uuid: its fed twin}, number of closed copies with several fed
    twins, left alone). `edges` are (source, target, edge_type)."""
    rows_by = {str(r["uuid"]): r for r in rows}
    inn: Dict[str, Set[str]] = defaultdict(set)
    etypes: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for s, t, et in edges:
        s, t = str(s), str(t)
        inn[t].add(s)
        etypes[(s, t)].add(str(et))
    leaves = {u: _leaves(r) for u, r in rows_by.items()}

    def releases(p: str, x: str) -> bool:
        lx = leaves.get(x, frozenset())
        return any("input" in etypes[(s, p)] and s != x and lx < leaves.get(s, frozenset())
                   for s in inn[p])

    closed: Set[str] = set()
    for x, r in rows_by.items():
        if r.get("node_kind") == "reaction" or not inn[x]:
            continue
        if all(etypes[(p, x)] == {"output"} for p in inn[x]) and all(releases(p, x) for p in inn[x]):
            closed.add(x)
    groups: Dict[tuple, List[str]] = defaultdict(list)
    for u, r in rows_by.items():
        groups[_identity(r)].append(u)
    joins: Dict[str, str] = {}
    ambiguous = 0
    for x in sorted(closed):
        fed = [w for w in groups[_identity(rows_by[x])] if w != x and w not in closed and inn[w]]
        if len(fed) == 1:
            joins[x] = fed[0]
        elif len(fed) > 1:
            ambiguous += 1
    return joins, ambiguous


def apply_joins(net: pd.DataFrame, mask: List[bool], joins: Mapping[str, str]
                ) -> Tuple[pd.DataFrame, List[bool]]:
    """Re-point every edge end in `joins` to its twin and drop exact duplicate
    edges (keeping the first, with its boundary flag). A self-loop the merge
    creates is kept: it is what the curated cycle says."""
    if not joins:
        return net, list(mask)
    out = net.copy()
    for col in ("source_id", "target_id"):
        out[col] = out[col].map(lambda u: joins.get(str(u), u) if pd.notna(u) else u)
    keep = ~out.astype(str).duplicated(keep="first")
    flags = [f for f, k in zip(mask, keep) if k]
    return out[keep.values].reset_index(drop=True), flags
