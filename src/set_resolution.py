"""Resolve an EntitySet to the leaf members it was split into.

The generator splits EntitySets into their member species, so a set has no
node of its own. Nothing recorded the link back, which silently cost the
benchmark 204 of 847 cases: every one of the 20 blocked readouts is a set,
and they are the canonical set-shaped readouts of the best-known pathways —
p-T,p-S-AKT, p-S9/21-GSK3, phospho-FOXO1/3/4/6, p-T,Y MAPK dimers, with 116
of the 204 in PIP3 alone.

**One hop is not enough.** 15 of those 20 resolve in a single hop, but the
other 5 are sets whose members are themselves sets — "p-T,Y MAPK monomers
and dimers" contains "p-T,Y MAPKs" and "p-T,Y MAPK dimers", neither of which
is a leaf. Recursion to leaves resolves all 20.

**On the cycle guard.** An earlier draft of the spec justified the visited
set by claiming Reactome's membership graph contains cycles. Measured on
Release97 that is false: zero self-membership, zero cycles at length 2, 3 or
4, across 5,440 nested sets, with maximum nesting depth 5. The guard stays
because a curation error would otherwise hang generation, not because
anything observed needs it. `MAX_DEPTH` is a second guard on the same risk.
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple, Optional, Set

from src.argument_parser import logger

# Observed maximum nesting depth is 5 on Release97. The bound is deliberately
# loose: exceeding it means curation changed shape, which should be reported
# rather than silently truncated.
MAX_DEPTH = 10


class SetLeaf(NamedTuple):
    """One leaf member of a set, and how far down it was found."""
    stable_id: str
    depth: int


class SetResolution(NamedTuple):
    """The result of resolving one set.

    ``truncated`` is the honest half: a caller must be able to tell a
    complete resolution from one that hit the depth bound, because combining
    member values over a partial set produces a plausible wrong number rather
    than a visible failure.
    """
    stable_id: str
    leaves: List[SetLeaf]
    truncated: bool
    max_depth_reached: int

    @property
    def leaf_ids(self) -> Set[str]:
        return {leaf.stable_id for leaf in self.leaves}


def resolve_set(
    stable_id: str,
    get_members,
    is_set,
    max_depth: int = MAX_DEPTH,
) -> SetResolution:
    """Expand ``stable_id`` through nested sets to its leaf members.

    ``get_members(stid) -> set[str]`` and ``is_set(stid) -> bool`` are passed
    in rather than imported so this is testable without Neo4j; production
    callers hand it ``neo4j_connector.get_set_members`` and a label check.

    A non-set input resolves to itself at depth 0, so callers do not need to
    branch on whether they hold a set.
    """
    if not is_set(stable_id):
        return SetResolution(stable_id, [SetLeaf(stable_id, 0)], False, 0)

    leaves: Dict[str, int] = {}
    visited: Set[str] = {stable_id}
    truncated = False
    deepest = 0
    frontier = [(stable_id, 0)]

    while frontier:
        current, depth = frontier.pop()
        if depth >= max_depth:
            # Report rather than silently stopping: a set deeper than the
            # bound means the graph changed shape.
            truncated = True
            logger.warning(
                f"set resolution hit the depth bound ({max_depth}) at {current} "
                f"while expanding {stable_id}; result is incomplete"
            )
            continue
        members = get_members(current) or set()
        if not members:
            # A set with no members is a leaf in practice; recording it keeps
            # the caller from seeing an empty resolution for a real entity.
            if current != stable_id:
                leaves.setdefault(current, depth)
                deepest = max(deepest, depth)
            continue
        for member in sorted(members):
            child_depth = depth + 1
            deepest = max(deepest, child_depth)
            if is_set(member):
                if member in visited:
                    # Unreachable on Release97; see the module docstring.
                    logger.warning(
                        f"cycle in set membership: {member} revisited while "
                        f"expanding {stable_id}"
                    )
                    continue
                visited.add(member)
                frontier.append((member, child_depth))
            else:
                # Keep the SHALLOWEST depth for a leaf reachable by several
                # routes; depth is "how far down this is", not "how we got here".
                if member not in leaves or child_depth < leaves[member]:
                    leaves[member] = child_depth

    resolved = [SetLeaf(sid, d) for sid, d in sorted(leaves.items())]
    if not resolved:
        logger.warning(f"set {stable_id} resolved to no leaves")
    return SetResolution(stable_id, resolved, truncated, deepest)


def make_neo4j_resolver(get_members_fn, get_labels_fn):
    """Bind ``resolve_set`` to Neo4j accessors.

    Kept separate so the recursion above stays free of the database.
    """
    def is_set(stable_id: str) -> bool:
        try:
            return "EntitySet" in (get_labels_fn(stable_id) or [])
        except Exception:
            # An unresolvable label means "not a set we can expand"; treating
            # it as a set would recurse into nothing and report a false empty.
            return False

    def resolve(stable_id: str, max_depth: int = MAX_DEPTH) -> SetResolution:
        return resolve_set(stable_id, get_members_fn, is_set, max_depth)

    return resolve
