"""Unit tests for recursive EntitySet resolution. No Neo4j."""
from src.set_resolution import MAX_DEPTH, resolve_set


def _fixture():
    """The real shape that broke one-hop resolution.

    R-HSA-5674340 "p-T,Y MAPK monomers and dimers" contains two sets, neither
    of which is a leaf. It is one of the 5 of 20 blocked readouts that a
    single hop could not resolve.
    """
    members = {
        "R-HSA-5674340": {"R-HSA-169289", "R-HSA-1268261"},
        "R-HSA-169289": {"MAPK1", "MAPK3"},
        "R-HSA-1268261": {"DIMER1", "DIMER2"},
        "R-HSA-202072": {"AKT1p", "AKT2p", "AKT3p"},
    }
    sets = set(members)
    return (lambda s: members.get(s, set()), lambda s: s in sets)


def test_single_hop_set_resolves_to_members():
    get, is_set = _fixture()
    r = resolve_set("R-HSA-202072", get, is_set)
    assert r.leaf_ids == {"AKT1p", "AKT2p", "AKT3p"}
    assert all(leaf.depth == 1 for leaf in r.leaves)
    assert not r.truncated


def test_nested_set_recurses_to_leaves():
    """One hop returns sets, not leaves — this is the 5-of-20 case."""
    get, is_set = _fixture()
    one_hop = get("R-HSA-5674340")
    assert all(is_set(m) for m in one_hop), "fixture must nest, or it tests nothing"

    r = resolve_set("R-HSA-5674340", get, is_set)
    assert r.leaf_ids == {"MAPK1", "MAPK3", "DIMER1", "DIMER2"}
    assert r.max_depth_reached == 2
    assert not r.truncated


def test_non_set_resolves_to_itself():
    get, is_set = _fixture()
    r = resolve_set("MAPK1", get, is_set)
    assert r.leaves == [("MAPK1", 0)]
    assert not r.truncated


def test_single_member_set_behaves_like_its_member():
    members = {"S": {"only"}}
    r = resolve_set("S", lambda s: members.get(s, set()), lambda s: s in members)
    assert r.leaf_ids == {"only"}


def test_diamond_keeps_the_shallowest_depth():
    """A leaf reachable two ways is reported once, at its shallowest depth.

    depth answers "how far down is this", not "which route did we take".
    """
    members = {"top": {"mid", "leaf"}, "mid": {"leaf"}}
    sets = {"top", "mid"}
    r = resolve_set("top", lambda s: members.get(s, set()), lambda s: s in sets)
    assert r.leaf_ids == {"leaf"}
    assert [leaf.depth for leaf in r.leaves] == [1]


def test_cycle_terminates_rather_than_hanging():
    """Release97 has no membership cycles; this proves the guard works anyway.

    Without the visited set this call does not return, so a curation error
    would hang generation rather than fail it.
    """
    members = {"A": {"B"}, "B": {"A", "leaf"}}
    sets = {"A", "B"}
    r = resolve_set("A", lambda s: members.get(s, set()), lambda s: s in sets)
    assert r.leaf_ids == {"leaf"}


def test_exceeding_the_depth_bound_is_reported_not_silently_truncated():
    """A partial resolution must be visibly partial.

    Combining member values over a silently-truncated set yields a plausible
    wrong number instead of a visible failure.
    """
    depth = MAX_DEPTH + 3
    members = {f"S{i}": {f"S{i+1}"} for i in range(depth)}
    members[f"S{depth}"] = {"leaf"}
    sets = set(members)
    r = resolve_set("S0", lambda s: members.get(s, set()), lambda s: s in sets)
    assert r.truncated, "hitting the depth bound must set truncated"
    assert "leaf" not in r.leaf_ids
