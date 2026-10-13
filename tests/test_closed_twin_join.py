"""deltasignal specs/053: a closed recycling copy is joined to its one fed twin."""
import pandas as pd

from src.closed_twin_join import apply_joins, find_joins


def row(uuid, kind="set_variant", entity="SPRED:NF1", leaves="NF1|SPRED2", chosen="SPRED2"):
    return {"uuid": uuid, "node_kind": kind, "diagram_entity_id": entity, "compartment": None,
            "member_leaves": leaves, "source_sets": "SPRED_set" if kind == "set_variant" else "",
            "chosen_members": chosen}


def rxn(uuid):
    return {"uuid": uuid, "node_kind": "reaction", "diagram_entity_id": uuid, "compartment": None,
            "member_leaves": "", "source_sets": "", "chosen_members": ""}


# The NF1 shape: F is made from NF1 by R_bind; G (same identity) is fed only by
# R_release, which releases it from RAS:GTP:GAP (a complex containing it).
NF1 = row("nf1", kind="simple_entity", entity="NF1", leaves="NF1", chosen="")
SPRED = row("spred", kind="simple_entity", entity="SPRED2", leaves="SPRED2", chosen="")
RAS = row("ras", kind="simple_entity", entity="RAS", leaves="RAS", chosen="")
F, G = row("f-fed"), row("g-closed")
RASGAP = row("rasgap", entity="RAS:GAP", leaves="NF1|RAS|SPRED2")
ROWS = [NF1, SPRED, RAS, F, G, RASGAP, rxn("R_bind"), rxn("R_gap"), rxn("R_release"), rxn("R_ub")]
EDGES = [
    ("nf1", "R_bind", "input"), ("spred", "R_bind", "input"), ("R_bind", "f-fed", "output"),
    ("f-fed", "R_ub", "input"),
    ("g-closed", "R_gap", "input"), ("ras", "R_gap", "input"), ("R_gap", "rasgap", "output"),
    ("rasgap", "R_release", "input"), ("rasgap", "R_release", "catalyst"),
    ("R_release", "g-closed", "output"), ("R_release", "ras", "output"),
]


def test_nf1_shape_joins_the_gap_copy_to_the_fed_copy():
    joins, ambiguous = find_joins(EDGES, ROWS)
    assert joins == {"g-closed": "f-fed"}
    assert ambiguous == 0


def test_two_fed_twins_are_left_alone_and_counted():
    rows = ROWS + [row("h-fed")]
    edges = EDGES + [("R_other", "h-fed", "output"), ("nf1", "R_other", "input")]
    joins, ambiguous = find_joins(edges, rows + [rxn("R_other")])
    assert joins == {}
    assert ambiguous == 1


def test_a_producer_that_does_not_release_it_is_not_closed():
    # R_release's input no longer contains G's leaves: a real synthesis, not a release
    rows = [r if r["uuid"] != "rasgap" else row("rasgap", entity="RAS:X", leaves="RAS|X") for r in ROWS]
    assert find_joins(EDGES, rows)[0] == {}


def test_a_non_output_in_edge_means_not_closed():
    edges = EDGES + [("nf1", "g-closed", "catalyst")]
    assert find_joins(edges, ROWS)[0] == {}


def test_different_identity_is_not_a_twin():
    rows = [r if r["uuid"] != "f-fed" else row("f-fed", chosen="SPRED1") for r in ROWS]
    assert find_joins(EDGES, rows)[0] == {}


def test_a_closed_twin_without_any_fed_copy_is_not_joined():
    edges = [e for e in EDGES if e[1] != "f-fed"]          # F loses its producer
    assert find_joins(edges, ROWS)[0] == {}


def test_apply_joins_repoints_drops_duplicates_and_keeps_flags_aligned():
    net = pd.DataFrame([
        {"source_id": "g-closed", "target_id": "R_gap", "pos_neg": "pos", "and_or": "and", "edge_type": "input"},
        {"source_id": "f-fed", "target_id": "R_gap", "pos_neg": "pos", "and_or": "and", "edge_type": "input"},
        {"source_id": "R_release", "target_id": "g-closed", "pos_neg": "pos", "and_or": "or", "edge_type": "output"},
        {"source_id": "a", "target_id": "b", "pos_neg": "pos", "and_or": "and", "edge_type": "assembly"},
    ])
    out, flags = apply_joins(net, [False, False, False, True], {"g-closed": "f-fed"})
    assert list(zip(out["source_id"], out["target_id"])) == [
        ("f-fed", "R_gap"), ("R_release", "f-fed"), ("a", "b")]
    assert flags == [False, False, True]
    assert "g-closed" not in set(out["source_id"]) | set(out["target_id"])


def test_no_joins_returns_the_network_unchanged():
    net = pd.DataFrame([{"source_id": "a", "target_id": "b", "edge_type": "input"}])
    out, flags = apply_joins(net, [True], {})
    assert out.equals(net) and flags == [True]
