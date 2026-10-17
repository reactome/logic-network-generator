"""Specialization bridges (deltasignal specs/067). Fixtures follow the traced
Signaling by PDGF case: the generic PDGF:Phospho-PDGF receptor dimer
(R-HSA-186811 = Phospho-PDGF receptor dimer set : Active PDGF dimers set) and
the specific Phospho PDGF beta receptor:PDGF chain B homodimer (R-HSA-389074 =
the beta homodimer : PDGF-BB), ordered by precedingEvent."""
import pandas as pd
import pytest

from src import specialization as sp
from src.env_flags import OTHER_FLAGS

# generic G = {receptor dimer set, PDGF dimer set}; specific S = {beta, BB}
COMPONENTS = {
    "G": ["setR", "setP"],
    "S": ["beta", "BB"],
    "S_ab": ["ab", "setAB_BB"],      # a narrower set in one slot
    "S_extra": ["beta", "BB", "PLCG"],
    "S_homo": ["beta", "beta"],       # stoichiometry: two receptors, no ligand
    "Other": ["beta", "GRB2"],
}
MEMBERS = {
    "setR": {"beta", "ab", "alpha"},
    "setP": {"BB", "AB", "AA"},
    "setAB_BB": {"AB", "BB"},
}


@pytest.fixture
def structure():
    return sp.Structure(COMPONENTS, MEMBERS)


def test_member_combination_specializes(structure):
    assert structure.specializes("S", "G")
    assert not structure.specializes("G", "S")


def test_narrower_set_specializes(structure):
    assert structure.specializes("S_ab", "G")


def test_component_count_and_stoichiometry_must_match(structure):
    assert not structure.specializes("S_extra", "G")
    assert not structure.specializes("S_homo", "G")


def test_a_component_outside_the_sets_does_not(structure):
    assert not structure.specializes("Other", "G")
    assert not structure.specializes("S", "S")


def test_nested_complex_slot(structure):
    s = sp.Structure({"G2": ["G", "X"], "S2": ["S", "X"], "G": COMPONENTS["G"],
                      "S": COMPONENTS["S"]}, MEMBERS)
    assert s.specializes("S2", "G2")


# reaction rows: (reaction, inputs, outputs, catalysts+regulators, preceding)
AUTOPHOS = ("autophos", {"ATP", "S0"}, {"S", "ADP"}, set(), set())
GRB2_BINDS = ("grb2", {"G", "GRB2"}, {"G:GRB2"}, set(), {"autophos"})


def test_linked_pair_is_found(structure):
    assert sp.linked_pairs([AUTOPHOS, GRB2_BINDS], structure) == [("S", "G")]


def test_unlinked_pair_is_not_bridged(structure):
    unlinked = ("grb2", {"G", "GRB2"}, {"G:GRB2"}, set(), set())
    assert sp.linked_pairs([AUTOPHOS, unlinked], structure) == []


def test_generic_produced_specific_consumed(structure):
    binds = ("binds", {"setP", "setR"}, {"G"}, set(), set())
    autophos = ("autophos", {"S", "ATP"}, {"pS"}, {"S"}, {"binds"})
    assert sp.linked_pairs([binds, autophos], structure) == [("G", "S")]


def test_regulator_consumption_counts(structure):
    regulated = ("r", {"Y"}, {"Z"}, {"G"}, {"autophos"})
    assert sp.linked_pairs([AUTOPHOS, regulated], structure) == [("S", "G")]


def rows():
    return [
        {"uuid": "s1", "node_kind": "simple_complex", "diagram_entity_id": "S",
         "member_leaves": "BB|beta"},
        {"uuid": "g_bb_beta", "node_kind": "set_variant", "diagram_entity_id": "G",
         "member_leaves": "BB|beta"},
        {"uuid": "g_aa_alpha", "node_kind": "set_variant", "diagram_entity_id": "G",
         "member_leaves": "AA|alpha"},
        {"uuid": "g_pool", "node_kind": "variant_pool", "diagram_entity_id": "G",
         "member_leaves": "BB|beta"},
    ]


def test_copies_join_only_on_identical_leaves():
    assert sp.copy_bridges([("S", "G")], rows()) == [("s1", "g_bb_beta")]


def network():
    return pd.DataFrame([
        # boundary decomposition of the generic copy (no producer in pathway)
        {"source_id": "beta", "target_id": "g_bb_beta", "pos_neg": "pos", "and_or": "and",
         "edge_type": "assembly", "stoichiometry": 1, "edge_reaction_id": None},
        {"source_id": "BB", "target_id": "g_bb_beta", "pos_neg": "pos", "and_or": "and",
         "edge_type": "assembly", "stoichiometry": 1, "edge_reaction_id": None},
        # the same decomposition into a copy that is not bridged
        {"source_id": "alpha", "target_id": "g_aa_alpha", "pos_neg": "pos", "and_or": "and",
         "edge_type": "assembly", "stoichiometry": 1, "edge_reaction_id": None},
        # a curated edge into the bridged copy
        {"source_id": "g_bb_beta", "target_id": "rx", "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-1"},
        {"source_id": "rx0", "target_id": "s1", "pos_neg": "pos", "and_or": "or",
         "edge_type": "output", "stoichiometry": 1, "edge_reaction_id": "R-0"},
    ]).astype({"stoichiometry": "Int64"})


MASK = [True, True, True, False, False]


def test_off_changes_nothing():
    net, mask = sp.apply_specialization_bridges(network(), MASK, [("s1", "g_bb_beta")], "off")
    assert net.equals(network()) and mask == MASK


def test_add_emits_a_curated_or_edge():
    net, mask = sp.apply_specialization_bridges(network(), MASK, [("s1", "g_bb_beta")], "add")
    assert len(net) == 6 and mask == MASK + [False]
    e = net.iloc[-1]
    assert (e.source_id, e.target_id, e.pos_neg, e.and_or, e.edge_type) == \
        ("s1", "g_bb_beta", "pos", "or", "specialization")
    assert str(net["stoichiometry"].dtype) == "Int64"


def test_replace_drops_boundary_assembly_into_bridged_copies_only():
    net, mask = sp.apply_specialization_bridges(network(), MASK, [("s1", "g_bb_beta")], "replace")
    into = net[net.target_id == "g_bb_beta"]
    assert list(into.edge_type) == ["specialization"]
    assert (net.target_id == "g_aa_alpha").sum() == 1      # other copy keeps its decomposition
    assert len(net) == len(mask) == 4
    assert mask == [True, False, False, False]


def test_replace_keeps_curated_assembly():
    curated = network()
    m = [False] + MASK[1:]       # the first assembly edge is curated, not boundary
    net, _ = sp.apply_specialization_bridges(curated, m, [("s1", "g_bb_beta")], "replace")
    assert ((net.source_id == "beta") & (net.target_id == "g_bb_beta")).sum() == 1


def test_existing_bridge_is_not_duplicated():
    net, mask = sp.apply_specialization_bridges(network(), MASK, [("s1", "g_bb_beta")], "add")
    net2, mask2 = sp.apply_specialization_bridges(net, mask, [("s1", "g_bb_beta")], "add")
    assert len(net2) == len(net) and mask2 == mask


def test_mode_is_validated(monkeypatch):
    monkeypatch.delenv("LNG_SPECIALIZATION_BRIDGE", raising=False)
    assert sp.specialization_bridge_mode() == "add"     # default since specs/067 adoption
    monkeypatch.setenv("LNG_SPECIALIZATION_BRIDGE", "on")
    with pytest.raises(ValueError):
        sp.specialization_bridge_mode()
    from src.logic_network_generator import _reject_removed_env
    with pytest.raises(ValueError):
        _reject_removed_env()


def test_flag_is_registered_and_fingerprinted():
    import src.pathway_generator as pg
    assert "LNG_SPECIALIZATION_BRIDGE" in OTHER_FLAGS
    assert "LNG_SPECIALIZATION_BRIDGE" in pg._FINGERPRINTED_ENV


def _result():
    from src.logic_network_generator import PathwayResult
    return PathwayResult(network(), {}, pd.DataFrame(), pd.DataFrame({"uid": [], "reactome_id": []}),
                         {}, boundary_mask=list(MASK))


def test_pipeline_step_off_returns_the_result_untouched(monkeypatch):
    import src.pathway_generator as pg
    monkeypatch.setenv("LNG_SPECIALIZATION_BRIDGE", "off")
    monkeypatch.setattr(pg, "specialization_bridges",
                        lambda *a: pytest.fail("off must not read Neo4j"))
    r = _result()
    assert pg.apply_specialization_step(r, "R-HSA-1") is r


def test_pipeline_step_applies_the_bridges(monkeypatch):
    import src.pathway_generator as pg
    import src.logic_network_generator as lng
    monkeypatch.delenv("LNG_SPECIALIZATION_BRIDGE", raising=False)      # default add
    monkeypatch.setattr(lng, "node_rows", lambda *a: rows())
    monkeypatch.setattr(pg, "specialization_bridges",
                        lambda pid, node_rows: ([("S", "G")], sp.copy_bridges([("S", "G")], node_rows)))
    out = pg.apply_specialization_step(_result(), "R-HSA-1")
    assert len(out.logic_network) == 6 and out.boundary_mask == MASK + [False]
    assert out.logic_network.iloc[-1].edge_type == "specialization"
