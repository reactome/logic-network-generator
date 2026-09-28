"""Interconversion pools (deltasignal specs/039): each protein's forms joined by
curated forward/reverse reactions, found at NODE level and oriented so the base
is the resting (unmodified) form. Shape: RAS. A = RAS:GDP, B = RAS:GTP;
F = "GEF exchange" consumes GTP (a donor), R = "GAP hydrolysis"."""
import pandas as pd

import src.logic_network_generator as m

# _uuid_to_stable_id_map detects the mapping direction by uuid SHAPE
U = {k: f"00000000-0000-0000-0000-{i:012d}" for i, k in enumerate(["a", "a2", "b", "vf", "vr"], 1)}

A, B, F, R = "R-HSA-GDP", "R-HSA-GTP", "R-HSA-F", "R-HSA-R"
PAIRS = [(F, R, A, B), (R, F, B, A)]


def e(s, t, et):
    return {"source_id": U[s], "target_id": U[t], "edge_type": et}


def net(loop=True):
    # a -> f -> b -> r -> a (the node-level loop); with loop=False, R's output is a
    # different copy of A (a2), which the uuids keep apart
    a_out = "a" if loop else "a2"
    edges = pd.DataFrame([e("a", "vf", "input"), e("vf", "b", "output"),
                          e("b", "vr", "input"), e("vr", a_out, "output")])
    rmap = pd.DataFrame({"uid": [U["vf"], U["vr"]], "reactome_id": [F, R]})
    umap = {U["a"]: A, U["a2"]: A, U["b"]: B}
    return edges, rmap, umap


def test_a_node_level_loop_is_one_pool_with_the_resting_form_as_base():
    edges, rmap, umap = net()
    forms, trans = m.find_pools(edges, rmap, umap, PAIRS, {A: (0, 2), B: (0, 2)}, {F})
    assert sorted((u, base) for _, u, _, base in forms) == sorted([(U["a"], True), (U["b"], False)])
    assert sorted((f, t, rx) for _, f, t, rx, _ in trans) == sorted([(U["a"], U["b"], U["vf"]), (U["b"], U["a"], U["vr"])])
    assert {p for p, *_ in forms} == {"pool1"}


def test_copies_the_uuids_separated_are_not_pooled():
    edges, rmap, umap = net(loop=False)
    forms, trans = m.find_pools(edges, rmap, umap, PAIRS, {}, {F})
    assert forms == [] and trans == []


def test_orientation_by_residues_then_donor_then_components():
    edges, rmap, umap = net()
    # residues decide first, whatever the donor says
    forms, _ = m.find_pools(edges, rmap, umap, PAIRS, {A: (2, 1), B: (0, 1)}, {F})
    assert dict((u, b) for _, u, _, b in forms)[U["b"]] is True
    # no residue difference: the donor sits on the REVERSE reaction, so GDP is modified
    forms, _ = m.find_pools(edges, rmap, umap, PAIRS, {A: (0, 2), B: (0, 2)}, {R})
    assert dict((u, b) for _, u, _, b in forms) == {U["a"]: False, U["b"]: True}
    # no residue or donor difference: the form with more components is modified
    forms, _ = m.find_pools(edges, rmap, umap, PAIRS, {A: (0, 1), B: (0, 3)}, set())
    assert dict((u, b) for _, u, _, b in forms)[U["a"]] is True


def test_undecided_pool_is_counted_as_a_tie():
    edges, rmap, umap = net()
    st = {}
    forms, _ = m.find_pools(edges, rmap, umap, PAIRS, {A: (0, 1), B: (0, 1)}, set(), st)
    assert st["ties"] == 1 and len(forms) == 2
    st = {}
    m.find_pools(edges, rmap, umap, PAIRS, {A: (0, 1), B: (0, 1)}, {F}, st)
    assert st["ties"] == 0


def test_a_reaction_converting_two_pairs_is_dropped():
    # vf also converts a2 -> b: its node would be written by two fluxes
    edges, rmap, umap = net()
    edges = pd.concat([edges, pd.DataFrame([e("a2", "vf", "input"), e("vr", "a2", "output")])], ignore_index=True)
    st = {}
    forms, trans = m.find_pools(edges, rmap, umap, PAIRS, {}, {F}, st)
    assert forms == [] and trans == [] and st["multi_use"] >= 1


def test_no_pairs_or_empty_network_give_no_pools():
    edges, rmap, umap = net()
    assert m.find_pools(edges, rmap, umap, [], {}, set()) == ([], [])
    assert m.find_pools(pd.DataFrame(columns=["source_id", "target_id", "edge_type"]), rmap, umap, PAIRS, {}, set()) == ([], [])
