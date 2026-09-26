"""LNG_BOUNDARY_HIERARCHY=1 (deltasignal specs/030): a root complex is decomposed one
hasComponent level at a time. A component the network already has is JOINED; a nested
complex it lacks is BUILT; the rest fall back to terminal leaves.

Shape: IFN alpha/beta. ISGF3 [cytosol] is produced upstream and consumed by nothing; the
root complex ISGF3:KPNA1:KPNB1 (produced by nothing) is translocated. Reactome has no
binding reaction between them, so flat decomposition (straight to STAT1/STAT2/IRF9/KPNA1/
KPNB1 leaves) left every perturbation upstream of ISGF3 severed at that point."""
import pytest
import src.neo4j_connector as nc
import src.logic_network_generator as lng
from src.logic_network_generator import _emit_boundary_decomposition_edges

K, N, I = "R-HSA-K", "R-HSA-N", "R-HSA-I"            # K = ISGF3:KPNA1:KPNB1, N = ISGF3:KPNA1, I = ISGF3
A1, B1, S1, S2, IR = "R-HSA-KPNA1", "R-HSA-KPNB1", "R-HSA-STAT1", "R-HSA-STAT2", "R-HSA-IRF9"
LABELS = {K: ["Complex"], N: ["Complex"], I: ["Complex"],
          **{x: ["EntityWithAccessionedSequence"] for x in (A1, B1, S1, S2, IR)}}
COMPONENTS = {K: {N: 1, B1: 1}, N: {I: 1, A1: 1}, I: {S1: 1, S2: 1, IR: 1}}
LEAVES = {K: {S1, S2, IR, A1, B1}, N: {S1, S2, IR, A1}, I: {S1, S2, IR}}


@pytest.fixture
def stub(monkeypatch):
    monkeypatch.setattr(nc, "get_labels", lambda s: LABELS.get(s, []))
    monkeypatch.setattr(lng, "get_labels", lambda s: LABELS.get(s, []), raising=False)
    monkeypatch.setattr(nc, "get_complex_components", lambda s: COMPONENTS.get(s, {}))
    monkeypatch.setattr(lng, "get_terminal_components", lambda s: LEAVES.get(s, {s}), raising=False)
    monkeypatch.delenv("LNG_COMPOSITION_EDGES", raising=False)


def network(isgf3_downstream_of_root=False):
    # u_X -> r0 -> u_I (ISGF3 produced, dead end);  u_K (root) -> r1 -> u_Z
    data = [
        {"source_id": "u_X", "target_id": "r0", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r0", "target_id": "u_I", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
        {"source_id": "u_K", "target_id": "r1", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        {"source_id": "r1", "target_id": "u_Z", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
    ]
    if isgf3_downstream_of_root:     # ISGF3 is produced FROM the root complex: reusing it would weld a cycle
        data[0] = {"source_id": "u_Z", "target_id": "r0", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1}
    return data, {"u_X": "R-HSA-X", "u_I": I, "u_K": K, "u_Z": "R-HSA-Z"}


def assembly(data):
    return [(e["source_id"], e["target_id"]) for e in data if e["edge_type"] == "assembly"]


def test_flat_mode_is_leaves_as_before(stub, monkeypatch):
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "0")
    data, r2u = network()
    _emit_boundary_decomposition_edges(data, r2u)
    edges = assembly(data)
    assert {t for _, t in edges} == {"u_K"}
    assert {r2u[s] for s, _ in edges} == {S1, S2, IR, A1, B1}
    assert ("u_I", "u_K") not in edges and all(s != "u_I" for s, _ in edges)


def test_hierarchy_builds_the_nested_complex_and_joins_the_produced_species(stub, monkeypatch):
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")
    data, r2u = network()
    _emit_boundary_decomposition_edges(data, r2u)
    edges = assembly(data)
    nested = [u for u, s in r2u.items() if s == N]
    assert len(nested) == 1                                    # ISGF3:KPNA1 built once
    n = nested[0]
    assert ("u_I", n) in edges                                 # produced ISGF3 joins it ...
    assert {r2u[s] for s, t in edges if t == n} == {I, A1}     # ... with KPNA1
    assert {r2u[s] for s, t in edges if t == "u_K"} == {N, B1} # and it joins the root with KPNB1
    assert not any(t == "u_I" for _, t in edges)               # ISGF3 is not descended into


def test_a_produced_copy_downstream_of_the_root_is_not_reused(stub, monkeypatch):
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")
    data, r2u = network(isgf3_downstream_of_root=True)
    _emit_boundary_decomposition_edges(data, r2u)
    edges = assembly(data)
    assert all(s != "u_I" for s, _ in edges)                   # no weld through the root's own output
    built_I = [u for u, s in r2u.items() if s == I and u != "u_I"]
    assert len(built_I) == 1                                   # a separate ISGF3 node is built instead
    assert {r2u[s] for s, t in edges if t == built_I[0]} == {S1, S2, IR}


def test_hierarchy_is_the_default(stub, monkeypatch):
    monkeypatch.delenv("LNG_BOUNDARY_HIERARCHY", raising=False)
    data, r2u = network()
    _emit_boundary_decomposition_edges(data, r2u)
    assert ("u_I", next(u for u, s in r2u.items() if s == N)) in assembly(data)
