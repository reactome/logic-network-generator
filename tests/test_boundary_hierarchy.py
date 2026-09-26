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


# --- review of PR #97 -------------------------------------------------------

def test_a_root_copy_is_preferred_over_a_produced_copy(stub, monkeypatch):
    # PDGF shape: two copies of KPNB1. u_Bprod is produced by an unrelated
    # reaction; u_Broot is a root (what the root-pinning benchmark perturbs).
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")
    data, r2u = network()
    data += [{"source_id": "u_Y", "target_id": "r9", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1},
             {"source_id": "r9", "target_id": "u_Bprod", "pos_neg": "pos", "and_or": "or", "edge_type": "output", "stoichiometry": 1},
             {"source_id": "u_Broot", "target_id": "r8", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1}]
    r2u.update({"u_Y": "R-HSA-Y", "u_Bprod": B1, "u_Broot": B1})
    _emit_boundary_decomposition_edges(data, r2u)
    into_root = {s for s, t in assembly(data) if t == "u_K"}
    assert "u_Broot" in into_root and "u_Bprod" not in into_root


def test_two_joins_that_close_a_cycle_together_are_not_both_made(monkeypatch):
    # Root R1 contains produced P2; root R2 contains produced P1. R1 -> rx -> P1 and
    # R2 -> ry -> P2. Joining P2 -> R1 alone is acyclic, and so is P1 -> R2 alone,
    # but both together close R1 -> P1 -> R2 -> P2 -> R1. A pre-loop snapshot
    # allows both; the updated test must refuse the second.
    R1, R2, P1, P2 = "R-HSA-R1", "R-HSA-R2", "R-HSA-P1", "R-HSA-P2"
    labels = {R1: ["Complex"], R2: ["Complex"], P1: ["Complex"], P2: ["Complex"]}
    comps = {R1: {P2: 1}, R2: {P1: 1}, P1: {}, P2: {}}
    monkeypatch.setattr(nc, "get_labels", lambda s: labels.get(s, []))
    monkeypatch.setattr(lng, "get_labels", lambda s: labels.get(s, []), raising=False)
    monkeypatch.setattr(nc, "get_complex_components", lambda s: comps.get(s, {}))
    monkeypatch.setattr(lng, "get_terminal_components", lambda s: {s}, raising=False)
    monkeypatch.delenv("LNG_COMPOSITION_EDGES", raising=False)
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")
    e = lambda a, b, t: {"source_id": a, "target_id": b, "pos_neg": "pos", "and_or": "and", "edge_type": t, "stoichiometry": 1}
    data = [e("u_R1", "rx", "input"), e("rx", "u_P1", "output"), e("u_R2", "ry", "input"), e("ry", "u_P2", "output")]
    r2u = {"u_R1": R1, "u_R2": R2, "u_P1": P1, "u_P2": P2}
    _emit_boundary_decomposition_edges(data, r2u)
    joins = {(s, t) for s, t in assembly(data)}
    assert not {("u_P2", "u_R1"), ("u_P1", "u_R2")} <= joins      # never both


def test_a_nested_complex_is_built_per_root(stub, monkeypatch):
    # Two root occurrences of K each get their own ISGF3:KPNA1 node.
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")
    data, r2u = network()
    data.append({"source_id": "u_K2", "target_id": "r2", "pos_neg": "pos", "and_or": "and", "edge_type": "input", "stoichiometry": 1})
    r2u["u_K2"] = K
    _emit_boundary_decomposition_edges(data, r2u)
    assert len([u for u, s in r2u.items() if s == N]) == 2


# --- second review of PR #97 ------------------------------------------------

def _fixture(monkeypatch, labels, comps):
    monkeypatch.setattr(nc, "get_labels", lambda s: labels.get(s, []))
    monkeypatch.setattr(lng, "get_labels", lambda s: labels.get(s, []), raising=False)
    monkeypatch.setattr(nc, "get_complex_components", lambda s: comps.get(s, {}))
    monkeypatch.setattr(lng, "get_terminal_components",
                        lambda s: set(comps[s]) if comps.get(s) else {s}, raising=False)
    monkeypatch.delenv("LNG_COMPOSITION_EDGES", raising=False)
    monkeypatch.setenv("LNG_BOUNDARY_HIERARCHY", "1")


def _e(a, b, t):
    return {"source_id": a, "target_id": b, "pos_neg": "pos", "and_or": "and", "edge_type": t, "stoichiometry": 1}


@pytest.mark.parametrize("k1, k2", [("u_z", "u_a"), ("u_a", "u_z")])
def test_copies_of_one_root_are_ordered_by_insertion_not_uuid_string(monkeypatch, k1, k2):
    # Two root copies of K = {A, B}: K1 -> r -> a_prod, K2 -> s -> b_prod. Whichever
    # copy is processed first takes the real join (the other is then downstream of
    # it), so the order must not follow the uuid STRING, which is new every run.
    K_, A_, B_ = "R-HSA-K", "R-HSA-A", "R-HSA-B"
    _fixture(monkeypatch, {K_: ["Complex"]}, {K_: {A_: 1, B_: 1}})
    data = [_e(k1, "r", "input"), _e("r", "a_prod", "output"),
            _e(k2, "s", "input"), _e("s", "b_prod", "output")]
    r2u = {k1: K_, k2: K_, "a_prod": A_, "b_prod": B_}   # k1 inserted first
    _emit_boundary_decomposition_edges(data, r2u)
    joins = set(assembly(data))
    assert ("b_prod", k1) in joins and ("a_prod", k2) not in joins


@pytest.mark.parametrize("depleted", [False, True])
def test_a_root_copy_stays_a_root_after_it_is_decomposed(monkeypatch, depleted):
    # S has a root copy (itself a root complex, decomposed FIRST because its stid
    # sorts first, which gives it incoming joins) and a produced copy. Root R
    # contains S and must still prefer the root copy. A depletion edge into the
    # root copy does not make it produced either.
    S_, R_, L_ = "R-HSA-AAA", "R-HSA-ZZZ", "R-HSA-L"
    _fixture(monkeypatch, {S_: ["Complex"], R_: ["Complex"]}, {S_: {L_: 1}, R_: {S_: 1}})
    data = [_e("s_root", "r1", "input"), _e("u_Y", "r2", "input"), _e("r2", "s_prod", "output"),
            _e("u_R", "r3", "input")]
    if depleted:
        data.append(_e("u_W", "s_root", "depletion"))
    # s_prod is registered first, so falling back to list order would pick it.
    r2u = {"s_prod": S_, "s_root": S_, "u_R": R_, "u_Y": "R-HSA-Y", "u_W": "R-HSA-W"}
    _emit_boundary_decomposition_edges(data, r2u)
    into_r = {s for s, t in assembly(data) if t == "u_R"}
    assert into_r == {"s_root"}
