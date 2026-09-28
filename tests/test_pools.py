"""Modification-cycle pools (deltasignal specs/039 amendment 2): a protein's
forms joined by R-steps, found at NODE level; states are the least-bound form
of each modification signature, the enzyme complexes on the way are
intermediates, and a transition is a path state -> intermediates* -> state.
Shapes: RAS (GEF / intrinsic / GAP bind-hydrolyse-release), a six-form
kinase/phosphatase ring, an enzyme's own binding loop."""
import pandas as pd

import src.logic_network_generator as m

# _uuid_to_stable_id_map detects the mapping direction by uuid SHAPE
NAMES = ["g", "t", "tx", "gap", "gef", "vgef", "vhyd", "vbind", "vrel", "a2",
         "s", "se", "sxe", "sx", "sxp", "sp", "e", "p", "v1", "v2", "v3", "v4", "v5", "v6",
         "x", "y", "xy", "va", "vb", "i1", "i2", "i3", "i4", "i5", "i6", "v7", "v8", "vback",
         "xsy", "xs", "ys", "vc", "vd", "ve", "vf"]
U = {k: f"00000000-0000-0000-0000-{i:012d}" for i, k in enumerate(NAMES, 1)}
RAS, GAP_P, GEF_P, SUB, KIN, PHOS = 1, 2, 3, 10, 11, 12   # reference entities


def e(s, t, et, st=1):
    return {"source_id": U[s], "target_id": U[t], "edge_type": et, "stoichiometry": st}


def prof(proteins, sig=None, slots=1, mods=0, comps=0, fixed=None):
    return {"proteins": set(proteins), "fixed": set(fixed if fixed is not None else proteins),
            "sig": {r: frozenset(sig or ()) for r in proteins}, "mods": mods, "comps": comps, "slots": slots}


def network(edges, reactions, umap):
    return (pd.DataFrame(edges),
            pd.DataFrame({"uid": [U[k] for k in reactions], "reactome_id": [f"R-HSA-{k}" for k in reactions]}),
            {U[k]: s for k, s in umap.items()})


# --- RAS: GDP <-> GTP by GEF (enzyme), intrinsic exchange / hydrolysis (self), GAP bind + release ----------

def ras_network(hydrolysis_catalyst="t"):
    edges = [e("g", "vgef", "input"), e("vgef", "t", "output"), e("gef", "vgef", "catalyst"),
             e("t", "vhyd", "input"), e("vhyd", "g", "output"), e(hydrolysis_catalyst, "vhyd", "catalyst"),
             e("t", "vbind", "input"), e("gap", "vbind", "input"), e("vbind", "tx", "output"),
             e("tx", "vrel", "input"), e("vrel", "g", "output"), e("vrel", "gap", "output")]
    net, rmap, umap = network(edges, ["vgef", "vhyd", "vbind", "vrel"],
                              {"g": "R-HSA-G", "t": "R-HSA-T", "tx": "R-HSA-TX", "gap": "R-HSA-GAP", "gef": "R-HSA-GEF"})
    profiles = {"R-HSA-G": prof([RAS], [("mol", "GDP")], slots=1, comps=2),
                "R-HSA-T": prof([RAS], [("mol", "GTP")], slots=1, comps=2),
                "R-HSA-TX": {"proteins": {RAS, GAP_P}, "fixed": {RAS, GAP_P},
                             "sig": {RAS: frozenset([("mol", "GTP")]), GAP_P: frozenset()}, "mods": 0, "comps": 3, "slots": 2},
                "R-HSA-GAP": prof([GAP_P]), "R-HSA-GEF": prof([GEF_P])}
    participants = [("R-HSA-vgef", "input", "R-HSA-G", 1), ("R-HSA-vgef", "output", "R-HSA-T", 1),
                    ("R-HSA-vgef", "catalyst", "R-HSA-GEF", 1),
                    ("R-HSA-vhyd", "input", "R-HSA-T", 1), ("R-HSA-vhyd", "output", "R-HSA-G", 1),
                    ("R-HSA-vhyd", "catalyst", "R-HSA-T", 1),
                    ("R-HSA-vbind", "input", "R-HSA-T", 1), ("R-HSA-vbind", "input", "R-HSA-GAP", 1),
                    ("R-HSA-vbind", "output", "R-HSA-TX", 1),
                    ("R-HSA-vrel", "input", "R-HSA-TX", 1), ("R-HSA-vrel", "output", "R-HSA-G", 1),
                    ("R-HSA-vrel", "output", "R-HSA-GAP", 1)]
    return net, rmap, umap, profiles, participants


def test_r_steps_flag_self_catalysis_as_non_enzyme():
    _, _, _, profiles, participants = ras_network()
    steps = m.r_steps(participants, profiles)
    by = {(rx, r): (a, b, enz) for rx, r, a, b, enz in steps}
    assert by[("R-HSA-vgef", RAS)] == ("R-HSA-G", "R-HSA-T", True)        # GEF catalyst: enzyme
    assert by[("R-HSA-vhyd", RAS)] == ("R-HSA-T", "R-HSA-G", False)       # catalysed by RAS:GTP itself
    assert by[("R-HSA-vbind", RAS)] == ("R-HSA-T", "R-HSA-TX", True)      # joining input: the GAP
    assert by[("R-HSA-vrel", RAS)] == ("R-HSA-TX", "R-HSA-G", False)      # release: nothing joins
    assert by[("R-HSA-vbind", GAP_P)] == ("R-HSA-GAP", "R-HSA-TX", True)  # for the GAP, RAS:GTP joins
    assert ("R-HSA-vrel", GAP_P) in by


def test_ras_pool_has_gdp_base_gtp_state_and_the_gap_complex_as_intermediate():
    net, rmap, umap, profiles, participants = ras_network()
    st = {}
    forms, trans, carriers = m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles,
                                          {"R-HSA-vgef"}, {}, st)
    assert [(u, role, base) for _, u, _, role, base in forms] == [
        (U["g"], "state", True), (U["t"], "state", False), (U["tx"], "intermediate", False)]
    paths = {}
    for _, pid, step, a, b, rx, rs, enz in trans:
        paths.setdefault(pid, []).append((step, a, b, rx, rs, enz))
    shapes = sorted(tuple((a, b, rs, enz) for _, a, b, _, rs, enz in sorted(p)) for p in paths.values())
    assert shapes == sorted([
        ((U["g"], U["t"], "R-HSA-vgef", True),),
        ((U["t"], U["g"], "R-HSA-vhyd", False),),
        ((U["t"], U["tx"], "R-HSA-vbind", True), (U["tx"], U["g"], "R-HSA-vrel", False)),
    ])
    assert carriers == [("pool1", U["gap"], U["vrel"])]
    assert st["pools"] == 1 and st["states"] == 2 and st["intermediates"] == 1 and st["paths"] == 3
    assert st["multi_step_pools"] == 1 and st["carriers"] == 1 and st["ties"] == 0
    assert st["carrier_loops"] == 1            # the GAP's own GAP -> RAS:GTP:GAP -> GAP loop
    assert st["autocat_source"] == 1 and st["autocat_product"] == 0 and st["autocat_other"] == 0


def test_a_step_whose_node_loop_is_broken_is_not_a_transition():
    net, rmap, umap, profiles, participants = ras_network()
    # the release reaction returns a DIFFERENT copy of RAS:GDP
    net.loc[(net["source_id"] == U["vrel"]) & (net["target_id"] == U["g"]), "target_id"] = U["a2"]
    umap[U["a2"]] = "R-HSA-G"
    st = {}
    forms, trans, _ = m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    assert {u for _, u, _, _, _ in forms} == {U["g"], U["t"]} and st["paths"] == 2 and st["intermediates"] == 0


def test_a_form_node_that_is_a_member_of_the_steps_set_maps_to_it():
    net, rmap, umap, profiles, participants = ras_network()
    participants = [(rx, role, "R-HSA-GSET" if pe == "R-HSA-G" else pe, s) for rx, role, pe, s in participants]
    profiles["R-HSA-GSET"] = profiles["R-HSA-G"]
    steps = m.r_steps(participants, profiles)
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, {})[0] == []
    forms, _, _ = m.find_pools(net, rmap, umap, steps, profiles, set(), {"R-HSA-GSET": {"R-HSA-G"}}, {})
    assert len(forms) == 3


# --- six-form ring S -> S:E -> S*:E -> S* -> S*:P -> S:P -> S -----------------------------------------------

def ring_network():
    edges = [e("s", "v1", "input"), e("e", "v1", "input"), e("v1", "se", "output"),
             e("se", "v2", "input"), e("v2", "sxe", "output"),
             e("sxe", "v3", "input"), e("v3", "sx", "output"), e("v3", "e", "output"),
             e("sx", "v4", "input"), e("p", "v4", "input"), e("v4", "sxp", "output"),
             e("sxp", "v5", "input"), e("v5", "sp", "output"),
             e("sp", "v6", "input"), e("v6", "s", "output"), e("v6", "p", "output")]
    net, rmap, umap = network(edges, ["v1", "v2", "v3", "v4", "v5", "v6"],
                              {"s": "R-HSA-S", "se": "R-HSA-SE", "sxe": "R-HSA-SXE", "sx": "R-HSA-SX",
                               "sxp": "R-HSA-SXP", "sp": "R-HSA-SP", "e": "R-HSA-E", "p": "R-HSA-P"})
    ph = [("res", "phospho-S at 1")]
    profiles = {"R-HSA-S": prof([SUB]), "R-HSA-SX": prof([SUB], ph, mods=1),
                "R-HSA-E": prof([KIN]), "R-HSA-P": prof([PHOS]),
                "R-HSA-SE": {"proteins": {SUB, KIN}, "fixed": {SUB, KIN}, "sig": {SUB: frozenset(), KIN: frozenset()},
                             "mods": 0, "comps": 2, "slots": 2},
                "R-HSA-SXE": {"proteins": {SUB, KIN}, "fixed": {SUB, KIN}, "sig": {SUB: frozenset(ph), KIN: frozenset()},
                              "mods": 1, "comps": 2, "slots": 2},
                "R-HSA-SXP": {"proteins": {SUB, PHOS}, "fixed": {SUB, PHOS}, "sig": {SUB: frozenset(ph), PHOS: frozenset()},
                              "mods": 1, "comps": 2, "slots": 2},
                "R-HSA-SP": {"proteins": {SUB, PHOS}, "fixed": {SUB, PHOS}, "sig": {SUB: frozenset(), PHOS: frozenset()},
                             "mods": 0, "comps": 2, "slots": 2}}
    participants = [("R-HSA-v1", "input", "R-HSA-S", 1), ("R-HSA-v1", "input", "R-HSA-E", 1), ("R-HSA-v1", "output", "R-HSA-SE", 1),
                    ("R-HSA-v2", "input", "R-HSA-SE", 1), ("R-HSA-v2", "output", "R-HSA-SXE", 1),
                    ("R-HSA-v3", "input", "R-HSA-SXE", 1), ("R-HSA-v3", "output", "R-HSA-SX", 1), ("R-HSA-v3", "output", "R-HSA-E", 1),
                    ("R-HSA-v4", "input", "R-HSA-SX", 1), ("R-HSA-v4", "input", "R-HSA-P", 1), ("R-HSA-v4", "output", "R-HSA-SXP", 1),
                    ("R-HSA-v5", "input", "R-HSA-SXP", 1), ("R-HSA-v5", "output", "R-HSA-SP", 1),
                    ("R-HSA-v6", "input", "R-HSA-SP", 1), ("R-HSA-v6", "output", "R-HSA-S", 1), ("R-HSA-v6", "output", "R-HSA-P", 1)]
    return net, rmap, umap, profiles, participants


def test_six_form_ring_gives_two_states_four_intermediates_two_paths_two_carriers():
    net, rmap, umap, profiles, participants = ring_network()
    st = {}
    forms, trans, carriers = m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    roles = {u: (role, base) for _, u, _, role, base in forms}
    assert roles == {U["s"]: ("state", True), U["sx"]: ("state", False),
                     U["se"]: ("intermediate", False), U["sxe"]: ("intermediate", False),
                     U["sxp"]: ("intermediate", False), U["sp"]: ("intermediate", False)}
    paths = {}
    for _, pid, step, a, b, rx, rs, enz in trans:
        paths.setdefault(pid, []).append((step, a, b, rx, enz))
    assert len(paths) == 2 and st["paths"] == 2 and st["multi_step_pools"] == 1
    for p in paths.values():
        assert [s for s, *_ in sorted(p)] == [1, 2, 3]
        assert any(enz for *_, enz in p)                         # the binding step joins the enzyme
        assert [enz for *_, enz in sorted(p)] == [True, False, False]
    assert {(p[0][1], p[-1][2]) for p in (sorted(v) for v in paths.values())} == {(U["s"], U["sx"]), (U["sx"], U["s"])}
    assert sorted(carriers) == sorted([("pool1", U["e"], U["v3"]), ("pool1", U["p"], U["v6"])])
    assert st["carriers"] == 2 and st["carrier_loops"] == 2 and st["ties"] == 0
    assert st["nonregenerating_paths"] == 0      # the kinase and the phosphatase both come back out


def test_an_enzymes_own_binding_loop_is_not_a_pool():
    # E -> E:S -> E: one signature for E, so a carrier loop, not a pool
    edges = [e("e", "v1", "input"), e("v1", "se", "output"), e("se", "v3", "input"), e("v3", "e", "output")]
    net, rmap, umap = network(edges, ["v1", "v3"], {"e": "R-HSA-E", "se": "R-HSA-SE"})
    profiles = {"R-HSA-E": prof([KIN]), "R-HSA-SE": prof([KIN, SUB], slots=2, comps=2)}
    steps = [("R-HSA-v1", KIN, "R-HSA-E", "R-HSA-SE", True), ("R-HSA-v3", KIN, "R-HSA-SE", "R-HSA-E", False)]
    st = {}
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["carrier_loops"] == 1 and st["pools"] == 0


# --- shared nodes, path cap, orientation --------------------------------------------------------------------

def test_mutual_modification_inside_a_complex_is_not_an_enzyme_cycle_for_either():
    # X + Y -> X:Y -> X*:Y* -> X* + Y*; X* -> X and Y* -> Y directly. Before amendment
    # 5 both pools claimed the complexes and were dropped as shared; now each path
    # consumes the partner and returns it transformed, so neither regenerates.
    edges = [e("x", "va", "input"), e("y", "va", "input"), e("va", "xy", "output"),
             e("xy", "vb", "input"), e("vb", "xsy", "output"),
             e("xsy", "vc", "input"), e("vc", "xs", "output"), e("vc", "ys", "output"),
             e("xs", "vd", "input"), e("vd", "x", "output"), e("ys", "ve", "input"), e("ve", "y", "output")]
    net, rmap, umap = network(edges, ["va", "vb", "vc", "vd", "ve"],
                              {"x": "R-HSA-X", "y": "R-HSA-Y", "xy": "R-HSA-XY", "xsy": "R-HSA-XSY",
                               "xs": "R-HSA-XS", "ys": "R-HSA-YS"})
    p1, p2 = [("res", "p1")], [("res", "p2")]
    profiles = {"R-HSA-X": prof([1]), "R-HSA-Y": prof([2]), "R-HSA-XS": prof([1], p1, mods=1), "R-HSA-YS": prof([2], p2, mods=1),
                "R-HSA-XY": {"proteins": {1, 2}, "fixed": {1, 2}, "mods": 0, "comps": 2, "slots": 2,
                             "sig": {1: frozenset(), 2: frozenset()}},
                "R-HSA-XSY": {"proteins": {1, 2}, "fixed": {1, 2}, "mods": 2, "comps": 2, "slots": 2,
                              "sig": {1: frozenset(p1), 2: frozenset(p2)}}}
    steps = [("R-HSA-va", 1, "R-HSA-X", "R-HSA-XY", True), ("R-HSA-vb", 1, "R-HSA-XY", "R-HSA-XSY", False),
             ("R-HSA-vc", 1, "R-HSA-XSY", "R-HSA-XS", False), ("R-HSA-vd", 1, "R-HSA-XS", "R-HSA-X", True),
             ("R-HSA-va", 2, "R-HSA-Y", "R-HSA-XY", True), ("R-HSA-vb", 2, "R-HSA-XY", "R-HSA-XSY", False),
             ("R-HSA-vc", 2, "R-HSA-XSY", "R-HSA-YS", False), ("R-HSA-ve", 2, "R-HSA-YS", "R-HSA-Y", True)]
    st = {}
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["nonregenerating_paths"] == 2 and st["no_path_pools"] == 2 and st["shared_nodes_dropped"] == 0


def test_a_node_claimed_by_two_pools_is_dropped_and_counted():
    # X:Y <-> X*:Y* by vb / vc (steps of X's pool and of Y's); X*:Y* <-> X**:Y* by va / vd,
    # steps of X's pool only (free Y is a second Y-containing input of va, so it is
    # no Y-step); X*:Y* <-> X*:Y** by ve / vf, Y's pool only. Two different pools
    # claim X:Y, X*:Y*, vb and vc: all four are removed and both pools collapse.
    edges = [e("xy", "vb", "input"), e("vb", "xsy", "output"), e("xsy", "vc", "input"), e("vc", "xy", "output"),
             e("xsy", "va", "input"), e("y", "va", "input"), e("va", "xs", "output"), e("va", "y", "output"),
             e("xs", "vd", "input"), e("vd", "xsy", "output"),
             e("xsy", "ve", "input"), e("x", "ve", "input"), e("ve", "ys", "output"), e("ve", "x", "output"),
             e("ys", "vf", "input"), e("vf", "xsy", "output")]
    net, rmap, umap = network(edges, ["va", "vb", "vc", "vd", "ve", "vf"],
                              {"x": "R-HSA-X", "y": "R-HSA-Y", "xy": "R-HSA-XY", "xsy": "R-HSA-XSY",
                               "xs": "R-HSA-XXY", "ys": "R-HSA-XYY"})
    p1, p2, pp1, pp2 = [("res", "p1")], [("res", "p2")], [("res", "p1"), ("res", "q1")], [("res", "p2"), ("res", "q2")]

    def dimer(s1, s2, mods):
        return {"proteins": {1, 2}, "fixed": {1, 2}, "mods": mods, "comps": 2, "slots": 2,
                "sig": {1: frozenset(s1), 2: frozenset(s2)}}
    profiles = {"R-HSA-X": prof([1]), "R-HSA-Y": prof([2]), "R-HSA-XY": dimer([], [], 0),
                "R-HSA-XSY": dimer(p1, p2, 2), "R-HSA-XXY": dimer(pp1, p2, 3), "R-HSA-XYY": dimer(p1, pp2, 3)}
    participants = [("R-HSA-vb", "input", "R-HSA-XY", 1), ("R-HSA-vb", "output", "R-HSA-XSY", 1),
                    ("R-HSA-vc", "input", "R-HSA-XSY", 1), ("R-HSA-vc", "output", "R-HSA-XY", 1),
                    ("R-HSA-va", "input", "R-HSA-XSY", 1), ("R-HSA-va", "input", "R-HSA-Y", 1),
                    ("R-HSA-va", "output", "R-HSA-XXY", 1), ("R-HSA-va", "output", "R-HSA-Y", 1),
                    ("R-HSA-vd", "input", "R-HSA-XXY", 1), ("R-HSA-vd", "output", "R-HSA-XSY", 1),
                    ("R-HSA-ve", "input", "R-HSA-XSY", 1), ("R-HSA-ve", "input", "R-HSA-X", 1),
                    ("R-HSA-ve", "output", "R-HSA-XYY", 1), ("R-HSA-ve", "output", "R-HSA-X", 1),
                    ("R-HSA-vf", "input", "R-HSA-XYY", 1), ("R-HSA-vf", "output", "R-HSA-XSY", 1)]
    steps = m.r_steps(participants, profiles)
    # the returns vd and vf are steps for both proteins, but only one direction of
    # each side branch is, so the branch is outside the other protein's cycle
    assert {(rx, r) for rx, r, *_ in steps} == {("R-HSA-vb", 1), ("R-HSA-vc", 1), ("R-HSA-va", 1), ("R-HSA-vd", 1),
                                                ("R-HSA-vb", 2), ("R-HSA-vc", 2), ("R-HSA-ve", 2), ("R-HSA-vf", 2),
                                                ("R-HSA-vd", 2), ("R-HSA-vf", 1)}
    for keep in (1, 2):
        st = {}
        m.find_pools(net, rmap, umap, [t for t in steps if t[1] == keep], profiles, set(), {}, st)
        assert st["pools"] == 1 and st["states"] == 3
    st = {}
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["shared_nodes_dropped"] == 4 and st["pools"] == 0


def test_forms_differing_only_by_a_bound_partner_are_not_a_pool():
    # Activated FGFR4 <-> FGFR4:PLCG1 <-> FGFR4:p-PLCG1 (R-HSA-5654743): the PLCG1
    # phosphorylation gives FGFR4 no new signature; the only core-only form is
    # FGFR4 itself, so this is a carrier loop, not a pool (amendment 3)
    edges = [e("x", "va", "input"), e("y", "va", "input"), e("va", "xy", "output"),
             e("xy", "vb", "input"), e("vb", "xsy", "output"),
             e("xsy", "vc", "input"), e("vc", "x", "output"), e("vc", "ys", "output")]
    net, rmap, umap = network(edges, ["va", "vb", "vc"],
                              {"x": "R-HSA-X", "y": "R-HSA-Y", "xy": "R-HSA-XY", "xsy": "R-HSA-XSY", "ys": "R-HSA-YS"})
    p2 = [("res", "p2")]
    profiles = {"R-HSA-X": prof([1]), "R-HSA-Y": prof([2]), "R-HSA-YS": prof([2], p2, mods=1),
                "R-HSA-XY": {"proteins": {1, 2}, "fixed": {1, 2}, "mods": 0, "comps": 2, "slots": 2,
                             "sig": {1: frozenset(), 2: frozenset()}},
                "R-HSA-XSY": {"proteins": {1, 2}, "fixed": {1, 2}, "mods": 1, "comps": 2, "slots": 2,
                              "sig": {1: frozenset(), 2: frozenset(p2)}}}
    steps = [("R-HSA-va", 1, "R-HSA-X", "R-HSA-XY", True), ("R-HSA-vb", 1, "R-HSA-XY", "R-HSA-XSY", False),
             ("R-HSA-vc", 1, "R-HSA-XSY", "R-HSA-X", False)]
    st = {}
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["carrier_loops"] == 1
    # and the same shape where the RECEPTOR is what gets phosphorylated, with the
    # partner released unchanged (X:Y -> X*:Y -> X* + Y; X* -> X), IS a pool
    profiles["R-HSA-XSY"]["sig"] = {1: frozenset(p2), 2: frozenset()}
    profiles["R-HSA-XS"] = prof([1], p2, mods=1)
    umap[U["ys"]] = "R-HSA-XS"
    edges = [e("x", "va", "input"), e("y", "va", "input"), e("va", "xy", "output"),
             e("xy", "vb", "input"), e("vb", "xsy", "output"),
             e("xsy", "vc", "input"), e("vc", "ys", "output"), e("vc", "y", "output"),
             e("ys", "vd", "input"), e("vd", "x", "output")]
    net, rmap, _ = network(edges, ["va", "vb", "vc", "vd"], {})
    steps = [("R-HSA-va", 1, "R-HSA-X", "R-HSA-XY", True), ("R-HSA-vb", 1, "R-HSA-XY", "R-HSA-XSY", False),
             ("R-HSA-vc", 1, "R-HSA-XSY", "R-HSA-XS", False), ("R-HSA-vd", 1, "R-HSA-XS", "R-HSA-X", True)]
    st = {}
    forms, _, _ = m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st)
    assert st["pools"] == 1 and {u: r for _, u, _, r, _ in forms} == {
        U["x"]: "state", U["ys"]: "state", U["xy"]: "intermediate", U["xsy"]: "intermediate"}


def test_a_carrier_that_another_pool_manages_is_excluded_and_counted():
    # the ring's kinase E is also a form of its own pool (E <-> E* by an outside
    # kinase and phosphatase): E is not a carrier of the substrate's pool
    net, rmap, umap, profiles, participants = ring_network()
    edges = [e("e", "v7", "input"), e("v7", "i1", "output"), e("i1", "v8", "input"), e("v8", "e", "output")]
    net = pd.concat([net, pd.DataFrame(edges)], ignore_index=True)
    rmap = pd.concat([rmap, pd.DataFrame({"uid": [U["v7"], U["v8"]], "reactome_id": ["R-HSA-v7", "R-HSA-v8"]})],
                     ignore_index=True)
    umap[U["i1"]] = "R-HSA-ES"
    profiles["R-HSA-ES"] = prof([KIN], [("res", "pE")], mods=1)
    participants += [("R-HSA-v7", "input", "R-HSA-E", 1), ("R-HSA-v7", "output", "R-HSA-ES", 1),
                     ("R-HSA-v7", "catalyst", "R-HSA-P", 1),
                     ("R-HSA-v8", "input", "R-HSA-ES", 1), ("R-HSA-v8", "output", "R-HSA-E", 1),
                     ("R-HSA-v8", "catalyst", "R-HSA-P", 1)]
    st = {}
    forms, _, carriers = m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    assert st["pools"] == 2 and st["shared_nodes_dropped"] == 0
    assert carriers == [(next(p for p, u, *_ in forms if u == U["s"]), U["p"], U["v6"])]
    assert st["carrier_conflicts"] == 1 and st["carriers"] == 1


def test_co_travelling_proteins_of_one_set_are_one_pool():
    # p21 RAS:GDP <-> RAS:GTP: the H/K/NRAS graphs are identical, so one pool, counted as merged
    net, rmap, umap, profiles, participants = ras_network()
    for s in ("R-HSA-G", "R-HSA-T"):
        profiles[s] = prof([RAS, 4], profiles[s]["sig"][RAS], comps=2)
    profiles["R-HSA-TX"]["proteins"].add(4)
    profiles["R-HSA-TX"]["sig"][4] = profiles["R-HSA-TX"]["sig"][RAS]
    st = {}
    forms, _, _ = m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    assert st["pools"] == 1 and st["merged_proteins"] == 1 and st["shared_nodes_dropped"] == 0 and len(forms) == 3


def chain_network(n_intermediates):
    inter = ["i1", "i2", "i3", "i4", "i5", "i6"][:n_intermediates]
    seq = ["s"] + inter + ["sx"]
    rxs = ["v1", "v2", "v3", "v4", "v5", "v6", "v7"][:len(seq) - 1]
    edges, umap, steps = [], {"s": "R-HSA-S", "sx": "R-HSA-SX"}, []
    profiles = {"R-HSA-S": prof([SUB]), "R-HSA-SX": prof([SUB], [("res", "p")], mods=1)}
    for k, rx in enumerate(rxs):
        a, b = seq[k], seq[k + 1]
        edges += [e(a, rx, "input"), e(rx, b, "output")]
        for x in (a, b):
            if x in inter:
                umap[x] = f"R-HSA-{x.upper()}"
                profiles[f"R-HSA-{x.upper()}"] = prof([SUB], slots=2, comps=2)
        steps.append((f"R-HSA-{rx}", SUB, umap[a], umap[b], True))
    edges += [e("sx", "vback", "input"), e("vback", "s", "output")]
    steps.append(("R-HSA-vback", SUB, "R-HSA-SX", "R-HSA-S", True))
    net, rmap, umap = network(edges, rxs + ["vback"], umap)
    return net, rmap, umap, steps, profiles


def test_paths_over_six_steps_are_dropped_and_counted():
    net, rmap, umap, steps, profiles = chain_network(5)     # 6 steps: kept
    st = {}
    m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st)
    assert st["paths"] == 2 and st["long_paths_dropped"] == 0 and st["intermediates"] == 5
    net, rmap, umap, steps, profiles = chain_network(6)     # 7 steps: dropped, and with it the pool,
    st = {}                                                  # which is left with no way back to S*
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["pools"] == 0 and st["long_paths_dropped"] == 1 and st["no_path_pools"] == 1
    # a second reaction for one of its steps is a second dropped PATH, not a prefix
    net = pd.concat([net, pd.DataFrame([e("i3", "v8", "input"), e("v8", "i4", "output")])], ignore_index=True)
    rmap = pd.concat([rmap, pd.DataFrame({"uid": [U["v8"]], "reactome_id": ["R-HSA-v8"]})], ignore_index=True)
    steps.append(("R-HSA-v8", SUB, "R-HSA-I3", "R-HSA-I4", True))
    st = {}
    m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st)
    assert st["pools"] == 0 and st["long_paths_dropped"] == 2


def test_orientation_by_residues_then_donor_then_components_and_ties_are_counted():
    net, rmap, umap, profiles, participants = ras_network()
    steps = m.r_steps(participants, profiles)

    def base_of(profiles, donors, st=None):
        forms, _, _ = m.find_pools(net, rmap, umap, steps, profiles, donors, {}, st if st is not None else {})
        return next(u for _, u, _, role, base in forms if base)
    # no residue difference: the donor is consumed on the way to GTP, so GDP is base
    assert base_of(profiles, {"R-HSA-vgef"}) == U["g"]
    assert base_of(profiles, {"R-HSA-vhyd"}) == U["t"]
    # residues decide first, whatever the donor says
    profiles["R-HSA-G"]["mods"] = 2
    assert base_of(profiles, {"R-HSA-vgef"}) == U["t"]
    profiles["R-HSA-G"]["mods"] = 0
    # no residue or donor difference: the form with more components is modified
    profiles["R-HSA-T"]["comps"] = 3
    assert base_of(profiles, set()) == U["g"]
    profiles["R-HSA-T"]["comps"] = 2
    # nothing decides: smaller stId, counted as a tie
    st = {}
    assert base_of(profiles, set(), st) == U["g"] and st["ties"] == 1


def test_no_steps_or_empty_network_give_no_pools():
    net, rmap, umap, profiles, participants = ras_network()
    assert m.find_pools(net, rmap, umap, [], profiles, set()) == ([], [], [])
    empty = pd.DataFrame(columns=["source_id", "target_id", "edge_type"])
    assert m.find_pools(empty, rmap, umap, m.r_steps(participants, profiles), profiles, set()) == ([], [], [])


def test_a_stoichiometry_two_input_is_not_a_step():
    _, _, _, profiles, participants = ras_network()
    participants = [(rx, role, pe, 2 if (rx, role) == ("R-HSA-vgef", "input") else s) for rx, role, pe, s in participants]
    assert ("R-HSA-vgef", RAS) not in {(rx, r) for rx, r, *_ in m.r_steps(participants, profiles)}


def test_copies_of_one_reaction_collapse_into_one_path_with_one_row_per_copy():
    # RAS with the GAP binding curated for two GAP variants (two reaction NODES of
    # one stId) and both release copies, plus an intrinsic exchange beside the GEF:
    # copies of one stId are one step (k rows); different reactions stay apart.
    net, rmap, umap, profiles, participants = ras_network()
    extra = [e("t", "v7", "input"), e("a2", "v7", "input"), e("v7", "tx", "output"),
             e("tx", "v8", "input"), e("v8", "g", "output"), e("v8", "a2", "output"),
             e("g", "vback", "input"), e("vback", "t", "output")]
    net = pd.concat([net, pd.DataFrame(extra)], ignore_index=True)
    rmap = pd.concat([rmap, pd.DataFrame({"uid": [U["v7"], U["v8"], U["vback"]],
                                          "reactome_id": ["R-HSA-vbind", "R-HSA-vrel", "R-HSA-vint"]})], ignore_index=True)
    umap[U["a2"]] = "R-HSA-GAP2"
    profiles["R-HSA-GAP2"] = prof([GAP_P])
    participants += [("R-HSA-vint", "input", "R-HSA-G", 1), ("R-HSA-vint", "output", "R-HSA-T", 1)]
    st = {}
    _, trans, carriers = m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    rows = {}
    for _, pid, step, a, b, rx, rs, enz in trans:
        rows.setdefault((pid, step), []).append((a, b, rx, rs, enz))
    by_shape = {}
    for (pid, step), rs_ in rows.items():
        assert len({(a, b, rs, enz) for a, b, _, rs, enz in rs_}) == 1     # one source, target, stId, flag per step
        a, b, _, rs, enz = rs_[0]
        by_shape.setdefault(pid, []).append((step, a, b, rs, enz, sorted(rx for _, _, rx, _, _ in rs_)))
    shapes = sorted(tuple(x[1:] for x in sorted(v)) for v in by_shape.values())
    assert shapes == sorted([
        ((U["g"], U["t"], "R-HSA-vgef", True, [U["vgef"]]),),
        ((U["g"], U["t"], "R-HSA-vint", False, [U["vback"]]),),           # not merged into the GEF step
        ((U["t"], U["g"], "R-HSA-vhyd", False, [U["vhyd"]]),),
        ((U["t"], U["tx"], "R-HSA-vbind", True, sorted([U["vbind"], U["v7"]])),
         (U["tx"], U["g"], "R-HSA-vrel", False, sorted([U["vrel"], U["v8"]]))),
    ])
    assert st["paths"] == 4 and st["copy_rows"] == 7 and st["multi_step_pools"] == 1
    assert sorted(carriers) == sorted([("pool1", U["gap"], U["vrel"]), ("pool1", U["a2"], U["v8"])])


# --- amendment 5: the enzyme must be regenerated ---------------------------------------------------------

def test_a_machine_path_whose_co_input_is_not_returned_is_dropped():
    # S + M -> S:M -> S* + M' (the partner leaves transformed); S* -> S directly.
    # The multi-step path is not an enzyme cycle, and S* then has no way back.
    edges = [e("s", "v1", "input"), e("e", "v1", "input"), e("v1", "se", "output"),
             e("se", "v2", "input"), e("v2", "sx", "output"), e("v2", "p", "output"),
             e("sx", "v3", "input"), e("v3", "s", "output")]
    net, rmap, umap = network(edges, ["v1", "v2", "v3"],
                              {"s": "R-HSA-S", "e": "R-HSA-M", "se": "R-HSA-SM", "sx": "R-HSA-SX", "p": "R-HSA-M2"})
    ph = [("res", "p")]
    profiles = {"R-HSA-S": prof([SUB]), "R-HSA-SX": prof([SUB], ph, mods=1), "R-HSA-M": prof([KIN]),
                "R-HSA-M2": prof([KIN], [("res", "m")], mods=1),
                "R-HSA-SM": {"proteins": {SUB, KIN}, "fixed": {SUB, KIN}, "sig": {SUB: frozenset(), KIN: frozenset()},
                             "mods": 0, "comps": 2, "slots": 2}}
    steps = [("R-HSA-v1", SUB, "R-HSA-S", "R-HSA-SM", True), ("R-HSA-v2", SUB, "R-HSA-SM", "R-HSA-SX", False),
             ("R-HSA-v3", SUB, "R-HSA-SX", "R-HSA-S", True)]
    st = {}
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["nonregenerating_paths"] == 1 and st["no_path_pools"] == 1
    # returned unchanged, the same shape is a pool
    net.loc[(net["source_id"] == U["v2"]) & (net["target_id"] == U["p"]), "target_id"] = U["e"]
    st = {}
    forms, _, carriers = m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st)
    assert st["pools"] == 1 and st["nonregenerating_paths"] == 0 and carriers == [("pool1", U["e"], U["v2"])]


def test_an_hdr_like_machine_is_not_a_pool():
    # RPA + DSB machine -> machine:RPA -> machine':RPA (the machine grows: BRCA2
    # joins) -> p-RPA + machine'' (with RAD51); p-RPA -> RPA by a phosphatase.
    # The machine that entered never comes back out: RPA is a passenger.
    edges = [e("s", "v1", "input"), e("i1", "v1", "input"), e("v1", "i2", "output"),
             e("i2", "v2", "input"), e("i3", "v2", "input"), e("v2", "i4", "output"),
             e("i4", "v3", "input"), e("v3", "sx", "output"), e("v3", "i5", "output"),
             e("sx", "v4", "input"), e("v4", "s", "output"), e("p", "v4", "catalyst")]
    net, rmap, umap = network(edges, ["v1", "v2", "v3", "v4"],
                              {"s": "R-HSA-RPA", "sx": "R-HSA-PRPA", "i1": "R-HSA-DSB", "i2": "R-HSA-DSB-RPA",
                               "i3": "R-HSA-BRCA2", "i4": "R-HSA-DSB-RPA-BRCA2", "i5": "R-HSA-DSB-RAD51", "p": "R-HSA-P"})
    ph = [("res", "pS")]

    def machine(proteins, sig_r):
        return {"proteins": set(proteins), "fixed": set(proteins), "mods": 0, "comps": len(proteins), "slots": len(proteins),
                "sig": {r: frozenset(sig_r if r == SUB else ()) for r in proteins}}
    profiles = {"R-HSA-RPA": prof([SUB]), "R-HSA-PRPA": prof([SUB], ph, mods=1), "R-HSA-P": prof([PHOS]),
                "R-HSA-DSB": machine([20, 21], []), "R-HSA-BRCA2": prof([22]),
                "R-HSA-DSB-RPA": machine([20, 21, SUB], []), "R-HSA-DSB-RPA-BRCA2": machine([20, 21, 22, SUB], []),
                "R-HSA-DSB-RAD51": machine([20, 21, 22, 23], [])}
    steps = [("R-HSA-v1", SUB, "R-HSA-RPA", "R-HSA-DSB-RPA", True),
             ("R-HSA-v2", SUB, "R-HSA-DSB-RPA", "R-HSA-DSB-RPA-BRCA2", True),
             ("R-HSA-v3", SUB, "R-HSA-DSB-RPA-BRCA2", "R-HSA-PRPA", False),
             ("R-HSA-v4", SUB, "R-HSA-PRPA", "R-HSA-RPA", True)]
    st = {}
    assert m.find_pools(net, rmap, umap, steps, profiles, set(), {}, st) == ([], [], [])
    assert st["nonregenerating_paths"] == 1 and st["no_path_pools"] == 1 and st["pools"] == 0


def test_small_molecules_are_exempt_from_regeneration():
    # S + E -> S:E; S:E + ATP -> S*:E + ADP; S*:E -> S* + E; S* -> S. ATP is a
    # cofactor consumed and not returned; the kinase is, so the path is kept.
    net, rmap, umap, profiles, participants = ring_network()
    net = pd.concat([net, pd.DataFrame([e("i1", "v2", "input"), e("v2", "i2", "output")])], ignore_index=True)
    umap[U["i1"]], umap[U["i2"]] = "R-ALL-ATP", "R-ALL-ADP"
    profiles["R-ALL-ATP"] = {"proteins": set(), "fixed": set(), "sig": {}, "mods": 0, "comps": 0, "slots": 0, "small": True}
    profiles["R-ALL-ADP"] = dict(profiles["R-ALL-ATP"])
    st = {}
    m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    assert st["pools"] == 1 and st["nonregenerating_paths"] == 0 and st["paths"] == 2
    profiles["R-ALL-ATP"]["small"] = False      # were ATP a partner, the path would fail
    st = {}
    m.find_pools(net, rmap, umap, m.r_steps(participants, profiles), profiles, set(), {}, st)
    assert st["nonregenerating_paths"] == 1 and st["pools"] == 0
