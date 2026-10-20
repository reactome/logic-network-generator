"""Recycle pools (deltasignal specs/081): a modified protein's return to the complex
it starts in. Shape: p130 (RBL2) starts in two copies of RBL2:E2F4:DP (cx1, cx2),
built by assembly from one starting copy s; a kinase moves it to p-RBL2 (m); PP2A
makes free RBL2 (f), which leads nowhere. The pool: states cx1, cx2 (base) and m,
paths cx -> m through each copy's kinase reaction and m -> cx through PP2A, and a
recycle edge f -> s."""
import pandas as pd

import src.logic_network_generator as m

NAMES = ["s", "e2f", "cx1", "cx2", "kin", "vk1", "vk2", "pm", "pp", "vp", "f", "w", "vw", "y", "vy", "q"]
U = {k: f"00000000-0000-0000-0000-{i:012d}" for i, k in enumerate(NAMES, 1)}
RBL2, E2F, KIN, PP, W = 1, 2, 3, 4, 5


def e(s, t, et):
    return {"source_id": U[s], "target_id": U[t], "edge_type": et, "stoichiometry": 1}


def prof(proteins, sig=()):
    return {"proteins": set(proteins), "fixed": set(proteins), "sig": {r: frozenset(sig) for r in proteins},
            "mods": 0, "comps": 0, "slots": 1}


def p130(extra=(), m_sig=(("res", "p-S672"),)):
    edges = [e("s", "cx1", "assembly"), e("s", "cx2", "assembly"), e("e2f", "cx1", "assembly"), e("e2f", "cx2", "assembly"),
             e("cx1", "vk1", "input"), e("kin", "vk1", "catalyst"), e("vk1", "pm", "output"),
             e("cx2", "vk2", "input"), e("kin", "vk2", "catalyst"), e("vk2", "pm", "output"),
             e("pm", "vp", "input"), e("pp", "vp", "catalyst"), e("vp", "f", "output")] + list(extra)
    net = pd.DataFrame(edges)
    rmap = pd.DataFrame({"uid": [U[k] for k in ("vk1", "vk2", "vp", "vw", "vy")],
                         "reactome_id": ["R-HSA-K", "R-HSA-K", "R-HSA-P", "R-HSA-W", "R-HSA-Y"]})
    umap = {U["s"]: "R-HSA-RBL2", U["f"]: "R-HSA-RBL2", U["e2f"]: "R-HSA-E2F",
            U["cx1"]: "R-HSA-CX1", U["cx2"]: "R-HSA-CX2", U["pm"]: "R-HSA-PM",
            U["kin"]: "R-HSA-KIN", U["pp"]: "R-HSA-PP", U["w"]: "R-HSA-W", U["y"]: "R-HSA-Y", U["q"]: "R-HSA-Q"}
    profiles = {"R-HSA-RBL2": prof([RBL2]), "R-HSA-E2F": prof([E2F]),
                "R-HSA-CX1": prof([RBL2, E2F]), "R-HSA-CX2": prof([RBL2, E2F]),
                "R-HSA-PM": prof([RBL2], m_sig), "R-HSA-KIN": prof([KIN]), "R-HSA-PP": prof([PP]),
                "R-HSA-W": prof([W]), "R-HSA-Y": prof([RBL2, W]), "R-HSA-Q": prof([RBL2])}
    return net, rmap, umap, profiles


def test_p130_returns_to_both_copies_of_its_starting_complex():
    st = {}
    forms, trans, edges = m.find_recycle_pools(*p130(), stats=st)
    assert [(u, role, base) for _, u, _, role, base in forms] == [
        (U["cx1"], "state", True), (U["cx2"], "state", True), (U["pm"], "state", False)]
    assert {f[0] for f in forms} == {"recycle1"}
    rows = sorted((t[1], t[2], t[3], t[4], t[5], t[7]) for t in trans)
    assert rows == sorted([("recycle1_f1", 1, U["cx1"], U["pm"], U["vk1"], True),
                           ("recycle1_r1", 1, U["pm"], U["cx1"], U["vp"], True),
                           ("recycle1_f2", 1, U["cx2"], U["pm"], U["vk2"], True),
                           ("recycle1_r2", 1, U["pm"], U["cx2"], U["vp"], True)])
    assert edges == [(U["f"], U["s"])]
    assert st["recycle_pools"] == 1 and st["recycle_states"] == 3 and st["recycle_paths"] == 4
    assert st["recycle_edges"] == 1 and st["recycle_claimed"] == 0


def test_a_partner_released_unchanged_is_not_a_recycle_pool():
    # the same route with no modification of RBL2 on the way: a carrier release
    assert m.find_recycle_pools(*p130(m_sig=())) == ([], [], [])


def test_a_free_copy_that_leads_somewhere_is_not_a_candidate():
    net, rmap, umap, profiles = p130(extra=[e("f", "vw", "input"), e("vw", "w", "output")])
    assert m.find_recycle_pools(net, rmap, umap, profiles) == ([], [], [])


def test_a_candidate_on_a_claimed_node_is_skipped_and_counted():
    st = {}
    assert m.find_recycle_pools(*p130(), claimed={U["pm"]}, stats=st) == ([], [], [])
    assert st["recycle_claimed"] == 1


def test_a_two_protein_free_copy_is_not_a_candidate():
    # a complex released at the end (RBL2:W) is not a single protein's return
    net, rmap, umap, profiles = p130()
    umap = dict(umap); umap[U["f"]] = "R-HSA-Y"
    assert m.find_recycle_pools(net, rmap, umap, profiles) == ([], [], [])


def test_the_flag_is_off_by_default():
    from src.env_flags import BOOL_FLAGS
    assert BOOL_FLAGS["LNG_RECYCLE_POOLS"] is False
