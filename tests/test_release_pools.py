"""Release pools (deltasignal specs/086): a protein bound, modified and released by a
complex, whose released copy carries on. Shape: beta-catenin (s) binds the destruction
complex (dc) into c1; a kinase phosphorylates it (c2); a phosphatase dephosphorylates it
(c3); the release (v4) makes p, which goes on to the nucleus (vt -> n). The pool: c3 (base)
and s, intermediates c1 and c2; paths s -> c1 -> c2 -> c3 and c3 -> s; p is the carrier
with release v4; edges s -> p (released_twin) and p -> s (recycle)."""
import pandas as pd

import src.logic_network_generator as m

NAMES = ["s", "dc", "v1", "c1", "kin", "v2", "c2", "pp", "v3", "c3", "v4", "p", "vt", "n"]
U = {k: f"00000000-0000-0000-0000-{i:012d}" for i, k in enumerate(NAMES, 1)}
CTNNB1, DC, KIN, PP = 1, 2, 3, 4


def e(s, t, et):
    return {"source_id": U[s], "target_id": U[t], "edge_type": et, "stoichiometry": 1}


def prof(proteins, sig=()):
    return {"proteins": set(proteins), "fixed": set(proteins), "sig": {r: frozenset(sig) for r in proteins},
            "mods": 0, "comps": 0, "slots": 1}


def bcat(modified=True, carries_on=True):
    edges = [e("s", "v1", "input"), e("dc", "v1", "input"), e("v1", "c1", "output"),
             e("c1", "v2", "input"), e("kin", "v2", "catalyst"), e("v2", "c2", "output"),
             e("c2", "v3", "input"), e("pp", "v3", "catalyst"), e("v3", "c3", "output"),
             e("c3", "v4", "input"), e("v4", "p", "output")]
    if carries_on:
        edges += [e("p", "vt", "input"), e("vt", "n", "output")]
    net = pd.DataFrame(edges)
    rx = ["v1", "v2", "v3", "v4", "vt"]
    rmap = pd.DataFrame({"uid": [U[k] for k in rx], "reactome_id": [f"R-HSA-{k.upper()}" for k in rx]})
    umap = {U["s"]: "R-HSA-CTNNB1", U["p"]: "R-HSA-CTNNB1", U["dc"]: "R-HSA-DC", U["c1"]: "R-HSA-C1",
            U["c2"]: "R-HSA-C2", U["c3"]: "R-HSA-C3", U["kin"]: "R-HSA-KIN", U["pp"]: "R-HSA-PP", U["n"]: "R-HSA-N"}
    phos = (("res", "p-S45"),) if modified else ()
    profiles = {"R-HSA-CTNNB1": prof([CTNNB1]), "R-HSA-DC": prof([DC]),
                "R-HSA-C1": prof([CTNNB1, DC]), "R-HSA-C2": prof([CTNNB1, DC], phos),
                "R-HSA-C3": prof([CTNNB1, DC]), "R-HSA-KIN": prof([KIN]), "R-HSA-PP": prof([PP]),
                "R-HSA-N": prof([CTNNB1])}
    return net, rmap, umap, profiles


def test_beta_catenin_is_a_free_bound_pool_whose_released_copy_reads_the_free_state():
    st = {}
    forms, trans, carriers, edges = m.find_release_pools(*bcat(), stats=st)
    assert [(u, role, base) for _, u, _, role, base in forms] == [
        (U["c3"], "state", True), (U["s"], "state", False),
        (U["c1"], "intermediate", False), (U["c2"], "intermediate", False)]
    rows = sorted((t[1], t[2], t[3], t[4], t[5]) for t in trans)
    assert rows == sorted([("release1_f1", 1, U["s"], U["c1"], U["v1"]), ("release1_f1", 2, U["c1"], U["c2"], U["v2"]),
                           ("release1_f1", 3, U["c2"], U["c3"], U["v3"]), ("release1_r1", 1, U["c3"], U["s"], U["v4"])])
    assert [(c[1], c[2]) for c in carriers] == [(U["p"], U["v4"])]
    assert edges == [(U["s"], U["p"], "released_twin"), (U["p"], U["s"], "recycle")]
    assert st["release_pools"] == 1 and st["release_paths"] == 2 and st["release_claimed"] == 0


def test_a_partner_released_unchanged_is_not_a_release_pool():
    assert m.find_release_pools(*bcat(modified=False)) == ([], [], [], [])


def test_a_released_copy_that_leads_nowhere_is_left_to_recycle_pools():
    assert m.find_release_pools(*bcat(carries_on=False)) == ([], [], [], [])


def test_a_candidate_on_a_claimed_node_is_skipped_and_counted():
    st = {}
    assert m.find_release_pools(*bcat(), claimed={U["c2"]}, stats=st) == ([], [], [], [])
    assert st["release_claimed"] == 1


def test_the_flag_is_on_by_default():
    # default since deltasignal specs/086 was adopted (Adam, 2026-10-21)
    from src.env_flags import BOOL_FLAGS
    assert BOOL_FLAGS["LNG_RELEASE_POOLS"] is True


def test_every_variant_chain_is_in_the_pool():
    # amendment 1: a second variant of the complex (dc2 -> c1b -> c2b -> c3b) also releases p
    net, rmap, umap, profiles = bcat()
    extra = ["dc2", "w1", "c1b", "w2", "c2b", "w3", "c3b", "w4"]
    for i, k in enumerate(extra, 100):
        U[k] = f"00000000-0000-0000-0000-{i:012d}"
    rows = [e("s", "w1", "input"), e("dc2", "w1", "input"), e("w1", "c1b", "output"),
            e("c1b", "w2", "input"), e("kin", "w2", "catalyst"), e("w2", "c2b", "output"),
            e("c2b", "w3", "input"), e("pp", "w3", "catalyst"), e("w3", "c3b", "output"),
            e("c3b", "w4", "input"), e("w4", "p", "output")]
    net = pd.concat([net, pd.DataFrame(rows)], ignore_index=True)
    rmap = pd.concat([rmap, pd.DataFrame({"uid": [U[k] for k in ("w1", "w2", "w3", "w4")],
                                          "reactome_id": ["R-HSA-V1", "R-HSA-V2", "R-HSA-V3", "R-HSA-V4"]})],
                     ignore_index=True)
    umap = dict(umap); umap.update({U["dc2"]: "R-HSA-DC", U["c1b"]: "R-HSA-C1", U["c2b"]: "R-HSA-C2", U["c3b"]: "R-HSA-C3"})
    st = {}
    forms, trans, carriers, edges = m.find_release_pools(net, rmap, umap, profiles, stats=st)
    assert sorted(u for _, u, _, role, base in forms if base) == sorted([U["c3"], U["c3b"]])
    assert sorted(c[2] for c in carriers) == sorted([U["v4"], U["w4"]])
    assert {t[1] for t in trans} == {"release1_f1", "release1_r1", "release1_f2", "release1_r2"}
    assert st["release_pools"] == 1 and st["release_paths"] == 4
