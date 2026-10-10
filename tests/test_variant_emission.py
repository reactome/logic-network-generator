"""specs/046 variant emitter: copies and output binding (no database)."""
import pytest

import src.variant_emission as ve
import src.variant_keys as vk

LABELS = {"CX": ["Complex"], "SET": ["EntitySet"], "OSET": ["EntitySet"], "GDP": ["SimpleEntity"],
          "DIMER": ["Complex"], "CAT": ["Complex"], "TIE": ["EntitySet"], "OTIE": ["EntitySet"]}
COMP = {"CX": ["GDP", "SET"], "DIMER": ["SET", "SET"], "CAT": ["SET", "GDP"]}
MEM = {"SET": ["RHOA", "RAC1"], "OSET": ["pRAC1", "pRHOA"], "TIE": ["A1", "A2"], "OTIE": ["B1", "B2"]}
SIG = {"RHOA": ("RHOA",), "pRHOA": ("RHOA",), "RAC1": ("RAC1",), "pRAC1": ("RAC1",),
       "A1": ("A",), "A2": ("A",), "B1": ("A",), "B2": ("A",)}


@pytest.fixture(autouse=True)
def stub(monkeypatch):
    vk._lookups.update(labels=lambda s: LABELS.get(s, ["EWAS"]), components=lambda s: COMP.get(s, []),
                       members=lambda s: MEM.get(s, []), atomic_sets=lambda: set(), max_variants=lambda: 512)
    vk.reset_caches(); ve.STATS.clear(); ve._sig_cache.clear()
    monkeypatch.setattr(ve, "_reference_signature", lambda m: SIG.get(m, (m,)))
    yield
    vk._lookups.clear(); vk.reset_caches()


def parts(**kw):
    p = {"in": {}, "out": {}, "cat": [], "pos": [], "neg": []}
    p.update(kw); return p


def test_catalyst_equal_to_input_is_one_choice_per_copy():
    choices, over = ve._choices(parts(**{"in": {"CX": 1}, "cat": ["CX"]}), 512)
    assert not over and len(choices) == 2
    for sig in choices:
        assert vk.vkey("CX", sig) == f"CX::variant::SET={sig['SET']}"


def test_dimer_of_one_set_gives_homo_copies_only():
    choices, _ = ve._choices(parts(**{"in": {"DIMER": 1}}), 512)
    assert sorted(c["SET"] for c in choices) == ["RAC1", "RHOA"]


def test_output_only_set_is_bound_by_identity_not_multiplied():
    choices, _ = ve._choices(parts(**{"in": {"SET": 1}, "out": {"OSET": 1}}), 512)
    assert len(choices) == 2
    assert {(c["SET"], c["OSET"]) for c in choices} == {("RHOA", "pRHOA"), ("RAC1", "pRAC1")}
    assert ve.STATS["output_slots_bound"] == 1


def test_a_signature_tie_is_paired_by_rank_and_counted():
    choices, _ = ve._choices(parts(**{"in": {"TIE": 1}, "out": {"OTIE": 1}}), 512)
    assert {(c["TIE"], c["OTIE"]) for c in choices} == {("A1", "B1"), ("A2", "B2")}
    assert ve.STATS["slot_binding_ties"] == 1


def test_negative_regulator_shares_the_copy_choice():
    choices, _ = ve._choices(parts(**{"in": {"CX": 1}, "neg": ["CAT"]}), 512)
    assert len(choices) == 2          # D4: expanded, but the shared set is one variable


def test_over_the_cap_is_reported():
    choices, over = ve._choices(parts(**{"in": {"CX": 1}, "out": {"OSET": 1}}), 1)
    assert over


def test_a_bound_member_complex_opens_its_own_slots():
    # output set bound to an input set; the mapped member is a complex that
    # holds a set of its own, which must be fanned out (IL-3 R-HSA-879909)
    LABELS.update({"ISET": ["EntitySet"], "OSET2": ["EntitySet"], "MCX": ["Complex"], "INNER": ["EntitySet"]})
    MEM.update({"ISET": ["RHOA"], "OSET2": ["MCX"], "INNER": ["X1", "X2"]})
    COMP["MCX"] = ["INNER"]
    SIG["MCX"] = ("RHOA",)
    choices, over = ve._choices(parts(**{"in": {"ISET": 1}, "out": {"OSET2": 1}}), 512)
    assert not over and sorted(c["INNER"] for c in choices) == ["X1", "X2"]
    for c in choices:
        assert vk.vkey("OSET2", c).startswith("MCX::variant::INNER=")


# --- D6 cap fallback (specs/046 research: the IFN alpha/beta cap seam) ---

CAP_LABELS = {"REG": ["Complex"], "RS1": ["EntitySet"], "RS2": ["EntitySet"], "OUT": ["EntitySet"],
              "BIG": ["EntitySet"], "IC": ["Complex"]}
CAP_COMP = {"REG": ["RS1", "RS2"], "IC": ["GENE", "OUT"]}
CAP_MEM = {"RS1": ["r1", "r2", "r3"], "RS2": ["s1", "s2", "s3"], "OUT": ["o1", "o2", "o3", "o4"],
           "BIG": [f"b{i}" for i in range(12)]}


@pytest.fixture
def cap_world():
    vk._lookups.update(labels=lambda s: CAP_LABELS.get(s, ["EWAS"]),
                       components=lambda s: CAP_COMP.get(s, []),
                       members=lambda s: CAP_MEM.get(s, []))
    vk.reset_caches()


def test_step1_pools_a_regulator_whose_choice_reaches_no_output(cap_world):
    # 9 regulator variants x 4 free output members = 36 copies, over a cap of 10.
    p = parts(**{"in": {"GENE": 1}, "out": {"OUT": 1}, "pos": ["REG"]})
    assert ve._choices(p, 10)[1]
    choices, pooled, step = ve.capped_fallback(p, 10)
    assert step == 1 and pooled == {"REG"}
    assert sorted(c["OUT"] for c in choices) == ["o1", "o2", "o3", "o4"]


def test_step1_keeps_a_participant_holding_an_output_slot(cap_world):
    # IC holds OUT's slot (reaches the output) and is kept; REG is pooled.
    p = parts(**{"in": {"IC": 1}, "out": {"OUT": 1}, "neg": ["REG"]})
    choices, pooled, step = ve.capped_fallback(p, 5)
    assert step == 1 and pooled == {"REG"} and len(choices) == 4
    assert all(vk.vkey("IC", c) == f"IC::variant::OUT={c['OUT']}" for c in choices)


def test_step2_gives_one_copy_per_output_variant(cap_world):
    # MIX holds OUT's slot AND RS1: kept by step 1 (4 x 3 = 12 > 5), pooled by step 2.
    CAP_LABELS["MIX"] = ["Complex"]; CAP_COMP["MIX"] = ["OUT", "RS1"]
    try:
        vk.reset_caches()
        p = parts(**{"in": {"MIX": 1}, "out": {"OUT": 1}})
        choices, pooled, step = ve.capped_fallback(p, 5)
        assert step == 2 and pooled == {"MIX"} and len(choices) == 4
    finally:
        del CAP_LABELS["MIX"], CAP_COMP["MIX"]


def test_step3_single_copy_when_outputs_alone_exceed_the_cap(cap_world):
    p = parts(**{"in": {"GENE": 1}, "out": {"BIG": 1}, "cat": ["REG"]})
    choices, pooled, step = ve.capped_fallback(p, 10)
    assert step == 3 and choices == [None] and pooled == {"REG"}


def test_emitted_copy_reads_the_pool_reference(cap_world, monkeypatch):
    monkeypatch.setattr(vk, "_max_variants", lambda: 10)
    monkeypatch.setattr(ve, "fetch_participants", lambda g, r: {
        "R1": parts(**{"in": {"GENE": 1}, "out": {"OUT": 1}, "pos": ["REG"]})})
    rid_map, vr, cat, neg, pos = ve.build_variant_reactions(None, ["R1"])
    assert len(rid_map) == 4 and set(pos["entity_id"]) == {"REG::pool"}
    assert {o for (_, outs, _, _) in vr.values() for o in outs} == {
        f"o{i}" for i in range(1, 5)}       # a bare set's key is its member's key
    assert ve.STATS["capped_step1"] == 1


def test_step1_does_not_pool_a_participant_sharing_a_slot_with_a_kept_one(cap_world):
    # review of vn6: catalyst A = {S, T} shares S with input B = {S, OUT}.
    # Pooling A would let each copy (fixed S for B) read every S variant of A.
    CAP_LABELS.update(A=["Complex"], B=["Complex"], S=["EntitySet"], T=["EntitySet"])
    CAP_COMP.update(A=["S", "T"], B=["S", "OUT"])
    CAP_MEM.update(S=["s1", "s2", "s3"], T=["t1", "t2", "t3"])
    try:
        vk.reset_caches()
        p = parts(**{"in": {"B": 1}, "out": {"OUT": 1}, "cat": ["A"]})
        choices, pooled, step = ve.capped_fallback(p, 20)
        assert "A" not in pooled or step != 1
    finally:
        for d, ks in ((CAP_LABELS, "ABST"), (CAP_COMP, "AB"), (CAP_MEM, "ST")):
            for k in ks:
                d.pop(k, None)
