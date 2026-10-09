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
