"""specs/046 variant keys (no database: lookups injected)."""
import pytest

import src.variant_keys as vk

# GPVI-like: CX = GDP + SET{RHOA, RAC1}; FGFR2-like: OUTER{INNER{C1, C2}} with
# C1 = HS + FGFSET{F7, F10} (a member complex that itself holds a set)
LABELS = {"CX": ["Complex"], "SET": ["EntitySet", "CandidateSet"], "GDP": ["SimpleEntity"],
          "RHOA": ["EWAS"], "RAC1": ["EWAS"], "OUTER": ["EntitySet", "DefinedSet"],
          "INNER": ["EntitySet", "DefinedSet"], "C1": ["Complex"], "C2": ["Complex"],
          "HS": ["SimpleEntity"], "FGFSET": ["EntitySet", "DefinedSet"], "F7": ["EWAS"],
          "F10": ["EWAS"], "R2": ["EWAS"], "UB": ["EntitySet"], "DIMER": ["Complex"],
          "CAT": ["Complex"]}
COMP = {"CX": ["GDP", "SET"], "C1": ["HS", "FGFSET"], "C2": ["R2"], "DIMER": ["SET", "SET"],
        "CAT": ["SET", "GDP"]}
MEM = {"SET": ["RHOA", "RAC1"], "OUTER": ["INNER"], "INNER": ["C1", "C2"], "FGFSET": ["F7", "F10"],
       "UB": ["x", "y"]}


@pytest.fixture(autouse=True)
def stub():
    vk._lookups.update(labels=lambda s: LABELS.get(s, []), components=lambda s: COMP.get(s, []),
                       members=lambda s: MEM.get(s, []), atomic_sets=lambda: {"UB"},
                       max_variants=lambda: 512)
    vk.reset_caches()
    yield
    vk._lookups.clear(); vk.reset_caches()


def test_complex_key_names_slot_and_member():
    assert vk.vkey("CX", {"SET": "RAC1"}) == "CX::variant::SET=RAC1"
    assert vk.parse_variant_key("CX::variant::SET=RAC1") == ("CX", {"SET": "RAC1"})
    assert vk.variant_leaves("CX::variant::SET=RAC1") == {"GDP", "RAC1"}


def test_bare_set_dissolves_to_member_key_and_nested_sets_flatten():
    assert vk.flat_members("OUTER") == ["C1", "C2"]
    # outer set, member complex with its own set: hoisted slot
    s = {"OUTER": "C1", "FGFSET": "F7"}
    assert vk.vkey("OUTER", s) == vk.vkey("C1", s) == "C1::variant::FGFSET=F7"


def test_producer_of_inner_set_and_consumer_of_outer_set_name_the_same_node():
    # F10: R-HSA-190408 outputs INNER, R-HSA-5654404 inputs OUTER
    assert vk.vkey("INNER", {"INNER": "C1", "FGFSET": "F10"}) == \
        vk.vkey("OUTER", {"OUTER": "C1", "FGFSET": "F10"})


def test_modifier_sets_and_set_free_entities_are_plain():
    assert vk.vkey("UB", {}) == "UB"
    assert vk.vkey("GDP", {}) == "GDP"
    assert vk.vkey("C2", {}) == "C2"


def test_capped_complex_is_its_plain_stid(monkeypatch):
    vk._lookups["max_variants"] = lambda: 1
    vk.reset_caches()
    assert vk.is_capped("CX")
    assert vk.vkey("CX", {}) == "CX"


def test_reaction_choices_share_one_variable_per_set():
    # D2: SET appears as input (CX) and as catalyst (CAT): one choice each copy
    choices = list(vk.reaction_choices(["CX", "CAT"]))
    assert choices == [{"SET": "RAC1"}, {"SET": "RHOA"}]
    # a stoichiometry-2 set gives homo variants only (D3)
    assert list(vk.reaction_choices(["DIMER"])) == [{"SET": "RAC1"}, {"SET": "RHOA"}]


def test_member_complex_opens_its_slots_recursively():
    got = list(vk.reaction_choices(["OUTER"]))
    assert {"OUTER": "C1", "FGFSET": "F10"} in got and {"OUTER": "C1", "FGFSET": "F7"} in got
    assert {"OUTER": "C2"} in got and len(got) == 3


def test_choices_respect_the_limit():
    assert len(list(vk.reaction_choices(["OUTER"], limit=2))) == 2


def test_keys_do_not_depend_on_lookup_order():
    a = vk.vkey("CX", {"SET": "RHOA"})
    MEM["SET"] = ["RAC1", "RHOA"][::-1]
    vk.reset_caches()
    assert vk.vkey("CX", {"SET": "RHOA"}) == a


def test_malformed_key_is_an_error():
    with pytest.raises(ValueError):
        vk.parse_variant_key("CX::variant::SETRAC1")
