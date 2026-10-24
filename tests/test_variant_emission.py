"""specs/046 variant emitter: copies and output binding (no database)."""
import pytest

import src.variant_emission as ve
import src.variant_keys as vk

REAL_SIGNATURE = ve._reference_signature      # the autouse stub replaces it per test

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


def test_output_with_a_bound_and_a_shared_slot_gives_no_duplicate_copies():
    # deltasignal week review 2026-10-15 (RAF R-HSA-5674373): the output complex
    # has one slot bound to an input slot (OSET <- SET) and one it shares with
    # the input (TIE). The bound slot was enumerated, then overwritten:
    # |OSET| identical copies of each real one.
    LABELS.update({"ICX": ["Complex"], "OCX": ["Complex"]})
    COMP.update({"ICX": ["SET", "TIE"], "OCX": ["OSET", "TIE"]})
    try:
        choices, over = ve._choices(parts(**{"in": {"ICX": 1}, "out": {"OCX": 1}}), 512)
        keys = [frozenset(c.items()) for c in choices]
        assert not over and len(keys) == len(set(keys)) == 4
        assert {(c["SET"], c["OSET"]) for c in choices} == {("RHOA", "pRHOA"), ("RAC1", "pRAC1")}
    finally:
        for k in ("ICX", "OCX"):
            LABELS.pop(k, None); COMP.pop(k, None)


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


# R-HSA-5675376 with its real ids, whose sorted order is what crossed the pairing.
PM, PD, M, D = "R-HSA-109853", "R-HSA-109855", "R-HSA-59282", "R-HSA-5675354"


def _dusp_lookups(monkeypatch):
    """p-T,Y MAPK monomers and dimers -> MAPK monomers and dimers. PM = p-MAPK1,
    PD = p-MAPK1 dimer (PM x2); M = MAPK1, D = MAPK1 dimer (M x2)."""
    import src.neo4j_connector as nc
    import src.reaction_generator as rg
    labels = {"IN": ["EntitySet"], "OUT": ["EntitySet"], PD: ["Complex"], D: ["Complex"]}
    vk._lookups.update(labels=lambda s: labels.get(s, ["EWAS"]), components=lambda s: [],
                       members=lambda s: {"IN": [PD, PM], "OUT": [D, M]}.get(s, []))
    vk.reset_caches()
    monkeypatch.setattr(ve, "_reference_signature", REAL_SIGNATURE)
    monkeypatch.setattr(nc, "get_reference_entity_id", lambda x: {PM: "MAPK1", M: "MAPK1"}.get(x))
    monkeypatch.setattr(nc, "get_complex_components", lambda x: {PD: {PM: 2}, D: {M: 2}}.get(x, {}))
    monkeypatch.setattr(rg, "get_terminal_components", lambda x: {PD: {PM}, D: {M}}.get(x, {x}))


def test_homodimer_binds_to_homodimer_not_monomer(monkeypatch):
    """specs/048: the signature keeps a complex's leaf multiset. Without it the
    dimer reads as its monomer and the tie is broken by sorted rank, which
    paired p-MAPK1 with the MAPK1 DIMER."""
    _dusp_lookups(monkeypatch)
    monkeypatch.setenv("LNG_BIND_STOICH", "1")
    assert ve._reference_signature(PD) == ("MAPK1", "MAPK1")
    assert ve._reference_signature(PM) == ("MAPK1",)
    (slot, mapping), = ve.bind_output_slots({"IN"}, {"OUT"}).values()
    assert slot == "IN" and mapping == {PM: M, PD: D}
    assert ve.STATS.get("slot_binding_ties", 0) == 0


def test_homodimer_binding_off_reproduces_the_crossed_pairing(monkeypatch):
    _dusp_lookups(monkeypatch)
    monkeypatch.setenv("LNG_BIND_STOICH", "0")
    (slot, mapping), = ve.bind_output_slots({"IN"}, {"OUT"}).values()
    assert mapping == {PM: D, PD: M}
    assert ve.STATS["slot_binding_ties"] == 1


def test_leaf_multiset_multiplies_stoichiometry_through_nesting(monkeypatch):
    """A 2x(homodimer) complex reads as four leaves, not two: stoichiometry
    multiplies at every level, and a set inside a complex contributes its
    alternatives once each, times its slot."""
    import src.neo4j_connector as nc
    import src.reaction_generator as rg
    labels = {"TETRA": ["Complex"], "DIM": ["Complex"], "MIX": ["Complex"], "S": ["EntitySet"]}
    vk._lookups.update(labels=lambda s: labels.get(s, ["EWAS"]))
    vk.reset_caches()
    comps = {"TETRA": {"DIM": 2}, "DIM": {"M": 2}, "MIX": {"S": 3, "M": 1}}
    monkeypatch.setattr(nc, "get_complex_components", lambda x: comps.get(x, {}))
    monkeypatch.setattr(rg, "get_terminal_components", lambda x: {"S": {"A", "B"}}.get(x, {x}))
    assert ve._leaf_multiset("TETRA") == ["M"] * 4
    assert ve._leaf_multiset("DIM") == ["M", "M"]
    assert ve._leaf_multiset("MIX") == sorted(["A", "B"] * 3 + ["M"])


def test_signature_cache_is_keyed_on_the_flag(monkeypatch):
    """A signature cached under one LNG_BIND_STOICH value is not reused under
    the other (the review of PR #110: the cache was keyed on the member only)."""
    _dusp_lookups(monkeypatch)
    monkeypatch.setenv("LNG_BIND_STOICH", "0")
    assert ve._reference_signature(PD) == ("MAPK1",)
    monkeypatch.setenv("LNG_BIND_STOICH", "1")
    assert ve._reference_signature(PD) == ("MAPK1", "MAPK1")



def test_priming_a_pathway_clears_set_member_signatures(monkeypatch):
    # deltasignal specs/065: a signature cached from an earlier pathway's
    # prefetch must not survive into the next pathway of a catalog build.
    import src.reaction_generator as rg
    monkeypatch.setattr(rg, "reaction_ids_from_connections", lambda rc: ["R1"])
    monkeypatch.setattr(rg, "prefetch_entity_data", lambda ids: None)
    ve._sig_cache[("M1", True)] = ("stale",)
    rg.prime_entity_caches(None)
    assert ve._sig_cache == {}


def test_priming_a_pathway_clears_derived_structure_caches(monkeypatch):
    # specs/065 amendment 1: every memo computed from the prefetch-limited
    # getters is per pathway.
    import src.reaction_generator as rg
    import src.logic_network_generator as lng
    monkeypatch.setattr(rg, "reaction_ids_from_connections", lambda rc: ["R1"])
    monkeypatch.setattr(rg, "prefetch_entity_data", lambda ids: None)
    caches = (lng._variant_leafsets_cache, lng._matching_leaves_cache,
              lng._handoff_leaf_cache, lng._sets_chosen_cache,
              vk._memo, vk._variants_cache, vk._struct)
    for c in caches:
        c["stale"] = None
    rg.prime_entity_caches(None)
    assert all("stale" not in c for c in caches)


def _mdm_world(monkeypatch):
    """R-HSA-6804741, MDM2 forms homo- or heterodimers: {MDM2, MDM4} + MDM2 x2 ->
    {MDM2 dimer, MDM2:MDM4}. The output set is not bound to the input set (other
    signatures), so the copies are the cross product."""
    import src.neo4j_connector as nc
    labels = {"INSET": ["EntitySet"], "OUTSET": ["EntitySet"], "HOMO": ["Complex"], "HET": ["Complex"]}
    vk._lookups.update(labels=lambda s: labels.get(s, ["EWAS"]),
                       components=lambda s: {"HOMO": ["MDM2"], "HET": ["MDM2", "MDM4"]}.get(s, []),
                       members=lambda s: {"INSET": ["MDM2", "MDM4"], "OUTSET": ["HOMO", "HET"]}.get(s, []))
    vk.reset_caches()
    monkeypatch.setattr(ve, "_reference_signature", lambda m: {"HOMO": ("MDM2", "MDM2"), "HET": ("MDM2", "MDM4")}.get(m, (m,)))
    monkeypatch.setattr(nc, "get_reference_entity_id", lambda x: x if x in ("MDM2", "MDM4") else None)
    return parts(**{"in": {"INSET": 1, "MDM2": 2}, "out": {"OUTSET": 1}})


def test_unbalanced_copies_are_dropped(monkeypatch):
    """specs/101: MDM2 alone must not make MDM2:MDM4, and MDM4 must not vanish into the MDM2 dimer."""
    p = _mdm_world(monkeypatch)
    found, over = ve._choices(p, 100)
    assert not over and len(found) == 4                  # the cross product
    kept = ve.balance_choices(p, found)
    pairs = sorted((c["INSET"], c["OUTSET"]) for c in kept)
    assert pairs == [("MDM2", "HOMO"), ("MDM4", "HET")]
    assert ve.STATS["balance_dropped"] == 2 and ve.STATS["balance_reactions"] == 1


def test_balance_keeps_an_input_with_no_balanced_partner(monkeypatch):
    """An input member no output member balances (incomplete curation) keeps a route: its
    best-balanced copy (MDMX -> MDM2 dimer loses MDMX, d 1; MDMX -> MDM2:MDM4 also makes MDM4, d 2)."""
    p = _mdm_world(monkeypatch)
    vk._lookups["members"] = lambda s: {"INSET": ["MDM2", "MDM4", "MDMX"], "OUTSET": ["HOMO", "HET"]}.get(s, [])
    vk.reset_caches()
    import src.neo4j_connector as nc
    monkeypatch.setattr(nc, "get_reference_entity_id", lambda x: x if x in ("MDM2", "MDM4", "MDMX") else None)
    kept = ve.balance_choices(p, ve._choices(p, 100)[0])
    pairs = sorted((c["INSET"], c["OUTSET"]) for c in kept)
    assert pairs == [("MDM2", "HOMO"), ("MDM4", "HET"), ("MDMX", "HOMO")]


def test_balance_flag_off_by_default():
    from src.env_flags import env_flag
    assert env_flag("LNG_BALANCE_CHOICES") is False
