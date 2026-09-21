"""Variant-node sharing: the same entity in the same role across the variants of
ONE Reactome reaction shares one UUID; across different reactions it does not.

Unconditional since deltasignal specs/020 (nodes -34.6%, edges -21.5%, no
pathway gains cyclic nodes, held-out net zero on both ground-truth axes). The
LNG_SHARE_VARIANT_NODES flag that gated it is gone, and setting it is an error.

Motivation (deltasignal specs/016): HDR's BCDX2 complex existed as 33 node copies
-- one per variant reaction -- and no solver rule over copies is right: min over
them measured -83 held-out, max -117. One node is the fix.
"""
import pytest

from src.logic_network_generator import _register_phase1


def phase1(vr_entities, vr_to_reaction, roots=frozenset(), terms=frozenset()):
    reg = {}
    _register_phase1(vr_entities, reg, set(roots), {}, set(terms), {},
                     vr_to_reaction=vr_to_reaction)
    return reg


# Reaction R1 has two variants (v1, v2) that both produce BCDX2 ("B") and consume "A";
# reaction R2 has one variant (v3) that also produces "B".
VR = {"v1": (["A", "S1"], ["B", "P1"]), "v2": (["A", "S2"], ["B", "P2"]), "v3": (["C"], ["B"])}
V2R = {"v1": "R1", "v2": "R1", "v3": "R2"}


def test_removed_flag_is_an_error_not_a_noop(monkeypatch):
    # Sharing is unconditional. A stale LNG_SHARE_VARIANT_NODES must fail loudly
    # rather than let a run silently measure the default under the old name.
    from src.logic_network_generator import create_pathway_logic_network
    import inspect
    src = inspect.getsource(create_pathway_logic_network)
    assert "LNG_SHARE_VARIANT_NODES was removed" in src
    assert 'os.environ.get("LNG_SHARE_VARIANT_NODES"' not in src


def test_sharing_is_on_with_no_flag_set(monkeypatch):
    monkeypatch.delenv("LNG_SHARE_VARIANT_NODES", raising=False)
    reg = phase1(VR, V2R)
    assert reg[("B", "v1", "output")] == reg[("B", "v2", "output")]
    assert reg[("A", "v1", "input")] == reg[("A", "v2", "input")]


def test_on_same_entity_same_reaction_same_role_shares_one_uuid():
    reg = phase1(VR, V2R)
    assert reg[("B", "v1", "output")] == reg[("B", "v2", "output")]
    assert reg[("A", "v1", "input")] == reg[("A", "v2", "input")]


def test_on_different_reactions_stay_distinct():
    reg = phase1(VR, V2R)
    assert reg[("B", "v1", "output")] != reg[("B", "v3", "output")]   # positional across reactions, as designed


def test_on_roles_are_not_conflated():
    vr = {"v1": (["B"], ["B"])}                       # B both consumed and produced by one variant
    reg = phase1(vr, {"v1": "R1"})
    assert reg[("B", "v1", "input")] != reg[("B", "v1", "output")]


def test_on_variant_specific_members_are_not_conflated():
    reg = phase1(VR, V2R)
    assert reg[("S1", "v1", "input")] != reg[("S2", "v2", "input")]  # different entities
    assert reg[("P1", "v1", "output")] != reg[("P2", "v2", "output")]


def test_boundary_entities_are_left_to_their_own_caches():
    # "A" is a root input: the boundary cache already shares it per stId; sharing
    # must not create a second uuid for it.
    reg = phase1(VR, V2R, roots={"A"})
    assert reg[("A", "v1", "input")] == reg[("A", "v2", "input")]
    assert len({reg[("A", v, "input")] for v in ("v1", "v2")}) == 1


def test_terminal_output_shared_across_reactions_by_its_own_cache_not_the_variant_cache():
    # "B" is an output in R1 (v1, v2) AND in R2 (v3), and a terminal output. The
    # boundary cache shares it per stId across ALL three; an implementation that
    # ignored the boundary exclusion and shared via the variant cache would give
    # v3 a different uuid from v1/v2. All three must be one uuid.
    reg = phase1(VR, V2R, terms={"B"})
    assert len({reg[("B", v, "output")] for v in ("v1", "v2", "v3")}) == 1


def test_missing_reaction_id_read_back_as_nan_string_is_not_a_shared_key():
    # reactome_id read from a cached CSV with no dtype turns NaN into "nan";
    # every variant of every reaction with a missing id must NOT collapse onto
    # one uuid through the key (eid, "nan", role).
    reg = phase1(VR, {"v1": "nan", "v2": "nan", "v3": "nan"})
    assert reg[("B", "v1", "output")] != reg[("B", "v2", "output")]
    assert reg[("B", "v1", "output")] != reg[("B", "v3", "output")]


def test_partial_reaction_map_is_reported():
    from src.logic_network_generator import _register_phase1
    reg = {}
    stats = _register_phase1(VR, reg, set(), {}, set(), {},
                             vr_to_reaction={"v1": "R1", "v3": "R2"})
    assert stats["unmapped"] == 1
    # the unmapped variant is registered per-variant, the mapped ones still share among themselves
    assert reg[("B", "v1", "output")] != reg[("B", "v2", "output")]


def test_unknown_reaction_falls_back_to_per_variant():
    reg = phase1(VR, {})                  # no vr -> reaction map
    assert reg[("B", "v1", "output")] != reg[("B", "v2", "output")]
