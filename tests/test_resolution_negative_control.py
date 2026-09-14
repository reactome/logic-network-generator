"""The completeness checks must FAIL on a corrupted mapping.

This file exists because of a specific past failure: `_decomposed_ids`
reported 11 of 11 passing while masking 18 dropped catalysts, because the
check had been relaxed until it could no longer fail. A completeness check
that always passes is worse than no check, since it converts an unknown into
a false assurance.

So every check here is exercised twice — once against a correct mapping,
where it must be silent, and once against a mapping broken in exactly the way
that check exists to catch, where it must speak up.
"""
import pytest

from src.resolution_validation import (
    check_exclusions_have_reasons, check_forward_completeness,
    check_glyph_pairing, check_referential_integrity, check_release_recorded,
    check_reverse_completeness, run_all, summarise,
)


def _row(stable_id, uuid, relation="self", **kw):
    row = {"stable_id": stable_id, "uuid": uuid, "relation": relation,
           "depth": 0, "role": "", "reaction_stid": "", "glyph_id": "",
           "diagram_stid": "", "release": "97"}
    row.update(kw)
    return row


@pytest.fixture
def good():
    """A small but complete mapping: two entity nodes and one set over them."""
    network = {"u1", "u2"}
    resolution = [
        _row("R-HSA-1", "u1"),
        _row("R-HSA-2", "u2"),
        _row("R-HSA-SET", "u1", "set_member", depth=1),
        _row("R-HSA-SET", "u2", "set_member", depth=1),
    ]
    exclusions = [{"stable_id": "R-HSA-UB", "release": "97",
                   "reason": "atomic modifier set, deliberately not expanded"}]
    participating = {"R-HSA-1", "R-HSA-2", "R-HSA-SET", "R-HSA-UB"}
    return network, resolution, exclusions, participating


def test_a_correct_mapping_passes(good):
    network, resolution, exclusions, participating = good
    assert run_all(network, resolution, exclusions, participating, "97") == []


# --- each check, broken in the way it exists to catch -----------------------

def test_reverse_completeness_fails_when_a_node_maps_to_nothing(good):
    network, resolution, _, _ = good
    network = network | {"u3_orphan"}
    violations = check_reverse_completeness(network, resolution)
    assert [v.subject for v in violations] == ["u3_orphan"]


def test_referential_integrity_fails_on_a_row_naming_a_missing_node(good):
    """This is issue #67's exact shape: the row named the undecomposed parent."""
    network, resolution, _, _ = good
    resolution = resolution + [_row("R-HSA-9", "parent_uuid_not_in_network")]
    violations = check_referential_integrity(network, resolution)
    assert len(violations) == 1
    assert violations[0].subject == "parent_uuid_not_in_network"


def test_forward_completeness_fails_on_an_undeclared_absence(good):
    _, resolution, exclusions, participating = good
    participating = participating | {"R-HSA-NEVER-GENERATED"}
    violations = check_forward_completeness(participating, resolution, exclusions)
    assert [v.subject for v in violations] == ["R-HSA-NEVER-GENERATED"]


def test_an_excluded_entity_is_not_a_forward_violation(good):
    """Declaring an absence is the point — it must not also be reported."""
    _, resolution, exclusions, participating = good
    assert check_forward_completeness(participating, resolution, exclusions) == []


def test_exclusion_without_a_reason_fails():
    violations = check_exclusions_have_reasons([{"stable_id": "R-HSA-X", "reason": ""}])
    assert len(violations) == 1


def test_vague_exclusion_reason_fails():
    """"missing" does not distinguish a design decision from a bug."""
    violations = check_exclusions_have_reasons([{"stable_id": "R-HSA-X", "reason": "missing"}])
    assert len(violations) == 1
    assert "vague" in violations[0].detail


def test_glyph_id_without_its_diagram_fails(good):
    _, resolution, _, _ = good
    resolution = resolution + [_row("R-HSA-3", "u1", glyph_id="535")]
    violations = check_glyph_pairing(resolution)
    assert len(violations) == 1


def test_diagram_without_its_glyph_also_fails(good):
    _, resolution, _, _ = good
    resolution = resolution + [_row("R-HSA-3", "u1", diagram_stid="R-HSA-1257604")]
    assert len(check_glyph_pairing(resolution)) == 1


def test_both_glyph_fields_together_pass(good):
    _, resolution, _, _ = good
    resolution = resolution + [_row("R-HSA-3", "u1", glyph_id="535",
                                    diagram_stid="R-HSA-1257604")]
    assert check_glyph_pairing(resolution) == []


def test_missing_release_fails(good):
    _, resolution, _, _ = good
    resolution = resolution + [_row("R-HSA-3", "u1", release="")]
    assert len(check_release_recorded(resolution)) == 1


def test_mixed_releases_in_one_table_fails(good):
    _, resolution, _, _ = good
    resolution = resolution + [_row("R-HSA-3", "u1", release="96")]
    violations = check_release_recorded(resolution)
    assert any("mixed releases" in v.detail for v in violations)


def test_release_mismatch_against_the_networks_fails(good):
    _, resolution, _, _ = good
    assert len(check_release_recorded(resolution, expected_release="96")) == 1


def test_run_all_reports_every_broken_check_at_once(good):
    """One corrupted mapping, several independent failures — none masked."""
    network, resolution, exclusions, participating = good
    network = network | {"u_orphan"}
    resolution = resolution + [
        _row("R-HSA-9", "not_in_network"),
        _row("R-HSA-3", "u1", glyph_id="535"),
    ]
    exclusions = exclusions + [{"stable_id": "R-HSA-Y", "reason": ""}]
    participating = participating | {"R-HSA-NEVER"}
    counts = summarise(run_all(network, resolution, exclusions, participating, "97"))
    assert counts == {
        "reverse_completeness": 1,
        "referential_integrity": 1,
        "exclusion_reason": 1,
        "glyph_pairing": 1,
        "forward_completeness": 1,
    }
