"""Tests for best_reaction_match: Hungarian pairing of decomposed inputs and outputs.

This module is the heart of input-to-output pairing within a single reaction.
Subtle bugs here (especially in the non-square padding and drop logic) would
silently mislabel virtual reactions, so these tests are intentionally low-level.

Decomposed-uid-mapping shape used throughout: each `uid` maps to one or more
component_id_or_reference_entity_id values; the matcher counts overlap between
the components of an input combination and the components of an output
combination.
"""

import pandas as pd

from src.best_reaction_match import (
    create_raw_counts_matrix,
    find_best_reaction_match,
)


def _mapping(rows):
    """Build a decomposed_uid_mapping DataFrame from a list of (uid, comp) pairs."""
    return pd.DataFrame(
        rows, columns=["uid", "component_id_or_reference_entity_id"]
    )


class TestCreateRawCountsMatrix:
    """Counts matrix M[i,j] = |components(input_i) ∩ components(output_j)|."""

    def test_no_overlap_yields_zeros(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("out1", "X"), ("out1", "Y"),
        ])
        m = create_raw_counts_matrix(["in1"], ["out1"], mapping)
        assert m.shape == (1, 1)
        assert m[0, 0] == 0

    def test_full_overlap_yields_count(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"), ("in1", "C"),
            ("out1", "A"), ("out1", "B"), ("out1", "C"),
        ])
        m = create_raw_counts_matrix(["in1"], ["out1"], mapping)
        assert m[0, 0] == 3

    def test_partial_overlap(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("out1", "A"), ("out1", "X"),
        ])
        m = create_raw_counts_matrix(["in1"], ["out1"], mapping)
        assert m[0, 0] == 1

    def test_multiple_inputs_and_outputs(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("in2", "C"), ("in2", "D"),
            ("out1", "A"), ("out1", "B"),
            ("out2", "C"), ("out2", "D"),
        ])
        m = create_raw_counts_matrix(["in1", "in2"], ["out1", "out2"], mapping)
        assert m[0, 0] == 2
        assert m[0, 1] == 0
        assert m[1, 0] == 0
        assert m[1, 1] == 2


class TestSquareHungarianMatching:
    """When input/output counts are equal, every combination is paired."""

    def test_optimal_pairing_picked(self):
        # Optimal: in1↔out1 (2 overlap), in2↔out2 (2 overlap), total 4.
        # Suboptimal (in1↔out2, in2↔out1) would give 0 overlap.
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("in2", "C"), ("in2", "D"),
            ("out1", "A"), ("out1", "B"),
            ("out2", "C"), ("out2", "D"),
        ])
        matches, counts = find_best_reaction_match(
            ["in1", "in2"], ["out1", "out2"], mapping
        )
        assert set(matches) == {("in1", "out1"), ("in2", "out2")}
        assert sorted(counts) == [2, 2]

    def test_zero_overlap_still_returns_pair(self):
        mapping = _mapping([
            ("in1", "A"),
            ("out1", "X"),
        ])
        matches, counts = find_best_reaction_match(["in1"], ["out1"], mapping)
        assert matches == [("in1", "out1")]
        assert counts == [0]


class TestNonSquarePaddingDrops:
    """Extra inputs (or outputs) get assigned to padding columns and dropped.

    This is the subtle correctness path: padding the cost matrix with zeros
    lets Hungarian run on a square problem, but the resulting pairs that hit
    the fake (padded) rows/columns must NOT be returned to the caller.
    """

    def test_more_inputs_than_outputs_drops_extras(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("in2", "X"),
            ("out1", "A"), ("out1", "B"),
        ])
        matches, _ = find_best_reaction_match(
            ["in1", "in2"], ["out1"], mapping
        )
        assert len(matches) == 1, "extra input must be dropped, not paired with a fake output"
        assert matches[0] == ("in1", "out1")

    def test_more_outputs_than_inputs_drops_extras(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("out1", "A"), ("out1", "B"),
            ("out2", "X"),
        ])
        matches, _ = find_best_reaction_match(
            ["in1"], ["out1", "out2"], mapping
        )
        assert len(matches) == 1, "extra output must be dropped, not paired with a fake input"
        assert matches[0] == ("in1", "out1")

    def test_non_square_picks_best_match_then_drops(self):
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"), ("in1", "C"),
            ("in2", "X"), ("in2", "Y"),
            ("in3", "Z"),
            ("out1", "A"), ("out1", "B"), ("out1", "C"),
            ("out2", "X"), ("out2", "Y"),
        ])
        matches, counts = find_best_reaction_match(
            ["in1", "in2", "in3"], ["out1", "out2"], mapping
        )
        assert len(matches) == 2
        assert set(matches) == {("in1", "out1"), ("in2", "out2")}
        assert sorted(counts) == [2, 3]


class TestEdgeCases:
    def test_empty_inputs_returns_no_matches(self):
        mapping = _mapping([("out1", "A")])
        matches, counts = find_best_reaction_match([], ["out1"], mapping)
        assert matches == []
        assert counts == []

    def test_empty_outputs_returns_no_matches(self):
        mapping = _mapping([("in1", "A")])
        matches, counts = find_best_reaction_match(["in1"], [], mapping)
        assert matches == []
        assert counts == []
