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


class TestNonSquareSurplusFansOut:
    """Surplus inputs/outputs get paired with their best counterpart, not dropped.

    EntitySet expansion routinely produces mismatched input vs output
    combination counts; cleavage reactions produce one input and several
    fragment outputs with zero refEntity overlap. In every case the
    surplus side represents real biology and should appear as virtual
    reactions, not be discarded. See docs/DESIGN_DECISIONS.md.
    """

    def test_more_inputs_than_outputs_pairs_each_input(self):
        # 2 inputs, 1 output. Hungarian pairs in1↔out1 (both share A,B).
        # in2 has zero overlap with out1 — but it's still a valid alternative
        # input alternative, so it pairs with out1 too.
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("in2", "X"),
            ("out1", "A"), ("out1", "B"),
        ])
        matches, _ = find_best_reaction_match(
            ["in1", "in2"], ["out1"], mapping
        )
        assert len(matches) == 2, "every input alternative must produce a pair"
        assert set(matches) == {("in1", "out1"), ("in2", "out1")}

    def test_more_outputs_than_inputs_pairs_each_output(self):
        # 1 input, 2 outputs. Hungarian pairs in1↔out1.
        # out2 has zero overlap (cleavage-style) — pair it with in1 too.
        mapping = _mapping([
            ("in1", "A"), ("in1", "B"),
            ("out1", "A"), ("out1", "B"),
            ("out2", "X"),
        ])
        matches, _ = find_best_reaction_match(
            ["in1"], ["out1", "out2"], mapping
        )
        assert len(matches) == 2, "every output alternative must produce a pair"
        assert set(matches) == {("in1", "out1"), ("in1", "out2")}

    def test_non_square_picks_best_match_and_pairs_surplus(self):
        # 3 inputs, 2 outputs. Hungarian picks the best 1-to-1 (in1↔out1, in2↔out2).
        # in3 has zero overlap with anything; argmax picks output 0 (out1).
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
        assert len(matches) == 3, "all 3 input alternatives must show up"
        assert ("in1", "out1") in matches
        assert ("in2", "out2") in matches
        # in3 pairs with whichever output is its argmax (zero-overlap → first)
        assert any(m[0] == "in3" for m in matches)

    def test_cleavage_one_input_three_outputs_each_pair(self):
        # Cleavage: input molecule X is cleaved into 3 fragments F1, F2, F3.
        # The fragments don't share refEntity with the input, so all overlaps
        # are zero — but every fragment must still be linked to the input.
        mapping = _mapping([
            ("in1", "X"),
            ("out1", "F1"),
            ("out2", "F2"),
            ("out3", "F3"),
        ])
        matches, _ = find_best_reaction_match(
            ["in1"], ["out1", "out2", "out3"], mapping
        )
        assert len(matches) == 3, "every cleavage product must be linked to the input"
        assert {m[1] for m in matches} == {"out1", "out2", "out3"}
        assert all(m[0] == "in1" for m in matches), \
            "all fragments come from the same input"


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
