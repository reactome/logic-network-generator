"""Hungarian assignment of decomposed inputs to decomposed outputs within one reaction.

Each input combination of a reaction is paired with one output combination based
on how many physical components they share. When the number of input combinations
differs from the number of output combinations, the cost matrix is padded with
zeros so Hungarian can run on a square problem; pairs that hit padding rows or
columns are filtered out before returning, so extras are dropped rather than
falsely matched to non-existent partners.

See tests/test_best_reaction_match.py for the contract this module promises.
"""

from typing import Dict, List, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment  # type: ignore

from src.argument_parser import logger


def _build_uid_to_components(
    uids: Sequence[str], decomposed_uid_mapping: pd.DataFrame
) -> Dict[str, Set[str]]:
    """One DataFrame scan instead of two-per-pair: pre-index uid → component set."""
    if not uids:
        return {}
    rows = decomposed_uid_mapping[decomposed_uid_mapping["uid"].isin(uids)]
    return {
        uid: set(group["component_id_or_reference_entity_id"])
        for uid, group in rows.groupby("uid")
    }


def create_raw_counts_matrix(
    input_reactions: Sequence[str],
    output_reactions: Sequence[str],
    decomposed_uid_mapping: pd.DataFrame,
) -> np.ndarray:
    """M[i,j] = |components(input_i) ∩ components(output_j)|."""
    inputs = list(input_reactions)
    outputs = list(output_reactions)

    in_components = _build_uid_to_components(inputs, decomposed_uid_mapping)
    out_components = _build_uid_to_components(outputs, decomposed_uid_mapping)

    counts = np.zeros((len(inputs), len(outputs)))
    for i, in_uid in enumerate(inputs):
        in_set = in_components.get(in_uid, set())
        for j, out_uid in enumerate(outputs):
            counts[i, j] = len(in_set & out_components.get(out_uid, set()))
    return counts


def find_best_match_both_decomposed_reactions(
    input_reactions: Sequence[str],
    output_reactions: Sequence[str],
    decomposed_uid_mapping: pd.DataFrame,
    reaction_id: str = None,
) -> Tuple[List[Tuple[str, str]], List[int]]:
    """Run Hungarian on the components-overlap matrix and return matched pairs."""
    inputs = list(input_reactions)
    outputs = list(output_reactions)

    if not inputs or not outputs:
        return [], []

    counts = create_raw_counts_matrix(inputs, outputs, decomposed_uid_mapping)
    num_rows, num_cols = counts.shape

    if num_rows != num_cols:
        unmatched_count = abs(num_rows - num_cols)
        side = "inputs" if num_rows > num_cols else "outputs"
        logger.warning(
            f"Reaction {reaction_id}: Hungarian matching dimension mismatch - "
            f"{num_rows} input combinations vs {num_cols} output combinations; "
            f"{unmatched_count} {side} will be unmatched"
        )
        max_dim = max(num_rows, num_cols)
        padded_counts = np.zeros((max_dim, max_dim))
        padded_counts[:num_rows, :num_cols] = counts
    else:
        padded_counts = counts

    # Negate to turn a maximum-overlap problem into a minimum-cost problem.
    # Padding rows/cols have cost 0; real pairs with overlap have cost
    # -count (negative), so Hungarian prefers real assignments and only
    # uses padding when there's nothing real left.
    row_indices, col_indices = linear_sum_assignment(-padded_counts)

    matched_pairs = [
        (i, j)
        for i, j in zip(row_indices, col_indices)
        if i < num_rows and j < num_cols  # filter assignments that hit padding
    ]
    matches = [(inputs[i], outputs[j]) for i, j in matched_pairs]
    counts_for_matches = [int(counts[i, j]) for i, j in matched_pairs]
    return matches, counts_for_matches


def find_best_reaction_match(
    input_reactions: Sequence[str],
    output_reactions: Sequence[str],
    decomposed_uid_mapping: pd.DataFrame,
    reaction_id: str = None,
) -> Tuple[List[Tuple[str, str]], List[int]]:
    """Public entry point — same signature, kept for callers.

    Returns (matched_pairs, match_counts). Both are empty if either input
    is empty.
    """
    return find_best_match_both_decomposed_reactions(
        input_reactions, output_reactions, decomposed_uid_mapping,
        reaction_id=reaction_id,
    )
