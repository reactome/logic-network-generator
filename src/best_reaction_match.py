"""Hungarian assignment of decomposed inputs to decomposed outputs within one reaction.

Each input combination of a reaction is paired with one output combination based
on how many physical components they share. When the number of input combinations
differs from the number of output combinations, the cost matrix is padded with
zeros so Hungarian can run on a square problem; pairs that hit padding rows or
columns are filtered out before returning, so extras are dropped rather than
falsely matched to non-existent partners.

See tests/test_best_reaction_match.py for the contract this module promises.
"""

from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

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
        str(uid): set(group["component_id_or_reference_entity_id"])
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
    input_reactions: Iterable[str],
    output_reactions: Iterable[str],
    decomposed_uid_mapping: pd.DataFrame,
    reaction_id: Optional[str] = None,
) -> Tuple[List[Tuple[str, str]], List[int]]:
    """Run Hungarian on the components-overlap matrix and return matched pairs.

    When the number of input combinations differs from the number of output
    combinations (typical when EntitySet expansion produces unequal cartesian
    products on the two sides, or for cleavage reactions where one input
    yields several fragment outputs), Hungarian alone would drop the surplus.
    We instead:

      1. Run Hungarian to find an optimal 1-to-1 matching across the square
         padded matrix.
      2. For each surplus input (assigned to a padding column by Hungarian),
         pair it with the output that maximises overlap.
      3. For each surplus output (assigned to a padding row), pair it with
         the input that maximises overlap.

    Both directions are symmetric: every input alternative and every output
    alternative shows up in the resulting list of pairs. Cleavage products
    (zero refEntity overlap with the input) and EntitySet alternatives both
    end up represented in the network. See docs/DESIGN_DECISIONS.md.
    """
    inputs = list(input_reactions)
    outputs = list(output_reactions)

    if not inputs or not outputs:
        return [], []

    counts = create_raw_counts_matrix(inputs, outputs, decomposed_uid_mapping)
    num_rows, num_cols = counts.shape

    if num_rows != num_cols:
        unmatched_count = abs(num_rows - num_cols)
        side = "inputs" if num_rows > num_cols else "outputs"
        logger.debug(
            f"Reaction {reaction_id}: matrix mismatch - "
            f"{num_rows} input combinations vs {num_cols} output combinations; "
            f"{unmatched_count} surplus {side} will be paired with their best counterpart"
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

    # Pair every surplus input with its best output (and vice versa) instead
    # of dropping. EntitySet alternatives and cleavage products are both
    # legitimate biological paths the curator-guide model expects to see.
    if num_rows > num_cols:
        matched_input_indices = {i for i, _ in matched_pairs}
        for surplus_i in range(num_rows):
            if surplus_i not in matched_input_indices and num_cols > 0:
                best_j = int(np.argmax(counts[surplus_i]))
                matched_pairs.append((surplus_i, best_j))
    elif num_cols > num_rows:
        matched_output_indices = {j for _, j in matched_pairs}
        for surplus_j in range(num_cols):
            if surplus_j not in matched_output_indices and num_rows > 0:
                best_i = int(np.argmax(counts[:, surplus_j]))
                matched_pairs.append((best_i, surplus_j))

    matches = [(inputs[i], outputs[j]) for i, j in matched_pairs]
    counts_for_matches = [int(counts[i, j]) for i, j in matched_pairs]
    return matches, counts_for_matches


def find_best_reaction_match(
    input_reactions: Iterable[str],
    output_reactions: Iterable[str],
    decomposed_uid_mapping: pd.DataFrame,
    reaction_id: Optional[str] = None,
) -> Tuple[List[Tuple[str, str]], List[int]]:
    """Public entry point — same signature, kept for callers.

    Returns (matched_pairs, match_counts). Both are empty if either input
    is empty.
    """
    return find_best_match_both_decomposed_reactions(
        input_reactions, output_reactions, decomposed_uid_mapping,
        reaction_id=reaction_id,
    )
