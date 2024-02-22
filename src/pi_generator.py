import uuid
from typing import Any, Dict

import pandas as pd
from pandas import DataFrame

from src.argument_parser import logger


def create_reaction_id_map(reactome_ids, decomposed_uid_mapping):
    reaction_id_map_column_types = {
        'uid': str,
        'reactome_id': pd.Int64Dtype(),
        'input_hash': str,
        'output_hash': str,
    }
    reaction_id_map = (
        pd.DataFrame(columns=reaction_id_map_column_types.keys())
        .astype(reaction_id_map_column_types)
    )

    rows = []
    for reactome_id in reactome_ids:
        associated_hashes = decomposed_uid_mapping[decomposed_uid_mapping['reactome_id'] == reactome_id]['uid'].unique().tolist()
        row = {
                "uid": str(uuid.uuid4()),
                "reactome_id": int(reactome_id),
                "input_hash": associated_hashes[0],
                "output_hash": associated_hashes[1],
                }
        rows.append(row)

    return pd.DataFrame(rows)


def create_pathway_pi(
    decomposed_uid_mapping: DataFrame,
    reaction_connections: DataFrame,
) -> DataFrame:
    """
    Create pathway_pi DataFrame based on decomposed_uid_mapping, reaction_connections, and best_matches.

    Args:
        decomposed_uid_mapping (DataFrame): DataFrame containing decomposed UID mapping.
        reaction_connections (DataFrame): DataFrame containing reaction connections.

    Returns:
        DataFrame: The created pathway_pi DataFrame.
    """
    logger.debug("Adding reaction pairs to pathway_pi")

    columns: Dict[str, pd.Series] = {
        "source_id": pd.Series(dtype="str"),
        "target_id": pd.Series(dtype="str"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
    }
    pathway_pi: DataFrame = pd.DataFrame(columns)

    print("reaction_connetions")
    print(reaction_connections)

    reaction_ids = pd.unique(
        reaction_connections[["preceding_reaction_id", "following_reaction_id"]]
        .stack()  # Stack the columns to convert them into a single series
        .dropna()  # Drop NaN values
    )

    reaction_id_map = create_reaction_id_map(reaction_ids, decomposed_uid_mapping)
    print("reaction_id_map")
    print(reaction_id_map)

    for reaction_id in reaction_ids:

        rows = decomposed_uid_mapping[
            decomposed_uid_mapping["reactome_id"] == reaction_id
        ]
        print(rows)

    return pathway_pi
