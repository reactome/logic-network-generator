import uuid
from typing import Any, Dict

import pandas as pd
from pandas import DataFrame

from src.argument_parser import logger


def create_reaction_id_map(reactome_ids, decomposed_uid_mapping):
    columns: Dict[str, pd.Series] = {
        "uid": pd.Series(dtype="str"),
        "reactome_id": pd.Series(dtype="int"),
        "input_hash": pd.Series(dtype="str"),
        "output_hash": pd.Series(dtype="str"),
    }
    reaction_id_map: DataFrame = pd.DataFrame(columns)
    print("reactome_ids")
    print(reactome_ids)

    print(decomposed_uid_mapping)
    rows = []
    for reactome_id in reactome_ids:
        print("reactome_id")
        print(reactome_id)
        print(type(reactome_id))
        associated_hashes = decomposed_uid_mapping[decomposed_uid_mapping['reactome_id'] == str(reactome_id)]['uid'].unique().tolist()
        print("associated_hashes")
        print(associated_hashes)
        row = {
                "uid": str(uuid.uuid4()),
                "reactome_id": reactome_id,
                "input_hash": associated_hashes[0],
                "output_hash": associated_hashes[1],
                }
        rows.append(row)

    return pd.DataFrame(rows)


def create_pathway_pi(
    decomposed_uid_mapping: DataFrame,
    reaction_connections: DataFrame,
    best_matches: Any,
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
            decomposed_uid_mapping["reactome_id"] == str(reaction_id)
        ]
        print(rows)

    return pathway_pi
