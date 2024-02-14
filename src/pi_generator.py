from typing import Any, Dict

import pandas as pd
from pandas import DataFrame

from src.argument_parser import logger


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
        best_matches (Any): Data representing best matches.

    Returns:
        DataFrame: The created pathway_pi DataFrame.
    """
    logger.debug("Adding reaction pairs to pathway_pi")

    columns: Dict[str, pd.Series] = {
        "parent_id": pd.Series(dtype="str"),
        "child_id": pd.Series(dtype="str"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
    }
    pathway_pi: DataFrame = pd.DataFrame(columns)

    print(reaction_connections)

    for idx, reaction_connection in reaction_connections.iterrows():
        parent_reaction_id = reaction_connection["parent_reaction_id"]
        child_reaction_id = reaction_connection["child_reaction_id"]

        print("parent_reaction_id")
        print(parent_reaction_id)
        print("child_reaction_id")
        print(child_reaction_id)
        print("fdsfsdfsf")
        print(decomposed_uid_mapping)
        # Assuming parent_reaction_id is a variable containing the value to filter by
        parent_rows = decomposed_uid_mapping[
            decomposed_uid_mapping["reactome_id"] == str(parent_reaction_id)
        ]
        child_rows = decomposed_uid_mapping[
            decomposed_uid_mapping["reactome_id"] == str(child_reaction_id)
        ]
        print(parent_rows)
        print(child_rows)
        print(best_matches)

    return pathway_pi
