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
        "source_id": pd.Series(dtype="str"),
        "target_id": pd.Series(dtype="str"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
    }
    pathway_pi: DataFrame = pd.DataFrame(columns)

    print(reaction_connections)

    for idx, reaction_connection in reaction_connections.iterrows():
        preceding_reaction_id = reaction_connection["preceding_reaction_id"]
        following_reaction_id = reaction_connection["following_reaction_id"]

        print("preceding_reaction_id")
        print(preceding_reaction_id)
        print("following_reaction_id")
        print(following_reaction_id)
        print("fdsfsdfsf")
        print(decomposed_uid_mapping)
        # Assuming preceding_reaction_id is a variable containing the value to filter by
        preceding_rows = decomposed_uid_mapping[
            decomposed_uid_mapping["reactome_id"] == str(preceding_reaction_id)
        ]
        following_rows = decomposed_uid_mapping[
            decomposed_uid_mapping["reactome_id"] == str(following_reaction_id)
        ]
        print(preceding_rows)
        print(following_rows)
        print(best_matches)

    return pathway_pi
