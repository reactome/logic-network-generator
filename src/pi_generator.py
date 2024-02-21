import sqlite3

import pandas as pd
from pandas import DataFrame

from src.argument_parser import logger
from src.reaction_generator import decomposed_uid_mapping

print(decomposed_uid_mapping.head())


def create_pathway_pi(
    pathway_id,
    decomposed_uid_mapping: DataFrame,
    reaction_connections: DataFrame,
) -> DataFrame:
    """
    Create pathway_pi DataFrame based on decomposed_uid_mapping, reaction_connections.

    Args:
        decomposed_uid_mapping (DataFrame): DataFrame containing decomposed UID mapping.
        reaction_connections (DataFrame): DataFrame containing reaction connections.

    Returns:
        DataFrame: The created pathway_pi DataFrame.
    """
    logger.debug("Adding reaction pairs to pathway_pi")

    conn = sqlite3.connect("pathway_data.db")

    # Convert DataFrames to SQLite tables
    decomposed_uid_mapping.to_sql(
        "decomposed_uid_mapping_table", conn, if_exists="replace", index=False
    )
    reaction_connections.to_sql(
        "reaction_connections_table", conn, if_exists="replace", index=False
    )

    # Query SQLite tables to create pathway_pi DataFrame
    query = """
    SELECT *
    FROM decomposed_uid_mapping_table d
    JOIN reaction_connections_table r ON d.reactome_id = r.parent_reaction_id
    """
    pathway_pi = pd.read_sql_query(query, conn)

    conn.close()

    return pathway_pi
    # Assuming you have your decomposed_uid_mapping and reaction_connections DataFrames ready
    pathway_pi_result = create_pathway_pi(decomposed_uid_mapping, reaction_connections)
    print(pathway_pi_result)
