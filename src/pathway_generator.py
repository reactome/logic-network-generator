import os

import pandas as pd

from src.argument_parser import logger
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types
from src.logic_network_generator import create_pathway_logic_network
from src.neo4j_connector import get_reaction_connections
from src.reaction_generator import get_decomposed_uid_mapping


def generate_pathway_file(
    pathway_id: str, taxon_id: str, pathway_name: str, decompose: bool = False
) -> None:
    logger.debug(f"Generating {pathway_id} {pathway_name}")
    print("pathway_id")
    print(pathway_id)

    # Define filenames for caching
    reaction_connections_file = f"reaction_connections_{pathway_id}.csv"
    decomposed_uid_mapping_file = f"decomposed_uid_mapping_{pathway_id}.csv"
    best_matches_file = f"best_matches_{pathway_id}.csv"

    if os.path.exists(reaction_connections_file):
        reaction_connections = pd.read_csv(reaction_connections_file)
    else:
        reaction_connections = get_reaction_connections(pathway_id)
        reaction_connections.to_csv(reaction_connections_file, index=False)

    number_of_reaction_connections: int = -1
    if number_of_reaction_connections > 0:
        reaction_connections = reaction_connections.iloc[
            :number_of_reaction_connections
        ]

    if os.path.exists(decomposed_uid_mapping_file) & os.path.exists(best_matches_file):
        decomposed_uid_mapping = pd.read_csv(
            decomposed_uid_mapping_file, dtype=decomposed_uid_mapping_column_types
        )
        best_matches = pd.read_csv(best_matches_file)
    else:
        [decomposed_uid_mapping, best_matches_list] = get_decomposed_uid_mapping(
            pathway_id, reaction_connections
        )
        best_matches = pd.DataFrame(
            best_matches_list, columns=["incomming", "outgoing"]
        )
        decomposed_uid_mapping.to_csv(decomposed_uid_mapping_file, index=False)
        best_matches.to_csv(best_matches_file, index=False)

    create_pathway_logic_network(
        decomposed_uid_mapping, reaction_connections, best_matches
    )
