from src.argument_parser import logger
from src.neo4j_connector import get_reaction_connections
from src.pi_generator import create_pathway_pi
from src.reaction_generator import get_decomposed_uid_mapping


def generate_pathway_file(
    pathway_id: str, taxon_id: str, pathway_name: str, decompose: bool = False
) -> None:
    logger.debug(f"Generating {pathway_id} {pathway_name}")
    print("pathway_id")
    print(pathway_id)
    reaction_connections = get_reaction_connections(pathway_id)

    number_of_reaction_connections: int = 2
    if number_of_reaction_connections > 0:
        reaction_connections = reaction_connections.iloc[
            :number_of_reaction_connections
        ]

    [decomposed_uid_mapping, best_matches] = get_decomposed_uid_mapping(
        pathway_id, reaction_connections
    )
    create_pathway_pi(decomposed_uid_mapping, reaction_connections, best_matches)
    exit()
