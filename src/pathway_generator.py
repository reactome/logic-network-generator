from src.argument_parser import logger
from src.neo4j_connector import get_reaction_connections
from src.pi_generator import create_pathway_pi
from src.reaction_generator import get_decomposed_uid_mapping


def generate_pathway_file(
    pathway_id: str, taxon_id: str, pathway_name: str, decompose: bool = False
) -> None:
    logger.debug(f"Generating {pathway_id} {pathway_name}")

    reaction_connections = get_reaction_connections(pathway_id)
    [decomposed_uid_mapping, best_matches] = get_decomposed_uid_mapping(
        pathway_id, reaction_connections.iloc[2].to_frame().T
    )
    print("dfdfd")
    print(decomposed_uid_mapping)
    create_pathway_pi(
        decomposed_uid_mapping, reaction_connections.iloc[2].to_frame().T, best_matches
    )
    exit()
