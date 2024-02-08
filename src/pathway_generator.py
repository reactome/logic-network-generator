from src.argument_parser import logger
from src.reaction_generator import get_decomposed_uid_mapping
from src.pi_generator import create_pathway_pi
from src.neo4j_connector import get_reaction_connections


def generate_pathway_file(pathway_id, taxon_id, pathway_name, decompose=False):
    logger.debug(f"Generating {pathway_id} {pathway_name}")

    reaction_connections = get_reaction_connections(pathway_id)
    [decomposed_uid_mapping, best_matches] = get_decomposed_uid_mapping(pathway_id,
                                                                        reaction_connections.iloc[2].to_frame().T)
    create_pathway_pi(decomposed_uid_mapping, reaction_connections.iloc[2].to_frame().T, best_matches)
    exit()
