import pprint

from src.argument_parser import logger
from src.reaction_generator import get_reactions
from src.pi_generator import create_pathway_pi_df
from src.neo4j_connector import get_reaction_connections

pp = pprint.PrettyPrinter(indent=4)


def generate_pathway_file(pathway_id, taxon_id, pathway_name, decompose=False):
    logger.debug(f"Generating {pathway_id} {pathway_name}")

    reaction_connections_df = get_reaction_connections(pathway_id)
    reaction_df = get_reactions(pathway_id, reaction_connections_df)
    pathway_pi_df = create_pathway_pi_df(reaction_df, reaction_connections_df)
    pathway_pi_df.to_csv('pathway_pi_' + pathway_id + '.csv', index=False)
