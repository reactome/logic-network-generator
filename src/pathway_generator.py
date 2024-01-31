import pprint

from src.argument_parser import logger
from src.reaction_generator import get_reactions_df
from src.pi_generator import create_pathway_pi_df

pp = pprint.PrettyPrinter(indent=4)


def generate_pathway_file(pathway_id, taxon_id, pathway_name, decompose=False):
    logger.debug(f"Generating {pathway_id} {pathway_name}")
    [reaction_inputs_and_outputs_df, reaction_connections_df] = get_reactions_df(pathway_id)
    pathway_pi_df = create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df)
    pathway_pi_df.to_csv('pathway_pi_' + pathway_id + '.csv', index=False)
    exit()
