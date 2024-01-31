import pandas as pd
import pprint

from src.argument_parser import logger

pp = pprint.PrettyPrinter(indent=4)


def create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df):
    logger.debug("Adding reaction pairs to pathway_pi_df")

    columns = {
        "parent_id": pd.Series(dtype='Int64'),
        "parent_reaction_id": pd.Series(dtype='Int64'),
        "parent_decomposed_reaction_id": pd.Series(dtype='str'),
        "child_id": pd.Series(dtype='Int64'),
        "child_reaction_id": pd.Series(dtype='Int64'),
        "child_decomposed_reaction_id": pd.Series(dtype='str'),
        "common_ids": pd.Series(dtype='str'),  # Common IDs between inputs and outputs
        "unmatched_inputs": pd.Series(dtype='str'),  # Unmatched input IDs
        "unmatched_outputs": pd.Series(dtype='str')  # Unmatched output IDs
    }
    pathway_pi_df = pd.DataFrame(columns)

    for idx, reaction_connection in reaction_connections_df.iterrows():
        parent_reaction_id = reaction_connection['parent_reaction_id']
        child_reaction_id = reaction_connection['child_reaction_id']

        parent_inputs = reaction_inputs_and_outputs_df[
            (reaction_inputs_and_outputs_df['reaction_id'] == parent_reaction_id)
            & (reaction_inputs_and_outputs_df['input_or_output'] == 'input')
        ]['decomposed_entity_id'].values

        child_outputs = reaction_inputs_and_outputs_df[
            (reaction_inputs_and_outputs_df['reaction_id'] == child_reaction_id)
            & (reaction_inputs_and_outputs_df['input_or_output'] == 'output')
        ]['decomposed_entity_id'].values

        common_ids = set(parent_inputs) & set(child_outputs)
        unmatched_inputs = set(parent_inputs) - common_ids
        unmatched_outputs = set(child_outputs) - common_ids

        row = {
            "parent_id": reaction_connection['parent_reaction_id'],
            "parent_reaction_id": parent_reaction_id,
            "parent_decomposed_reaction_id": reaction_connection['parent_decomposed_reaction_id'],
            "child_id": reaction_connection['child_reaction_id'],
            "child_reaction_id": child_reaction_id,
            "child_decomposed_reaction_id": reaction_connection['child_decomposed_reaction_id'],
            "common_ids": '-'.join(map(str, sorted(list(common_ids)))),
            "unmatched_inputs": '-'.join(map(str, sorted(list(unmatched_inputs)))),
            "unmatched_outputs": '-'.join(map(str, sorted(list(unmatched_outputs))))
        }

        pathway_pi_df = pathway_pi_df.append(row, ignore_index=True)

    return pathway_pi_df
