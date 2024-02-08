import pandas as pd

from src.argument_parser import logger

def create_pathway_pi(decomposed_uid_mapping, reaction_connections, best_matches):
    logger.debug("Adding reaction pairs to pathway_pi")

    columns = {
        "parent_id": pd.Series(dtype='str'),
        "child_id": pd.Series(dtype='str'),
        "pos_neg": pd.Series(dtype='str'),
        "and_or": pd.Series(dtype='str'),
    }
    pathway_pi = pd.DataFrame(columns)

    print(reaction_connections)

    for idx, reaction_connection in reaction_connections.iterrows():
        parent_reaction_id = reaction_connection['parent_reaction_id']
        child_reaction_id = reaction_connection['child_reaction_id']

        print("parent_reaction_id")
        print(parent_reaction_id)
        print("child_reaction_id")
        print(child_reaction_id)

        # Assuming parent_reaction_id is a variable containing the value to filter by
        parent_rows = decomposed_uid_mapping[decomposed_uid_mapping['reactome_id'] == parent_reaction_id]
        child_rows = decomposed_uid_mapping[decomposed_uid_mapping['reactome_id'] == child_reaction_id]
        print(parent_rows)
        print(child_rows)
        print(best_matches)


    return pathway_pi
