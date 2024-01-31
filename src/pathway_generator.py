import os
import itertools
import pandas as pd
import pprint
import uuid

from src.argument_parser import logger
from src.neo4j_connector import get_reaction_connections
from src.neo4j_connector import get_labels
from src.neo4j_connector import get_complex_components
from src.neo4j_connector import get_set_members
from src.neo4j_connector import get_reaction_input_output_ids
from src.neo4j_connector import get_reference_entities


pp = pprint.PrettyPrinter(indent=4)

decomposed_entity_uid_mapping = pd.DataFrame(columns=['uid', 'components', 'complex_id'])


def get_components_from_list(broken_apart_members):
    global decomposed_entity_uid_mapping

    # Initialize an empty list to store components
    components = []

    # Iterate over members in the list
    for member in broken_apart_members:
        # Check if the member is a UID in the DataFrame
        if member in decomposed_entity_uid_mapping['uid'].values:
            # If yes, append the components from the DataFrame
            member_components = decomposed_entity_uid_mapping.loc[
                decomposed_entity_uid_mapping['uid'] == member, 'components'].iloc[0]
            components += get_components_from_list(member_components)
        else:
            # If not, it's a leaf component, so append it
            components.append(member)

    return components


def break_apart_entity(entity_id):
    global decomposed_entity_uid_mapping
    labels = get_labels(entity_id)

    if "EntitySet" in labels:
        member_ids = get_set_members(entity_id)
        broken_apart_members = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            for member in members:
                broken_apart_members.append(member)

        logger.debug(f"Debugging: break_apart_entity - entity_id: {entity_id}")
        logger.debug(f"Debugging: break_apart_entity - labels: {labels}")
        logger.debug(f"Debugging: break_apart_entity - broken_apart_members: {broken_apart_members}")

        return broken_apart_members

    elif "Complex" in labels:
        member_ids = get_complex_components(entity_id)
        broken_apart_members = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            for member in members:
                broken_apart_members.append(member)

        logger.debug(f"Debugging: break_apart_entity - entity_id: {entity_id}")
        logger.debug(f"Debugging: break_apart_entity - labels: {labels}")
        logger.debug(f"Debugging: break_apart_entity - broken_apart_members: {broken_apart_members}")

        if set(broken_apart_members) == set(member_ids):
            return [[entity_id]]
        else:
            uid = str(uuid.uuid4())
            logger.debug(
                f"Generated UID {uid} for entity with different broken_apart_members: {entity_id}")

            components = get_components_from_list(broken_apart_members)
            decomposed_entity_uid_mapping = decomposed_entity_uid_mapping.append({
                'uid': uid,
                'components': components,
                'complex_id': entity_id  # Assuming entity_id is the complex Reactome ID
            }, ignore_index=True)

            return [[uid]]

    elif any(entity_label in labels for entity_label in [
            "ChemicalDrug",
            "Drug",
            "EntityWithAccessionedSequence",
            "GenomeEncodedEntity",
            "OtherEntity",
            "Polymer",
            "SimpleEntity"]):
        return [[entity_id]]
    else:
        logger.error(f"Not handling labels correctly for: {entity_id}")
        exit(1)


def add_outputs_for_reaction():
    logger.debug("Adding output_reactions")


def add_reaction_pair(pathway_pi_df, reaction_pair):
    logger.debug("Adding reaction pair")
    logger.debug(reaction_pair)
    exit()
    add_outputs_for_reaction(reaction_pair["parent_reaction_id"], )


def generate_combinations(entity_ids):
    decomposed_entities = []
    for entity_id in entity_ids:
        decomposed_entities.append(break_apart_entity(entity_id))
    return list(itertools.product(*decomposed_entities))


def create_entity_combinations_dict(reactions_entities):
    entity_combinations = {}
    for entities in reactions_entities:
        uid = str(uuid.uuid4())
        components = []
        for entity in entities:
            components += get_components_from_list(entity)
        entity_combinations[uid] = components
    return entity_combinations


def create_rows(reaction_id, decomposed_combinations, input_or_output):
    rows = []
    for entities in decomposed_combinations():
        for entity in entities:
            row = {
                "reaction_id": reaction_id,
                "decomposed_reaction_id": str(uuid.uuid4()),
                "input_or_output": input_or_output,
                "decomposed_entity_id": "-".join(map(str, sorted(list(entity['reactome_id'])))),
                "reactome_id": entity['reactome_id'],
            }
            rows.append(row)
    return rows


def match_input_to_output(input_combination_key, input_combination_key_parts, output_combinations):
    best_match_count = 0
    output_entities = []

    for output_combination_key, output_combination_value in output_combinations.items():
        output_combination_key_parts = output_combination_key.split("-")
        elements_in_common = len(
            set(output_combination_key_parts) & set(input_combination_key_parts))

        if elements_in_common > best_match_count:
            output_entities = output_combination_value
            best_match_count = elements_in_common

    logger.debug(
        f"Debugging: match_input_to_output - input_combination_key: {input_combination_key}")
    logger.debug(
        f"Debugging: match_input_to_output - input_combination_key_parts: {input_combination_key_parts}")
    logger.debug(
        f"Debugging: match_input_to_output - best_match_count: {best_match_count}")
    logger.debug(
        f"Debugging: match_input_to_output - output_entities: {output_entities}")

    return output_entities


def matching_input_and_output_decomposed_reactions(reaction_id, input_combinations, output_combinations):
    best_match_stats = {
        'num_inputs': None,
        'num_outputs': None,
        'num_matches': 0,
        'match_percentage': 0.0
    }

    match_stats_list = []

    for input_combination_key, input_entities in input_combinations.items():
        for output_combination_key, output_entities in output_combinations.items():
            # Compare input_entities and output_entities to see how well they match
            common_ids = set(input_entities) & set(output_entities)
            num_matches = len(common_ids)
            num_inputs = len(input_entities)
            num_outputs = len(output_entities)

            match_percentage = num_matches / max(num_inputs, num_outputs) * 100 \
                if max(num_inputs, num_outputs) > 0 else 0.0

            # Create a table with the number of inputs, number of outputs, and number of matches
            match_stats = {
                'input_combination_key': input_combination_key,
                'output_combination_key': output_combination_key,
                'num_inputs': num_inputs,
                'num_outputs': num_outputs,
                'num_matches': num_matches,
                'match_percentage': match_percentage
                }

            match_stats_list.append(match_stats)

            # Update best_match_stats if the current match is better
            if num_matches > best_match_stats['num_matches']:
                best_match_stats = {
                    'input_combination_key': input_combination_key,
                    'output_combination_key': output_combination_key,
                    'num_inputs': num_inputs,
                    'num_outputs': num_outputs,
                    'num_matches': num_matches,
                    'match_percentage': match_stats['match_percentage']
                }

    # Create a DataFrame of match statistics
    match_stats_df = pd.DataFrame(match_stats_list)

    # Now you can use match_stats_df for further analysis or export to a file
    match_stats_df.to_csv(
        f'match_stats_{reaction_id}.csv', index=False)

    return best_match_stats


def decompose_unmatched_entities_with_references(unmatched_entities, neo4j_connector):
    decomposed_entities = []
    reference_entities = []

    for entity_id in unmatched_entities:
        decomposed_entities.extend(break_apart_entity(entity_id))

        # Query Neo4j for reference entities
        reference_df = get_reference_entities(entity_id)
        reference_entities.extend(reference_df['reference_entity_id'].tolist())

    return decomposed_entities, reference_entities


def get_reaction_inputs_and_outputs(reaction_ids):
    logger.debug("Creating reaction inputs and outputs dataframe")
    rows = []

    for reaction_id in reaction_ids:
        logger.debug(reaction_id)
        input_ids = get_reaction_input_output_ids(
            reaction_id, "input")

        broken_apart_input_id_set = [
            break_apart_entity(input_id) for input_id in input_ids]
        iterproduct_inputs = generate_combinations(
            broken_apart_input_id_set)
        input_combinations = create_entity_combinations_dict(
            iterproduct_inputs)

        output_ids = get_reaction_input_output_ids(
            reaction_id, "output")
        broken_apart_output_id_set = [
            break_apart_entity(output_id) for output_id in output_ids]
        iterproduct_outputs = generate_combinations(
            broken_apart_output_id_set)
        output_combinations = create_entity_combinations_dict(
            iterproduct_outputs)

        reaction_rows = matching_input_and_output_decomposed_reactions(
            reaction_id, input_combinations, output_combinations)
        rows.append(reaction_rows)
    return pd.DataFrame.from_records(rows)


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


def decompose_unmatched_entities(unmatched_entities):
    decomposed_entities = []
    for entity_id in unmatched_entities:
        decomposed_entities.extend(break_apart_entity(entity_id))
    return decomposed_entities


def generate_pathway_file(pathway_id, taxon_id, pathway_name, decompose=False):
    logger.debug(f"Generating {pathway_id} {pathway_name}")
    reaction_connections_df = get_reaction_connections(pathway_id)
    reaction_ids = pd.unique(reaction_connections_df[['parent_reaction_id', 'child_reaction_id']].values.ravel('K'))
    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]  # removing NA value from list

    reaction_inputs_and_outputs_filename = 'reaction_inputs_and_outputs_df_' + pathway_id + '.tsv'
    if os.path.isfile(reaction_inputs_and_outputs_filename):
        reaction_inputs_and_outputs_df = pd.read_table(reaction_inputs_and_outputs_filename, delimiter="\t")

        reaction_inputs_and_outputs_df = get_reaction_inputs_and_outputs(reaction_ids)
        reaction_inputs_and_outputs_df.to_csv(reaction_inputs_and_outputs_filename, sep="\t")

    pathway_pi_df = create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df)
    pathway_pi_df.to_csv('pathway_pi_' + pathway_id + '.csv', index=False)
    exit()
