import os
import itertools
import pandas as pd
import pprint
import uuid

from src.argument_parser import logger
from src.neo4j_connector import get_labels
from src.neo4j_connector import get_complex_components
from src.neo4j_connector import get_set_members
from src.neo4j_connector import get_reaction_input_output_ids
from src.neo4j_connector import get_reference_entities
from src.best_reaction_match import find_best_reaction_match

pp = pprint.PrettyPrinter(indent=4)

decomposed_uid_mapping = pd.DataFrame(columns=['uid', 'component_id', 'input_or_output_id', 'reactome_id'])


def is_valid_uuid(value):
    try:
        uuid_obj = uuid.UUID(str(value), version=4)
        return str(uuid_obj) == value
    except ValueError:
        return False


def get_broken_apart_ids(broken_apart_members, reactome_id):
    if any(isinstance(member, set) for member in broken_apart_members):
        for i in range(len(broken_apart_members)):
            if not isinstance(broken_apart_members[i], set):
                broken_apart_members[i] = {broken_apart_members[i]}
        iterproduct_components = list(itertools.product(*broken_apart_members))

        return get_uid_for_iterproduct_components(iterproduct_components, reactome_id)
    else:
        return reactome_id


def get_or_assign_uid(input_or_output_ids):
    global decomposed_uid_mapping

    uid_to_input_or_output = {}
    for index, row in decomposed_uid_mapping.iterrows():
        uid = row['uid']
        input_or_output_id = row['input_or_output_id']

        if uid in uid_to_input_or_output:
            uid_to_input_or_output[uid].add(input_or_output_id)
        else:
            uid_to_input_or_output[uid] = {input_or_output_id}

    matching_uid = None
    for uid, input_or_output_set in uid_to_input_or_output.items():
        if input_or_output_set == input_or_output_ids:
            matching_uid = uid
            break

    return matching_uid if matching_uid else str(uuid.uuid4())


def get_uid_for_iterproduct_components(iterproduct_components, reactome_id):
    global decomposed_uid_mapping
    uids = set()
    for component in iterproduct_components:
        component_to_input_or_output = dict()
        for item in component:
            if is_valid_uuid(item):
                selected_rows = decomposed_uid_mapping.loc[decomposed_uid_mapping['uid'] == item]
                for index, selected_row in selected_rows.iterrows():
                    component_id = selected_row['component_id']
                    component_to_input_or_output[component_id] = item
            else:
                component_to_input_or_output[item] = item

        uid = get_or_assign_uid(iterproduct_components)

        rows = []
        for component_id, input_or_output_id in component_to_input_or_output.items():
            row = {
                'uid': uid,
                'component_id': component_id,
                'input_or_output_id': input_or_output_id,
                'reactome_id': reactome_id}
            rows.append(row)

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])
        uids.add(uid)

    return uids


def break_apart_entity(entity_id):
    global decomposed_uid_mapping
    labels = get_labels(entity_id)

    if "EntitySet" in labels:
        member_ids = get_set_members(entity_id)
        broken_apart_members = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            if isinstance(members, set):
                for member in members:
                    broken_apart_members.append(member)
            else:
                broken_apart_members.append(members)

        logger.debug(f"Debugging: break_apart_entity - entity_id: {entity_id}")
        logger.debug(f"Debugging: break_apart_entity - labels: {labels}")
        logger.debug(f"Debugging: break_apart_entity - broken_apart_members: {broken_apart_members}")

        return set(broken_apart_members)

    elif "Complex" in labels:
        member_ids = get_complex_components(entity_id)
        broken_apart_members = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            broken_apart_members.append(members)

        logger.debug(f"Debugging: break_apart_entity - entity_id: {entity_id}")
        logger.debug(f"Debugging: break_apart_entity - labels: {labels}")
        logger.debug(f"Debugging: break_apart_entity - broken_apart_members: {broken_apart_members}")

        return get_broken_apart_ids(broken_apart_members, entity_id)

    elif any(entity_label in labels for entity_label in [
            "ChemicalDrug",
            "Drug",
            "EntityWithAccessionedSequence",
            "GenomeEncodedEntity",
            "OtherEntity",
            "Polymer",
            "SimpleEntity"]):

        return entity_id

    else:
        logger.error(f"Not handling labels correctly for: {entity_id}")
        exit(1)


def generate_combinations(entity_ids):
    decomposed_entities = []
    for entity_id in entity_ids:
        decomposed_entities.append(break_apart_entity(entity_id))
    return decomposed_entities


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


def match_input_and_output_decomposed_reactions(reaction_id, input_reactions, output_reactions):
    global decomposed_uid_mapping

    [best_matches, match_counts] = find_best_reaction_match(input_reactions, output_reactions, decomposed_uid_mapping)
    print("best_matches, match_counts")
    print(best_matches, match_counts)
    # create stats and return with best matches

    exit()

    return best_matches


def decompose_unmatched_entities_with_references(unmatched_entities, neo4j_connector):
    decomposed_entities = []
    reference_entities = []

    for entity_id in unmatched_entities:
        decomposed_entities.extend(break_apart_entity(entity_id))

        # Query Neo4j for reference entities
        reference_df = get_reference_entities(entity_id)
        reference_entities.extend(reference_df['reference_entity_id'].tolist())

    return decomposed_entities, reference_entities


def create_reaction_df(reaction_ids):
    logger.debug("Creating reaction inputs and outputs dataframe")
    rows = []

    for reaction_id in reaction_ids:
        logger.debug(reaction_id)
        input_ids = get_reaction_input_output_ids(
            reaction_id, "input")

        broken_apart_input_id_set = [
            break_apart_entity(input_id) for input_id in input_ids]

        input_combinations = get_broken_apart_ids(
            broken_apart_input_id_set, reaction_id)

        output_ids = get_reaction_input_output_ids(
            reaction_id, "output")
        broken_apart_output_id_set = [
            break_apart_entity(output_id) for output_id in output_ids]
        output_combinations = get_broken_apart_ids(
            broken_apart_output_id_set, reaction_id)

        reaction_rows = match_input_and_output_decomposed_reactions(
            reaction_id, input_combinations, output_combinations)
        rows.append(reaction_rows)
    return pd.DataFrame.from_records(rows)


def decompose_unmatched_entities(unmatched_entities):
    decomposed_entities = []
    for entity_id in unmatched_entities:
        decomposed_entities.extend(break_apart_entity(entity_id))
    return decomposed_entities


def get_reactions(pathway_id, reaction_connections_df):
    reaction_ids = pd.unique(reaction_connections_df[['parent_reaction_id', 'child_reaction_id']].values.ravel('K'))
    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]  # removing NA value from list

    reaction_df = None

    reaction_filename = 'reaction_df_' + str(pathway_id) + '.tsv'
    if os.path.isfile(reaction_filename):
        reaction_df = pd.read_table(reaction_filename, delimiter="\t")
    else:
        reaction_df = create_reaction_df(reaction_ids)
        reaction_df.to_csv(reaction_filename, sep="\t")

    return reaction_df
