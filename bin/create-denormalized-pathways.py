#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import logging
import itertools
from py2neo import Graph
import pandas as pd
import numpy as np
import pprint
import uuid


import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.argument_parser import parse_args, configure_logging

pp = pprint.PrettyPrinter(indent=4)

decomposed_entity_uid_mapping = pd.DataFrame(columns=['uid', 'components', 'complex_id'])

uri = "bolt://localhost:7687"
graph = Graph(uri, auth=('neo4j', 'test'))

# Define a logger
logger = logging.getLogger(__name__)


def get_reaction_connections(pathway_id):
    query = """
       MATCH (pathway:Pathway)-[:hasEvent*]->(r1:ReactionLikeEvent)
           WHERE pathway.dbId = %s
       OPTIONAL MATCH (r1:ReactionLikeEvent)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway:Pathway)
           WHERE pathway.dbId = %s
       RETURN r1.dbId AS parent_reaction_id, r2.dbId AS child_reaction_id
    """ % (pathway_id, pathway_id)

    try:
        df = pd.DataFrame(graph.run(query).data())
        df = df.astype({'parent_reaction_id': 'Int64', 'child_reaction_id': 'Int64'})
        return df
    except Exception:
        logger.error("Error in get_reaction_connections", exc_info=True)
        raise


def get_all_pathways():
    query = """
        MATCH (p:Pathway)
        WHERE p.speciesName='Homo sapiens'
        RETURN
            p.stId AS id,
            p.name[0] AS name
        LIMIT 10
        """

    try:
        return graph.run(query).data()
    except Exception:
        logger.error("Error in get_all_pathways", exc_info=True)
        raise


def get_labels(entity_id):
    query_get_labels_template = """
       MATCH (e)
          WHERE e.dbId = %s
       RETURN labels(e) AS labels
       """
    query = query_get_labels_template % entity_id

    try:
        return graph.run(query).data()[0]["labels"]
    except Exception:
        logger.error("Error in get_labels", exc_info=True)
        raise


def get_complex_components(entity_id):
    query_get_components_template = """
       MATCH (entity)-[:hasComponent]->(component)
           WHERE entity.dbId = %s
       RETURN collect(component.dbId) AS component_ids
       """
    query = query_get_components_template % entity_id

    try:
        return set(graph.run(query).data()[0]["component_ids"])
    except Exception:
        logger.error("Error in get_complex_components", exc_info=True)
        raise


def get_set_members(entity_id):
    query_get_members_template = """
        MATCH (entity)-[:hasCandidate|hasMember]->(member)
            WHERE entity.dbId = %s
        RETURN collect(member.dbId) as member_ids
        """
    query = query_get_members_template % entity_id

    try:
        return set(graph.run(query).data()[0]["member_ids"])
    except Exception:
        logger.error("Error in get_set_members", exc_info=True)
        raise


# def decompose_sets(entity_ids):
# decomposed_entities = []
# for entity_id in entity_ids:
# decomposed_entities.append(break_apart_entity(entity_id))
# return list(itertools.product(*decomposed_entities))


def get_reactions(pathway_id, taxon_id):
    query_reaction_template = """
        MATCH (reaction)<-[:hasEvent*]-(pathway:Pathway)-[:species]->(species:Species)
             WHERE (reaction:Reaction OR reaction:ReactionLikeEvent)
                   AND pathway.dbId=%s AND species.taxId="%s"
        RETURN COLLECT(reaction.dbId) AS reaction_ids
    """
    query = query_reaction_template % (pathway_id, taxon_id)

    try:
        return graph.run(query).data()[0]["reaction_ids"]
    except Exception:
        logger.error("Error in get_reactions", exc_info=True)
        raise


def get_reaction_input_output_ids(reaction_id, input_or_output):
    query_template = """
       MATCH (reaction)-[:%s]-(io)
           WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) AND reaction.dbId=%s
       RETURN COLLECT(io.dbId) AS io_ids
    """
    relation_type = "input" if input_or_output == "input" else "output"
    query = query_template % (relation_type, reaction_id)

    try:
        return set(graph.run(query).data()[0]["io_ids"])
    except Exception:
        logger.error("Error in get_reaction_input_output_ids", exc_info=True)
        raise


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

            # Create a table with the number of inputs, number of outputs, and number of matches
            match_stats = {
                'input_combination_key': input_combination_key,
                'output_combination_key': output_combination_key,
                'num_inputs': num_inputs,
                'num_outputs': num_outputs,
                'num_matches': num_matches,
                'match_percentage': num_matches / max(num_inputs, num_outputs) * 100 if max(num_inputs, num_outputs) > 0 else 0.0
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

def get_reference_entities(entity_id):
    query = """
        MATCH (e:Entity {entity_id: %s})-[:HAS_REFERENCE_ENTITY]->(ref:Entity)
        RETURN ref.entity_id AS reference_entity_id
    """ % entity_id

    try:
        df = pd.DataFrame(graph.run(query).data())
        df = df.astype({'reference_entity_id': 'str'})
        return df
    except Exception:
        logger.error("Error in get_reference_entities", exc_info=True)
        raise
        
        
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
            (reaction_inputs_and_outputs_df['reaction_id'] == parent_reaction_id) & (reaction_inputs_and_outputs_df['input_or_output'] == 'input')
        ]['decomposed_entity_id'].values

        parent_outputs = reaction_inputs_and_outputs_df[
            (reaction_inputs_and_outputs_df['reaction_id'] == parent_reaction_id) & (reaction_inputs_and_outputs_df['input_or_output'] == 'output')
        ]['decomposed_entity_id'].values

        child_inputs = reaction_inputs_and_outputs_df[
            (reaction_inputs_and_outputs_df['reaction_id'] == child_reaction_id) & (reaction_inputs_and_outputs_df['input_or_output'] == 'input')
        ]['decomposed_entity_id'].values

        child_outputs = reaction_inputs_and_outputs_df[
            (reaction_inputs_and_outputs_df['reaction_id'] == child_reaction_id) & (reaction_inputs_and_outputs_df['input_or_output'] == 'output')
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


def decompose_unmatched_ids(pathway_pi_df):
    for idx, row in pathway_pi_df.iterrows():
        unmatched_inputs = row['unmatched_inputs'].split('-') if row['unmatched_inputs'] else []
        unmatched_outputs = row['unmatched_outputs'].split('-') if row['unmatched_outputs'] else []

        decomposed_unmatched_inputs = decompose_unmatched_entities(unmatched_inputs)
        decomposed_unmatched_outputs = decompose_unmatched_entities(unmatched_outputs)
        common_reference_entities = set(decomposed_unmatched_inputs) & set(decomposed_unmatched_outputs)

        pathway_pi_df.at[idx, 'decomposed_unmatched_inputs'] = '-'.join(map(str, sorted(list(decomposed_unmatched_inputs))))
        pathway_pi_df.at[idx, 'decomposed_unmatched_outputs'] = '-'.join(map(str, sorted(list(decomposed_unmatched_outputs))))
        pathway_pi_df.at[idx, 'common_reference_entities'] = '-'.join(map(str, sorted(list(common_reference_entities))))

    return pathway_pi_df


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
    exit()


def main():
    # parse command line arguments
    args = parse_args()

    # configure logging based on debug flag
    configure_logging(args.debug, args.verbose)

    taxon_id = "9606"

    if args.input_file:
        # Read pathways from the input file
        try:
            pathways_df = pd.read_csv(args.input_file, sep='\t')
            pathways = dict(zip(pathways_df['ID'], pathways_df['PathwayName']))
        except Exception as e:
            logger.error(f"Error reading input file: {e}")
            return
    else:
        logger.error("Input file (--input_file) is required.")
        return

    # create a .tsv file for pathways list
    pathways_list_df = pd.DataFrame(list(pathways.items()), columns=['ID', 'PathwayName'])
    pathways_list_df.to_csv(args.output, sep='\t', index=False)

    for pathway_id, pathway_name in pathways.items():
        generate_pathway_file(pathway_id, taxon_id, pathway_name, decompose=args.decompose)

        reaction_connections_df = get_reaction_connections(pathway_id)
        reaction_ids = reaction_connections_df['reaction_id'].unique()
        reaction_inputs_and_outputs_df = get_reaction_inputs_and_outputs(reaction_ids)

        pathway_pi_df = create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df)
        pathway_pi_df = decompose_unmatched_ids(pathway_pi_df)

        print(f"Unmatched IDs for Pathway {pathway_id}:")
        for idx, row in pathway_pi_df.iterrows():
            unmatched_inputs = row['unmatched_inputs'].split('-') if row['unmatched_inputs'] else []
            unmatched_outputs = row['unmatched_outputs'].split('-') if row['unmatched_outputs'] else []

            print(f"Pathway {pathway_id} - Unmatched Inputs: {unmatched_inputs}")
            print(f"Pathway {pathway_id} - Unmatched Outputs: {unmatched_outputs}")

        # Now you can use pathway_pi_df for further analysis or export to a file
        pathway_pi_df.to_csv('pathway_pi_df.csv', index=False)


if __name__ == "__main__":
    main()
