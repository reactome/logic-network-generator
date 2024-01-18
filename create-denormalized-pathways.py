#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import logging
import itertools
from py2neo import Graph
import pandas as pd
import numpy as np
import pprint

pp = pprint.PrettyPrinter(indent=4)

# Configure the logging settings
logging.basicConfig(filename='debug_log.txt', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

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


def break_apart_entity(entity_id):
    labels = get_labels(entity_id)
    if "Complex" in labels or "EntitySet" in labels:
        member_ids = get_set_members(entity_id) if "EntitySet" in labels else get_complex_components(entity_id)
        broken_apart_members = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            for member in members:
                broken_apart_members.append(member)
        logger.debug(f"Debugging: break_apart_entity - entity_id: {entity_id}")
        logger.debug(f"Debugging: break_apart_entity - labels: {labels}")
        logger.debug(f"Debugging: break_apart_entity - broken_apart_members: {broken_apart_members}")

        return broken_apart_members if broken_apart_members else [[entity_id]]
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
        logger.error("Labels not handled")
        logger.error(f"Labels: {labels}")
        logger.error(f"For entity: {entity_id}")
        raise ValueError("Labels not handled for entity")


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


def create_entity_combinations_dict(combinations):
    entity_combinations = {}
    for combination in combinations:
        key = "-".join(map(str, sorted(list(np.concatenate(combination)))))
        entity_combinations[key] = combination
    return entity_combinations


def create_rows(reaction_id, decomposed_combinations, input_or_output):
    rows = []
    for key, entities in decomposed_combinations.items():
        for entity in entities:
            row = {
                "reaction_id": reaction_id,
                "decomposed_reaction_id": key,
                "input_or_output": input_or_output,
                "decomposed_entity_id": "-".join(map(str, sorted(list(entity))))
            }
            rows.append(row)
    return rows


def match_input_to_output(input_combination_key, input_combination_key_parts, output_combinations):
    best_match_count = 0
    output_entities = []

    for output_combination_key, output_combination_value in output_combinations.items():
        output_combination_key_parts = output_combination_key.split("-")
        elements_in_common = len(set(output_combination_key_parts) & set(input_combination_key_parts))

        if elements_in_common > best_match_count:
            output_entities = output_combination_value
            best_match_count = elements_in_common

    logger.debug(f"Debugging: match_input_to_output - input_combination_key: {input_combination_key}")
    logger.debug(f"Debugging: match_input_to_output - input_combination_key_parts: {input_combination_key_parts}")
    logger.debug(f"Debugging: match_input_to_output - best_match_count: {best_match_count}")
    logger.debug(f"Debugging: match_input_to_output - output_entities: {output_entities}")

    return output_entities


def get_reaction_inputs_and_outputs(reaction_ids):
    logger.debug("Creating reaction inputs and outputs dataframe")
    rows = []

    for reaction_id in reaction_ids:
        logger.debug(reaction_id)
        input_ids = get_reaction_input_output_ids(reaction_id, "input")

        broken_apart_input_id_set = [break_apart_entity(input_id) for input_id in input_ids]
        iterproduct_inputs = generate_combinations(broken_apart_input_id_set)
        input_combinations = create_entity_combinations_dict(iterproduct_inputs)

        output_ids = get_reaction_input_output_ids(reaction_id, "output")
        broken_apart_output_id_set = [break_apart_entity(output_id) for output_id in output_ids]
        iterproduct_outputs = generate_combinations(broken_apart_output_id_set)
        output_combinations = create_entity_combinations_dict(iterproduct_outputs)

        for input_combination_key, input_entities in input_combinations.items():
            rows.extend(create_rows(reaction_id, {"input": input_entities}, "input"))

            if len(output_combinations) == 1:
                output_entities = list(output_combinations.values())[0]
            elif input_combination_key in output_combinations:
                output_entities = output_combinations[input_combination_key]
            else:
                input_combination_key_parts = input_combination_key.split("-")
                output_entities = match_input_to_output(input_combination_key, input_combination_key_parts,
                                                        output_combinations)

            for output_entity in output_entities:
                rows.extend(create_rows(reaction_id, {"output": output_entity}, "output"))

    return pd.DataFrame.from_records(rows)


def create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df):
    logger.debug("Adding reaction pairs to pathway_pi_df")

    columns = {"parent_id": pd.Series(dtype='Int64'),
               "parent_reaction_id": pd.Series(dtype='Int64'),
               "parent_decomposed_reaction_id": pd.Series(dtype='str'),
               "child_id": pd.Series(dtype='Int64'),
               "child_reaction_id": pd.Series(dtype='Int64'),
               "child_decomposed_reaction_id": pd.Series(dtype='str')}
    pathway_pi_df = pd.DataFrame(columns)
    logger.debug(f"Debugging: create_pathway_pi_df - reaction_connections_df: {reaction_connections_df}")
    for idx, reaction_connection in reaction_connections_df.iterrows():
        logger.debug("reaction_connection")
        logger.debug(reaction_connection)
        exit()


def generate_pathway_file(pathway_id, taxon_id, pathway_name):
    logger.debug(f"Generating {pathway_id} {pathway_name}")
    reaction_connections_df = get_reaction_connections(pathway_id)
    reaction_ids = pd.unique(reaction_connections_df[['parent_reaction_id', 'child_reaction_id']].values.ravel('K'))
    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]  # removing NA value from list

    reaction_inputs_and_outputs_filename = 'reaction_inputs_and_outputs_df_' + pathway_id + '.tsv'
    if os.path.isfile(reaction_inputs_and_outputs_filename):
        reaction_inputs_and_outputs_df = pd.read_table(reaction_inputs_and_outputs_filename, delimiter="\t")
    else:
        reaction_inputs_and_outputs_df = get_reaction_inputs_and_outputs(reaction_ids)
        reaction_inputs_and_outputs_df.to_csv(reaction_inputs_and_outputs_filename, sep="\t")
    return
    pathway_pi_df = create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df)
    exit()


def main():
    taxon_id = "9606"
    # pathways = get_all_pathways()
    pathways = {"69620": "Cell_Cycle_Checkpoints",
                "5693567": "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA_",
                "453274": "Mitotic_G2-G2_M_phases",
                "68875": "Mitotic_Prophase",
                "453279": "Mitotic_G1-G1_S_phases",
                "1257604": "PIP3_activates_AKT_signaling",
                "5673001": "RAF_MAP_kinase_cascade",
                "1227986": "Signaling_by_ERBB2",
                "195721": "Signaling_by_WNT",
                "69242": "S_Phase",
                "3700989": "Transcriptional_Regulation_by_TP53"
                }

    for pathway_id, pathway_name in pathways.items():
        generate_pathway_file(pathway_id, taxon_id, pathway_name)


if __name__ == "__main__":
    main()
