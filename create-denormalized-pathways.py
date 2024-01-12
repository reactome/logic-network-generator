#!/usr/bin/python

import os
from py2neo import Graph
import pandas as pd
import numpy as np
import pprint
import itertools
pp = pprint.PrettyPrinter(indent=4)

uri = "bolt://localhost:7687"
graph = Graph(uri, auth=('neo4j', 'test'))


def get_reaction_connections(pathway_id):
    query = """
       MATCH (pathway:Pathway)-[:hasEvent*]->(r1:ReactionLikeEvent)
           WHERE pathway.dbId = %s
       OPTIONAL MATCH (r1:ReactionLikeEvent)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway:Pathway)
           WHERE pathway.dbId = %s
       RETURN r1.dbId AS parent_reaction_id, r2.dbId AS child_reaction_id
    """ %  (pathway_id, pathway_id)
    df = pd.DataFrame(graph.run(query).data())
    df = df.astype({'parent_reaction_id': 'Int64',
                    'child_reaction_id': 'Int64'})
    return df


def get_all_pathways():
    query = """
        MATCH (p:Pathway)
        WHERE p.speciesName='Homo sapiens'
        RETURN
            p.stId AS id,
            p.name[0] AS name
        LIMIT 10
        """
    return graph.run(query).data();

def get_labels(entity_id):
    query_get_labels_template = """
       MATCH (e)
          WHERE e.dbId = %s
       RETURN labels(e) AS labels
       """
    query = query_get_labels_template % entity_id
    return graph.run(query).data()[0]["labels"]


def get_complex_components(entity_id):
    query_get_components_template = """
       MATCH (entity)-[:hasComponent]->(component)
           WHERE entity.dbId = %s
       RETURN collect(component.dbId) AS component_ids
       """
    query = query_get_components_template % entity_id

    return set(graph.run(query).data()[0]["component_ids"])

def get_set_members(entity_id):
    query_get_members_template = """
        MATCH (entity)-[:hasCandidate|hasMember]->(member)
            where entity.dbId = %s
        RETURN collect(member.dbId) as member_ids
        """
    query = query_get_members_template % entity_id

    return set(graph.run(query).data()[0]["member_ids"])


def get_reactions(pathway_id, taxon_id):
    query_reaction_template = """
        MATCH (reaction)<-[:hasEvent*]-(pathway:Pathway)-[:species]->(species:Species)
             WHERE (reaction:Reaction OR reaction:ReactionLikeEvent)
                   AND pathway.dbId=%s AND species.taxId="%s"
        RETURN COLLECT(reaction.dbId) AS reaction_ids
    """
    query = query_reaction_template % (pathway_id, taxon_id)

    return graph.run(query).data()[0]["reaction_ids"]


def get_reaction_input_ids(reaction_id):
    query_reaction_inputs_template = """
       MATCH (reaction)-[:input]-(input)
           WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) AND reaction.dbId=%s
       RETURN COLLECT(input.dbId) AS inputs
    """
    query = query_reaction_inputs_template % (reaction_id)

    return set(graph.run(query).data()[0]["inputs"])


def get_reaction_output_ids(reaction_id):
    query_reaction_outputs_template = """
       MATCH (reaction)-[:output]-(output)
            WHERE (reaction:Reaction or reaction:ReactionLikeEvent) AND reaction.dbId=%s
       RETURN COLLECT(output.dbId) AS outputs
    """
    query = query_reaction_outputs_template % (reaction_id)

    return set(graph.run(query).data()[0]["outputs"])

def add_entity_to_members(members, entity_id):
    broken_apart_members = []
    for member in members:
        member.append(entity_id)
        broken_apart_members.append(member)

    return broken_apart_members

def break_apart_complex(entity_id):
    component_ids = get_complex_components(entity_id)

    broken_apart_members = []
    for component_id in component_ids:
        members = break_apart_entity_based_on_sets(component_id)
        if not broken_apart_members:
            broken_apart_members = members
        else:
            current_members = broken_apart_members
            broken_apart_members = []
            for current_member in current_members:
                for member in members:
                    new_member_combination = current_member + member
                    broken_apart_members.append(new_member_combination)
    #pp.pprint(component_ids)
    #pp.pprint(broken_apart_members)
    #exit()
    if len(broken_apart_members) == 1:
       return [[entity_id]]
    else:
       return broken_apart_members

def break_apart_set(entity_id):
    member_ids = get_set_members(entity_id)
    broken_apart_members = []
    for member_id in member_ids:
        members = break_apart_entity_based_on_sets(member_id)
        for member in members:
            broken_apart_members.append(member)

    return broken_apart_members

def break_apart_entity_based_on_sets(entity_id):
    labels = get_labels(entity_id)
    if "Complex" in labels:
        return break_apart_complex(entity_id)
    elif "EntitySet" in labels:
        return break_apart_set(entity_id)
    elif "EntityWithAccessionedSequence" in labels:
        return [[entity_id]]
    elif "SimpleEntity" in labels:
        return [[entity_id]]
    elif "OtherEntity" in labels:
        return [[entity_id]]
    elif "GenomeEncodedEntity" in labels:
        return [[entity_id]]
    elif "Polymer" in labels:
        return [[entity_id]]
    elif "ChemicalDrug" in labels:
        return [[entity_id]]
    elif "Drug" in labels:
        return [[entity_id]]
    else:
      print("labels not handled")
      print(labels)
      print("for entity")
      print(entity_id)
      exit()

def add_outputs_for_reaction():
    print("adding output_reactions")


def add_reaction_pair(pathway_pi_df, reaction_pair):
    print("adding reaction pair")
    print(reaction_pair)
    exit()
    add_outputs_for_reaction(reaction_pair["parent_reaction_id"], )

def get_reaction_inputs_and_outputs(reaction_ids):

    print("Create reaction inputs and outputs dataframe")
    rows = []
    for reaction_id in reaction_ids:
        print(reaction_id)
        input_ids = get_reaction_input_ids(reaction_id)

        broken_apart_input_id_set = []
        for input_id in input_ids:
            broken_apart_input_id_set.append(break_apart_entity_based_on_sets(input_id))
        #pp.pprint(broken_apart_input_id_set)
        iterproduct_inputs = list(itertools.product(*broken_apart_input_id_set))
        #pp.pprint(iterproduct_inputs)
        input_combinations = {}
        for inputs in iterproduct_inputs:
            input_combinations["-".join(map(str, sorted(list(np.concatenate(inputs)))))] = inputs
        output_ids = get_reaction_output_ids(reaction_id)
        broken_apart_output_id_set = []
        for output_id in output_ids:
            broken_apart_output_id = break_apart_entity_based_on_sets(output_id)
            broken_apart_output_id_set.append(broken_apart_output_id)

        iterproduct_outputs = list(itertools.product(*broken_apart_output_id_set))
        output_combinations = {}
        for outputs in iterproduct_outputs:
            output_combinations["-".join(map(str, sorted(list(np.concatenate(outputs)))))] = outputs

        for input_combination_key, input_entities in input_combinations.items():
            for input_entity in input_entities:
                row = {"reaction_id": reaction_id,
                       "decomposed_reaction_id": input_combination_key,
                       "input_or_output": "input",
                       "decomposed_entity_id": "-".join(map(str, sorted(list(input_entity))))
                       }
                rows.append(row)
            if len(output_combinations.keys()) == 1:
                output_entities = list(output_combinations.values())[0]
            elif input_combination_key in output_combinations:
                output_entities = output_combinations[input_combination_key]
            else:
                input_combination_key_parts = input_combination_key.split("-")
                best_match_count = 0
                for output_combination_key, output_combination_value in output_combinations.items():
                    output_combination_key_parts = output_combination_key.split("-")
                    elements_in_common = len(set(output_combination_key_parts) & set(input_combination_key_parts))
                    if elements_in_common > best_match_count:
                        output_entities = output_combination_value
                        best_match_count = elements_in_common
                if best_match_count == 0:
                   print("cant match inputs with outputs for reaction")
                   print(reaction_id)
                   print(output_combinations)
                   print(input_combinations)
                   print(input_combination_key)
                   continue
            for output_entity in output_entities:
                row = {"reaction_id": reaction_id,
                       "decomposed_reaction_id": input_combination_key,
                       "input_or_output": "output",
                       "decomposed_entity_id": "-".join(map(str, sorted(list(output_entity))))
                       }
                rows.append(row)
    return pd.DataFrame.from_records(rows)


def create_pathway_pi_df(reaction_inputs_and_outputs_df, reaction_connections_df):
    print("Add reaction pairs to pathway_pi_df")

    columns = {"parent_id": pd.Series(dtype='Int64'),
               "parent_reaction_id": pd.Series(dtype='Int64'),
               "parent_decomposed_reaction_id": pd.Series(dtype='str'),
               "child_id": pd.Series(dtype='Int64'),
               "child_reaction_id": pd.Series(dtype='Int64'),
               "child_decomposed_reaction_id": pd.Series(dtype='str')}
    pathway_pi_df = pd.DataFrame(columns)

    for idx, reaction_connection in reaction_connections_df.iterrows():
        print("reaction_connection")
        print(reaction_connection)
        exit()



def generate_pathway_file(pathway_id, taxon_id, pathway_name):
    print("Generating " + pathway_id + " " + pathway_name)
    reaction_connections_df =  get_reaction_connections(pathway_id)
    reaction_ids = pd.unique(reaction_connections_df[['parent_reaction_id', 'child_reaction_id']].values.ravel('K'))
    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]#removing NA value from list

    reaction_inputs_and_outputs_filename = 'reaction_inputs_and_outputs_df_' +  pathway_id + '.tsv'
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
    #pathways = get_all_pathways();
    pathways = { "69620": "Cell_Cycle_Checkpoints",
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


main()
