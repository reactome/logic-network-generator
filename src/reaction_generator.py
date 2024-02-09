import os
import itertools
import pandas as pd
import pprint
import uuid
from typing import List, Dict, Any, Tuple

from src.argument_parser import logger
from src.neo4j_connector import get_reaction_connections
from src.neo4j_connector import get_labels
from src.neo4j_connector import get_complex_components
from src.neo4j_connector import get_set_members
from src.neo4j_connector import get_reaction_input_output_ids
from src.neo4j_connector import get_reference_entities

pp = pprint.PrettyPrinter(indent=4)

decomposed_entity_uid_mapping = pd.DataFrame(columns=['uid', 'components', 'complex_id'])


def get_components_from_list(broken_apart_members: List[str]) -> List[str]:
    global decomposed_entity_uid_mapping
    components = []
    for member in broken_apart_members:
        if member in decomposed_entity_uid_mapping['uid'].values:
            member_components = decomposed_entity_uid_mapping.loc[
                decomposed_entity_uid_mapping['uid'] == member, 'components'].iloc[0]
            components += get_components_from_list(member_components)
        else:
            components.append(member)
    return components


def break_apart_entity(entity_id: str) -> List[List[str]]:
    global decomposed_entity_uid_mapping
    labels = get_labels(entity_id)
    if "EntitySet" in labels:
        member_ids = get_set_members(entity_id)
        broken_apart_members = []
        for member_id in member_ids:
            members = break_apart_entity(member_id)
            for member in members:
                broken_apart_members.append(member)
        return broken_apart_members
    elif "Complex" in labels:
        member_ids = get_complex_components(entity_id)
        broken_apart_members = []
        for member_id in member_ids:
            members = break_apart_entity(member_id)
            for member in members:
                broken_apart_members.append(member)
        if any(isinstance(member, list) for member in broken_apart_members):
            return [[entity_id]]
        else:
            uid = str(uuid.uuid4())
            components = get_components_from_list(broken_apart_members)
            decomposed_entity_uid_mapping = decomposed_entity_uid_mapping.append({
                'uid': uid,
                'components': components,
                'complex_id': entity_id
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


def generate_combinations(entity_ids: List[List[str]]) -> List[Tuple[str]]:
    decomposed_entities = []
    for entity_id in entity_ids:
        decomposed_entities.append(break_apart_entity(entity_id))
    return list(itertools.product(*decomposed_entities))


def create_entity_combinations_dict(reactions_entities: List[List[str]]) -> Dict[str, List[str]]:
    entity_combinations = {}
    for entities in reactions_entities:
        uid = str(uuid.uuid4())
        components = []
        for entity in entities:
            components += get_components_from_list(entity)
        entity_combinations[uid] = components
    return entity_combinations


def create_rows(reaction_id: str, decomposed_combinations: List[Dict[str, Any]], input_or_output: str) -> List[Dict[str, Any]]:
    rows = []
    for entities in decomposed_combinations:
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


def match_input_to_output(input_combination_key: str, input_combination_key_parts: List[str], output_combinations: Dict[str, Any]) -> List[Any]:
    best_match_count = 0
    output_entities = []
    for output_combination_key, output_combination_value in output_combinations.items():
        output_combination_key_parts = output_combination_key.split("-")
        elements_in_common = len(
            set(output_combination_key_parts) & set(input_combination_key_parts))
        if elements_in_common > best_match_count:
            output_entities = output_combination_value
            best_match_count = elements_in_common
    return output_entities


def matching_input_and_output_decomposed_reactions(reaction_id: str, input_combinations: Dict[str, Any], output_combinations: Dict[str, Any]) -> Dict[str, Any]:
    best_match_stats = {
        'num_inputs': None,
        'num_outputs': None,
        'num_matches': 0,
        'match_percentage': 0.0
    }
    match_stats_list = []
    for input_combination_key, input_entities in input_combinations.items():
        for output_combination_key, output_entities in output_combinations.items():
            common_ids = set(input_entities) & set(output_entities)
            num_matches = len(common_ids)
            num_inputs = len(input_entities)
            num_outputs = len(output_entities)
            match_percentage = num_matches / max(num_inputs, num_outputs) * 100 \
                if max(num_inputs, num_outputs) > 0 else 0.0
            match_stats = {
                'input_combination_key': input_combination_key,
                'output_combination_key': output_combination_key,
                'num_inputs': num_inputs,
                'num_outputs': num_outputs,
                'num_matches': num_matches,
                'match_percentage': match_percentage
                }
            match_stats_list.append(match_stats)
            if num_matches > best_match_stats['num_matches']:
                best_match_stats = {
                    'input_combination_key': input_combination_key,
                    'output_combination_key': output_combination_key,
                    'num_inputs': num_inputs,
                    'num_outputs': num_outputs,
                    'num_matches': num_matches,
                    'match_percentage': match_stats['match_percentage']
                }
    match_stats_df = pd.DataFrame(match_stats_list)
    match_stats_df.to_csv(
        f'match_stats_{reaction_id}.csv', index=False)
    return best_match_stats


def decompose_unmatched_entities_with_references(unmatched_entities: List[str], neo4j_connector: Any) -> Tuple[List[str], List[str]]:
    decomposed_entities = []
    reference_entities = []
    for entity_id in unmatched_entities:
        decomposed_entities.extend(break_apart_entity(entity_id))
        reference_df = get_reference_entities(entity_id)
        reference_entities.extend(reference_df['reference_entity_id'].tolist())
    return decomposed_entities, reference_entities


def get_reaction_inputs_and_outputs(reaction_ids: List[str]) -> pd.DataFrame:
    rows = []
    for reaction_id in reaction_ids:
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


def decompose_unmatched_entities(unmatched_entities: List[str]) -> List[str]:
    decomposed_entities = []
    for entity_id in unmatched_entities:
        decomposed_entities.extend(break_apart_entity(entity_id))
    return decomposed_entities


def get_reactions_df(pathway_id: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    reaction_connections_df = get_reaction_connections(pathway_id)
    reaction_ids = pd.unique(reaction_connections_df[['parent_reaction_id', 'child_reaction_id']].values.ravel('K'))
    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]
    reaction_inputs_and_outputs_df = None
    reaction_inputs_and_outputs_filename = 'reaction_inputs_and_outputs_df_' + str(pathway_id) + '.tsv'
    if os.path.isfile(reaction_inputs_and_outputs_filename):
        reaction_inputs_and_outputs_df = pd.read_table(reaction_inputs_and_outputs_filename, delimiter="\t")
    else:
        reaction_inputs_and_outputs_df = get_reaction_inputs_and_outputs(reaction_ids)
        reaction_inputs_and_outputs_df.to_csv(reaction_inputs_and_outputs_filename, sep="\t")
    return reaction_inputs_and_outputs_df, reaction_connections_df

