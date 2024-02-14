import itertools
import uuid
from typing import Any, Dict, List, Set, Tuple, Union

import pandas as pd

from src.argument_parser import logger
from src.best_reaction_match import find_best_reaction_match
from src.neo4j_connector import (get_complex_components, get_labels,
                                 get_reaction_input_output_ids,
                                 get_set_members)

# Define types
UID = str
ComponentID = str
InputOutputID = str
ReactomeID = str
DataFrameRow = Dict[str, Any]

column_types = {
    "uid": str,
    "reactome_id": int,
    "component_id": str,
    "input_or_output_uid": str,
    "input_or_output_reactome_id": int,
}

decomposed_uid_mapping = pd.DataFrame(columns=list(column_types.keys()))


def is_valid_uuid(value: Any) -> bool:
    """Check if the given value is a valid UUID."""
    try:
        uuid_obj = uuid.UUID(str(value), version=4)
        return str(uuid_obj) == value
    except ValueError:
        return False


def get_broken_apart_ids(
    broken_apart_members: list[set[str]], reactome_id: ReactomeID
) -> Set[UID]:
    """Get broken apart IDs."""
    global decomposed_uid_mapping

    uids: Set[UID]
    if any(isinstance(member, set) for member in broken_apart_members):
        new_broken_apart_members = []
        for member in broken_apart_members:
            if isinstance(member, set):
                new_broken_apart_members.append(member)
            else:
                new_broken_apart_members.append({member})

        iterproduct_components = list(itertools.product(*new_broken_apart_members))
        iterproduct_components_as_sets = [set(map(str, item)) for item in iterproduct_components]
        uids = get_uids_for_iterproduct_components(iterproduct_components_as_sets, reactome_id)
    else:
        uid = str(uuid.uuid4())
        rows: List[DataFrameRow] = []
        row: DataFrameRow
        for member in broken_apart_members:
            if is_valid_uuid(member):
                component_ids = decomposed_uid_mapping.loc[
                    decomposed_uid_mapping["uid"] == member, "component_id"
                ].tolist()
                for component_id in component_ids:
                    row = {
                        "uid": uid,
                        "component_id": component_id,
                        "reactome_id": reactome_id,
                        "input_or_output_uid": member,
                        "input_or_output_reactome_id": None,
                    }
                    rows.append(row)
            else:
                row = {
                    "uid": uid,
                    "component_id": member,
                    "reactome_id": reactome_id,
                    "input_or_output_uid": None,
                    "input_or_output_reactome_id": member,
                }
                rows.append(row)
        uids = {uid}
        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])

    return uids


def get_or_assign_uid(input_or_output_ids: Set[InputOutputID]) -> UID:
    """Get or assign UID."""
    global decomposed_uid_mapping

    uid_to_input_or_output: Dict[UID, Set[InputOutputID]] = {}
    for index, row in decomposed_uid_mapping.iterrows():
        uid = row["uid"]
        input_or_output_uid = row["input_or_output_uid"]
        input_or_output_reactome_id = row["input_or_output_reactome_id"]

        if uid in uid_to_input_or_output:
            if input_or_output_uid:
                uid_to_input_or_output[uid].add(input_or_output_uid)
            if input_or_output_reactome_id:
                uid_to_input_or_output[uid].add(input_or_output_reactome_id)
        else:
            if input_or_output_uid:
                uid_to_input_or_output[uid] = {input_or_output_uid}
                if input_or_output_reactome_id:
                    uid_to_input_or_output[uid].add(input_or_output_reactome_id)
            elif input_or_output_reactome_id:
                uid_to_input_or_output[uid] = {input_or_output_reactome_id}

    matching_uid: UID = ""
    for uid, input_or_output_set in uid_to_input_or_output.items():
        if input_or_output_set == input_or_output_ids:
            matching_uid = uid
            break

    return matching_uid if matching_uid else str(uuid.uuid4())


def get_uids_for_iterproduct_components(
    iterproduct_components: List[Set[ComponentID]], reactome_id: ReactomeID
) -> Set[UID]:
    """Get UID for iterproduct components."""
    global decomposed_uid_mapping

    uids: Set[UID] = set()
    for component in iterproduct_components:
        component_to_input_or_output: Dict[ComponentID, InputOutputID] = {}
        for item in component:
            if is_valid_uuid(item):
                selected_rows = decomposed_uid_mapping.loc[
                    decomposed_uid_mapping["uid"] == item
                ]
                for index, selected_row in selected_rows.iterrows():
                    component_id = selected_row["component_id"]
                    component_to_input_or_output[component_id] = item
            else:
                component_to_input_or_output[item] = item

        uid = get_or_assign_uid(set(component_to_input_or_output.values()))

        rows: List[DataFrameRow] = []
        for component_id, input_or_output_id in component_to_input_or_output.items():
            input_or_output_uid = input_or_output_id if is_valid_uuid(input_or_output_id) else None
            input_or_output_reactome_id = input_or_output_id if not is_valid_uuid(input_or_output_id) else None
            row: DataFrameRow = {
                "uid": uid,
                "component_id": component_id,
                "input_or_output_uid": input_or_output_uid,
                "input_or_output_reactome_id": input_or_output_reactome_id,
                "reactome_id": reactome_id,
            }
            rows.append(row)

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])
        uids.add(uid)

    return uids


def break_apart_entity(entity_id: int) -> Set[str]:
    """Break apart entity."""
    global decomposed_uid_mapping


    labels = get_labels(entity_id)
    if "EntitySet" in labels:
        member_ids = get_set_members(entity_id)

        member_list: List[str] = []
        for member_id in member_ids:
            members = break_apart_entity(member_id)

            if isinstance(members, set):
                member_list.extend(members)
            else:
                member_list.extend(set(members))

        return set(member_list)

    elif "Complex" in labels:
        broken_apart_members: List[Set[str]] = []
        member_ids = get_complex_components(entity_id)

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            broken_apart_members.append(members)

        return get_broken_apart_ids(broken_apart_members, str(entity_id))

    elif any(
        entity_label in labels
        for entity_label in [
            "ChemicalDrug",
            "Drug",
            "EntityWithAccessionedSequence",
            "GenomeEncodedEntity",
            "OtherEntity",
            "Polymer",
            "SimpleEntity",
        ]
    ):

        return {str(entity_id)}

    else:
        logger.error(f"Not handling labels correctly for: {entity_id}")
        exit(1)


def decompose_by_reactions(reaction_ids: List[int]) -> List[Any]:
    """Decompose by reactions."""
    global decomposed_uid_mapping

    logger.debug("Decomposing reactions")

    all_best_matches = []

    for reaction_id in reaction_ids:
        input_ids = get_reaction_input_output_ids(reaction_id, "input")
        broken_apart_input_id = [break_apart_entity(input_id) for input_id in input_ids]
        input_combinations = get_broken_apart_ids(
            broken_apart_input_id, str(reaction_id)
        )

        output_ids = get_reaction_input_output_ids(reaction_id, "output")
        broken_apart_output_id = [
            break_apart_entity(output_id) for output_id in output_ids
        ]
        output_combinations = get_broken_apart_ids(
            broken_apart_output_id, str(reaction_id)
        )

        [best_matches, _] = find_best_reaction_match(
            input_combinations, output_combinations, decomposed_uid_mapping
        )

        all_best_matches += best_matches

    return all_best_matches


def get_decomposed_uid_mapping(
    pathway_id: str, reaction_connections: pd.DataFrame
) -> Tuple[pd.DataFrame, List[Any]]:
    """Get decomposed UID mapping."""
    global decomposed_uid_mapping

    decomposed_uid_mapping.drop(decomposed_uid_mapping.index, inplace=True)

    reaction_ids = pd.unique(
        reaction_connections[["parent_reaction_id", "child_reaction_id"]].values.ravel(
            "K"
        )
    )

    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]  # removing NA value from list
    reaction_ids = reaction_ids.astype(int).tolist()  # converting to integer
    best_matches = decompose_by_reactions(list(reaction_ids))

    return (decomposed_uid_mapping, best_matches)
