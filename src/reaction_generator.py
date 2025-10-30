import hashlib
import itertools
import uuid
import warnings
from typing import Any, Dict, List, Set, Tuple, Union

import pandas as pd

from src.argument_parser import logger
from src.best_reaction_match import find_best_reaction_match
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types
from src.neo4j_connector import (
    contains_reference_gene_product_molecule_or_isoform,
    get_complex_components,
    get_labels,
    get_reaction_input_output_ids,
    get_reference_entity_id,
    get_set_members,
)

warnings.filterwarnings(
    "ignore",
    message="The behavior of DataFrame concatenation with empty or all-NA entries is deprecated.",
    category=FutureWarning,
)

# Define types
UID = str
ComponentID = str
InputOutputID = str
ReactomeID = str
DataFrameRow = Dict[str, Any]

decomposed_uid_mapping = pd.DataFrame(
    columns=decomposed_uid_mapping_column_types.keys()
).astype(  # type: ignore
    decomposed_uid_mapping_column_types
)

reference_entity_dict: Dict[str, str] = {}


def get_component_id_or_reference_entity_id(reactome_id: int) -> Union[str, int]:
    """Get the reference entity ID for a Reactome ID, with caching.

    Args:
        reactome_id: Reactome database ID for the entity

    Returns:
        Reference entity ID (string) if it exists, otherwise the reactome_id (int)
    """
    global reference_entity_dict

    if reactome_id in reference_entity_dict:
        component_id = reference_entity_dict[reactome_id]
        return component_id

    reference_entity_id = get_reference_entity_id(reactome_id)

    reference_entity_dict[reactome_id] = (
        reference_entity_id if reference_entity_id else reactome_id
    )

    return reference_entity_dict[reactome_id]


def is_valid_uuid(identifier: Any) -> bool:
    """Check if the given value is a valid UUID."""
    return True if len(identifier) == 64 else False


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
        iterproduct_components_as_sets = [
            set(map(str, item)) for item in iterproduct_components
        ]
        uids = get_uids_for_iterproduct_components(
            iterproduct_components_as_sets, reactome_id
        )
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
                        "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(
                            component_id
                        ),
                        "input_or_output_uid": member,
                        "input_or_output_reactome_id": None,
                    }
                    rows.append(row)
            else:
                row = {
                    "uid": uid,
                    "component_id": member,
                    "reactome_id": reactome_id,
                    "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(
                        component_id
                    ),
                    "input_or_output_uid": None,
                    "input_or_output_reactome_id": member,
                }
                rows.append(row)
        uids = {uid}
        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])

    return uids


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

        uid = hashlib.sha256(
            str(sorted(component_to_input_or_output.values())).encode()
        ).hexdigest()

        rows: List[DataFrameRow] = []
        for component_id, input_or_output_id in component_to_input_or_output.items():
            input_or_output_uid = (
                input_or_output_id if is_valid_uuid(input_or_output_id) else None
            )
            input_or_output_reactome_id = (
                input_or_output_id if not is_valid_uuid(input_or_output_id) else None
            )
            row: DataFrameRow = {
                "uid": uid,
                "component_id": component_id,
                "reactome_id": reactome_id,
                "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(
                    component_id
                ),
                "input_or_output_uid": input_or_output_uid,
                "input_or_output_reactome_id": input_or_output_reactome_id,
            }
            rows.append(row)

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])
        uids.add(uid)

    return uids


def break_apart_entity(entity_id: int) -> Set[str]:
    """Break apart entity."""
    global decomposed_uid_mapping

    labels = get_labels(entity_id)
    if "EntitySet" in labels or "Complex" in labels:
        filtered_rows = decomposed_uid_mapping[
            decomposed_uid_mapping["reactome_id"] == entity_id
        ]
        if not filtered_rows.empty:
            input_output_uid_values = filtered_rows["input_or_output_uid"]
            input_output_reactome_id_values = filtered_rows[
                "input_or_output_reactome_id"
            ]
            return set(input_output_uid_values.dropna()) | set(
                input_output_reactome_id_values.dropna()
            )

    if "EntitySet" in labels:
        if entity_id == 68524:  # ubiquitin
            return set([str(entity_id)])

        contains_thing = contains_reference_gene_product_molecule_or_isoform(entity_id)
        if contains_thing:
            return set([str(entity_id)])
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
        reaction_connections[
            ["preceding_reaction_id", "following_reaction_id"]
        ].values.ravel("K")
    )

    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]  # removing NA value from list
    reaction_ids = reaction_ids.astype(int).tolist()  # converting to integer
    best_matches = decompose_by_reactions(list(reaction_ids))

    return (decomposed_uid_mapping, best_matches)
