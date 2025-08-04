import hashlib
import itertools
import uuid
import warnings
from typing import Any, Dict, List, Set, Tuple, Optional

import pandas as pd

from src.argument_parser import logger
from src.best_reaction_match import find_best_reaction_match
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types
from src.neo4j_connector import (
    contains_reference_gene_product_molecule_or_isoform,
    get_complex_components_with_stoichiometry,  # MODIFIED IMPORT
    get_labels,
    get_reaction_input_output_ids_with_stoichiometry,  # MODIFIED IMPORT
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

# Enhanced functional equivalence groups configuration
FUNCTIONAL_EQUIVALENT_GROUPS = {
    # Ubiquitin - classic case of functionally equivalent isoforms
    "68524": {
        "group_name": "ubiquitin_isoforms",
        "strategy": "use_canonical",
        "canonical_representative": "68524",
        "description": "Ubiquitin isoforms - all serve identical protein degradation signaling function",
        "biological_rationale": "Multiple ubiquitin genes produce identical/near-identical proteins"
    },
}

# MODIFIED: Add stoichiometry columns to the schema
enhanced_decomposed_uid_mapping_column_types = {
    **decomposed_uid_mapping_column_types,
    "stoichiometry": pd.Int64Dtype(),
    "entity_stoichiometry": pd.Int64Dtype(),
    "component_stoichiometry": pd.Int64Dtype(),
}

decomposed_uid_mapping = pd.DataFrame(
    columns=enhanced_decomposed_uid_mapping_column_types.keys()
).astype(enhanced_decomposed_uid_mapping_column_types)

reference_entity_dict: Dict[str, str] = {}


def get_component_id_or_reference_entity_id(reactome_id):
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


def is_functionally_equivalent_entity(entity_id: str) -> Optional[Dict[str, Any]]:
    """Check if entity is part of a functional equivalence group."""
    return FUNCTIONAL_EQUIVALENT_GROUPS.get(entity_id)


def apply_functional_equivalence_strategy(
    entity_id: str, 
    functional_group: Dict[str, Any],
    original_result: Set[str]
) -> Set[str]:
    """Apply functional equivalence strategy to entity decomposition."""
    strategy = functional_group.get("strategy", "use_canonical")
    
    if strategy == "use_canonical":
        canonical_rep = functional_group.get("canonical_representative", entity_id)
        return {canonical_rep}
    else:
        return original_result


def get_broken_apart_ids(
    broken_apart_members: list[set[str]], 
    reactome_id: ReactomeID,
    stoichiometry_map: Dict[str, int] = None  # NEW PARAMETER
) -> Set[UID]:
    """Get broken apart IDs with stoichiometry support."""
    global decomposed_uid_mapping

    if stoichiometry_map is None:
        stoichiometry_map = {}

    uids: Set[UID]
    if any(isinstance(member, set) for member in broken_apart_members):
        new_broken_apart_members = []
        for member in broken_apart_members:
            if isinstance(member, set):
                new_broken_apart_members.append(member)
            else:
                new_broken_apart_members.append({member})

        iterproduct_components = list(itertools.product(*new_broken_apart_members))
        total_combinations = len(iterproduct_components)
        
        if total_combinations > 1000:
            logger.warning(f"Large combination set for reaction {reactome_id}: {total_combinations} combinations")

        iterproduct_components_as_sets = [
            set(map(str, item)) for item in iterproduct_components
        ]
        uids = get_uids_for_iterproduct_components(
            iterproduct_components_as_sets, reactome_id, stoichiometry_map
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
                    # NEW: Add stoichiometry information
                    stoichiometry = stoichiometry_map.get(component_id, 1)
                    
                    row = {
                        "uid": uid,
                        "component_id": component_id,
                        "reactome_id": reactome_id,
                        "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(
                            component_id
                        ),
                        "input_or_output_uid": member,
                        "input_or_output_reactome_id": None,
                        "stoichiometry": stoichiometry,  # NEW
                        "entity_stoichiometry": stoichiometry,  # NEW
                        "component_stoichiometry": 1,  # NEW
                    }
                    rows.append(row)
            else:
                # NEW: Add stoichiometry information
                stoichiometry = stoichiometry_map.get(member, 1)
                
                row = {
                    "uid": uid,
                    "component_id": member,
                    "reactome_id": reactome_id,
                    "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(
                        member
                    ),
                    "input_or_output_uid": None,
                    "input_or_output_reactome_id": member,
                    "stoichiometry": stoichiometry,  # NEW
                    "entity_stoichiometry": stoichiometry,  # NEW
                    "component_stoichiometry": 1,  # NEW
                }
                rows.append(row)
        uids = {uid}

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])

    return uids


def get_uids_for_iterproduct_components(
    iterproduct_components: List[Set[ComponentID]], 
    reactome_id: ReactomeID,
    stoichiometry_map: Dict[str, int] = None  # NEW PARAMETER
) -> Set[UID]:
    """Get UID for iterproduct components with stoichiometry support."""
    global decomposed_uid_mapping

    if stoichiometry_map is None:
        stoichiometry_map = {}

    uids: Set[UID] = set()
    for component in iterproduct_components:
        component_to_input_or_output: Dict[ComponentID, InputOutputID] = {}
        component_stoichiometries: Dict[ComponentID, int] = {}  # NEW
        
        for item in component:
            if is_valid_uuid(item):
                selected_rows = decomposed_uid_mapping.loc[
                    decomposed_uid_mapping["uid"] == item
                ]
                for index, selected_row in selected_rows.iterrows():
                    component_id = selected_row["component_id"]
                    component_to_input_or_output[component_id] = item
                    # NEW: Get stoichiometry from existing data
                    component_stoichiometries[component_id] = selected_row.get("stoichiometry", 1)
            else:
                component_to_input_or_output[item] = item
                # NEW: Get stoichiometry from the map
                component_stoichiometries[item] = stoichiometry_map.get(item, 1)

        # MODIFIED: Include stoichiometry in UID calculation
        uid_data = [(comp, io, component_stoichiometries.get(comp, 1)) 
                   for comp, io in sorted(component_to_input_or_output.items())]
        uid = hashlib.sha256(str(uid_data).encode()).hexdigest()

        rows: List[DataFrameRow] = []
        for component_id, input_or_output_id in component_to_input_or_output.items():
            input_or_output_uid = (
                input_or_output_id if is_valid_uuid(input_or_output_id) else None
            )
            input_or_output_reactome_id = (
                input_or_output_id if not is_valid_uuid(input_or_output_id) else None
            )
            
            # NEW: Get stoichiometry information
            stoichiometry = component_stoichiometries.get(component_id, 1)
            
            row: DataFrameRow = {
                "uid": uid,
                "component_id": component_id,
                "reactome_id": reactome_id,
                "component_id_or_reference_entity_id": get_component_id_or_reference_entity_id(
                    component_id
                ),
                "input_or_output_uid": input_or_output_uid,
                "input_or_output_reactome_id": input_or_output_reactome_id,
                "stoichiometry": stoichiometry,  # NEW
                "entity_stoichiometry": stoichiometry,  # NEW
                "component_stoichiometry": stoichiometry,  # NEW
            }
            rows.append(row)

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])
        uids.add(uid)

    return uids


def break_apart_entity(entity_id: int) -> Set[str]:
    """Break apart entity with functional equivalence support."""
    global decomposed_uid_mapping

    entity_str = str(entity_id)

    # Check for functional equivalence rules first
    functional_group = is_functionally_equivalent_entity(entity_str)
    if functional_group:
        original_result = break_apart_entity_standard(entity_id)
        return apply_functional_equivalence_strategy(entity_str, functional_group, original_result)

    # Use original decomposition logic
    return break_apart_entity_standard(entity_id)


def break_apart_entity_standard(entity_id: int) -> Set[str]:
    """Standard break_apart_entity logic with stoichiometry support."""
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
        
        # MODIFIED: Get components with stoichiometry
        component_data = get_complex_components_with_stoichiometry(entity_id)
        print(f"Complex {entity_id} component data: {component_data}")
        
        component_stoichiometry = {}
        member_ids = []
        
        for component_info in component_data:
            component_id = component_info['component_id']
            stoichiometry = component_info['stoichiometry']
            component_stoichiometry[str(component_id)] = stoichiometry
            member_ids.append(component_id)

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            broken_apart_members.append(members)

        return get_broken_apart_ids(
            broken_apart_members, str(entity_id), component_stoichiometry
        )

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
    """Decompose by reactions with stoichiometry support."""
    global decomposed_uid_mapping

    logger.debug("Decomposing reactions with stoichiometry")

    all_best_matches = []
    for reaction_id in reaction_ids:
        print(f"Processing reaction {reaction_id}")
        
        # MODIFIED: Get inputs with stoichiometry
        input_data = get_reaction_input_output_ids_with_stoichiometry(reaction_id, "input")
        input_ids = [item['entity_id'] for item in input_data]
        input_stoichiometry = {str(item['entity_id']): item['stoichiometry'] for item in input_data}
        
        print(f"Reaction {reaction_id} inputs: {input_ids}")
        print(f"Input stoichiometry: {input_stoichiometry}")
        
        broken_apart_input_id = [break_apart_entity(input_id) for input_id in input_ids]
        input_combinations = get_broken_apart_ids(
            broken_apart_input_id, str(reaction_id), input_stoichiometry
        )

        # MODIFIED: Get outputs with stoichiometry
        output_data = get_reaction_input_output_ids_with_stoichiometry(reaction_id, "output")
        output_ids = [item['entity_id'] for item in output_data]
        output_stoichiometry = {str(item['entity_id']): item['stoichiometry'] for item in output_data}
        
        print(f"Reaction {reaction_id} outputs: {output_ids}")
        print(f"Output stoichiometry: {output_stoichiometry}")
        
        broken_apart_output_id = [
            break_apart_entity(output_id) for output_id in output_ids
        ]
        output_combinations = get_broken_apart_ids(
            broken_apart_output_id, str(reaction_id), output_stoichiometry
        )

        [best_matches, _] = find_best_reaction_match(
            input_combinations, output_combinations, decomposed_uid_mapping
        )

        all_best_matches += best_matches

    return all_best_matches


def get_decomposed_uid_mapping(
    pathway_id: str, reaction_connections: pd.DataFrame
) -> Tuple[pd.DataFrame, List[Any]]:
    """Get decomposed UID mapping with stoichiometry support."""
    global decomposed_uid_mapping

    # MODIFIED: Use enhanced column types
    decomposed_uid_mapping = pd.DataFrame(
        columns=enhanced_decomposed_uid_mapping_column_types.keys()
    ).astype(enhanced_decomposed_uid_mapping_column_types)

    reaction_ids = pd.unique(
        reaction_connections[
            ["preceding_reaction_id", "following_reaction_id"]
        ].values.ravel("K")
    )

    reaction_ids = reaction_ids[~pd.isna(reaction_ids)]  # removing NA value from list
    reaction_ids = reaction_ids.astype(int).tolist()  # converting to integer
    
    print(f"Processing {len(reaction_ids)} reactions for pathway {pathway_id}")
    
    best_matches = decompose_by_reactions(list(reaction_ids))
    
    print(f"Final decomposed_uid_mapping shape: {decomposed_uid_mapping.shape}")
    print(f"Stoichiometry column summary:")
    if 'stoichiometry' in decomposed_uid_mapping.columns:
        print(decomposed_uid_mapping['stoichiometry'].value_counts())

    return (decomposed_uid_mapping, best_matches)
