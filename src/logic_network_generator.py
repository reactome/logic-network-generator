import uuid
from typing import Dict, List, Any, NamedTuple, Optional, Set

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore
from src.argument_parser import logger
from src.reaction_generator import _complex_contains_entity_set, _UBIQUITIN_ENTITY_SET_IDS

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


class PathwayResult(NamedTuple):
    """Result of pathway logic network generation.

    Attributes:
        logic_network: DataFrame containing the pathway logic network edges
        uuid_mapping: Dictionary mapping Reactome IDs to UUIDs
        catalyst_regulator_map: DataFrame containing catalyst and regulator information
        reaction_id_map: DataFrame mapping reaction UUIDs to Reactome reaction IDs
    """
    logic_network: pd.DataFrame
    uuid_mapping: Dict[str, str]
    catalyst_regulator_map: pd.DataFrame
    reaction_id_map: pd.DataFrame


def _get_reactome_id_from_hash(decomposed_uid_mapping: pd.DataFrame, hash_value: str) -> str:
    """Extract reactome_id (stable ID) for a given hash from decomposed_uid_mapping."""
    return decomposed_uid_mapping.loc[
        decomposed_uid_mapping["uid"] == hash_value, "reactome_id"
    ].values[0]


def create_reaction_id_map(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_ids: List[str],
    best_matches: pd.DataFrame
) -> pd.DataFrame:
    """Create a mapping between reaction UIDs, Reactome IDs, and input/output hashes.

    This function creates "virtual reactions" from best_matches, which pairs input
    and output combinations within biological reactions. Each best_match represents
    one possible transformation within a reaction.

    Why Virtual Reactions?
        A biological reaction in Reactome might have:
        - Multiple inputs (e.g., ATP, Water)
        - Multiple outputs (e.g., ADP, Phosphate)

        After decomposition (breaking down complexes and sets), we need to pair
        specific input combinations with specific output combinations. The Hungarian
        algorithm (used to create best_matches) optimally pairs these combinations.

        Each pairing becomes a "virtual reaction" with:
        - A unique UID (UUID v4)
        - The original Reactome reaction ID
        - An input_hash (identifying the input combination)
        - An output_hash (identifying the output combination)

    UID Strategy:
        - Each virtual reaction gets a NEW unique UID (UUID v4)
        - This UID is distinct from the original Reactome reaction ID
        - The UID is used to track transformations through the logic network
        - The Reactome ID preserves the link to the original biological reaction

    Example:
        Biological Reaction (Reactome ID: 12345):
            Inputs: Complex(A,B), ATP
            Outputs: Complex(A,B,P), ADP

        After decomposition and best matching:
            Virtual Reaction 1 (UID: uuid-1, Reactome ID: 12345):
                input_hash: "hash-of-A,B,ATP"
                output_hash: "hash-of-A,B,P,ADP"

        This virtual reaction can then be used to create entity→reaction→entity edges:
            A→VR1, B→VR1, ATP→VR1 (inputs), VR1→A, VR1→B, VR1→P, VR1→ADP (outputs)

    Args:
        decomposed_uid_mapping: Maps hashes to decomposed physical entities
        reaction_ids: List of Reactome reaction IDs (currently unused in function)
        best_matches: DataFrame with 'incomming' and 'outgoing' hash columns
                     Each row represents an optimal input/output pairing

    Returns:
        DataFrame with columns:
            - uid: Unique identifier for this virtual reaction (UUID v4 string)
            - reactome_id: Original Reactome reaction ID
            - input_hash: Hash identifying the input combination
            - output_hash: Hash identifying the output combination

    Note:
        The function assumes best_matches comes from Hungarian algorithm optimal
        pairing, ensuring each input combination maps to exactly one output combination.
    """
    
    reaction_id_map_column_types = {
        "uid": str,
        "reactome_id": str,
        "input_hash": str,
        "output_hash": str,
    }

    rows = []
    for _, match in best_matches.iterrows():
        incomming_hash = match["incomming"]
        outgoing_hash = match["outgoing"]
        reactome_id = _get_reactome_id_from_hash(decomposed_uid_mapping, incomming_hash)

        row = {
            "uid": str(uuid.uuid4()),
            "reactome_id": reactome_id,
            "input_hash": incomming_hash,
            "output_hash": outgoing_hash,
        }
        rows.append(row)

    reaction_id_map = pd.DataFrame(rows).astype(reaction_id_map_column_types)
    
    return reaction_id_map


def _execute_regulator_query(
    graph: Graph,
    query: str,
    reaction_uuid: str,
    function_name: str
) -> List[Dict[str, Any]]:
    """Execute a regulator query and return processed results."""
    try:
        result = graph.run(query)
        regulators = []

        for record in result:
            regulator_uuid = str(uuid.uuid4())
            regulators.append({
                "reaction": record.get("reaction"),
                "PhysicalEntity": record.get("PhysicalEntity"),  # Keep stId from query
                "edge_type": "regulator",
                "uuid": regulator_uuid,
                "reaction_uuid": reaction_uuid,
            })

        return regulators

    except Exception as e:
        logger.error(f"Error in {function_name}", exc_info=True)
        raise e


def get_catalysts_for_reaction(reaction_id_map: DataFrame, graph: Graph) -> DataFrame:
    """Get catalysts for reactions using Neo4j graph queries."""
    catalyst_list = []
    
    for _, row in reaction_id_map.iterrows():
        reaction_id = row["reactome_id"]
        reaction_uuid = row["uid"]
        
        query = (
            f"MATCH (reaction:ReactionLikeEvent{{stId: '{reaction_id}'}})-[:catalystActivity]->(catalystActivity:CatalystActivity)-[:physicalEntity]->(catalyst:PhysicalEntity) "
            f"RETURN reaction.stId AS reaction_id, catalyst.stId AS catalyst_id, 'catalyst' AS edge_type"
        )
        
        try:
            data = graph.run(query).data()
            for item in data:
                item["uuid"] = str(uuid.uuid4())
                item["reaction_uuid"] = reaction_uuid
            catalyst_list.extend(data)
            
        except Exception as e:
            logger.error("Error in get_catalysts_for_reaction", exc_info=True)
            raise e
    
    return pd.DataFrame(
        catalyst_list,
        columns=["reaction_id", "catalyst_id", "edge_type", "uuid", "reaction_uuid"],
    )


def get_positive_regulators_for_reaction(
    reaction_id_mapping: DataFrame, 
    graph: Graph
) -> DataFrame:
    """Get positive regulators for reactions using Neo4j graph queries."""
    regulators_list = []
    
    for _, row in reaction_id_mapping.iterrows():
        reaction_id = row["reactome_id"]
        reaction_uuid = row["uid"]
        
        if pd.isna(reaction_uuid):
            logger.error(f"No UUID found for reaction ID {reaction_id}")
            continue
        
        query = (
            f"MATCH (reaction)-[:regulatedBy]->(regulator:PositiveRegulation)-[:regulator]->(pe:PhysicalEntity) "
            f"WHERE reaction.stId = '{reaction_id}' "
            "RETURN reaction.stId as reaction, pe.stId as PhysicalEntity"
        )

        regulators = _execute_regulator_query(
            graph, query, reaction_uuid, "get_positive_regulators_for_reaction"
        )
        regulators_list.extend(regulators)
    
    return pd.DataFrame(
        regulators_list,
        columns=["reaction", "PhysicalEntity", "edge_type", "uuid", "reaction_uuid"],
        index=None,
    )


def get_negative_regulators_for_reaction(
    reaction_id_mapping: DataFrame, 
    graph: Graph
) -> DataFrame:
    """Get negative regulators for reactions using Neo4j graph queries."""
    regulators_list = []
    
    for _, row in reaction_id_mapping.iterrows():
        reaction_id = row["reactome_id"]
        reaction_uuid = row["uid"]
        
        if pd.isna(reaction_uuid):
            logger.error(f"No UUID found for reaction ID {reaction_id}")
            continue
        
        query = (
            f"MATCH (reaction)-[:regulatedBy]->(regulator:NegativeRegulation)-[:regulator]->(pe:PhysicalEntity) "
            f"WHERE reaction.stId = '{reaction_id}' "
            "RETURN reaction.stId as reaction, pe.stId as PhysicalEntity"
        )

        regulators = _execute_regulator_query(
            graph, query, reaction_uuid, "get_negative_regulators_for_reaction"
        )
        regulators_list.extend(regulators)
    
    return pd.DataFrame(
        regulators_list,
        columns=["reaction", "PhysicalEntity", "edge_type", "uuid", "reaction_uuid"],
        index=None,
    )


def _get_non_null_values(df: pd.DataFrame, column: str) -> List[Any]:
    """Extract non-null values from a DataFrame column."""
    return [value for value in df[column].tolist() if pd.notna(value)]


def _get_hash_for_reaction(reaction_id_map: pd.DataFrame, uid: str, hash_type: str) -> str:
    """Get input_hash or output_hash for a given reaction UID."""
    return reaction_id_map.loc[
        reaction_id_map["uid"] == uid, hash_type
    ].iloc[0]


def _extract_uid_and_reactome_values(decomposed_uid_mapping: pd.DataFrame, hash_value: str) -> tuple:
    """Extract UID and Reactome ID values for a given hash."""
    filtered_rows = decomposed_uid_mapping[decomposed_uid_mapping["uid"] == hash_value]

    uid_values = _get_non_null_values(filtered_rows, "input_or_output_uid")
    reactome_id_values = _get_non_null_values(filtered_rows, "input_or_output_reactome_id")

    return uid_values, reactome_id_values


def _build_uid_index(decomposed_uid_mapping: pd.DataFrame) -> Dict[str, tuple]:
    """Build a lookup index from decomposed_uid_mapping for fast UID resolution.

    Returns a dict mapping each uid to (list_of_nested_uids, list_of_terminal_reactome_ids, stoich_map).
    stoich_map maps reference IDs (nested UIDs or terminal Reactome IDs) to their stoichiometry.
    """
    index: Dict[str, tuple] = {}
    for uid_val, group in decomposed_uid_mapping.groupby("uid"):
        nested_uids = _get_non_null_values(group, "input_or_output_uid")
        terminal_ids = _get_non_null_values(group, "input_or_output_reactome_id")
        stoich_map: Dict[str, int] = {}
        for _, row in group.iterrows():
            stoich = row.get("stoichiometry")
            if pd.isna(stoich):
                stoich = 1
            else:
                stoich = int(stoich)
            if pd.notna(row.get("input_or_output_uid")):
                stoich_map[row["input_or_output_uid"]] = stoich
            if pd.notna(row.get("input_or_output_reactome_id")):
                stoich_map[row["input_or_output_reactome_id"]] = stoich
        index[uid_val] = (nested_uids, terminal_ids, stoich_map)
    return index


def _resolve_to_terminal_reactome_ids(
    uid_index: Dict[str, tuple],
    hash_value: str,
    visited: set = None
) -> Dict[str, int]:
    """Recursively resolve a hash to its terminal Reactome IDs with stoichiometry.

    With full EntitySet decomposition, the decomposed_uid_mapping contains nested UIDs:
    a hash may point to other UIDs (input_or_output_uid) rather than terminal Reactome IDs
    (input_or_output_reactome_id). This function follows the UID chain to find the actual
    terminal entity IDs, multiplying stoichiometry through each level.

    Args:
        uid_index: Pre-built lookup index from _build_uid_index
        hash_value: The hash/UID to resolve
        visited: Set of already-visited hashes (cycle detection)

    Returns:
        Dict mapping terminal Reactome ID → cumulative stoichiometry
    """
    if visited is None:
        visited = set()
    if hash_value in visited:
        return {}
    visited.add(hash_value)

    entry = uid_index.get(hash_value)
    if entry is None:
        return {}

    nested_uids, terminal_ids, stoich_map = entry
    result: Dict[str, int] = {}

    for tid in terminal_ids:
        stoich = stoich_map.get(tid, 1)
        result[tid] = result.get(tid, 0) + stoich

    for nested_uid in nested_uids:
        parent_stoich = stoich_map.get(nested_uid, 1)
        sub_results = _resolve_to_terminal_reactome_ids(uid_index, nested_uid, visited)
        for tid, sub_stoich in sub_results.items():
            result[tid] = result.get(tid, 0) + parent_stoich * sub_stoich

    return result


def _get_or_create_entity_uuid(
    entity_dbId: str,
    source_reaction_uuid: str,
    target_reaction_uuid: str,
    entity_uuid_registry: Dict[tuple, str]
) -> str:
    """
    Get or create UUID for entity based on its position in the pathway.

    Uses union-find logic to ensure entities in the same connected component
    get the same UUID, while entities at different pathway positions get different UUIDs.

    Args:
        entity_dbId: Reactome database ID of the entity
        source_reaction_uuid: UUID of reaction that outputs this entity
        target_reaction_uuid: UUID of reaction that receives this entity as input
        entity_uuid_registry: Registry mapping (entity_dbId, reaction_uuid, role) -> entity_uuid

    Returns:
        UUID for this entity at this position
    """
    # Create keys for this connection
    target_key = (entity_dbId, target_reaction_uuid, "input")
    source_key = (entity_dbId, source_reaction_uuid, "output")

    target_uuid = entity_uuid_registry.get(target_key)
    source_uuid = entity_uuid_registry.get(source_key)

    if target_uuid and source_uuid and target_uuid == source_uuid:
        # Already registered with same UUID (shouldn't happen but handle gracefully)
        logger.debug(f"Entity {entity_dbId} already has same UUID at both positions")
        return target_uuid
    elif target_uuid and source_uuid:
        # Entity has different UUIDs at source and target - merge them
        # Keep target_uuid, update all source_uuid references to target_uuid
        merge_count = 0
        for key, uuid_val in list(entity_uuid_registry.items()):
            if uuid_val == source_uuid:
                entity_uuid_registry[key] = target_uuid
                merge_count += 1
        logger.debug(
            f"Merged UUIDs for entity {entity_dbId}: "
            f"{source_uuid[:8]}... -> {target_uuid[:8]}... ({merge_count} position entries merged)"
        )
        return target_uuid
    elif target_uuid:
        # Entity already has UUID at target - share it with source
        entity_uuid_registry[source_key] = target_uuid
        logger.debug(f"Entity {entity_dbId} sharing UUID {target_uuid[:8]}... from target to source")
        return target_uuid
    elif source_uuid:
        # Entity already has UUID at source - share it with target
        entity_uuid_registry[target_key] = source_uuid
        logger.debug(f"Entity {entity_dbId} sharing UUID {source_uuid[:8]}... from source to target")
        return source_uuid
    else:
        # New position - create new UUID
        new_uuid = str(uuid.uuid4())
        entity_uuid_registry[target_key] = new_uuid
        entity_uuid_registry[source_key] = new_uuid
        logger.debug(f"Created new UUID {new_uuid[:8]}... for entity {entity_dbId}")
        return new_uuid


def _assign_uuids(
    reactome_ids: List[str],
    source_reaction_uuid: str,
    target_reaction_uuid: str,
    entity_uuid_registry: Dict[tuple, str]
) -> List[str]:
    """
    Assign position-aware UUIDs to entities based on their connections.

    Args:
        reactome_ids: List of entity Reactome database IDs
        source_reaction_uuid: UUID of reaction that outputs these entities
        target_reaction_uuid: UUID of reaction that receives these entities as inputs
        entity_uuid_registry: Registry for tracking entity UUIDs by position

    Returns:
        List of UUIDs for the entities
    """
    return [
        _get_or_create_entity_uuid(
            entity_dbId, source_reaction_uuid, target_reaction_uuid, entity_uuid_registry
        )
        for entity_dbId in reactome_ids
    ]


def _register_entity_uuid(
    entity_dbId: str,
    reaction_uuid: str,
    role: str,
    entity_uuid_registry: Dict[tuple, str],
    boundary_eids: Optional[Set[str]] = None,
    boundary_cache: Optional[Dict[str, str]] = None,
) -> str:
    """Register an entity with a single role key, creating a new UUID if needed.

    Unlike _get_or_create_entity_uuid which creates both input and output keys,
    this only creates the specified role key. Used in Phase 1 to avoid spurious
    cross-role entries.

    When boundary_eids and boundary_cache are provided, entities in boundary_eids
    share a single UUID across all their appearances (via the cache). This ensures
    root inputs and terminal outputs get one UUID per stId within their role.

    Args:
        entity_dbId: Reactome database ID of the entity
        reaction_uuid: UUID of the reaction
        role: "input" or "output"
        entity_uuid_registry: Registry mapping (entity_dbId, reaction_uuid, role) -> UUID
        boundary_eids: Optional set of entity IDs that are boundary entities
        boundary_cache: Optional cache mapping entity_dbId -> shared UUID for boundary entities

    Returns:
        UUID for this entity at this position
    """
    key = (entity_dbId, reaction_uuid, role)
    if key not in entity_uuid_registry:
        if boundary_eids and boundary_cache is not None and entity_dbId in boundary_eids:
            if entity_dbId not in boundary_cache:
                boundary_cache[entity_dbId] = str(uuid.uuid4())
            entity_uuid_registry[key] = boundary_cache[entity_dbId]
        else:
            entity_uuid_registry[key] = str(uuid.uuid4())
    return entity_uuid_registry[key]


def _build_entity_producer_count(vr_entities: Dict[str, tuple]) -> Dict[str, int]:
    """Count how many VRs produce each entity as output.

    Used to determine OR logic on output edges: entities produced by
    multiple VRs get and_or="or" (either source can provide it).
    """
    count: Dict[str, int] = {}
    for vr_uid, (input_ids, output_ids, *_) in vr_entities.items():
        for eid in output_ids:
            count[eid] = count.get(eid, 0) + 1
    return count


def _build_reactome_to_vr_map(reaction_id_map: pd.DataFrame) -> Dict[str, List[str]]:
    """Build mapping from original Reactome reaction stable ID to list of virtual reaction UIDs.

    A single Reactome reaction can produce multiple virtual reactions (one per
    input/output pairing from the Hungarian algorithm).

    Args:
        reaction_id_map: DataFrame with 'reactome_id' and 'uid' columns

    Returns:
        Dict mapping reactome_id (stId) -> list of VR UIDs
    """
    reactome_to_vr: Dict[str, List[str]] = {}
    for _, row in reaction_id_map.iterrows():
        reactome_id = row["reactome_id"]
        vr_uid = row["uid"]
        reactome_to_vr.setdefault(reactome_id, []).append(vr_uid)
    return reactome_to_vr


def _resolve_vr_entities(
    reaction_id_map: pd.DataFrame,
    uid_index: Dict[str, tuple]
) -> Dict[str, tuple]:
    """Resolve each virtual reaction's input/output hashes to terminal Reactome IDs.

    Caches the resolution so Phase 2 and Phase 3 don't re-resolve.

    Args:
        reaction_id_map: DataFrame with 'uid', 'input_hash', 'output_hash' columns
        uid_index: Pre-built lookup index from _build_uid_index

    Returns:
        Dict mapping vr_uid -> (input_reactome_ids, output_reactome_ids,
                                input_stoich_map, output_stoich_map)
        where stoich maps are Dict[str, int] mapping entity_id → stoichiometry
    """
    vr_entities: Dict[str, tuple] = {}
    for _, row in reaction_id_map.iterrows():
        vr_uid = row["uid"]
        input_stoich = _resolve_to_terminal_reactome_ids(uid_index, row["input_hash"])
        output_stoich = _resolve_to_terminal_reactome_ids(uid_index, row["output_hash"])
        input_ids = list(input_stoich.keys())
        output_ids = list(output_stoich.keys())
        vr_entities[vr_uid] = (input_ids, output_ids, input_stoich, output_stoich)
    return vr_entities


def _decompose_regulator_entity(entity_id: str) -> List[tuple]:
    """Decompose a catalyst/regulator entity to terminal members.

    Returns list of (terminal_id, logic_type, stoichiometry) tuples.
    Complex members -> "and" (all needed), stoichiometry multiplied through.
    EntitySet members -> "or" (any suffices), stoichiometry preserved from sub-components.
    Simple entities -> returned as-is with "and" and stoichiometry 1.
    """
    from src.neo4j_connector import get_labels, get_complex_components, get_set_members

    labels = get_labels(entity_id)

    if "Complex" in labels:
        # Only decompose complexes that contain EntitySets (consistent with break_apart_entity)
        if not _complex_contains_entity_set(entity_id):
            return [(entity_id, "and", 1)]
        components = get_complex_components(entity_id)  # Dict[str, int]
        result = []
        for member_id, stoich in components.items():
            sub_results = _decompose_regulator_entity(member_id)
            for mid, logic, sub_stoich in sub_results:
                result.append((mid, logic, stoich * sub_stoich))
        return result if result else [(entity_id, "and", 1)]

    elif "EntitySet" in labels or "DefinedSet" in labels or "CandidateSet" in labels:
        # Skip ubiquitin EntitySets (consistent with break_apart_entity)
        if entity_id in _UBIQUITIN_ENTITY_SET_IDS:
            return [(entity_id, "or", 1)]
        members = get_set_members(entity_id)
        result = []
        for member_id in members:
            sub_results = _decompose_regulator_entity(member_id)
            # EntitySet members are OR alternatives — override logic_type
            result.extend((mid, "or", sub_stoich) for mid, _, sub_stoich in sub_results)
        return result if result else [(entity_id, "or", 1)]

    else:
        return [(entity_id, "and", 1)]


def append_regulators(
    catalyst_map: pd.DataFrame,
    negative_regulator_map: pd.DataFrame,
    positive_regulator_map: pd.DataFrame,
    pathway_logic_network_data: List[Dict[str, Any]],
    reactome_id_to_uuid: Dict[str, str],
    entity_uuid_registry: Optional[Dict[tuple, str]] = None,
) -> None:
    """Append regulatory relationships to the pathway network.

    Decomposes Complex/EntitySet catalysts and regulators to their terminal
    members so that perturbation of individual subunits can be traced through
    the network.

    When entity_uuid_registry is provided, reuses existing UUIDs for entities
    that already appear in the pathway (e.g., a protein that is both an input
    and a catalyst). This prevents the same protein from appearing as two
    disconnected nodes.
    """
    # Build reverse lookup: stId → first existing UUID from the registry
    stid_to_existing_uuid: Dict[str, str] = {}
    if entity_uuid_registry:
        for (entity_dbId, _reaction_uuid, _role), entity_uuid in entity_uuid_registry.items():
            if entity_dbId not in stid_to_existing_uuid:
                stid_to_existing_uuid[entity_dbId] = entity_uuid

    regulator_configs = [
        (catalyst_map, "pos", "catalyst", "catalyst_id"),
        (negative_regulator_map, "neg", "regulator", "PhysicalEntity"),
        (positive_regulator_map, "pos", "regulator", "PhysicalEntity"),
    ]

    for map_df, pos_neg, edge_type, entity_col in regulator_configs:
        for _, row in map_df.iterrows():
            entity_id = row.get(entity_col)
            if pd.isna(entity_id):
                entity_id = row.get("uuid")
            entity_id = str(entity_id)

            terminal_members = _decompose_regulator_entity(entity_id)

            for member_id, member_logic, member_stoich in terminal_members:
                # Reuse existing UUID if this entity already appears in the pathway
                if member_id in stid_to_existing_uuid:
                    member_uuid = stid_to_existing_uuid[member_id]
                else:
                    member_uuid = str(uuid.uuid4())
                and_or = member_logic
                pathway_logic_network_data.append({
                    "source_id": member_uuid,
                    "target_id": row["reaction_uuid"],
                    "pos_neg": pos_neg,
                    "and_or": and_or,
                    "edge_type": edge_type,
                    "stoichiometry": member_stoich,
                })
                reactome_id_to_uuid[member_uuid] = member_id


def _calculate_reaction_statistics(reaction_connections: pd.DataFrame) -> None:
    """Calculate and print statistics about reactions without preceding events."""
    reactions_without_preceding_events = reaction_connections[
        ~reaction_connections["following_reaction_id"].isin(
            reaction_connections["preceding_reaction_id"]
        )
    ]
    
    num_reactions_without_preceding = len(reactions_without_preceding_events)
    num_total_reactions = len(reaction_connections)

    if num_total_reactions > 0:
        percentage_without_preceding = (num_reactions_without_preceding / num_total_reactions) * 100
        logger.info(
            f"Percentage of reactions without preceding events: {percentage_without_preceding:.1f}%"
        )


def _print_regulator_statistics(
    positive_regulator_map: pd.DataFrame,
    negative_regulator_map: pd.DataFrame,
    catalyst_map: pd.DataFrame
) -> None:
    """Log statistics about regulators and catalysts."""
    logger.info(
        f"Regulator statistics - "
        f"Positive: {len(positive_regulator_map)}, "
        f"Negative: {len(negative_regulator_map)}, "
        f"Catalysts: {len(catalyst_map)}"
    )


def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: Any,
) -> PathwayResult:
    """Create a pathway logic network from decomposed UID mappings and reaction connections.

    This function generates a logic network with position-aware UUIDs. Entities at different
    pathway positions get different UUIDs, while entities in the same connected component
    share UUIDs (via union-find algorithm). This minimizes self-loops while maintaining
    proper entity tracking.

    Args:
        decomposed_uid_mapping: DataFrame containing mappings from hashes to physical entities.
            Required columns: 'uid', 'reactome_id', 'input_or_output_reactome_id'
        reaction_connections: DataFrame containing connections between reactions.
            Required columns: 'preceding_reaction_id', 'following_reaction_id'
        best_matches: DataFrame containing pairings of input/output hashes.
            Required columns: 'incomming', 'outgoing'

    Returns:
        PathwayResult containing:
            - logic_network: DataFrame with edges between physical entities
            - uuid_mapping: Dict[str, str] mapping UUIDs to Reactome database IDs
            - catalyst_regulator_map: DataFrame with catalyst and regulator information
            - reaction_id_map: DataFrame mapping reaction UIDs to Reactome IDs

    Raises:
        ValueError: If input DataFrames are empty or missing required columns.

    Notes:
        - Uses entity_uuid_registry to track (entity_dbId, reaction_uuid, role) -> UUID mappings
        - Union-find algorithm merges UUIDs for entities in same connected component
        - See POSITION_AWARE_UUID_DESIGN.md for detailed design documentation
    """
    logger.debug("Adding reaction pairs to pathway_logic_network")

    # Validate inputs
    if decomposed_uid_mapping.empty:
        raise ValueError("decomposed_uid_mapping cannot be empty")

    required_mapping_cols = {'uid', 'reactome_id', 'input_or_output_reactome_id'}
    missing_cols = required_mapping_cols - set(decomposed_uid_mapping.columns)
    if missing_cols:
        raise ValueError(
            f"decomposed_uid_mapping is missing required columns: {missing_cols}. "
            f"Available columns: {list(decomposed_uid_mapping.columns)}"
        )

    if reaction_connections.empty:
        raise ValueError("reaction_connections cannot be empty")

    required_connection_cols = {'preceding_reaction_id', 'following_reaction_id'}
    missing_cols = required_connection_cols - set(reaction_connections.columns)
    if missing_cols:
        raise ValueError(
            f"reaction_connections is missing required columns: {missing_cols}. "
            f"Available columns: {list(reaction_connections.columns)}"
        )

    # best_matches can be a DataFrame or other iterable
    if isinstance(best_matches, pd.DataFrame):
        if best_matches.empty:
            raise ValueError("best_matches cannot be empty")

        required_match_cols = {'incomming', 'outgoing'}
        missing_cols = required_match_cols - set(best_matches.columns)
        if missing_cols:
            raise ValueError(
                f"best_matches is missing required columns: {missing_cols}. "
                f"Available columns: {list(best_matches.columns)}"
            )

    logger.info(
        f"Input validation passed: {len(decomposed_uid_mapping)} mappings, "
        f"{len(reaction_connections)} connections, "
        f"{len(best_matches)} matches"
    )

    # Initialize data structures
    columns = {
        "source_id": pd.Series(dtype="Int64"),
        "target_id": pd.Series(dtype="Int64"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
        "edge_type": pd.Series(dtype="str"),
        "stoichiometry": pd.Series(dtype="Int64"),
    }
    pathway_logic_network_data: List[Dict[str, Any]] = []
    
    # Extract unique reaction IDs
    reaction_ids = pd.unique(
        reaction_connections[["preceding_reaction_id", "following_reaction_id"]]
        .stack()
        .dropna()
    )
    
    # Calculate and print statistics
    _calculate_reaction_statistics(reaction_connections)
    
    # Create mappings and connections
    reaction_id_map = create_reaction_id_map(decomposed_uid_mapping, reaction_ids, best_matches)
    catalyst_map = get_catalysts_for_reaction(reaction_id_map, graph)
    negative_regulator_map = get_negative_regulators_for_reaction(reaction_id_map, graph)
    positive_regulator_map = get_positive_regulators_for_reaction(reaction_id_map, graph)

    # Print regulator statistics
    _print_regulator_statistics(positive_regulator_map, negative_regulator_map, catalyst_map)

    # 3-Phase entity UUID assignment for inter-reaction connectivity
    entity_uuid_registry: Dict[tuple, str] = {}
    reactome_id_to_uuid: Dict[str, str] = {}

    # Pre-build index for fast UID resolution (O(1) lookups instead of O(N) DataFrame scans)
    uid_index = _build_uid_index(decomposed_uid_mapping)
    logger.debug(f"Built UID index with {len(uid_index)} entries")

    # Resolve VR entities and build reactome-to-VR map
    vr_entities = _resolve_vr_entities(reaction_id_map, uid_index)
    reactome_to_vr = _build_reactome_to_vr_map(reaction_id_map)

    logger.debug(f"Processing {len(vr_entities)} virtual reactions in 3 phases")

    # Pre-compute boundary entity sets for UUID caching.
    # Root inputs (never produced as output) and terminal outputs (never consumed
    # as input) should share one UUID per stId within their role.
    all_input_eids: Set[str] = set()
    all_output_eids: Set[str] = set()
    for vr_uid, (input_ids, output_ids, *_) in vr_entities.items():
        all_input_eids.update(input_ids)
        all_output_eids.update(output_ids)
    root_input_eids = all_input_eids - all_output_eids
    terminal_output_eids = all_output_eids - all_input_eids
    root_input_uuid_cache: Dict[str, str] = {}
    terminal_output_uuid_cache: Dict[str, str] = {}

    logger.debug(
        f"Boundary entities: {len(root_input_eids)} root inputs, "
        f"{len(terminal_output_eids)} terminal outputs"
    )

    # Phase 1: Register entities with correct role keys
    # Each entity gets a unique UUID per (entity, reaction, role) triple.
    # No cross-role keys are created (unlike the old self-loop approach).
    # Boundary entities (root inputs / terminal outputs) share one UUID per stId.
    for vr_uid, (input_ids, output_ids, *_) in vr_entities.items():
        for eid in input_ids:
            _register_entity_uuid(eid, vr_uid, "input", entity_uuid_registry,
                                  root_input_eids, root_input_uuid_cache)
        for eid in output_ids:
            _register_entity_uuid(eid, vr_uid, "output", entity_uuid_registry,
                                  terminal_output_eids, terminal_output_uuid_cache)

    logger.debug(f"Phase 1 complete: {len(entity_uuid_registry)} registry entries")

    # Phase 2: Merge UUIDs based on reaction topology
    # For each (preceding, following) connection, find shared entities
    # (preceding VR's outputs ∩ following VR's inputs) and merge their UUIDs.
    merge_count = 0
    for _, conn in reaction_connections.iterrows():
        if pd.isna(conn["preceding_reaction_id"]) or pd.isna(conn["following_reaction_id"]):
            continue
        preceding_rid = conn["preceding_reaction_id"]
        following_rid = conn["following_reaction_id"]

        preceding_vr_uids = reactome_to_vr.get(preceding_rid, [])
        following_vr_uids = reactome_to_vr.get(following_rid, [])

        for p_vr in preceding_vr_uids:
            p_outputs = set(vr_entities.get(p_vr, ([], [], {}, {}))[1])
            for f_vr in following_vr_uids:
                f_inputs = set(vr_entities.get(f_vr, ([], [], {}, {}))[0])
                shared = p_outputs & f_inputs
                for eid in shared:
                    _get_or_create_entity_uuid(
                        eid, p_vr, f_vr, entity_uuid_registry
                    )
                    merge_count += 1

    logger.debug(f"Phase 2 complete: {merge_count} merges performed")

    # Phase 3: Create edges using merged UUIDs
    # Look up the now-merged UUIDs from the registry and create
    # input→VR + VR→output edges.
    # Output edges get "or" when the entity is produced by multiple VRs.
    entity_producer_count = _build_entity_producer_count(vr_entities)

    for vr_uid, (input_ids, output_ids, input_stoich, output_stoich) in vr_entities.items():
        if not input_ids or not output_ids:
            continue

        for eid in input_ids:
            input_uuid = entity_uuid_registry[(eid, vr_uid, "input")]
            pathway_logic_network_data.append({
                "source_id": input_uuid,
                "target_id": vr_uid,
                "pos_neg": "pos",
                "and_or": "and",
                "edge_type": "input",
                "stoichiometry": input_stoich.get(eid, 1),
            })

        for eid in output_ids:
            output_uuid = entity_uuid_registry[(eid, vr_uid, "output")]
            and_or = "or" if entity_producer_count.get(eid, 0) > 1 else ""
            pathway_logic_network_data.append({
                "source_id": vr_uid,
                "target_id": output_uuid,
                "pos_neg": "pos",
                "and_or": and_or,
                "edge_type": "output",
                "stoichiometry": output_stoich.get(eid, 1),
            })

    # Log UUID registry statistics
    unique_uuids = set(entity_uuid_registry.values())
    unique_entities = set(key[0] for key in entity_uuid_registry.keys())
    logger.info(
        f"Position-aware UUID registry: {len(entity_uuid_registry)} position entries, "
        f"{len(unique_uuids)} unique UUIDs, {len(unique_entities)} unique entities"
    )

    # Build UUID -> stId mapping for export from the entity_uuid_registry
    for (entity_dbId, reaction_uuid, role), entity_uuid in entity_uuid_registry.items():
        reactome_id_to_uuid[entity_uuid] = entity_dbId

    # Pre-fetch decomposition data for catalyst/regulator entities
    cat_reg_entity_ids: Set[str] = set()
    for _, row in catalyst_map.iterrows():
        if pd.notna(row.get("catalyst_id")):
            cat_reg_entity_ids.add(str(row["catalyst_id"]))
    for _, row in pd.concat([negative_regulator_map, positive_regulator_map]).iterrows():
        if pd.notna(row.get("PhysicalEntity")):
            cat_reg_entity_ids.add(str(row["PhysicalEntity"]))

    if cat_reg_entity_ids:
        from src.neo4j_connector import prefetch_entity_decomposition_data
        prefetch_entity_decomposition_data(list(cat_reg_entity_ids))

    append_regulators(
        catalyst_map,
        negative_regulator_map,
        positive_regulator_map,
        pathway_logic_network_data,
        reactome_id_to_uuid,
        entity_uuid_registry=entity_uuid_registry,
    )
    
    # Create final DataFrame
    pathway_logic_network = pd.DataFrame(pathway_logic_network_data, columns=columns.keys())
    
    # Find root inputs and terminal outputs
    root_inputs = find_root_inputs(pathway_logic_network)
    terminal_outputs = find_terminal_outputs(pathway_logic_network)

    logger.info(
        f"Generated network with {len(pathway_logic_network)} edges, "
        f"{len(root_inputs)} root inputs, {len(terminal_outputs)} terminal outputs"
    )

    # Combine catalyst and regulator maps for export
    catalyst_regulator_uuid_map = pd.concat([
        catalyst_map,
        negative_regulator_map,
        positive_regulator_map
    ], ignore_index=True)

    return PathwayResult(
        logic_network=pathway_logic_network,
        uuid_mapping=reactome_id_to_uuid,
        catalyst_regulator_map=catalyst_regulator_uuid_map,
        reaction_id_map=reaction_id_map
    )

def find_root_inputs(pathway_logic_network: pd.DataFrame) -> List[Any]:
    """Find root input physical entities that are only sources, never targets.

    Args:
        pathway_logic_network: DataFrame with source_id and target_id columns

    Returns:
        List of physical entity IDs that appear as sources but never as targets
    """
    root_inputs = pathway_logic_network[
        (pathway_logic_network["source_id"].notnull())
        & (~pathway_logic_network["source_id"].isin(pathway_logic_network["target_id"]))
    ]["source_id"].tolist()
    return root_inputs


def find_terminal_outputs(pathway_logic_network: pd.DataFrame) -> List[Any]:
    """Find terminal output physical entities that are only targets, never sources.

    Args:
        pathway_logic_network: DataFrame with source_id and target_id columns

    Returns:
        List of physical entity IDs that appear as targets but never as sources
    """
    terminal_outputs = pathway_logic_network[
        ~pathway_logic_network["target_id"].isin(
            pathway_logic_network["source_id"].unique()
        )
    ]["target_id"].tolist()
    return terminal_outputs


def export_uuid_to_reactome_mapping(
    pathway_logic_network: pd.DataFrame,
    reaction_id_map: pd.DataFrame,
    reactome_id_to_uuid: Dict[str, str],
    catalyst_regulator_map: pd.DataFrame,
    output_file: str
) -> None:
    """Export mapping from UUIDs in logic network to Reactome stable IDs.

    Creates a simple two-column mapping file for all UUIDs that appear in the logic network.

    Args:
        pathway_logic_network: DataFrame with the logic network edges
        reaction_id_map: DataFrame with reaction UIDs and their Reactome IDs
        reactome_id_to_uuid: Dictionary mapping Reactome IDs to entity UUIDs
        catalyst_regulator_map: DataFrame with catalyst/regulator information
        output_file: Path to save the mapping CSV file

    Output CSV columns:
        - uuid: The UUID used in the logic network
        - stable_id: The Reactome stable ID (e.g., R-HSA-12345)
    """
    # Get all UUIDs from the logic network
    all_uuids: set[str] = set()
    all_uuids.update(pathway_logic_network['source_id'].dropna().unique())
    all_uuids.update(pathway_logic_network['target_id'].dropna().unique())

    # Create reverse mapping: UUID -> reactome_id
    uuid_to_reactome = {}

    # 1. Add entity UUIDs
    # With position-aware UUIDs, we iterate the other direction
    # The passed dict might be stId->UUID or UUID->stId, check first entry
    if reactome_id_to_uuid:
        sample_key = next(iter(reactome_id_to_uuid.keys()))
        # If key looks like a UUID (contains dashes), it's already uuid->stId
        if '-' in str(sample_key):
            # Already UUID -> stId mapping
            for entity_uuid, reactome_id in reactome_id_to_uuid.items():
                if entity_uuid in all_uuids:
                    uuid_to_reactome[entity_uuid] = str(reactome_id)
        else:
            # Old format: stId -> UUID mapping (may miss some UUIDs with position-awareness)
            for reactome_id, entity_uuid in reactome_id_to_uuid.items():
                if entity_uuid in all_uuids:
                    uuid_to_reactome[entity_uuid] = str(reactome_id)

    # 2. Add reaction UUIDs (from reaction_id_map)
    for _, row in reaction_id_map.iterrows():
        reaction_uuid = row['uid']
        if reaction_uuid in all_uuids:
            uuid_to_reactome[reaction_uuid] = str(row['reactome_id'])

    # 3. Add catalyst and regulator UUIDs (from catalyst_regulator_map)
    for _, row in catalyst_regulator_map.iterrows():
        cat_reg_uuid = row['uuid']
        if cat_reg_uuid in all_uuids:
            # Get the entity stId (catalyst_id or regulator PhysicalEntity)
            if 'catalyst_id' in row and pd.notna(row['catalyst_id']):
                entity_id = str(row['catalyst_id'])
            elif 'PhysicalEntity' in row and pd.notna(row['PhysicalEntity']):
                entity_id = str(row['PhysicalEntity'])
            else:
                continue  # Skip if we can't find the entity ID

            uuid_to_reactome[cat_reg_uuid] = entity_id

    # Create DataFrame and save
    mapping_rows = [{'uuid': uuid, 'stable_id': stable_id}
                    for uuid, stable_id in uuid_to_reactome.items()]

    mapping_df = pd.DataFrame(mapping_rows, columns=['uuid', 'stable_id'])
    mapping_df = mapping_df.sort_values('uuid')  # Sort for easier lookup

    mapping_df.to_csv(output_file, index=False)
    logger.info(f"Exported UUID to Reactome stable ID mapping with {len(mapping_df)} entries")
