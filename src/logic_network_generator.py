import uuid
from typing import Dict, List, Any

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore
from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


def _get_reactome_id_from_hash(decomposed_uid_mapping: pd.DataFrame, hash_value: str) -> int:
    """Extract reactome_id for a given hash from decomposed_uid_mapping."""
    return decomposed_uid_mapping.loc[
        decomposed_uid_mapping["uid"] == hash_value, "reactome_id"
    ].values[0]


def create_reaction_id_map(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_ids: List[int],
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

        This virtual reaction can then be used to create transformation edges:
            A→A, A→B, A→P, A→ADP, B→A, B→B, B→P, B→ADP, ATP→A, ATP→B, ATP→P, ATP→ADP

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
        "reactome_id": pd.Int64Dtype(),
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
            "reactome_id": int(reactome_id),
            "input_hash": incomming_hash,
            "output_hash": outgoing_hash,
        }
        rows.append(row)
    
    reaction_id_map = pd.DataFrame(rows).astype(reaction_id_map_column_types)
    
    return reaction_id_map


def create_uid_reaction_connections(
    reaction_id_map: pd.DataFrame, 
    best_matches: pd.DataFrame, 
    decomposed_uid_mapping: pd.DataFrame
) -> pd.DataFrame:
    """Create connections between reaction UIDs based on best matches."""
    
    reactome_id_to_uid_mapping = dict(
        zip(reaction_id_map["reactome_id"], reaction_id_map["uid"])
    )
    
    uid_reaction_connections_data = []
    
    for _, match in best_matches.iterrows():
        incomming_hash = match["incomming"]
        outgoing_hash = match["outgoing"]
        
        # Get reactome IDs for both hashes
        preceding_reaction_id = _get_reactome_id_from_hash(decomposed_uid_mapping, incomming_hash)
        following_reaction_id = _get_reactome_id_from_hash(decomposed_uid_mapping, outgoing_hash)
        
        # Get corresponding UIDs
        preceding_uid = reactome_id_to_uid_mapping.get(preceding_reaction_id)
        following_uid = reactome_id_to_uid_mapping.get(following_reaction_id)
        
        # Only add connection if both UIDs exist
        if preceding_uid is not None and following_uid is not None:
            uid_reaction_connections_data.append({
                "preceding_uid": preceding_uid, 
                "following_uid": following_uid
            })

    uid_reaction_connections = pd.DataFrame(
        uid_reaction_connections_data, columns=["preceding_uid", "following_uid"]
    )
    return uid_reaction_connections


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
                "reaction": reaction_uuid,
                "PhysicalEntity": regulator_uuid,
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
            f"MATCH (reaction:ReactionLikeEvent{{dbId: {reaction_id}}})-[:catalystActivity]->(catalystActivity:CatalystActivity)-[:physicalEntity]->(catalyst:PhysicalEntity) "
            f"RETURN reaction.dbId AS reaction_id, catalyst.dbId AS catalyst_id, 'catalyst' AS edge_type"
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
            f"WHERE reaction.dbId = {reaction_id} "
            "RETURN reaction.dbId as reaction, pe.dbId as PhysicalEntity"
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
            f"WHERE reaction.dbId = {reaction_id} "
            "RETURN reaction.dbId as reaction, pe.dbId as PhysicalEntity"
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


def _assign_uuids(reactome_ids: List[str], reactome_id_to_uuid: Dict[str, str]) -> List[str]:
    """Assign UUIDs to Reactome IDs, creating new ones if they don't exist."""
    return [
        reactome_id_to_uuid.setdefault(reactome_id, str(uuid.uuid4()))
        for reactome_id in reactome_ids
    ]


def _determine_edge_properties(num_preceding_reactions: int) -> tuple:
    """Determine AND/OR logic and edge type based on preceding reaction count.

    This function implements the user requirement for logic network semantics:
    - All inputs to reactions are AND relationships (required)
    - Multiple sources producing the same entity create OR relationships (alternatives)

    Logic Rules:
        1. Multiple sources (num_preceding > 1) → OR relationship
           - Multiple reactions can produce the same physical entity
           - Entity can come from ANY of the preceding reactions (alternative paths)
           - edge_type: "output" (entity is output of multiple reactions)

        2. Single source (num_preceding == 1) → AND relationship
           - Entity comes from exactly one source
           - Entity is REQUIRED from that source
           - edge_type: "input" (entity is required input)

    Examples:
        Scenario 1: Single pathway
            R1: Glucose → Glucose-6-P
            num_preceding = 1 → ("and", "input")
            Meaning: Glucose-6-P must come from R1

        Scenario 2: Multiple pathways converge
            R1: PathwayA → ATP
            R2: PathwayB → ATP
            R3: ATP → Energy

            For R3's perspective:
            - ATP can come from R1 OR R2
            - num_preceding = 2 → ("or", "output")
            - Edges: R1→ATP (OR), R2→ATP (OR)

            Then ATP→R3 would be AND (ATP is required input to R3)

        Scenario 3: Complex formation
            R1: ProteinA + ProteinB → Complex(A,B)
            Both inputs are required (AND)
            num_preceding = 1 → ("and", "input")

    Args:
        num_preceding_reactions: Number of reactions feeding into the current reaction.
                                For a given reaction, this counts how many preceding
                                reactions produce outputs consumed by current reaction.

    Returns:
        Tuple[str, str]: (and_or, edge_type)
            - and_or: "and" (required) or "or" (alternative)
            - edge_type: "input" (single source) or "output" (multiple sources)

    Note:
        This function doesn't directly handle regulator/catalyst logic, which is
        managed separately in append_regulators().
    """
    if num_preceding_reactions > 1:
        return "or", "output"
    else:
        return "and", "input"


def _add_pathway_connections(
    input_uuids: List[str], 
    output_uuids: List[str], 
    and_or: str, 
    edge_type: str,
    pathway_logic_network_data: List[Dict[str, Any]]
) -> None:
    """Add all input-output connections to the pathway network data."""
    for input_uuid in input_uuids:
        for output_uuid in output_uuids:
            pathway_logic_network_data.append({
                "source_id": input_uuid,
                "target_id": output_uuid,
                "pos_neg": "pos",
                "and_or": and_or,
                "edge_type": edge_type,
            })


def extract_inputs_and_outputs(
    reaction_uid: str,
    reaction_uids: List[str],
    uid_reaction_connections: pd.DataFrame,
    reaction_id_map: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame,
    reactome_id_to_uuid: Dict[str, str],
    pathway_logic_network_data: List[Dict[str, Any]],
) -> None:
    """Extract inputs and outputs for reactions and create transformation edges.

    IMPORTANT: This function creates edges representing biochemical transformations
    WITHIN each reaction, not connections BETWEEN reactions. Edges connect input
    physical entities (reactants) to output physical entities (products) using a
    cartesian product: every input connects to every output.

    Edge Semantics:
        Edges represent transformations within reactions:
        - Reaction: ATP + Water → ADP + Phosphate
        - Creates 4 edges: ATP→ADP, ATP→Phosphate, Water→ADP, Water→Phosphate

        Reactions connect IMPLICITLY through shared physical entities:
        - Reaction 1: A → B (creates edge: A is source, B is target)
        - Reaction 2: B → C (creates edge: B is source, C is target)
        - Result: Pathway flow A → B → C (B connects the reactions)

    AND/OR Logic Assignment:
        The function assigns AND/OR relationships based on how many preceding
        reactions feed into the current reaction:

        - Multiple sources (len(preceding_uids) > 1) → OR relationship
          Example: R1→EntityX (OR), R2→EntityX (OR)
          Meaning: Entity X can come from either R1 OR R2

        - Single source (len(preceding_uids) == 1) → AND relationship
          Example: R1→EntityX (AND)
          Meaning: Entity X must come from R1 (required input)

    Args:
        reaction_uid: Current reaction being processed (not actually used - iterates over all)
        reaction_uids: List of all reaction UIDs to process
        uid_reaction_connections: DataFrame with 'preceding_uid' and 'following_uid' columns
        reaction_id_map: Maps reaction UIDs to input/output hashes
        decomposed_uid_mapping: Maps hashes to physical entity Reactome IDs
        reactome_id_to_uuid: Cache mapping Reactome IDs to UUIDs (modified in-place)
        pathway_logic_network_data: Output list of edge dictionaries (modified in-place)

    Side Effects:
        - Modifies reactome_id_to_uuid by adding new UUID assignments
        - Appends edge dictionaries to pathway_logic_network_data

    Example:
        For a reaction with 2 inputs (A, B) and 2 outputs (C, D):
        - Creates 4 edges: A→C, A→D, B→C, B→D
        - Each edge has: source_id, target_id, pos_neg, and_or, edge_type
    """

    logger.debug(f"Processing {len(reaction_uids)} reaction UIDs")

    for idx, reaction_uid in enumerate(reaction_uids):
        # Extract input information
        input_hash = _get_hash_for_reaction(reaction_id_map, reaction_uid, "input_hash")
        input_uid_values, input_reactome_id_values = _extract_uid_and_reactome_values(
            decomposed_uid_mapping, input_hash
        )

        # Process preceding reactions (outputs)
        preceding_uids = uid_reaction_connections[
            uid_reaction_connections["following_uid"] == reaction_uid
        ]["preceding_uid"].tolist()

        for preceding_uid in preceding_uids:
            # Extract output information
            output_hash = _get_hash_for_reaction(reaction_id_map, preceding_uid, "output_hash")
            output_uid_values, output_reactome_id_values = _extract_uid_and_reactome_values(
                decomposed_uid_mapping, output_hash
            )

            # Assign UUIDs
            input_uuids = _assign_uuids(input_reactome_id_values, reactome_id_to_uuid)
            output_uuids = _assign_uuids(output_reactome_id_values, reactome_id_to_uuid)

            # Determine edge properties based on number of preceding reactions
            # If multiple preceding reactions produce outputs for this reaction → OR
            # If single source → AND
            and_or, edge_type = _determine_edge_properties(len(preceding_uids))

            # Add connections to pathway network
            _add_pathway_connections(
                input_uuids, output_uuids, and_or, edge_type, pathway_logic_network_data
            )


def append_regulators(
    catalyst_map: pd.DataFrame,
    negative_regulator_map: pd.DataFrame,
    positive_regulator_map: pd.DataFrame,
    pathway_logic_network_data: List[Dict[str, Any]],
    reactome_id_to_uuid: Dict[str, str],
    and_or: str,
    edge_type: str,
) -> None:
    """Append regulatory relationships to the pathway network."""
    
    regulator_configs = [
        (catalyst_map, "pos", "catalyst"),
        (negative_regulator_map, "neg", "regulator"),
        (positive_regulator_map, "pos", "regulator"),
    ]
    
    for map_df, pos_neg, edge_type_override in regulator_configs:
        for _, row in map_df.iterrows():
            pathway_logic_network_data.append({
                "source_id": row["uuid"],
                "target_id": row["reaction_uuid"],
                "pos_neg": pos_neg,
                "and_or": and_or,
                "edge_type": edge_type_override,
            })


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
) -> pd.DataFrame:
    """Create a pathway logic network from decomposed UID mappings and reaction connections.

    Args:
        decomposed_uid_mapping: DataFrame containing mappings from hashes to physical entities.
            Required columns: 'uid', 'reactome_id', 'input_or_output_reactome_id'
        reaction_connections: DataFrame containing connections between reactions.
            Required columns: 'preceding_reaction_id', 'following_reaction_id'
        best_matches: DataFrame containing pairings of input/output hashes.
            Required columns: 'incomming', 'outgoing'

    Returns:
        DataFrame representing the logic network with edges between physical entities.

    Raises:
        ValueError: If input DataFrames are empty or missing required columns.
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
    
    uid_reaction_connections = create_uid_reaction_connections(
        reaction_id_map, best_matches, decomposed_uid_mapping
    )
    
    reaction_uids = pd.unique(
        uid_reaction_connections[["preceding_uid", "following_uid"]].stack().dropna()
    )
    
    # Print regulator statistics
    _print_regulator_statistics(positive_regulator_map, negative_regulator_map, catalyst_map)
    
    # Process reactions and regulators
    reactome_id_to_uuid: Dict[str, str] = {}
    
    for reaction_uid in reaction_uids:
        extract_inputs_and_outputs(
            reaction_uid,
            reaction_uids,
            uid_reaction_connections,
            reaction_id_map,
            decomposed_uid_mapping,
            reactome_id_to_uuid,
            pathway_logic_network_data,
        )
    
    and_or = ""
    edge_type = ""
    append_regulators(
        catalyst_map,
        negative_regulator_map,
        positive_regulator_map,
        pathway_logic_network_data,
        reactome_id_to_uuid,
        and_or,
        edge_type,
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

    return pathway_logic_network

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
