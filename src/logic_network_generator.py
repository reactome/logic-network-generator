import uuid
from typing import Dict, List, Any

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore
from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


def get_reactome_id_from_hash(
    decomposed_uid_mapping: pd.DataFrame, hash_value: str
) -> int:
    """Extract reactome_id for a given hash from decomposed_uid_mapping."""
    return decomposed_uid_mapping.loc[
        decomposed_uid_mapping["uid"] == hash_value, "reactome_id"
    ].values[0]


def create_reaction_id_map(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_ids: List[int],
    best_matches: pd.DataFrame,
) -> pd.DataFrame:
    """Create a mapping between reaction UIDs, reactome IDs, and input/output hashes."""

    reaction_id_map_column_types = {
        "uid": str,
        "reactome_id": pd.Int64Dtype(),
        "input_hash": str,
        "output_hash": str,
    }

    print("Checking best_matches contents:")

    rows = []
    for _, match in best_matches.iterrows():
        incoming_hash = match["incoming"]
        outgoing_hash = match["outgoing"]
        reactome_id = get_reactome_id_from_hash(decomposed_uid_mapping, incoming_hash)

        row = {
            "uid": str(uuid.uuid4()),
            "reactome_id": int(reactome_id),
            "input_hash": incoming_hash,
            "output_hash": outgoing_hash,
        }
        print("row")
        print(row)
        rows.append(row)

    reaction_id_map = pd.DataFrame(rows).astype(reaction_id_map_column_types)

    return reaction_id_map


def create_uid_reaction_connections(
    reaction_id_map: pd.DataFrame,
    best_matches: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """Create connections between reaction UIDs based on best matches."""

    reactome_id_to_uid_mapping = dict(
        zip(reaction_id_map["reactome_id"], reaction_id_map["uid"])
    )

    uid_reaction_connections_data = []

    for _, match in best_matches.iterrows():
        incoming_hash = match["incoming"]
        outgoing_hash = match["outgoing"]

        # Get reactome IDs for both hashes
        preceding_reaction_id = get_reactome_id_from_hash(
            decomposed_uid_mapping, incoming_hash
        )
        following_reaction_id = get_reactome_id_from_hash(
            decomposed_uid_mapping, outgoing_hash
        )

        # Get corresponding UIDs
        preceding_uid = reactome_id_to_uid_mapping.get(preceding_reaction_id)
        following_uid = reactome_id_to_uid_mapping.get(following_reaction_id)

        # Only add connection if both UIDs exist
        if preceding_uid is not None and following_uid is not None:
            uid_reaction_connections_data.append(
                {"preceding_uid": preceding_uid, "following_uid": following_uid}
            )

    # Fix: Always create DataFrame with proper columns, even if empty
    if uid_reaction_connections_data:
        return pd.DataFrame(uid_reaction_connections_data)
    else:
        # Return empty DataFrame with correct column structure
        return pd.DataFrame(columns=["preceding_uid", "following_uid"])


def _execute_regulator_query(
    graph: Graph, query: str, reaction_uuid: str, function_name: str
) -> List[Dict[str, Any]]:
    """Execute a regulator query and return processed results."""
    try:
        result = graph.run(query)
        regulators = []

        for record in result:
            regulator_uuid = str(uuid.uuid4())
            regulators.append(
                {
                    "reaction": reaction_uuid,
                    "PhysicalEntity": regulator_uuid,
                    "edge_type": "regulator",
                    "uuid": regulator_uuid,
                    "reaction_uuid": reaction_uuid,
                }
            )

        return regulators

    except Exception as e:
        logger.error(f"Error in {function_name}", exc_info=True)
        raise e


def get_catalysts_for_reaction(reaction_id_map: DataFrame, graph: Graph) -> DataFrame:
    """Get catalysts for reactions with stoichiometry using Neo4j graph queries."""
    catalyst_list = []

    for _, row in reaction_id_map.iterrows():
        reaction_id = row["reactome_id"]
        reaction_uuid = row["uid"]

        query = (
            f"MATCH (reaction:ReactionLikeEvent{{dbId: {reaction_id}}})-[r:catalystActivity]->(catalystActivity:CatalystActivity)-[:physicalEntity]->(catalyst:PhysicalEntity) "
            f"RETURN reaction.dbId AS reaction_id, catalyst.dbId AS catalyst_id, "
            f"COALESCE(r.stoichiometry, 1) AS stoichiometry, 'catalyst' AS edge_type"
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
        columns=[
            "reaction_id",
            "catalyst_id",
            "stoichiometry",
            "edge_type",
            "uuid",
            "reaction_uuid",
        ],
    )


def get_positive_regulators_for_reaction(
    reaction_id_mapping: DataFrame, graph: Graph
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
    reaction_id_mapping: DataFrame, graph: Graph
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


def get_non_null_values(df: pd.DataFrame, column: str) -> List[Any]:
    """Extract non-null values from a DataFrame column."""
    return [value for value in df[column].tolist() if pd.notna(value)]


def get_hash_for_reaction(
    reaction_id_map: pd.DataFrame, uid: str, hash_type: str
) -> str:
    """Get input_hash or output_hash for a given reaction UID."""
    return reaction_id_map.loc[reaction_id_map["uid"] == uid, hash_type].iloc[0]


def extract_uid_and_reactome_values(
    decomposed_uid_mapping: pd.DataFrame, hash_value: str
) -> tuple:
    """Extract UID and Reactome ID values for a given hash."""
    filtered_rows = decomposed_uid_mapping[decomposed_uid_mapping["uid"] == hash_value]

    uid_values = get_non_null_values(filtered_rows, "input_or_output_uid")
    reactome_id_values = get_non_null_values(
        filtered_rows, "input_or_output_reactome_id"
    )

    return uid_values, reactome_id_values


def assign_uuids(
    reactome_ids: List[str], reactome_id_to_uuid: Dict[str, str]
) -> List[str]:
    """Assign UUIDs to Reactome IDs, creating new ones if they don't exist."""
    return [
        reactome_id_to_uuid.setdefault(reactome_id, str(uuid.uuid4()))
        for reactome_id in reactome_ids
    ]


def determine_edge_properties(input_uid_values: List[Any]) -> tuple:
    """Determine and_or and edge_type based on input UID values."""
    if input_uid_values:
        return "and", "input"
    else:
        return "or", "output"


def add_pathway_connections(
    input_uuids: List[str],
    output_uuids: List[str],
    and_or: str,
    edge_type: str,
    pathway_logic_network_data: List[Dict[str, Any]],
) -> None:
    """Add all input-output connections to the pathway network data."""
    for input_uuid in input_uuids:
        for output_uuid in output_uuids:
            pathway_logic_network_data.append(
                {
                    "source_id": input_uuid,
                    "target_id": output_uuid,
                    "pos_neg": "pos",
                    "and_or": and_or,
                    "edge_type": edge_type,
                }
            )


def extract_inputs_and_outputs(
    reaction_uid: str,
    reaction_uids: List[str],
    uid_reaction_connections: pd.DataFrame,
    reaction_id_map: pd.DataFrame,
    decomposed_uid_mapping: pd.DataFrame,
    reactome_id_to_uuid: Dict[str, str],
    pathway_logic_network_data: List[Dict[str, Any]],
) -> None:
    """Extract inputs and outputs with stoichiometry support."""

    for reaction_uid in reaction_uids:
        # Extract input information WITH stoichiometry
        input_hash = get_hash_for_reaction(reaction_id_map, reaction_uid, "input_hash")
        input_uid_values, input_reactome_id_values = extract_uid_and_reactome_values(
            decomposed_uid_mapping, input_hash
        )

        # NEW: Extract stoichiometry for inputs
        input_stoichiometry = extract_stoichiometry_for_entity(
            decomposed_uid_mapping, input_hash, input_reactome_id_values
        )

        # Process preceding reactions (outputs)
        preceding_uids = uid_reaction_connections[
            uid_reaction_connections["following_uid"] == reaction_uid
        ]["preceding_uid"].tolist()

        for preceding_uid in preceding_uids:
            # Extract output information WITH stoichiometry
            output_hash = get_hash_for_reaction(
                reaction_id_map, preceding_uid, "output_hash"
            )
            output_uid_values, output_reactome_id_values = (
                extract_uid_and_reactome_values(decomposed_uid_mapping, output_hash)
            )

            # NEW: Extract stoichiometry for outputs
            output_stoichiometry = extract_stoichiometry_for_entity(
                decomposed_uid_mapping, output_hash, output_reactome_id_values
            )

            # Assign UUIDs
            input_uuids = assign_uuids(input_reactome_id_values, reactome_id_to_uuid)
            output_uuids = assign_uuids(output_reactome_id_values, reactome_id_to_uuid)

            # Determine edge properties
            and_or, edge_type = determine_edge_properties(input_uid_values)

            # Add connections with stoichiometry
            add_pathway_connections_with_stoichiometry(
                input_uuids,
                output_uuids,
                input_stoichiometry,
                output_stoichiometry,
                and_or,
                edge_type,
                pathway_logic_network_data,
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
    """Append regulatory relationships with stoichiometry."""

    regulator_configs = [
        (catalyst_map, "pos", "catalyst"),
        (negative_regulator_map, "neg", "regulator"),
        (positive_regulator_map, "pos", "regulator"),
    ]

    for map_df, pos_neg, edge_type_override in regulator_configs:
        for _, row in map_df.iterrows():
            stoichiometry = row.get(
                "stoichiometry", 1
            )  # Get stoichiometry if available

            pathway_logic_network_data.append(
                {
                    "source_id": row["uuid"],
                    "target_id": row["reaction_uuid"],
                    "pos_neg": pos_neg,
                    "and_or": and_or,
                    "edge_type": edge_type_override,
                    "source_stoichiometry": stoichiometry,
                    "target_stoichiometry": 1,  # Reactions are typically 1:1 with regulators
                    "stoichiometry_ratio": 1.0,
                }
            )


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
        percentage_without_preceding = (
            num_reactions_without_preceding / num_total_reactions
        ) * 100
        print("Percentage of reactions without preceding events")
        print(percentage_without_preceding)


def _print_regulator_statistics(
    positive_regulator_map: pd.DataFrame,
    negative_regulator_map: pd.DataFrame,
    catalyst_map: pd.DataFrame,
) -> None:
    """Print statistics about regulators and catalysts."""
    print(
        f"Positive regulator count: {len(positive_regulator_map)}\n"
        f"Negative regulator count: {len(negative_regulator_map)}\n"
        f"Number of catalysts: {len(catalyst_map)}"
    )


def extract_stoichiometry_for_entity(
    decomposed_uid_mapping: pd.DataFrame,
    hash_value: str,
    entity_reactome_ids: List[str],
) -> Dict[str, int]:
    """Extract stoichiometry information for entities."""
    stoichiometry_map = {}

    for reactome_id in entity_reactome_ids:
        # Get stoichiometry from decomposed mapping
        filtered_rows = decomposed_uid_mapping[
            (decomposed_uid_mapping["uid"] == hash_value)
            & (decomposed_uid_mapping["input_or_output_reactome_id"] == reactome_id)
        ]

        if not filtered_rows.empty:
            stoichiometry = filtered_rows["stoichiometry"].iloc[0]
            stoichiometry_map[reactome_id] = (
                stoichiometry if pd.notna(stoichiometry) else 1
            )
        else:
            stoichiometry_map[reactome_id] = 1

    return stoichiometry_map


def add_pathway_connections_with_stoichiometry(
    input_uuids: List[str],
    output_uuids: List[str],
    input_stoichiometry: Dict[str, int],
    output_stoichiometry: Dict[str, int],
    and_or: str,
    edge_type: str,
    pathway_logic_network_data: List[Dict[str, Any]],
) -> None:
    """Add connections with stoichiometry information."""

    for input_uuid in input_uuids:
        for output_uuid in output_uuids:
            # Calculate stoichiometry ratio
            input_stoich = input_stoichiometry.get(input_uuid, 1)
            output_stoich = output_stoichiometry.get(output_uuid, 1)
            stoich_ratio = output_stoich / input_stoich if input_stoich > 0 else 1.0

            pathway_logic_network_data.append(
                {
                    "source_id": input_uuid,
                    "target_id": output_uuid,
                    "pos_neg": "pos",
                    "and_or": and_or,
                    "edge_type": edge_type,
                    "source_stoichiometry": input_stoich,
                    "target_stoichiometry": output_stoich,
                    "stoichiometry_ratio": stoich_ratio,
                }
            )


def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: Any,
) -> pd.DataFrame:
    """Create a pathway logic network with stoichiometry support."""

    # Enhanced column structure
    columns = {
        "source_id": pd.Series(dtype="Int64"),
        "target_id": pd.Series(dtype="Int64"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
        "edge_type": pd.Series(dtype="str"),
        "source_stoichiometry": pd.Series(dtype="Int64"),  # NEW
        "target_stoichiometry": pd.Series(dtype="Int64"),  # NEW
        "stoichiometry_ratio": pd.Series(dtype="float"),  # NEW
    }
    pathway_logic_network_data = []

    # Extract unique reaction IDs
    reaction_ids = pd.unique(
        reaction_connections[["preceding_reaction_id", "following_reaction_id"]]
        .stack()
        .dropna()
    )

    # Calculate and print statistics
    _calculate_reaction_statistics(reaction_connections)

    # Create mappings and connections
    reaction_id_map = create_reaction_id_map(
        decomposed_uid_mapping, reaction_ids, best_matches
    )
    catalyst_map = get_catalysts_for_reaction(reaction_id_map, graph)
    negative_regulator_map = get_negative_regulators_for_reaction(
        reaction_id_map, graph
    )
    positive_regulator_map = get_positive_regulators_for_reaction(
        reaction_id_map, graph
    )

    uid_reaction_connections = create_uid_reaction_connections(
        reaction_id_map, best_matches, decomposed_uid_mapping
    )

    reaction_uids = pd.unique(
        uid_reaction_connections[["preceding_uid", "following_uid"]].stack().dropna()
    )

    # Print regulator statistics
    _print_regulator_statistics(
        positive_regulator_map, negative_regulator_map, catalyst_map
    )

    # Process reactions and regulators
    reactome_id_to_uuid = {}

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
    pathway_logic_network = pd.DataFrame(
        pathway_logic_network_data, columns=columns.keys()
    )

    # Find root inputs and terminal outputs
    root_inputs = find_root_inputs(pathway_logic_network)
    terminal_outputs = find_terminal_outputs(pathway_logic_network)

    print(
        f"root_inputs: {root_inputs}\n"
        f"terminal_outputs: {terminal_outputs}\n"
        f"pathway_logic_network: {pathway_logic_network}"
    )

    return pathway_logic_network


def find_root_inputs(pathway_logic_network):
    root_inputs = pathway_logic_network[
        (pathway_logic_network["source_id"].notnull())
        & (~pathway_logic_network["source_id"].isin(pathway_logic_network["target_id"]))
    ]["source_id"].tolist()
    return root_inputs


def find_terminal_outputs(pathway_logic_network):
    terminal_outputs = pathway_logic_network[
        ~pathway_logic_network["target_id"].isin(
            pathway_logic_network["source_id"].unique()
        )
    ]["target_id"].tolist()
    return terminal_outputs
