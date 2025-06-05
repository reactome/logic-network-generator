import uuid

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore
from typing import List, Dict, Any
from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


def create_reaction_id_map(decomposed_uid_mapping, reaction_ids, best_matches):
    reaction_id_map_column_types = {
        "uid": str,
        "reactome_id": pd.Int64Dtype(),
        "input_hash": str,
        "output_hash": str,
    }
    reaction_id_map = pd.DataFrame(columns=reaction_id_map_column_types.keys()).astype(
        reaction_id_map_column_types
    )

    rows = []
    print("Checking best_matches contents:")

    for index, match in best_matches.iterrows():
        incomming_hash = match["incomming"]
        outgoing_hash = match["outgoing"]
        reactome_id = decomposed_uid_mapping.loc[
            decomposed_uid_mapping["uid"] == incomming_hash, "reactome_id"
        ].values[0]
        row = {
            "uid": str(uuid.uuid4()),
            "reactome_id": int(reactome_id),
            "input_hash": incomming_hash,
            "output_hash": outgoing_hash,
        }
        print("row")
        print(row)
        rows.append(row)

    new_rows_df = pd.DataFrame(rows)
    reaction_id_map = pd.concat([reaction_id_map, new_rows_df], ignore_index=True)

    return reaction_id_map


def create_uid_reaction_connections(
    reaction_id_map: pd.DataFrame, best_matches: pd.DataFrame, decomposed_uid_mapping
) -> pd.DataFrame:
    uid_reaction_connections_data = []
    reactome_id_to_uid_mapping = dict(
        zip(reaction_id_map["reactome_id"], reaction_id_map["uid"])
    )

    # Create uid_reaction_connections from best_matches
    for _, match in best_matches.iterrows():
        incomming_hash = match["incomming"]
        outgoing_hash = match["outgoing"]
        preceding_reaction_id = decomposed_uid_mapping.loc[
            decomposed_uid_mapping["uid"] == incomming_hash, "reactome_id"
        ].values[0]
        following_reaction_id = decomposed_uid_mapping.loc[
            decomposed_uid_mapping["uid"] == outgoing_hash, "reactome_id"
        ].values[0]
        preceding_uid = reactome_id_to_uid_mapping.get(preceding_reaction_id)
        following_uid = reactome_id_to_uid_mapping.get(following_reaction_id)
        if preceding_uid is not None and following_uid is not None:
            uid_reaction_connections_data.append(
                {"preceding_uid": preceding_uid, "following_uid": following_uid}
            )

    uid_reaction_connections = pd.DataFrame(uid_reaction_connections_data)
    return uid_reaction_connections


def get_catalysts_for_reaction(reaction_id_map: DataFrame, graph: Graph) -> DataFrame:
    catalyst_list = []

    for _, row in reaction_id_map.iterrows():
        reaction_id = row["reactome_id"]
        query = (
            f"MATCH (reaction:ReactionLikeEvent{{dbId: {reaction_id}}})-[:catalystActivity]->(catalystActivity:CatalystActivity)-[:physicalEntity]->(catalyst:PhysicalEntity) "
            f"RETURN reaction.dbId AS reaction_id, catalyst.dbId AS catalyst_id, 'catalyst' AS edge_type"
        )
        try:
            data = graph.run(query).data()
            for item in data:
                item["uuid"] = str(uuid.uuid4())  # Generate UUID for each entity
                # Map the reaction ID to the UUID
                item["reaction_uuid"] = row["uid"]
            catalyst_list.extend(data)
        except Exception as e:
            logger.error("Error in get_catalysts_for_reaction", exc_info=True)
            raise e

    return pd.DataFrame(
        catalyst_list,
        columns=["reaction_id", "catalyst_id", "edge_type", "uuid", "reaction_uuid"],
    )


def get_positive_regulators_for_reaction(
    reaction_id_mapping: DataFrame, graph: Graph
) -> DataFrame:
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
        try:
            result = graph.run(query)
            for record in result:
                regulator_uuid = str(uuid.uuid4())  # Generate UUID for each entity
                regulators_list.append(
                    {
                        "reaction": reaction_uuid,
                        "PhysicalEntity": regulator_uuid,
                        "edge_type": "regulator",
                        "uuid": regulator_uuid,
                        "reaction_uuid": reaction_uuid,
                    }
                )
        except Exception as e:
            logger.error("Error in get_positive_regulators_for_reaction", exc_info=True)
            raise e

    return pd.DataFrame(
        regulators_list,
        columns=["reaction", "PhysicalEntity", "edge_type", "uuid", "reaction_uuid"],
        index=None,
    )


def get_negative_regulators_for_reaction(
    reaction_id_mapping: DataFrame, graph: Graph
) -> DataFrame:
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
        try:
            result = graph.run(query)
            for record in result:
                regulator_uuid = str(uuid.uuid4())  # Generate UUID for each entity
                regulators_list.append(
                    {
                        "reaction": reaction_uuid,
                        "PhysicalEntity": regulator_uuid,
                        "edge_type": "regulator",
                        "uuid": regulator_uuid,
                        "reaction_uuid": reaction_uuid,
                    }
                )
        except Exception as e:
            logger.error("Error in get_negative_regulators_for_reaction", exc_info=True)
            raise e

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


def _determine_edge_properties(input_uid_values: List[Any]) -> tuple:
    """Determine and_or and edge_type based on input UID values."""
    if input_uid_values:
        return "and", "input"
    else:
        return "or", "output"


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
    """Extract inputs and outputs for reactions and add them to the pathway network."""
    
    for reaction_uid in reaction_uids:
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
            
            # Determine edge properties
            and_or, edge_type = _determine_edge_properties(input_uid_values)
            
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
        print("Percentage of reactions without preceding events")
        print(percentage_without_preceding)


def _print_regulator_statistics(
    positive_regulator_map: pd.DataFrame,
    negative_regulator_map: pd.DataFrame,
    catalyst_map: pd.DataFrame
) -> None:
    """Print statistics about regulators and catalysts."""
    print(
        f"Positive regulator count: {len(positive_regulator_map)}\n"
        f"Negative regulator count: {len(negative_regulator_map)}\n"
        f"Number of catalysts: {len(catalyst_map)}"
    )


def create_pathway_logic_network(
    decomposed_uid_mapping: pd.DataFrame,
    reaction_connections: pd.DataFrame,
    best_matches: Any,
) -> pd.DataFrame:
    """Create a pathway logic network from decomposed UID mappings and reaction connections."""
    logger.debug("Adding reaction pairs to pathway_logic_network")

    # Initialize data structures
    columns = {
        "source_id": pd.Series(dtype="Int64"),
        "target_id": pd.Series(dtype="Int64"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
        "edge_type": pd.Series(dtype="str"),
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
    
    # Note: Using empty strings to maintain original functionality
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
