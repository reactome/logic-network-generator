import uuid
from typing import Dict

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore

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


def create_pathway_logic_network(
    decomposed_uid_mapping: DataFrame, reaction_connections: DataFrame, best_matches
) -> DataFrame:
    """
    Create pathway_logic_network DataFrame based on decomposed_uid_mapping, reaction_connections, and best_matches.

    Args:
        decomposed_uid_mapping (DataFrame): DataFrame containing decomposed UID mapping.
        reaction_connections (DataFrame): DataFrame containing reaction connections.

    Returns:
        DataFrame: The created pathway_logic_network DataFrame.
    """
    logger.debug("Adding reaction pairs to pathway_logic_network")

    columns: Dict[str, pd.Series] = {
        "source_id": pd.Series(dtype="Int64"),
        "target_id": pd.Series(dtype="Int64"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
        "edge_type": pd.Series(dtype="str"),
    }
    pathway_logic_network_data = []
    reaction_ids = pd.unique(
        reaction_connections[["preceding_reaction_id", "following_reaction_id"]]
        .stack()  # Stack the columns to convert them into a single series
        .dropna()  # Drop NaN values
    )
    print("reaction_connections")
    print(reaction_connections)
    # Count the number of reactions without preceding events
    reactions_without_preceding_events = reaction_connections[
        ~reaction_connections["following_reaction_id"].isin(
            reaction_connections["preceding_reaction_id"]
        )
    ]
    num_reactions_without_preceding_events = len(reactions_without_preceding_events)
    num_total_reactions = len(reaction_connections)
    print(
        "Number of reactions without preceding events:",
        len(reactions_without_preceding_events),
    )
    print("Total number of reactions:", num_total_reactions)
    if num_total_reactions > 0:
        percentage_reactions_without_preceding_events = (
            num_reactions_without_preceding_events / num_total_reactions
        ) * 100
        print(
            f"Percentage of reactions without preceding events: {percentage_reactions_without_preceding_events:.2f}%"
        )
    else:
        print("Total number of reactions is zero, cannot calculate percentage.")

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
    print("reaction_id_map")
    print(reaction_id_map)
    uid_reaction_connections = create_uid_reaction_connections(
        reaction_id_map, best_matches, decomposed_uid_mapping
    )
    print("uid_RC")
    print(uid_reaction_connections)

    reaction_uids = pd.unique(
        uid_reaction_connections[["preceding_uid", "following_uid"]]
        .stack()  # Stack the columns to convert them into a single series
        .dropna()  # Drop NaN values
    )

    positive_regulator_count = len(positive_regulator_map)
    negative_regulator_count = len(negative_regulator_map)

    # Count the number of catalysts
    num_catalysts = len(catalyst_map)
    print("Number of catalysts:", num_catalysts)
    print("Number of positive regulators:", positive_regulator_count)
    print("Number of negative regulators:", negative_regulator_count)

    # Create a dictionary to store relationships between Reactome IDs and UUIDs
    reactome_id_to_uuid = {}

    for reaction_uid in reaction_uids:
        input_hash = reaction_id_map.loc[
            reaction_id_map["uid"] == reaction_uid, "input_hash"
        ].iloc[0]
        filtered_rows_input = decomposed_uid_mapping[
            decomposed_uid_mapping["uid"] == input_hash
        ]
        input_or_output_uid_values_input = filtered_rows_input[
            "input_or_output_uid"
        ].tolist()
        input_or_output_uid_values_input = [
            value for value in input_or_output_uid_values_input if pd.notna(value)
        ]
        input_or_output_reactome_id_values_input = filtered_rows_input[
            "input_or_output_reactome_id"
        ].tolist()
        input_or_output_reactome_id_values_input = [
            value
            for value in input_or_output_reactome_id_values_input
            if pd.notna(value)
        ]

        preceding_uids = uid_reaction_connections[
            uid_reaction_connections["following_uid"] == reaction_uid
        ]["preceding_uid"].tolist()
        for preceding_uid in preceding_uids:
            output_hash = reaction_id_map.loc[
                reaction_id_map["uid"] == preceding_uid, "output_hash"
            ].iloc[0]
            filtered_rows_output = decomposed_uid_mapping[
                decomposed_uid_mapping["uid"] == output_hash
            ]
            input_or_output_uid_values_output = filtered_rows_output[
                "input_or_output_uid"
            ].tolist()
            input_or_output_uid_values_output = [
                value for value in input_or_output_uid_values_output if pd.notna(value)
            ]
            input_or_output_reactome_id_values_output = filtered_rows_output[
                "input_or_output_reactome_id"
            ].tolist()
            input_or_output_reactome_id_values_output = [
                value
                for value in input_or_output_reactome_id_values_output
                if pd.notna(value)
            ]

            # Assign UUIDs to extracted inputs and outputs
            input_uuids = [
                reactome_id_to_uuid.setdefault(input_reactome_id, str(uuid.uuid4()))
                for input_reactome_id in input_or_output_reactome_id_values_input
            ]
            output_uuids = [
                reactome_id_to_uuid.setdefault(output_reactome_id, str(uuid.uuid4()))
                for output_reactome_id in input_or_output_reactome_id_values_output
            ]

            # Determine "and_or" value based on input or output
            if input_or_output_uid_values_input:
                and_or = "and"
                edge_type = "input"  # Set edge_type for inputs
            else:
                and_or = "or"
                edge_type = "output"  # Set edge_type for outputs

            # Append the input and output UUIDs to pathway_logic_network_data
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

    for _, row_output in filtered_rows_output.iterrows():
        pathway_logic_network_data.append(
            {
                "source_id": (
                    input_or_output_uid_values_input[0]
                    if input_or_output_uid_values_input
                    else None
                ),
                "target_id": row_output["input_or_output_uid"],
                "pos_neg": "pos",
                "and_or": and_or,
                "edge_type": edge_type,
            }
        )
    for map_df, pos_neg, edge_type in zip(
        [catalyst_map, negative_regulator_map, positive_regulator_map],
        ["pos", "neg", "pos"],
        ["catalyst", "regulator", "regulator"],
    ):
        for _, row in map_df.iterrows():
            pathway_logic_network_data.append(
                {
                    "source_id": row["uuid"],
                    "target_id": row["reaction_uuid"],
                    "pos_neg": pos_neg,
                    "and_or": and_or,
                    "edge_type": edge_type,
                }
            )

    print("Appended data:", pathway_logic_network_data[-1])
    pathway_logic_network = pd.DataFrame(
        pathway_logic_network_data, columns=columns.keys()
    )
    print("pathway_logic_network")
    print(pathway_logic_network)

    root_inputs = find_root_inputs(pathway_logic_network)
    print("root_inputs")
    print(root_inputs)
    terminal_outputs = find_terminal_outputs(pathway_logic_network)
    print("terminal_outputs")
    print(terminal_outputs)

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
