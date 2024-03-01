import uuid
from typing import Dict

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore

from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


def create_reaction_id_map(reactome_ids, decomposed_uid_mapping):
    reaction_id_map_column_types = {
        "uid": str,
        "reactome_id": pd.Int64Dtype(),
        "input_hash": str,
        "output_hash": str,
    }
    reaction_id_map = pd.DataFrame(columns=reaction_id_map_column_types.keys()).astype(
        reaction_id_map_column_types
    )

    # Dictionary to store assigned UUIDs for each Reactome ID within each reaction
    assigned_uids = {}

    rows = []
    for reactome_id in reactome_ids:
        associated_hashes = (
            decomposed_uid_mapping[
                decomposed_uid_mapping["reactome_id"] == reactome_id
            ]["uid"]
            .unique()
            .tolist()
        )

        # Check if the Reactome ID is already assigned a UUID within the same reaction
        if reactome_id in assigned_uids:
            uid = assigned_uids[reactome_id]
        else:
            uid = str(uuid.uuid4())
            assigned_uids[reactome_id] = uid

        row = {
            "uid": uid,
            "reactome_id": int(reactome_id),
            "input_hash": associated_hashes[0],
            "output_hash": associated_hashes[1],
        }
        rows.append(row)

    if rows:
        new_rows_df = pd.DataFrame(rows)
        reaction_id_map = pd.concat([reaction_id_map, new_rows_df], ignore_index=True)

    return reaction_id_map


def get_reactions_for_pathway(pathway_id: int, graph: Graph) -> pd.DataFrame:
    query = (
        f"MATCH (pathway:Pathway{{dbId: {pathway_id}}})-[:hasEvent]->(reaction:ReactionLikeEvent) "
        "RETURN reaction.dbId AS reaction_id"
    )
    try:
        result = graph.run(query)
        pathway_reactions_list = [
            {k: v for k, v in record.items()} for record in result
        ]
        return pd.DataFrame(pathway_reactions_list)
    except Exception as e:
        logger.error("Error in get_reactions_for_pathway", exc_info=True)
        raise e


# Usage example
pathway_id = 1257604
pathway_reactions_list = get_reactions_for_pathway(pathway_id, graph)
print("pathway_reactions_list")
print(pathway_reactions_list)


def get_catalysts_for_reaction(reaction_id_map: DataFrame, graph: Graph) -> DataFrame:
    catalyst_list = []

    for _, row in reaction_id_map.iterrows():
        reaction_id = row["reactome_id"]
        print("reaction_id")
        print(reaction_id)
        query = (
            f"MATCH (reaction:ReactionLikeEvent{{dbId: {reaction_id}}})-[:catalystActivity]->(catalystActivity:CatalystActivity)-[:physicalEntity]->(catalyst:PhysicalEntity) "
            f"RETURN reaction.dbId AS reaction_id, catalyst.dbId AS catalyst_id, 'catalyst' AS edge_type"
        )
        print("query")
        print(query)
        try:
            data = graph.run(query).data()
            print("data")
            print(data)
            catalyst_list.extend(data)
        except Exception as e:
            logger.error("Error in get_catalysts_for_reaction", exc_info=True)
            raise e

    return pd.DataFrame(
        catalyst_list, columns=["reaction_id", "catalyst_id", "edge_type"]
    )


def get_positive_regulators_for_reaction(
    reaction_id_mapping: DataFrame, graph: Graph
) -> DataFrame:
    regulators_list = []

    for _, row in reaction_id_mapping.iterrows():
        reaction_id = row["reactome_id"]
        query = (
            f"MATCH (reaction)-[:regulatedBy]->(regulator:PositiveRegulation)-[:regulator]->(pe:PhysicalEntity) "
            f"WHERE reaction.dbId = {reaction_id} "
            "RETURN reaction.dbId as reaction, pe.dbId as PhysicalEntity"
        )
        try:
            result = graph.run(query)
            regulators_list.extend(
                [{k: v for k, v in record.items()} for record in result]
            )
        except Exception as e:
            logger.error("Error in get_negative_regulators_for_reaction", exc_info=True)
            raise e

    for record in regulators_list:
        record["edge_type"] = "regulator"

    return pd.DataFrame(
        regulators_list, columns=["reaction", "PhysicalEntity"], index=None
    )


def get_negative_regulators_for_reaction(
    reaction_id_mapping: DataFrame, graph: Graph
) -> DataFrame:
    regulators_list = []

    for _, row in reaction_id_mapping.iterrows():
        reaction_id = row["reactome_id"]
        query = (
            f"MATCH (reaction)-[:regulatedBy]->(regulator:NegativeRegulation)-[:regulator]->(pe:PhysicalEntity) "
            f"WHERE reaction.dbId = {reaction_id} "
            "RETURN reaction.dbId as reaction, pe.dbId as PhysicalEntity"
        )
        try:
            result = graph.run(query)
            regulators_list.extend(
                [{k: v for k, v in record.items()} for record in result]
            )
        except Exception as e:
            logger.error("Error in get_negative_regulators_for_reaction", exc_info=True)
            raise e

    for record in regulators_list:
        record["edge_type"] = "regulator"

    return pd.DataFrame(
        regulators_list, columns=["reaction", "PhysicalEntity"], index=None
    )


def create_pathway_logic_network(
    decomposed_uid_mapping: DataFrame,
    reaction_connections: DataFrame,
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
        "source_id": pd.Series(dtype="str"),
        "target_id": pd.Series(dtype="str"),
        "pos_neg": pd.Series(dtype="str"),
        "and_or": pd.Series(dtype="str"),
        "edge_type": pd.Series(dtype="str"),
    }
    pathway_logic_network_data = []

    print("reaction_connections")
    print(reaction_connections)

    reaction_ids = pd.unique(
        reaction_connections[["preceding_reaction_id", "following_reaction_id"]]
        .stack()  # Stack the columns to convert them into a single series
        .dropna()  # Drop NaN values
    )

    assigned_uids = {}
    reaction_id_map = create_reaction_id_map(reaction_ids, decomposed_uid_mapping)
    catalyst_map = get_catalysts_for_reaction(reaction_id_map, graph)
    negative_regulator_map = get_negative_regulators_for_reaction(
        reaction_id_map, graph
    )
    positive_regulator_map = get_positive_regulators_for_reaction(
        reaction_id_map, graph
    )
    print("negative_regulator_map")
    print(negative_regulator_map)
    print("positive_regulator_map")
    print(positive_regulator_map)
    print("catalyst_map")
    print(catalyst_map)
    print("reaction_id_map")
    print(reaction_id_map)
    missing_reactions = pathway_reactions_list[
        ~pathway_reactions_list["reaction_id"].isin(reaction_id_map["reactome_id"])
    ]

    missing_reaction_ids = missing_reactions["reaction_id"].tolist()
    print("missing_reaction_ids")
    print(missing_reaction_ids)

    for reaction_id in reaction_ids:
        input_hash = reaction_id_map.loc[
            reaction_id_map["reactome_id"] == reaction_id, "input_hash"
        ].iloc[0]
        print("input_hash")
        print(input_hash)
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

        preceding_reaction_ids = reaction_connections[
            reaction_connections["following_reaction_id"] == reaction_id
        ]["preceding_reaction_id"].tolist()
        for preceding_reaction_id in preceding_reaction_ids:
            output_hash = reaction_id_map.loc[
                reaction_id_map["reactome_id"] == preceding_reaction_id, "output_hash"
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

            for reactome_id, uid in assigned_uids.items():
                reactome_id_tuple = tuple(
                    sorted(
                        [
                            input_or_output_reactome_id_values_input,
                            input_or_output_reactome_id_values_output,
                        ]
                    )
                )
                if reactome_id_tuple in assigned_uids:
                    input_uid = assigned_uids[reactome_id_tuple]
                    output_uid = assigned_uids[reactome_id_tuple]
                    break
                else:
                    input_uid = str(uuid.uuid4())
                    output_uid = str(uuid.uuid4())
                    assigned_uids[reactome_id_tuple] = input_uid
                    assigned_uids[reactome_id_tuple] = output_uid

            # Determine "and_or" value based on input or output
            if input_or_output_reactome_id_values_input:
                and_or = "and"
            else:
                and_or = "or"
            # Determine "edge_type" based on whether the entity is an input or output
            if input_or_output_reactome_id_values_input:
                edge_type = "input"
            elif input_or_output_reactome_id_values_output:
                edge_type = "output"

            pathway_logic_network_data.append(
                {
                    "source_id": (
                        input_or_output_reactome_id_values_input[0]
                        if input_or_output_reactome_id_values_input
                        else None
                    ),
                    "target_id": (
                        input_or_output_reactome_id_values_output[0]
                        if input_or_output_reactome_id_values_output
                        else None
                    ),
                    "pos_neg": "pos",  # replace with actual value
                    "and_or": and_or,
                    "edge_type": edge_type,
                }
            )

            catalyst_ids = catalyst_map[
                catalyst_map["reaction_id"] == preceding_reaction_id
            ]["catalyst_id"].tolist()
            for catalyst_id in catalyst_ids:
                pathway_logic_network_data.append(
                    {
                        "source_id": catalyst_id,
                        "target_id": preceding_reaction_id,
                        "pos_neg": "pos",  # replace with actual value
                        "and_or": and_or,
                        "edge_type": "catalyst",
                    }
                )

            for _, row in negative_regulator_map.iterrows():
                pathway_logic_network_data.append(
                    {
                        "source_id": row["PhysicalEntity"],
                        "target_id": row["reaction"],
                        "pos_neg": "neg",
                        "and_or": "and",
                        "edge_type": "regulator",
                    }
                )

            for _, row in positive_regulator_map.iterrows():
                pathway_logic_network_data.append(
                    {
                        "source_id": row["PhysicalEntity"],
                        "target_id": row["reaction"],
                        "pos_neg": "pos",
                        "and_or": "and",
                        "edge_type": "regulator",
                    }
                )
        # Add missing reactions to pathway_logic_network_data
    for reaction_id in missing_reaction_ids:
        pathway_logic_network_data.append(
            {
                "source_id": None,  # Set this to the appropriate value if available
                "target_id": reaction_id,
                "pos_neg": None,  # Set this to the appropriate value if needed
                "and_or": None,  # Set this to the appropriate value if needed
                "edge_type": None,  # Set this to the appropriate value if needed
            }
        )
    pathway_logic_network = pd.DataFrame(
        pathway_logic_network_data, columns=columns.keys()
    )
    print("pathway_logic_network")
    print(pathway_logic_network)
    return pathway_logic_network
