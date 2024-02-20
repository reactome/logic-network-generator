import itertools
import uuid
from typing import Any, Dict, List, Set, Tuple, Union
import sqlite3
import pandas as pd
from py2neo import Graph

from src.argument_parser import logger
from src.best_reaction_match import find_best_reaction_match
from src.neo4j_connector import (get_complex_components, get_labels,
                                 get_reaction_input_output_ids,
                                 get_set_members)


conn = sqlite3.connect('your_database.db')

def convert_df_to_sqlite_table(df: pd.DataFrame, table_name: str, schema: str):
    cursor = conn.cursor()
    cursor.execute(f"CREATE TABLE IF NOT EXISTS {table_name} ({schema})")
    df.to_sql(table_name, conn, if_exists='replace', index=False)
    conn.commit()

# Define your schema
table_schema = "uid TEXT, component_id TEXT, input_or_output_uid TEXT, input_or_output_reactome_id INT, reactome_id INT"

# Print debug information
print("Checking if table exists:")
cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='decomposed_uid_mapping';")
if cursor.fetchone():
    print("Table 'decomposed_uid_mapping' exists.")
else:
    print("Table 'decomposed_uid_mapping' does not exist.")


uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


# Define types
UID = str
ComponentID = str
InputOutputID = str
ReactomeID = str
DataFrameRow = Dict[str, Any]

# Dataframe type
DecomposedUIDDataFrame = pd.DataFrame
column_types = {
    "uid": str,
    "component_id": str,
    "input_or_output_uid": str,
    "input_or_output_reactome_id": int,
    "reactome_id": int,
}

decomposed_uid_mapping: DecomposedUIDDataFrame = pd.DataFrame(
    columns=["uid", "component_id", "input_or_output_id", "reactome_id"]
)

# Create new columns for 'input_or_output_uid' and 'input_or_output_reactome_id'
decomposed_uid_mapping["input_or_output_uid"] = decomposed_uid_mapping[
    "input_or_output_id"
].apply(lambda x: str(x) if not x.isdigit() else None)
decomposed_uid_mapping["input_or_output_reactome_id"] = decomposed_uid_mapping[
    "input_or_output_id"
].apply(lambda x: int(x) if x.isdigit() else None)

# Drop the original 'input_or_output_id' column
decomposed_uid_mapping.drop(columns=["input_or_output_id"], inplace=True)


def is_valid_uuid(value: Any) -> bool:
    """Check if the given value is a valid UUID."""
    try:
        uuid_obj = uuid.UUID(str(value), version=4)
        return str(uuid_obj) == value
    except ValueError:
        return False


def get_broken_apart_ids(
    broken_apart_members: List[Union[Set[str], str]], reactome_id: ReactomeID
) -> UID:
    """Get broken apart IDs."""
    global decomposed_uid_mapping

    uid: UID
    if any(isinstance(member, set) for member in broken_apart_members):
        new_broken_apart_members = []
        for member in broken_apart_members:
            if isinstance(member, set):
                new_broken_apart_members.append(member)
            else:
                new_broken_apart_members.append({member})

        iterproduct_components = list(itertools.product(*new_broken_apart_members))

        uid = get_uid_for_iterproduct_components(iterproduct_components, reactome_id)
    else:
        uid = str(uuid.uuid4())
        rows: List[DataFrameRow] = []
        row: DataFrameRow
        for broken_apart_member in broken_apart_members:
            if is_valid_uuid(broken_apart_member):
                component_ids = decomposed_uid_mapping.loc[
                    decomposed_uid_mapping["uid"] == broken_apart_member, "component_id"
                ].tolist()
                for component_id in component_ids:
                    row = {
                        "uid": uid,
                        "component_id": component_id,
                        "reactome_id": reactome_id,
                        "input_or_output_uid": broken_apart_member,
                        "input_or_output_reactome_id": None,
                    }
                    rows.append(row)
            else:
                broken_apart_member = next(iter(broken_apart_member))
                row = {
                    "uid": uid,
                    "component_id": broken_apart_member,
                    "reactome_id": reactome_id,
                    "input_or_output_uid": None,
                    "input_or_output_reactome_id": (
                        int(broken_apart_member)
                        if broken_apart_member.isdigit()
                        else None
                    ),
                }
                rows.append(row)

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])

    return uid


def get_or_assign_uid(input_or_output_ids: Set[InputOutputID]) -> UID:
    """Get or assign UID."""
    global decomposed_uid_mapping

    uid_to_input_or_output: Dict[UID, Set[InputOutputID]] = {}
    for index, row in decomposed_uid_mapping.iterrows():
        uid = row["uid"]
        input_or_output_id = row["input_or_output_id"]

        if uid in uid_to_input_or_output:
            uid_to_input_or_output[uid].add(input_or_output_id)
        else:
            uid_to_input_or_output[uid] = {input_or_output_id}

    matching_uid: UID = ""
    for uid, input_or_output_set in uid_to_input_or_output.items():
        if input_or_output_set == input_or_output_ids:
            matching_uid = uid
            break

    return matching_uid if matching_uid else str(uuid.uuid4())


def get_uid_for_iterproduct_components(
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
            row: DataFrameRow = {
                "uid": uid,
                "component_id": component_id,
                "input_or_output_id": input_or_output_id,
                "reactome_id": reactome_id,
            }
            rows.append(row)

        decomposed_uid_mapping = pd.concat([decomposed_uid_mapping, pd.DataFrame(rows)])
        uids.add(uid)

    return uids


def break_apart_entity(entity_id: int) -> Union[str, Set[str]]:
    """Break apart entity."""
    global decomposed_uid_mapping
    labels = get_labels(entity_id)

    if "EntitySet" in labels:
        member_ids = get_set_members(entity_id)
        broken_apart_members: List[Union[Set[str], str]] = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            if isinstance(members, set):
                broken_apart_members.extend(members)
            else:
                broken_apart_members.append(members)
        x = broken_apart_members
        print("x")
        print(x)
        return set(broken_apart_members)

    elif "Complex" in labels:
        member_ids = get_complex_components(entity_id)
        broken_apart_members = []

        for member_id in member_ids:
            members = break_apart_entity(member_id)
            broken_apart_members.extend(members)

        return get_broken_apart_ids(list(set(broken_apart_members)), str(entity_id))

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

        return str(entity_id)

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

# Use the modified function to convert DataFrame to SQLite table
convert_df_to_sqlite_table(decomposed_uid_mapping, 'decomposed_uid_mapping', table_schema)

print("Contents of decomposed_uid_mapping SQLite table:")
cursor = conn.execute("SELECT * FROM decomposed_uid_mapping")
for row in cursor.fetchall():
    print(row)

# Print debug information
print("Checking if table exists:")
cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='decomposed_uid_mapping';")
if cursor.fetchone():
    print("Table 'decomposed_uid_mapping' exists.")
else:
    print("Table 'decomposed_uid_mapping' does not exist.")



def get_decomposed_uid_mapping(
    pathway_id: str, reaction_connections: pd.DataFrame
) -> Tuple[DecomposedUIDDataFrame, List[Any]]:
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
