import uuid
from typing import Dict, List, Any, NamedTuple, Optional, Set

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore
from src.argument_parser import logger
from src.neo4j_connector import get_graph
from src.reaction_generator import (
    _complex_contains_entity_set,
    _UBIQUITIN_ENTITY_SET_IDS,
    get_terminal_components,
)


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
    best_matches: pd.DataFrame,
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
        # reaction_id was attached when this best_match was emitted by
        # decompose_by_reactions, avoiding the ambiguity of reverse-deriving
        # it from the input hash (the same hash can appear under multiple
        # reactome_ids in decomposed_uid_mapping).
        reactome_id = match["reactome_id"]

        row = {
            "uid": str(uuid.uuid4()),
            "reactome_id": reactome_id,
            "input_hash": incomming_hash,
            "output_hash": outgoing_hash,
        }
        rows.append(row)

    reaction_id_map = pd.DataFrame(rows).astype(reaction_id_map_column_types)
    
    return reaction_id_map


def _bulk_fetch_reaction_links(
    graph: Graph,
    reaction_ids: List[str],
    cypher_template: str,
) -> List[Dict[str, Any]]:
    """Run a single Cypher query that joins many reactions to their linked
    entities (catalysts or regulators).

    cypher_template must use $reaction_ids and return one row per
    (reaction_id, entity_id) pair. Replaces the previous N+1-query loop
    that hit Neo4j once per reaction.
    """
    if not reaction_ids:
        return []
    try:
        return graph.run(cypher_template, reaction_ids=reaction_ids).data()
    except Exception:
        logger.error("Error in _bulk_fetch_reaction_links", exc_info=True)
        raise


_CATALYST_CYPHER = (
    "UNWIND $reaction_ids AS rid "
    "MATCH (reaction:ReactionLikeEvent {stId: rid})"
    "-[:catalystActivity]->(:CatalystActivity)"
    "-[:physicalEntity]->(catalyst:PhysicalEntity) "
    "RETURN reaction.stId AS reaction_id, catalyst.stId AS entity_id"
)

_POS_REG_CYPHER = (
    "UNWIND $reaction_ids AS rid "
    "MATCH (reaction:ReactionLikeEvent {stId: rid})"
    "-[:regulatedBy]->(:PositiveRegulation)"
    "-[:regulator]->(pe:PhysicalEntity) "
    "RETURN reaction.stId AS reaction_id, pe.stId AS entity_id"
)

_NEG_REG_CYPHER = (
    "UNWIND $reaction_ids AS rid "
    "MATCH (reaction:ReactionLikeEvent {stId: rid})"
    "-[:regulatedBy]->(:NegativeRegulation)"
    "-[:regulator]->(pe:PhysicalEntity) "
    "RETURN reaction.stId AS reaction_id, pe.stId AS entity_id"
)


# Both catalyst and regulator DataFrames share this schema so they can be
# concatenated and consumed by a single column-name reader downstream.
_CAT_REG_COLUMNS = ["reaction_id", "entity_id", "edge_type", "uuid", "reaction_uuid"]


def _bulk_fetch_reaction_entity_links(
    reaction_id_map: DataFrame,
    graph: Graph,
    cypher: str,
    edge_type: str,
) -> DataFrame:
    """Shared body for catalyst and regulator fetching."""
    rids_to_uuids: Dict[str, List[str]] = {}
    for _, row in reaction_id_map.iterrows():
        if pd.isna(row["uid"]):
            logger.error(f"No UUID found for reaction ID {row['reactome_id']}")
            continue
        rids_to_uuids.setdefault(row["reactome_id"], []).append(row["uid"])

    rows = _bulk_fetch_reaction_links(graph, list(rids_to_uuids.keys()), cypher)

    out_rows = []
    for record in rows:
        for reaction_uuid in rids_to_uuids.get(record["reaction_id"], []):
            out_rows.append({
                "reaction_id": record["reaction_id"],
                "entity_id": record["entity_id"],
                "edge_type": edge_type,
                "uuid": str(uuid.uuid4()),
                "reaction_uuid": reaction_uuid,
            })

    return pd.DataFrame(out_rows, columns=_CAT_REG_COLUMNS)


def get_catalysts_for_reaction(reaction_id_map: DataFrame, graph: Graph) -> DataFrame:
    """Fetch catalysts for all reactions in one bulk Cypher query."""
    return _bulk_fetch_reaction_entity_links(
        reaction_id_map, graph, _CATALYST_CYPHER, edge_type="catalyst"
    )


def get_positive_regulators_for_reaction(
    reaction_id_map: DataFrame,
    graph: Graph,
) -> DataFrame:
    """Fetch positive regulators for all reactions in one bulk Cypher query."""
    return _bulk_fetch_reaction_entity_links(
        reaction_id_map, graph, _POS_REG_CYPHER, edge_type="regulator"
    )


def get_negative_regulators_for_reaction(
    reaction_id_map: DataFrame,
    graph: Graph,
) -> DataFrame:
    """Fetch negative regulators for all reactions in one bulk Cypher query."""
    return _bulk_fetch_reaction_entity_links(
        reaction_id_map, graph, _NEG_REG_CYPHER, edge_type="regulator"
    )


def _get_non_null_values(df: pd.DataFrame, column: str) -> List[Any]:
    """Extract non-null values from a DataFrame column."""
    return [value for value in df[column].tolist() if pd.notna(value)]


def _get_hash_for_reaction(reaction_id_map: pd.DataFrame, uid: str, hash_type: str) -> str:
    """Get input_hash or output_hash for a given reaction UID."""
    return reaction_id_map.loc[
        reaction_id_map["uid"] == uid, hash_type
    ].iloc[0]


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
            stoich_raw = row.get("stoichiometry")
            stoich = 1 if stoich_raw is None or pd.isna(stoich_raw) else int(stoich_raw)
            if pd.notna(row.get("input_or_output_uid")):
                stoich_map[row["input_or_output_uid"]] = stoich
            if pd.notna(row.get("input_or_output_reactome_id")):
                stoich_map[row["input_or_output_reactome_id"]] = stoich
        index[str(uid_val)] = (nested_uids, terminal_ids, stoich_map)
    return index


def _resolve_to_terminal_reactome_ids(
    uid_index: Dict[str, tuple],
    hash_value: str,
    visited: Optional[Set[str]] = None
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


def _uf_find(uid: str, unions: Dict[str, str]) -> str:
    """Walk the union-find chain to the root, with path compression."""
    if uid not in unions:
        return uid
    root = uid
    while root in unions:
        root = unions[root]
    cur = uid
    while cur in unions and unions[cur] != root:
        nxt = unions[cur]
        unions[cur] = root
        cur = nxt
    return root


def _canonicalize_registry(
    entity_uuid_registry: Dict[tuple, str],
    uuid_unions: Dict[str, str],
) -> None:
    """Rewrite every value in the registry to its union-find root.

    Called once after Phase 2's merges; turns the deferred unions into
    canonical UUIDs that downstream code can read directly.
    """
    if not uuid_unions:
        return
    for key, u in entity_uuid_registry.items():
        entity_uuid_registry[key] = _uf_find(u, uuid_unions)


def _get_or_create_entity_uuid(
    entity_dbId: str,
    source_reaction_uuid: str,
    target_reaction_uuid: str,
    entity_uuid_registry: Dict[tuple, str],
    uuid_unions: Optional[Dict[str, str]] = None,
) -> str:
    """Get or create UUID for entity based on its position in the pathway.

    Uses union-find to merge entities at connected positions in the pathway.
    Merges record source→target unions in `uuid_unions` rather than scanning
    the registry on every call — without that, repeated merges over a large
    registry are O(N²). The caller (`create_pathway_logic_network`) does a
    single canonicalization pass after Phase 2 finishes.

    If `uuid_unions` is omitted (e.g. tests calling this in isolation), a
    fresh map is allocated for the call. That keeps the previous semantics
    for single-shot use without forcing every call site to thread the map.
    """
    if uuid_unions is None:
        uuid_unions = {}

    target_key = (entity_dbId, target_reaction_uuid, "input")
    source_key = (entity_dbId, source_reaction_uuid, "output")

    target_uuid = entity_uuid_registry.get(target_key)
    if target_uuid is not None:
        target_uuid = _uf_find(target_uuid, uuid_unions)
    source_uuid = entity_uuid_registry.get(source_key)
    if source_uuid is not None:
        source_uuid = _uf_find(source_uuid, uuid_unions)

    if target_uuid and source_uuid and target_uuid == source_uuid:
        return target_uuid
    elif target_uuid and source_uuid:
        # Different roots — record source → target union (no registry scan)
        uuid_unions[source_uuid] = target_uuid
        return target_uuid
    elif target_uuid:
        entity_uuid_registry[source_key] = target_uuid
        return target_uuid
    elif source_uuid:
        entity_uuid_registry[target_key] = source_uuid
        return source_uuid
    else:
        new_uuid = str(uuid.uuid4())
        entity_uuid_registry[target_key] = new_uuid
        entity_uuid_registry[source_key] = new_uuid
        return new_uuid


def _assign_uuids(
    reactome_ids: List[str],
    source_reaction_uuid: str,
    target_reaction_uuid: str,
    entity_uuid_registry: Dict[tuple, str],
    uuid_unions: Optional[Dict[str, str]] = None,
) -> List[str]:
    """Assign position-aware UUIDs to entities based on their connections."""
    if uuid_unions is None:
        uuid_unions = {}
    return [
        _get_or_create_entity_uuid(
            entity_dbId, source_reaction_uuid, target_reaction_uuid,
            entity_uuid_registry, uuid_unions,
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
    """Decompose a catalyst/regulator entity to (terminal_id, stoichiometry) pairs.

    The decomposition rules mirror break_apart_entity (matching layer):
    Complex with EntitySet → cartesian over members; simple Complex
    returned intact; EntitySet → flat alternatives; ubiquitin sets are
    treated as atomic to avoid combinatorial explosion. AND/OR semantics
    for the resulting edges are decided by append_regulators based on
    pos_neg, not by within-entity decomposition shape.
    """
    from src.neo4j_connector import get_labels, get_complex_components, get_set_members

    labels = get_labels(entity_id)

    if "Complex" in labels:
        if not _complex_contains_entity_set(entity_id):
            return [(entity_id, 1)]
        components = get_complex_components(entity_id)  # Dict[str, int]
        result = []
        for member_id, stoich in components.items():
            for mid, sub_stoich in _decompose_regulator_entity(member_id):
                result.append((mid, stoich * sub_stoich))
        return result if result else [(entity_id, 1)]

    if "EntitySet" in labels or "DefinedSet" in labels or "CandidateSet" in labels:
        if entity_id in _UBIQUITIN_ENTITY_SET_IDS:
            return [(entity_id, 1)]
        members = get_set_members(entity_id)
        result = []
        for member_id in members:
            result.extend(_decompose_regulator_entity(member_id))
        return result if result else [(entity_id, 1)]

    return [(entity_id, 1)]


def _emit_boundary_decomposition_edges(
    pathway_logic_network_data: List[Dict[str, Any]],
    root_input_eids: Set[str],
    terminal_output_eids: Set[str],
    root_input_uuid_cache: Dict[str, str],
    terminal_output_uuid_cache: Dict[str, str],
    reactome_id_to_uuid: Dict[str, str],
) -> None:
    """Append synthetic edges that expose leaves of root/terminal complexes.

    For each root-input complex C with components {A, B, ...}, emit
    ``A → C``, ``B → C``, ... edges of edge_type='assembly'. For each
    terminal-output complex, emit ``C → A``, ``C → B``, ... of
    edge_type='dissociation'. Each leaf shares a single UUID across all
    boundary contexts so that perturbing a leaf at the assembly side
    propagates through any downstream dissociation that reads the same
    species.

    Intermediate complexes (those produced by some reaction AND consumed
    by another in this pathway) are intentionally NOT expanded — they're
    real biological species flowing between reactions, and the AB dimer
    is a different molecule from free A and free B. See
    docs/DESIGN_DECISIONS.md, "Two layers of decomposition."

    Simple-leaf root/terminal entities (proteins, small molecules) are
    skipped: they're already perturbable as themselves.

    A leaf reuses any UUID the entity already has elsewhere in the network
    (regular VR inputs/outputs, regulators, catalysts) so that perturbing a
    protein in one role propagates through every other role. Without this,
    boundary leaves would be disconnected duplicate nodes for the same
    biological entity.
    """
    from src.neo4j_connector import get_labels

    # Build stId → existing UUID lookup from everything assigned so far
    # (entity registry from VR phases, plus regulator/catalyst UUIDs added
    # by append_regulators). reactome_id_to_uuid is keyed by UUID, so invert.
    stid_to_existing_uuid: Dict[str, str] = {}
    for existing_uuid, stid in reactome_id_to_uuid.items():
        if stid not in stid_to_existing_uuid:
            stid_to_existing_uuid[stid] = existing_uuid

    leaf_uuid_registry: Dict[str, str] = {}

    def _leaf_uuid(leaf_stid: str) -> str:
        if leaf_stid in stid_to_existing_uuid:
            return stid_to_existing_uuid[leaf_stid]
        if leaf_stid not in leaf_uuid_registry:
            leaf_uuid_registry[leaf_stid] = str(uuid.uuid4())
            reactome_id_to_uuid[leaf_uuid_registry[leaf_stid]] = leaf_stid
        return leaf_uuid_registry[leaf_stid]

    def _is_complex(entity_id: str) -> bool:
        return "Complex" in get_labels(entity_id)

    assembly_count = 0
    for eid in root_input_eids:
        if not _is_complex(eid):
            continue
        complex_uuid = root_input_uuid_cache.get(eid)
        if not complex_uuid:
            continue
        leaves = get_terminal_components(eid)
        # If the only "leaf" is the complex itself, there's nothing to expose.
        if leaves == {str(eid)}:
            continue
        for leaf in leaves:
            pathway_logic_network_data.append({
                "source_id": _leaf_uuid(leaf),
                "target_id": complex_uuid,
                "pos_neg": "pos",
                "and_or": "and",
                "edge_type": "assembly",
                "stoichiometry": 1,
            })
            assembly_count += 1

    dissociation_count = 0
    for eid in terminal_output_eids:
        if not _is_complex(eid):
            continue
        complex_uuid = terminal_output_uuid_cache.get(eid)
        if not complex_uuid:
            continue
        leaves = get_terminal_components(eid)
        if leaves == {str(eid)}:
            continue
        for leaf in leaves:
            pathway_logic_network_data.append({
                "source_id": complex_uuid,
                "target_id": _leaf_uuid(leaf),
                "pos_neg": "pos",
                "and_or": "and",
                "edge_type": "dissociation",
                "stoichiometry": 1,
            })
            dissociation_count += 1

    if assembly_count or dissociation_count:
        logger.info(
            f"Boundary expansion: {assembly_count} assembly edges, "
            f"{dissociation_count} dissociation edges, "
            f"{len(leaf_uuid_registry)} unique boundary leaves"
        )


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
    # Build stId → UUID lookup. Seeded from entity_uuid_registry so that a
    # protein appearing both as a regular VR input and as a regulator shares
    # one UUID. We also UPDATE this map when minting fresh UUIDs below, so
    # the same regulator entity used across multiple regulator-rows (e.g.
    # MDM2 regulating R1 and R2) shares a single UUID across emissions
    # rather than getting disconnected duplicate nodes per emission.
    stid_to_existing_uuid: Dict[str, str] = {}
    if entity_uuid_registry:
        for (entity_dbId, _reaction_uuid, _role), entity_uuid in entity_uuid_registry.items():
            if entity_dbId not in stid_to_existing_uuid:
                stid_to_existing_uuid[entity_dbId] = entity_uuid

    regulator_configs = [
        (catalyst_map, "pos", "catalyst"),
        (negative_regulator_map, "neg", "regulator"),
        (positive_regulator_map, "pos", "regulator"),
    ]

    for map_df, pos_neg, edge_type in regulator_configs:
        for _, row in map_df.iterrows():
            entity_id = str(row["entity_id"])

            terminal_members = _decompose_regulator_entity(entity_id)

            # and_or expresses reaction-level requirement, not within-entity
            # decomposition logic. Anything that contributes to a reaction
            # proceeding (catalyst, positive regulator) is "and"; anything
            # that blocks it (negative regulator) is "or" because any one
            # blocker suffices. The Complex/EntitySet decomposition tree
            # is preserved in decomposed_uid_mapping.csv.
            and_or = "and" if pos_neg == "pos" else "or"

            for member_id, member_stoich in terminal_members:
                if member_id in stid_to_existing_uuid:
                    member_uuid = stid_to_existing_uuid[member_id]
                else:
                    member_uuid = str(uuid.uuid4())
                    stid_to_existing_uuid[member_id] = member_uuid
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
        - See docs/UUID_DESIGN.md for detailed design documentation
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
    
    # Calculate and print statistics
    _calculate_reaction_statistics(reaction_connections)
    
    # Create mappings and connections
    reaction_id_map = create_reaction_id_map(decomposed_uid_mapping, best_matches)
    graph = get_graph()
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
    # Merges accumulate in a union-find map; we canonicalize the registry
    # once after the loop instead of scanning every entry on every merge.
    uuid_unions: Dict[str, str] = {}
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
                        eid, p_vr, f_vr, entity_uuid_registry, uuid_unions,
                    )
                    merge_count += 1

    _canonicalize_registry(entity_uuid_registry, uuid_unions)
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
            # 'or' when this entity is produced by multiple VRs (any one
            # source can supply it). None for single-producer outputs:
            # there's no AND/OR relationship to express, the entity simply
            # IS produced. None > "" because "" silently leaks into CSV
            # as an inconsistent empty string distinct from NaN.
            and_or = "or" if entity_producer_count.get(eid, 0) > 1 else None
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
    for _, row in pd.concat(
        [catalyst_map, negative_regulator_map, positive_regulator_map]
    ).iterrows():
        if pd.notna(row.get("entity_id")):
            cat_reg_entity_ids.add(str(row["entity_id"]))

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

    # Boundary expansion: root-input and terminal-output complexes get
    # synthetic assembly / dissociation edges to their leaf components so
    # individual proteins are perturbable / readable at the network
    # boundary. Intermediate complexes are deliberately left intact —
    # they're the actual biological species flowing between reactions.
    # See docs/DESIGN_DECISIONS.md, "Two layers of decomposition."
    _emit_boundary_decomposition_edges(
        pathway_logic_network_data=pathway_logic_network_data,
        root_input_eids=root_input_eids,
        terminal_output_eids=terminal_output_eids,
        root_input_uuid_cache=root_input_uuid_cache,
        terminal_output_uuid_cache=terminal_output_uuid_cache,
        reactome_id_to_uuid=reactome_id_to_uuid,
    )

    # Create final DataFrame
    pathway_logic_network = pd.DataFrame(pathway_logic_network_data, columns=list(columns.keys()))
    
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
    """Return distinct UUIDs that appear as a source but never as a target."""
    sources = pathway_logic_network["source_id"].dropna()
    targets = set(pathway_logic_network["target_id"].dropna().unique())
    return sources[~sources.isin(targets)].unique().tolist()


def find_terminal_outputs(pathway_logic_network: pd.DataFrame) -> List[Any]:
    """Return distinct UUIDs that appear as a target but never as a source."""
    targets = pathway_logic_network["target_id"].dropna()
    sources = set(pathway_logic_network["source_id"].dropna().unique())
    return targets[~targets.isin(sources)].unique().tolist()


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
        if cat_reg_uuid in all_uuids and pd.notna(row.get('entity_id')):
            uuid_to_reactome[cat_reg_uuid] = str(row['entity_id'])

    # Create DataFrame and save
    mapping_rows = [{'uuid': uuid, 'stable_id': stable_id}
                    for uuid, stable_id in uuid_to_reactome.items()]

    mapping_df = pd.DataFrame(mapping_rows, columns=['uuid', 'stable_id'])
    mapping_df = mapping_df.sort_values('uuid')  # Sort for easier lookup

    mapping_df.to_csv(output_file, index=False)
    logger.info(f"Exported UUID to Reactome stable ID mapping with {len(mapping_df)} entries")
