import os
import uuid
from typing import Dict, List, Any, NamedTuple, Optional, Set, Tuple

import pandas as pd
from pandas import DataFrame
from py2neo import Graph  # type: ignore
from src.argument_parser import logger
from src.neo4j_connector import get_graph
from src.reaction_generator import (
    _complex_contains_entity_set,
    _UBIQUITIN_ENTITY_SET_IDS,
    get_terminal_components,
    MAX_VARIANTS,
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
    entity_uuid_registry: Dict[tuple, str] = {}


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
    return [value for value in df[column].tolist() if not _is_missing_value(value)]


def _is_missing_value(value: Any) -> bool:
    """Return True for real nulls and stringified nulls from cached CSV reloads."""
    if pd.isna(value):
        return True
    if isinstance(value, str) and value.strip() in {"None", "nan", "NaN", "<NA>"}:
        return True
    return False


def _is_missing_reference_value(value: Any) -> bool:
    """Return True when a UID/reference field should not become a graph ID."""
    if _is_missing_value(value):
        return True
    return isinstance(value, str) and value.strip() == ""


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
        nested_uids = [
            value for value in group["input_or_output_uid"].tolist()
            if not _is_missing_reference_value(value)
        ]
        terminal_ids = [
            value for value in group["input_or_output_reactome_id"].tolist()
            if not _is_missing_reference_value(value)
        ]
        stoich_map: Dict[str, int] = {}
        for _, row in group.iterrows():
            stoich_raw = row.get("stoichiometry")
            stoich = 1 if stoich_raw is None or pd.isna(stoich_raw) else int(stoich_raw)
            if not _is_missing_reference_value(row.get("input_or_output_uid")):
                stoich_map[row["input_or_output_uid"]] = stoich
            if not _is_missing_reference_value(row.get("input_or_output_reactome_id")):
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


def _parse_variant_members(variant_id: str) -> Set[str]:
    """Terminal member stIds encoded in a ``{parent}::variant::{m1_m2}`` id."""
    if "::variant::" not in variant_id:
        return set()
    tail = variant_id.split("::variant::", 1)[1]
    return {m for m in tail.split("_") if m}


_variant_leafsets_cache: Dict[str, List[frozenset]] = {}
# Complexes whose set-variant expansion exceeded MAX_VARIANTS (issue #40) and were
# therefore bundled into one opaque node. Emission uses this to name them by their
# plain stId rather than a giant ``::variant::<all-leaves>`` id.
_variant_capped: Set[str] = set()


def _complex_variant_leafsets(complex_id: str) -> List[frozenset]:
    """Enumerate a complex's set-variants as FULL terminal-leaf sets.

    One frozenset per variant, each the complete terminal membership of that
    variant (fixed components + one choice per internal EntitySet + recursively
    for nested complexes). Unlike :func:`_expand_complex_variants` this returns
    flat leaf sets (no nested ``::variant::`` ids), which is what the emission
    node id and the member-set matching need.
    """
    import itertools
    from src.neo4j_connector import get_labels, get_complex_components, get_set_members

    if complex_id in _variant_leafsets_cache:
        return _variant_leafsets_cache[complex_id]

    components = get_complex_components(complex_id)
    if not components:
        result = [frozenset(get_terminal_components(complex_id))]
        _variant_leafsets_cache[complex_id] = result
        return result

    per_component_choices: List[List[frozenset]] = []
    for member_id in components:
        labels = get_labels(member_id)
        if "Complex" in labels and _complex_contains_entity_set(member_id):
            per_component_choices.append(_complex_variant_leafsets(member_id))
        elif (
            any(lbl in labels for lbl in ("EntitySet", "DefinedSet", "CandidateSet"))
            and member_id not in _UBIQUITIN_ENTITY_SET_IDS
        ):
            choices = [frozenset(get_terminal_components(sm)) for sm in get_set_members(member_id)]
            per_component_choices.append(choices or [frozenset(get_terminal_components(member_id))])
        else:
            per_component_choices.append([frozenset(get_terminal_components(member_id))])

    cartesian_size = 1
    for choices in per_component_choices:
        cartesian_size *= max(1, len(choices))
    if MAX_VARIANTS > 0 and cartesian_size > MAX_VARIANTS:
        # Too many variants (issue #40): bundle the whole complex into one
        # opaque leaf-set rather than enumerating the product. Mirrors the
        # matching-layer cap in reaction_generator.get_broken_apart_ids so both
        # layers agree this complex is a single node.
        logger.warning(
            f"variant cap hit for complex {complex_id}: {cartesian_size} "
            f"variants > {MAX_VARIANTS}; bundling into one node"
        )
        _variant_capped.add(complex_id)
        result = [frozenset(get_terminal_components(complex_id))]
        _variant_leafsets_cache[complex_id] = result
        return result

    variants: List[frozenset] = []
    seen: Set[frozenset] = set()
    for combo in itertools.product(*per_component_choices):
        leaves = frozenset().union(*combo) if combo else frozenset()
        if leaves and leaves not in seen:
            seen.add(leaves)
            variants.append(leaves)
    result = variants or [frozenset(get_terminal_components(complex_id))]
    _variant_leafsets_cache[complex_id] = result
    return result


def _map_annotated_entity_to_nodes(entity_id: str, member_set: Set[str]) -> Set[str]:
    """Map one reaction-annotated input/output entity to its emission node(s).

    This is where "a complex is a single node, split only by its internal sets"
    is enforced — at *emission* time, using the reaction's real annotated
    entities (unambiguous), not the content-addressed matching hashes.

    * simple entity / simple complex → the entity's own stId (one bundled node)
    * complex-that-contains-a-set → the specific set-VARIANT node selected by
      this virtual reaction, named ``{complex}::variant::{sorted members}``.
      The variant is chosen by matching its members against ``member_set`` (the
      terminal members this VR actually resolved to), so we pick the right
      alternative rather than emitting all of them.
    * bare EntitySet → the terminal member(s) present in this VR (sets expand)
    """
    from src.neo4j_connector import get_labels
    labels = get_labels(entity_id)

    if "Complex" in labels:
        if not _complex_contains_entity_set(entity_id):
            return {str(entity_id)}  # simple complex → single bundled node
        # Set-variant node: OUTERMOST complex stId + a FLAT sorted list of the
        # variant's FULL terminal membership. We enumerate the complex's true
        # variants (each with complete membership) and pick the one this VR
        # selected — the variant whose leaves are all present in member_set
        # (largest such, to prefer the fullest match). Enumerating true variants
        # (rather than intersecting all-possible-leaves with member_set) avoids
        # emitting spurious partial variants missing a subunit. Id is kept flat
        # so the parent is always ``id.split("::variant::")[0]``.
        variant_leafsets = _complex_variant_leafsets(entity_id)
        if entity_id in _variant_capped:
            # Over the variant cap (issue #40): opaque bundle, named by stId.
            return {str(entity_id)}
        subset = [ls for ls in variant_leafsets if ls <= member_set]
        if subset:
            chosen = max(subset, key=len)
        else:
            chosen = max(variant_leafsets, key=lambda ls: len(ls & member_set),
                         default=frozenset())
            if not (chosen & member_set):
                return {str(entity_id)}  # no fit → fall back to plain complex
        return {f"{entity_id}::variant::{'_'.join(sorted(chosen))}"}

    if any(lbl in labels for lbl in ("EntitySet", "DefinedSet", "CandidateSet")):
        if entity_id in _UBIQUITIN_ENTITY_SET_IDS:
            return {str(entity_id)}
        present = get_terminal_components(entity_id) & member_set
        return present if present else {str(entity_id)}

    return {str(entity_id)}  # simple entity (protein / small molecule / …)


def _resolve_vr_entities(
    reaction_id_map: pd.DataFrame,
    uid_index: Dict[str, tuple]
) -> Dict[str, tuple]:
    """Resolve each virtual reaction's inputs/outputs to emission NODES.

    A node is a bundled complex (or set-variant of one), a bare-set member, or a
    free simple entity — NOT the complex's individual member proteins. The
    matching layer (content hashes) is used only to learn which terminal members
    this VR selected; the node *identities* come from the reaction's real
    annotated entities via :func:`_map_annotated_entity_to_nodes`, which keeps
    parent-complex provenance unambiguous (content hashes collide across
    complexes with identical member content).

    Caches the resolution so Phase 2 and Phase 3 don't re-resolve.

    Returns:
        Dict mapping vr_uid -> (input_node_ids, output_node_ids,
                                input_stoich_map, output_stoich_map)
    """
    from src.neo4j_connector import get_reaction_input_output_ids

    annotated_cache: Dict[tuple, Set[str]] = {}

    def _annotated(reaction_id: str, io: str) -> Set[str]:
        key = (reaction_id, io)
        if key not in annotated_cache:
            annotated_cache[key] = set(get_reaction_input_output_ids(reaction_id, io))
        return annotated_cache[key]

    vr_entities: Dict[str, tuple] = {}
    for _, row in reaction_id_map.iterrows():
        vr_uid = row["uid"]
        reaction_id = str(row["reactome_id"])
        input_members = set(_resolve_to_terminal_reactome_ids(uid_index, row["input_hash"]))
        output_members = set(_resolve_to_terminal_reactome_ids(uid_index, row["output_hash"]))

        input_ids: Set[str] = set()
        for e in _annotated(reaction_id, "input"):
            input_ids |= _map_annotated_entity_to_nodes(str(e), input_members)
        output_ids: Set[str] = set()
        for e in _annotated(reaction_id, "output"):
            output_ids |= _map_annotated_entity_to_nodes(str(e), output_members)

        vr_entities[vr_uid] = (
            list(input_ids), list(output_ids),
            {n: 1 for n in input_ids}, {n: 1 for n in output_ids},
        )
    return vr_entities


def _decompose_regulator_entity(
    entity_id: str,
    variant_decomposition: bool = False,
) -> List[tuple]:
    """Decompose a catalyst/regulator entity to (terminal_id, stoichiometry) pairs.

    Two decomposition modes:

    **subunit decomposition** (default, used for catalysts and positive
    regulators): Complex with EntitySet → cartesian over members down to
    terminal proteins; simple Complex returned intact; EntitySet → flat
    alternatives. Each terminal subunit becomes its own edge — biologically
    appropriate for catalysts where every subunit of the holoenzyme is
    required (AND).

    **variant decomposition** (``variant_decomposition=True``, used for
    negative regulators): Complex with EntitySet → one entity per cartesian
    variant of the EntitySet expansion, but the complex itself is preserved
    (NOT broken into individual proteins). A bare EntitySet expands to its
    member alternatives. This is what you want for inhibitors, because:
    the inhibitor complex acts as a single biological unit; breaking it
    into individual proteins (HSP90, CDC37, etc.) would make those
    bystander proteins act as standalone inhibitors, which spuriously
    crushes downstream signal whenever unrelated reactions produce them.

    Ubiquitin sets are treated as atomic in both modes to avoid
    combinatorial explosion. AND/OR semantics for the resulting edges are
    decided by ``append_regulators`` based on pos_neg, not by within-entity
    decomposition shape.
    """
    from src.neo4j_connector import get_labels, get_complex_components, get_set_members

    labels = get_labels(entity_id)

    if "Complex" in labels:
        if not _complex_contains_entity_set(entity_id):
            return [(entity_id, 1)]
        if variant_decomposition:
            return _expand_complex_variants(entity_id)
        components = get_complex_components(entity_id)  # Dict[str, int]
        result = []
        for member_id, stoich in components.items():
            for mid, sub_stoich in _decompose_regulator_entity(
                member_id, variant_decomposition=variant_decomposition
            ):
                result.append((mid, stoich * sub_stoich))
        return result if result else [(entity_id, 1)]

    if "EntitySet" in labels or "DefinedSet" in labels or "CandidateSet" in labels:
        if entity_id in _UBIQUITIN_ENTITY_SET_IDS:
            return [(entity_id, 1)]
        members = get_set_members(entity_id)
        result = []
        for member_id in members:
            result.extend(
                _decompose_regulator_entity(
                    member_id, variant_decomposition=variant_decomposition
                )
            )
        return result if result else [(entity_id, 1)]

    return [(entity_id, 1)]


def _expand_complex_variants(complex_id: str) -> List[tuple]:
    """Expand a Complex-with-EntitySet into its cartesian-product variants.

    Each variant is itself a complex (same proteins, one specific choice
    per internal EntitySet), kept as a single biological entity. Returns
    [(variant_id, 1), ...] where each variant_id is a deterministic
    synthetic ID of the form ``{parent_stid}::variant::{sorted_member_ids}``.
    The synthetic ID is content-addressed: the same combination of leaf
    members under the same parent complex always produces the same ID, so
    cross-pathway references to "the same variant" are consistent.

    For a simple Complex (no EntitySet inside) the input is returned
    unchanged — caller already short-circuits on this case, but we
    re-check here for safety.
    """
    import itertools
    from src.neo4j_connector import get_labels, get_complex_components, get_set_members

    if not _complex_contains_entity_set(complex_id):
        return [(complex_id, 1)]

    components = get_complex_components(complex_id)
    if not components:
        return [(complex_id, 1)]

    # For each component, collect the list of identities it can take in a
    # single variant. A simple member contributes [member_id]; an internal
    # EntitySet contributes [alt_1, alt_2, ...]; a nested Complex (rare)
    # contributes its own variant IDs.
    per_component_choices: List[List[str]] = []
    for member_id, _stoich in components.items():
        labels = get_labels(member_id)
        if "Complex" in labels:
            sub_variants = _expand_complex_variants(member_id)
            per_component_choices.append([vid for vid, _ in sub_variants])
        elif (
            ("EntitySet" in labels or "DefinedSet" in labels or "CandidateSet" in labels)
            and member_id not in _UBIQUITIN_ENTITY_SET_IDS
        ):
            alts: List[str] = []
            for set_member in get_set_members(member_id):
                set_member_labels = get_labels(set_member)
                if (
                    "Complex" in set_member_labels
                    and _complex_contains_entity_set(set_member)
                ):
                    alts.extend(vid for vid, _ in _expand_complex_variants(set_member))
                else:
                    alts.append(set_member)
            per_component_choices.append(alts if alts else [member_id])
        else:
            per_component_choices.append([member_id])

    variants: List[tuple] = []
    seen_variant_ids: set = set()
    for combo in itertools.product(*per_component_choices):
        combo_sorted = sorted(combo)
        variant_id = f"{complex_id}::variant::{'_'.join(combo_sorted)}"
        if variant_id in seen_variant_ids:
            continue
        seen_variant_ids.add(variant_id)
        variants.append((variant_id, 1))

    return variants if variants else [(complex_id, 1)]


_COFACTOR_STIDS: frozenset = frozenset({
    "R-ALL-113592",  # ATP
    "R-ALL-29358",   # ATP variant
    "R-ALL-113582",  # ADP
    "R-ALL-29370",   # ADP variant
    "R-ALL-29360",   # ADP variant
    "R-ALL-29356",   # H2O
    "R-ALL-29372",   # Pi
    "R-ALL-29390",   # Pi variant
    "R-ALL-29438",   # PPi
    "R-ALL-217093",  # NADP+
    "R-ALL-110114",  # NADPH
    "R-ALL-29986",   # NAD+
    "R-ALL-73473",   # NADH
})

# Ubiquitin entity stIds (human + cross-species variants). A reaction that
# takes one of these as INPUT is a ubiquitination reaction (Ub is consumed
# and attached to a target protein). Reactions whose OUTPUT is Ub are
# deubiquitinations (we don't emit depletion for those).
_UBIQUITIN_STIDS: frozenset = frozenset({
    "R-HSA-68524",   # Ub [nucleoplasm]
    "R-HSA-113595",  # Ub [cytosol]
    "R-HSA-9660007", # Ub [lysosomal lumen]
})


def _emit_substrate_depletion_edges(
    pathway_logic_network_data: List[Dict[str, Any]],
    reactome_id_to_uuid: Dict[str, str],
    catalyst_map: pd.DataFrame,
) -> None:
    """Emit catalyst→input "depletion" edges for PHOSPHATASE reactions only.

    The biology to capture: when a catalyst REMOVES its substrate from the
    available pool, the substrate's level should respond to the catalyst's
    activity. PTEN dephosphorylates PIP3 → PTEN-KO needs to BOOST PIP3 in
    the propagator (de-repression via divide-form H on the depletion edge).

    SELECTION CRITERION: the reaction's OUTPUTS include Pi (R-ALL-29372).
    Cleanly identifies dephosphorylation reactions (PTEN, PTPN12, PHLPP,
    PP1, PP2A, etc.) without name heuristics. The free substrate is the
    direct input, so catalyst→input depletion targets the right node.

    NOTE: ubiquitin-ligase reactions look superficially similar but
    Reactome models them more precisely — the "MDM2 ubiquitinates TP53"
    reaction takes the ASSEMBLED `p-MDM2:MDM4:TP53` complex as input, not
    free TP53. So a catalyst→input depletion edge points MDM2 → complex,
    not MDM2 → TP53, and the heuristic doesn't help (empirically -0.5pp
    on the experimental benchmark). Capturing ubiquitin-driven depletion
    would require modeling the multi-step bind→ubiquitinate→degrade chain.

    The deltasignal solver applies divide-form inhibition to depletion edges
    specifically (DS_DEPLETION_H_MAX caps the de-repression boost, default
    10). Regular regulator edges stay on devspec (no false boost).

    Cofactor inputs (ATP/ADP/H2O/Pi/etc.) are excluded.
    """
    # Build reaction_uuid → list of input edges + list of catalyst edges
    # + list of output edges (to detect phosphatases by Pi output).
    by_target_inputs: Dict[str, List[str]] = {}
    by_target_catalysts: Dict[str, List[str]] = {}
    by_source_outputs: Dict[str, List[str]] = {}  # reaction_uuid → output stids
    PI_STID = "R-ALL-29372"  # inorganic phosphate
    for edge in pathway_logic_network_data:
        if edge.get("pos_neg") != "pos":
            continue
        et = edge.get("edge_type", "")
        if et == "input":
            by_target_inputs.setdefault(edge["target_id"], []).append(edge["source_id"])
        elif et == "catalyst":
            by_target_catalysts.setdefault(edge["target_id"], []).append(edge["source_id"])
        elif et == "output":
            # source is the reaction, target is the output entity
            by_source_outputs.setdefault(edge["source_id"], []).append(edge["target_id"])

    # Identify phosphatase reactions: those whose outputs include Pi (R-ALL-29372).
    phosphatase_rxn_uuids = set()
    for rxn_uuid, output_uuids in by_source_outputs.items():
        output_stids = {reactome_id_to_uuid.get(u, "") for u in output_uuids}
        if PI_STID in output_stids:
            phosphatase_rxn_uuids.add(rxn_uuid)

    # Identify ubiquitin-ligase reactions TOPOLOGICALLY: those that take
    # ubiquitin (Ub) as an INPUT. Reactome models ubiquitination by
    # consuming free Ub and producing a ubiquitinated form of the target,
    # so this is a clean structural signal (no display-name string matching
    # required). For each such reaction we then identify the SUBSTRATE by
    # gene-level matching between input-complex members and output members:
    # the protein appearing in both (by referenceEntity gene) is the one
    # being modified, and its free-form network nodes are the correct
    # depletion targets — NOT the assembled complex that appears as the
    # literal reaction input.

    # Collect all Reactome reaction ids in this pathway from catalyst_map
    # (each catalyst row has a reaction_id and reaction_uuid).
    rxn_uuid_to_rstid: Dict[str, str] = {}
    if not catalyst_map.empty:
        for _, row in catalyst_map.iterrows():
            ru = str(row["reaction_uuid"]); rs = str(row["reaction_id"])
            if ru and rs:
                rxn_uuid_to_rstid[ru] = rs
    unique_rstids = list({rs for rs in rxn_uuid_to_rstid.values() if rs})

    # For each candidate reaction, ask neo4j: is Ub an input? Then collect
    # input/output protein gene names (decomposing complexes/sets to leaves).
    # Returns reaction_stid → (set of substrate genes, set of input-protein stids).
    ubiquitin_subst_by_rstid: Dict[str, Tuple[Set[str], Set[str]]] = {}
    if unique_rstids:
        from src.neo4j_connector import get_graph
        try:
            # Reactions with any Ub stId as input
            ub_rows = get_graph().run(
                "UNWIND $ids AS id "
                "MATCH (rle:ReactionLikeEvent {stId: id})-[:input]->(ub:PhysicalEntity) "
                "WHERE ub.stId IN $ub_stids RETURN DISTINCT id AS rxn",
                ids=unique_rstids,
                ub_stids=list(_UBIQUITIN_STIDS),
            ).data()
            ubiq_rstids = [r["rxn"] for r in ub_rows]
            # For each ubiquitination reaction, gather input-leaf proteins
            # (decomposed through hasComponent/hasMember/hasCandidate) and
            # their referenceEntity genes; do the same for outputs.
            in_rows = get_graph().run(
                "UNWIND $ids AS id "
                "MATCH (rle:ReactionLikeEvent {stId: id})-[:input]->(inp:PhysicalEntity) "
                "WHERE NOT inp.stId IN $ub_stids "
                "OPTIONAL MATCH (inp)-[:hasComponent|hasMember|hasCandidate*0..3]->(leaf:PhysicalEntity) "
                "OPTIONAL MATCH (leaf)-[:referenceEntity]->(re:ReferenceEntity) "
                "RETURN id AS rxn, leaf.stId AS leaf_stid, re.geneName AS genes",
                ids=ubiq_rstids,
                ub_stids=list(_UBIQUITIN_STIDS),
            ).data() if ubiq_rstids else []
            out_rows = get_graph().run(
                "UNWIND $ids AS id "
                "MATCH (rle:ReactionLikeEvent {stId: id})-[:output]->(o:PhysicalEntity) "
                "WHERE NOT o.stId IN $ub_stids "
                "OPTIONAL MATCH (o)-[:hasComponent|hasMember|hasCandidate*0..3]->(leaf:PhysicalEntity) "
                "OPTIONAL MATCH (leaf)-[:referenceEntity]->(re:ReferenceEntity) "
                "RETURN id AS rxn, re.geneName AS genes",
                ids=ubiq_rstids,
                ub_stids=list(_UBIQUITIN_STIDS),
            ).data() if ubiq_rstids else []
            in_by_rxn: Dict[str, List[Tuple[Optional[str], Optional[List[str]]]]] = {}
            for r in in_rows:
                in_by_rxn.setdefault(r["rxn"], []).append((r["leaf_stid"], r["genes"]))
            out_by_rxn: Dict[str, Set[str]] = {}
            for r in out_rows:
                g = r["genes"]
                if g:
                    for gn in g:
                        out_by_rxn.setdefault(r["rxn"], set()).add(gn)
            for rxn in ubiq_rstids:
                output_genes = out_by_rxn.get(rxn, set())
                if not output_genes: continue
                # The substrate: input-leaf proteins whose gene also appears
                # in the output (modified form). Collect their leaf stids.
                subst_stids: Set[str] = set()
                subst_genes: Set[str] = set()
                for leaf_stid, leaf_genes in in_by_rxn.get(rxn, []):
                    if not leaf_stid or not leaf_genes: continue
                    common = output_genes & set(leaf_genes)
                    if common:
                        subst_stids.add(leaf_stid)
                        subst_genes.update(common)
                if subst_stids:
                    ubiquitin_subst_by_rstid[rxn] = (subst_genes, subst_stids)
        except Exception as exc:
            logger.warning(f"Ubiquitin topology lookup failed: {exc}")

    # Build network stid → list of UUIDs index from pathway_logic_network_data,
    # so we can target the substrate's NETWORK NODES (not just the leaf
    # protein stId, which may not be a direct node).
    stid_to_uuids_in_net: Dict[str, List[str]] = {}
    seen_uuids: Set[str] = set()
    for edge in pathway_logic_network_data:
        for uid in (edge.get("source_id"), edge.get("target_id")):
            if not uid or uid in seen_uuids: continue
            seen_uuids.add(uid)
            sid = reactome_id_to_uuid.get(uid, "")
            if sid:
                stid_to_uuids_in_net.setdefault(sid, []).append(uid)

    # Emit depletion edges. For phosphatase reactions: catalyst → input
    # (free substrate IS the input). For ubiquitin reactions: catalyst →
    # all network UUIDs of the substrate protein stId (free form, not the
    # complex that's the literal reaction input).
    seen_edges: set = set()
    n_emitted_phos = 0
    n_emitted_ub = 0
    # Phosphatase pass — existing behavior.
    for rxn_uuid in phosphatase_rxn_uuids:
        catalyst_uuids = by_target_catalysts.get(rxn_uuid, [])
        inputs = by_target_inputs.get(rxn_uuid, [])
        if not catalyst_uuids or not inputs:
            continue
        for cat_uuid in catalyst_uuids:
            cat_stid = reactome_id_to_uuid.get(cat_uuid, "")
            for inp_uuid in inputs:
                if cat_uuid == inp_uuid:
                    continue
                inp_stid = reactome_id_to_uuid.get(inp_uuid, "")
                if inp_stid == cat_stid and inp_stid:
                    continue  # same biological entity at different positions
                if inp_stid in _COFACTOR_STIDS:
                    continue
                key = (cat_uuid, inp_uuid)
                if key in seen_edges:
                    continue
                seen_edges.add(key)
                pathway_logic_network_data.append({
                    "source_id": cat_uuid,
                    "target_id": inp_uuid,
                    "pos_neg": "neg",
                    "and_or": "and",
                    "edge_type": "depletion",
                    "stoichiometry": 1.0,
                })
                n_emitted_phos += 1

    # Ubiquitin pass — emit catalyst → free-substrate-UUIDs depletion edges.
    # For each VR uuid corresponding to a ubiquitin reaction, get its
    # catalysts, then for each substrate stid, look up all network UUIDs of
    # that substrate (free TP53 nodes, not the assembled complex) and emit.
    for ru, rs in rxn_uuid_to_rstid.items():
        info = ubiquitin_subst_by_rstid.get(rs)
        if info is None:
            continue
        _, subst_stids = info
        catalyst_uuids = by_target_catalysts.get(ru, [])
        if not catalyst_uuids:
            continue
        for cat_uuid in catalyst_uuids:
            cat_stid = reactome_id_to_uuid.get(cat_uuid, "")
            for subst_stid in subst_stids:
                if subst_stid == cat_stid: continue
                if subst_stid in _COFACTOR_STIDS: continue
                target_uuids = stid_to_uuids_in_net.get(subst_stid, [])
                for tgt_uuid in target_uuids:
                    if tgt_uuid == cat_uuid: continue
                    key = (cat_uuid, tgt_uuid)
                    if key in seen_edges: continue
                    seen_edges.add(key)
                    pathway_logic_network_data.append({
                        "source_id": cat_uuid,
                        "target_id": tgt_uuid,
                        "pos_neg": "neg",
                        "and_or": "and",
                        "edge_type": "depletion",
                        "stoichiometry": 1.0,
                    })
                    n_emitted_ub += 1
    logger.info(
        f"Emitted {n_emitted_phos + n_emitted_ub} substrate-depletion edges "
        f"({len(phosphatase_rxn_uuids)} phosphatase reactions, "
        f"{len(ubiquitin_subst_by_rstid)} ubiquitin-ligase reactions with "
        f"identified substrate; topology-based)"
    )


_handoff_leaf_cache: Dict[str, frozenset] = {}


def _node_leaves(node_id: str) -> frozenset:
    """Non-cofactor terminal leaf stIds contained in a node id.

    variant node → the leaves encoded after ``::variant::``; simple complex →
    its terminal components; anything else → itself. Cofactors and ubiquitin are
    excluded so they can't act as spurious connectivity carriers.
    """
    if node_id in _handoff_leaf_cache:
        return _handoff_leaf_cache[node_id]
    if "::variant::" in node_id:
        s = {m for m in node_id.split("::variant::")[-1].split("_") if m.startswith("R-")}
    else:
        try:
            from src.neo4j_connector import get_labels
            s = set(get_terminal_components(node_id)) if "Complex" in get_labels(node_id) else {node_id}
        except Exception:
            s = {node_id}
    leaves = frozenset(s - _COFACTOR_STIDS - _UBIQUITIN_STIDS)
    _handoff_leaf_cache[node_id] = leaves
    return leaves


def _emit_precedingevent_handoff_edges(
    pathway_logic_network_data: List[Dict[str, Any]],
    reaction_connections: pd.DataFrame,
    reactome_to_vr: Dict[str, List[str]],
    vr_entities: Dict[str, tuple],
    entity_uuid_registry: Dict[tuple, str],
) -> None:
    """Restore curator-intended connectivity dropped by complex bundling.

    A ``precedingEvent`` is the curator asserting that at least one entity flows
    from the preceding reaction to the following one. The generator realizes that
    only when the two reactions share a *whole* entity. But a molecule is often
    handed off as a *component* — bound in a complex on one side, free/other-
    complex on the other — so bundling makes them different nodes and the link is
    lost. This is NOT a heuristic: the ``precedingEvent`` edge is in Neo4j.

    Rule (per Adam): for each precedingEvent pair, do nothing if the reactions are
    ALREADY connected by a shared whole entity (the curator's asserted entity is
    represented). Only when they share no whole entity — yet the curator says
    they're connected — add ONE bridge, between the output/input nodes that share
    the most non-cofactor components (the most likely real carrier). "Skip already-
    connected pairs" + "one bridge per gap" keeps this bounded by the number of
    precedingEvent gaps (hundreds), not the variant×variant blow-up.
    """
    # Carrier specificity: a leaf that appears in many nodes is a promiscuous
    # subunit (e.g. RBL2, a pocket protein in many complexes). Bridging on such a
    # hub spreads a perturbation to readouts it doesn't affect (false positives).
    # Count how many distinct nodes each leaf appears in; only leaves appearing
    # in <= HUB_MAX nodes are allowed to act as the transferred carrier.
    HUB_MAX = int(os.environ.get("LNG_HANDOFF_HUB_MAX", "3"))
    leaf_nodes: Dict[str, set] = {}
    all_node_ids: Set[str] = set()
    for (ins_, outs_, _si, _so) in vr_entities.values():
        all_node_ids.update(ins_); all_node_ids.update(outs_)
    for nid in all_node_ids:
        for lf in _node_leaves(nid):
            leaf_nodes.setdefault(lf, set()).add(nid)
    specific = {lf for lf, ns in leaf_nodes.items() if len(ns) <= HUB_MAX}

    existing = {(e["source_id"], e["target_id"]) for e in pathway_logic_network_data}
    seen: Set[tuple] = set()
    n = 0
    for _, conn in reaction_connections.iterrows():
        pre = conn.get("preceding_reaction_id"); fol = conn.get("following_reaction_id")
        if pd.isna(pre) or pd.isna(fol):
            continue
        # Gather (node_id, uuid) for the preceding reaction's outputs and the
        # following reaction's inputs, across all their virtual reactions.
        outs = []
        for p_vr in reactome_to_vr.get(pre, []):
            for on in vr_entities.get(p_vr, ([], [], {}, {}))[1]:
                u = entity_uuid_registry.get((on, p_vr, "output"))
                if u:
                    outs.append((on, u))
        ins = []
        for f_vr in reactome_to_vr.get(fol, []):
            for inn in vr_entities.get(f_vr, ([], [], {}, {}))[0]:
                u = entity_uuid_registry.get((inn, f_vr, "input"))
                if u:
                    ins.append((inn, u))
        if not outs or not ins:
            continue
        # Already connected? (whole-entity match -> Phase 2 gave a shared UUID)
        if {u for _, u in outs} & {u for _, u in ins}:
            continue
        # Otherwise bridge the single best output/input node pair, scored by the
        # number of shared SPECIFIC (non-hub) carrier components. A pair sharing
        # only hub subunits scores 0 and is not bridged.
        best = None; best_n = 0
        for on, ou in outs:
            lon = _node_leaves(on) & specific
            if not lon:
                continue
            for inn, iu in ins:
                if ou == iu:
                    continue
                sh = len(lon & _node_leaves(inn))
                if sh > best_n:
                    best_n = sh; best = (ou, iu)
        if best and best_n > 0 and best[0] != best[1] and best not in existing and best not in seen:
            seen.add(best)
            pathway_logic_network_data.append({
                "source_id": best[0],
                "target_id": best[1],
                "pos_neg": "pos",
                "and_or": "or",
                "edge_type": "handoff",
                "stoichiometry": 1,
                "edge_reaction_id": None,
            })
            n += 1
    logger.info(
        f"Emitted {n} precedingEvent hand-off edges "
        f"(one bridge per otherwise-disconnected precedingEvent gap)"
    )


def _emit_boundary_decomposition_edges(
    pathway_logic_network_data: List[Dict[str, Any]],
    reactome_id_to_uuid: Dict[str, str],
) -> None:
    """Expose the members of every root-input and terminal-output complex.

    Boundary membership is decided **positionally, per network occurrence** —
    not by a global stId set difference. A node is a *root input* if it is a
    source but never a target (no reaction produces it); a *terminal output* if
    it is a target but never a source (no reaction consumes it). The same
    complex stId can appear at several positions: the occurrences that are root
    inputs get decomposed, the occurrence sitting intermediate between two
    reactions is left intact — exactly as a real species should be.

    (The previous implementation used ``all_input_eids - all_output_eids`` over
    the whole pathway, which removed a complex's stId entirely the moment it was
    produced *anywhere*, leaving its genuine root-input occurrences undecomposed.
    Curator perturbations target individual proteins, so those proteins must be
    addressable wherever they enter or leave the pathway.)

    The two boundary directions produce deliberately OPPOSITE kinds of node:

    * **Assembly** (root input): for complex C with members {A, B, ...} emit
      ``A → C``, ``B → C`` (``edge_type='assembly'``). Members are *shared*
      upstream perturbation handles — a member reuses any UUID it already has in
      the network (else one freshly-minted UUID shared across occurrences). So a
      complex appearing as a root input at 15 positions gives member A one node
      with 15 assembly edges; knocking out A then propagates into all of them
      (the complex requires A — a real forward dependency).

    * **Dissociation** (terminal output): for complex C emit ``C → A``, ``C → B``
      (``edge_type='dissociation'``) where each member is a FRESH, SEPARATE
      readout node — one per (terminal-complex occurrence, member), carrying the
      member's stId, value inherited from the complex, and with NO outgoing
      edges. These are downstream *sinks*, NOT shared with the member's
      functional/assembly node. A terminal output is by definition consumed by
      nothing, so its members carry no forward signal; sharing them with the
      functional protein would (wrongly) inject the complex's activity into every
      other reaction that protein touches. Keeping them separate lets you measure
      how affected each member is *at that location* (read the sink; aggregate
      across a member's sinks for an overall figure) without any cross-talk.
    """
    from src.neo4j_connector import get_labels

    # Positional roots / terminals from the current edge list.
    sources: Set[str] = set()
    targets: Set[str] = set()
    for edge in pathway_logic_network_data:
        sources.add(edge["source_id"])
        targets.add(edge["target_id"])
    root_uuids = sources - targets       # produced by no reaction in this pathway
    terminal_uuids = targets - sources   # consumed by no reaction in this pathway

    # stId → existing UUID, so a member reuses the node it already has elsewhere
    # (free protein, regulator, catalyst) rather than becoming a disconnected dup.
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
        # Synthetic variant IDs (from negative-regulator variant decomposition)
        # are not in Neo4j and shouldn't be further decomposed at the boundary.
        # Tolerate any other unknown stIds too (e.g., entities added by other
        # synthetic emissions) — if the lookup fails, assume it's not a complex
        # and skip decomposition rather than crashing.
        if "::variant::" in entity_id:
            return False
        try:
            return "Complex" in get_labels(entity_id)
        except IndexError:
            return False

    seen_edges: Set[tuple] = set()
    assembly_count = 0
    for complex_uuid in root_uuids:
        stid = reactome_id_to_uuid.get(complex_uuid) or ""
        if not stid or not _is_complex(stid):
            continue
        leaves = get_terminal_components(stid)
        if leaves == {str(stid)}:  # nothing below the complex to expose
            continue
        for leaf in leaves:
            leaf_uuid = _leaf_uuid(leaf)
            if (leaf_uuid, complex_uuid) in seen_edges:
                continue
            seen_edges.add((leaf_uuid, complex_uuid))
            pathway_logic_network_data.append({
                "source_id": leaf_uuid,
                "target_id": complex_uuid,
                "pos_neg": "pos",
                "and_or": "and",
                "edge_type": "assembly",
                "stoichiometry": 1,
            })
            assembly_count += 1

    dissociation_count = 0
    for complex_uuid in terminal_uuids:
        stid = reactome_id_to_uuid.get(complex_uuid) or ""
        if not stid or not _is_complex(stid):
            continue
        leaves = get_terminal_components(stid)
        if leaves == {str(stid)}:
            continue
        for leaf in leaves:
            # Fresh per-occurrence readout sink — NOT _leaf_uuid (which would
            # share the member's functional node and re-introduce cross-talk).
            readout_uuid = str(uuid.uuid4())
            reactome_id_to_uuid[readout_uuid] = leaf
            pathway_logic_network_data.append({
                "source_id": complex_uuid,
                "target_id": readout_uuid,
                "pos_neg": "pos",
                "and_or": "and",
                "edge_type": "dissociation",
                "stoichiometry": 1,
            })
            dissociation_count += 1

    if assembly_count or dissociation_count:
        logger.info(
            f"Boundary expansion (positional): {assembly_count} assembly edges "
            f"(shared member handles), {dissociation_count} dissociation edges "
            f"(separate readout sinks), {len(leaf_uuid_registry)} new assembly leaves"
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
        # Negative regulators use VARIANT decomposition: an inhibitor complex
        # with an internal EntitySet expands into one entity per cartesian
        # variant of that EntitySet, but each variant is kept as a single
        # complex — NOT broken down into individual subunits. Without this,
        # HSP90 / CDC37 / ERBIN of an ERBB2 inhibitor complex would each
        # become standalone inhibitor edges, and any unrelated reaction
        # producing those bystander proteins would spuriously crush
        # downstream signal.
        #
        # Catalysts and positive regulators keep SUBUNIT decomposition: each
        # holoenzyme subunit is biologically AND-required for catalysis, so
        # decomposing to terminal proteins is correct there.
        variant_decomposition = (pos_neg == "neg")

        for _, row in map_df.iterrows():
            entity_id = str(row["entity_id"])

            terminal_members = _decompose_regulator_entity(
                entity_id, variant_decomposition=variant_decomposition
            )

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
        "edge_reaction_id": pd.Series(dtype="str"),
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

    # vr_uid -> Reactome reaction stId, for edge provenance (edge_reaction_id).
    vr_to_reaction: Dict[str, str] = dict(
        zip(reaction_id_map["uid"].astype(str), reaction_id_map["reactome_id"].astype(str))
    )

    for vr_uid, (input_ids, output_ids, input_stoich, output_stoich) in vr_entities.items():
        if not input_ids or not output_ids:
            continue
        reaction_stid = vr_to_reaction.get(str(vr_uid))

        for eid in input_ids:
            input_uuid = entity_uuid_registry[(eid, vr_uid, "input")]
            pathway_logic_network_data.append({
                "source_id": input_uuid,
                "target_id": vr_uid,
                "pos_neg": "pos",
                "and_or": "and",
                "edge_type": "input",
                "stoichiometry": input_stoich.get(eid, 1),
                "edge_reaction_id": reaction_stid,
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
                "edge_reaction_id": reaction_stid,
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

    # Substrate-depletion edges: for each catalytic reaction, emit edges
    # `catalyst → input` with edge_type="depletion". These capture the
    # biology that a catalyst REMOVES its substrate (PTEN dephosphorylates
    # PIP3, MDM2 ubiquitinates TP53, etc.). The deltasignal solver applies
    # divide-form inhibition to these edges specifically, so catalyst-knockout
    # boosts the substrate via de-repression. Skips small-molecule cofactor
    # inputs (ATP/H2O/Pi/etc.) to avoid noise.
    _emit_substrate_depletion_edges(
        pathway_logic_network_data=pathway_logic_network_data,
        reactome_id_to_uuid=reactome_id_to_uuid,
        catalyst_map=catalyst_map,
    )

    # Boundary expansion: every root-input and terminal-output complex
    # occurrence gets synthetic assembly / dissociation edges to its member
    # proteins, so individual subunits are perturbable / readable wherever
    # they enter or leave the pathway. Decided positionally per occurrence
    # (a complex that is also intermediate elsewhere keeps that intermediate
    # node intact). See docs/DESIGN_DECISIONS.md, "Two layers of decomposition."
    _emit_boundary_decomposition_edges(
        pathway_logic_network_data=pathway_logic_network_data,
        reactome_id_to_uuid=reactome_id_to_uuid,
    )

    # Restore curator-intended connectivity that complex-bundling drops: two
    # precedingEvent-linked reactions that hand off a shared COMPONENT (bound in
    # a complex on one side, free/other-complex on the other) share no whole
    # entity, so they were left disconnected. Honoring the curated precedingEvent
    # by bridging on the shared component was tried three ways (naive / max-shared
    # / hub-guarded) and was net-NEGATIVE on the benchmark every time: the bridge
    # injects roughly as much spurious coupling as real signal it recovers (same
    # reason the member-exploded network tied set-variant at ~69%). Kept for
    # future work but OFF by default. See memory project_complex_as_node_result.
    if os.environ.get("LNG_HANDOFF_EDGES", "0") != "0":
        _emit_precedingevent_handoff_edges(
            pathway_logic_network_data=pathway_logic_network_data,
            reaction_connections=reaction_connections,
            reactome_to_vr=reactome_to_vr,
            vr_entities=vr_entities,
            entity_uuid_registry=entity_uuid_registry,
        )

    # Create final DataFrame
    pathway_logic_network = pd.DataFrame(pathway_logic_network_data, columns=list(columns.keys()))
    # Coerce stoichiometry to nullable Int64 — emission sites use a mix of
    # int (`1`) and float (`1.0`) literals, which makes pandas infer float64
    # for the column and serialize as `1.0` in the CSV. Force integer so the
    # column reads as a clean whole-number count (Reactome's stoichiometries
    # are all integers for the pathways we care about; if a fractional value
    # ever appears it will surface as a parse error rather than silent loss).
    if not pathway_logic_network.empty:
        pathway_logic_network["stoichiometry"] = (
            pathway_logic_network["stoichiometry"].astype("Int64")
        )
    
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
        reaction_id_map=reaction_id_map,
        entity_uuid_registry=entity_uuid_registry,
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


def export_entity_reaction_proxy_mapping(
    pathway_logic_network: pd.DataFrame,
    reaction_id_map: pd.DataFrame,
    reactome_id_to_uuid: Dict[str, str],
    pathway_id: str,
    output_file: str,
) -> None:
    """Map curated species absent from the UUID mapping to a reaction-flux proxy.

    Complexes/EntitySets that contain an EntitySet are expanded into virtual
    variants when the logic network is built (see docs/DESIGN_DECISIONS.md,
    "EntitySet expansion produces multiple virtual reactions"). The variants get
    their own UUIDs but the *parent's* stId is not preserved anywhere in
    ``stid_to_uuid_mapping.csv`` — so a consumer that knows a curated species by
    its Reactome stId (e.g. MP-BioPath key-outputs, which are predominantly
    Complexes) can't find it, even though the reaction that produces it *is* in
    the network.

    The right proxy for "is species S present?" is the flux through the reaction
    that produces it: if S is the output of reaction R, then R's node activity is
    a direct, tight readout of S's production. Mapping S to its *terminal
    components* instead would be far too lax — a hub protein like β-catenin
    appears in dozens of reaction contexts, so a max over its component UUIDs
    reads "active" almost everywhere and discriminates nothing.

    For each PhysicalEntity that participates in the pathway but is not directly
    present in the UUID mapping, we therefore record the UUIDs of the reactions
    that **output** it (the producing reactions). If a species has no producing
    reaction in the network we fall back to reactions that **consume** it as an
    input — its presence is still implied by that reaction's activity.

    The primary ``stid_to_uuid_mapping.csv`` is intentionally left untouched: it
    keeps its one-row-per-UUID identity contract, and this supplementary file
    carries the (many-to-many) species → proxy-reaction relationship.

    Output CSV columns:
        - entity_stable_id: a species stId absent from the UUID mapping
        - proxy_uuid: a reaction UUID present in the network whose flux proxies it
        - proxy_role: 'producing' (entity is the reaction's output) or
                      'consuming' (entity is the reaction's input)
    """
    from src.neo4j_connector import (
        get_pathway_participating_entities,
        get_pathway_entity_reactions,
    )

    network_uuids: Set[str] = set()
    network_uuids.update(pathway_logic_network['source_id'].dropna().unique())
    network_uuids.update(pathway_logic_network['target_id'].dropna().unique())

    # stable_id of every entity directly addressable in the mapping.
    present_stids: Set[str] = set()
    if reactome_id_to_uuid:
        sample_key = next(iter(reactome_id_to_uuid.keys()))
        uuid_keyed = '-' in str(sample_key)
        for k, v in reactome_id_to_uuid.items():
            entity_uuid, stid = (k, v) if uuid_keyed else (v, k)
            if entity_uuid in network_uuids:
                present_stids.add(str(stid))

    # reaction stId -> [reaction UUIDs present in the network]
    reaction_stid_to_uuids: Dict[str, List[str]] = {}
    for _, row in reaction_id_map.iterrows():
        ruuid = str(row['uid'])
        if ruuid in network_uuids:
            reaction_stid_to_uuids.setdefault(str(row['reactome_id']), []).append(ruuid)

    participating = get_pathway_participating_entities(pathway_id)
    missing = {e for e in participating if e not in present_stids}
    if not missing:
        pd.DataFrame(columns=['entity_stable_id', 'proxy_uuid', 'proxy_role']).to_csv(
            output_file, index=False)
        logger.info("Entity-reaction proxy mapping: no missing species to proxy")
        return

    # entity stId -> {'output': [reaction stIds], 'input': [reaction stIds]}
    entity_reactions = get_pathway_entity_reactions(pathway_id, list(missing))

    rows: List[Dict[str, str]] = []
    for entity_stid in missing:
        roles = entity_reactions.get(entity_stid, {})
        # Prefer producing reactions; fall back to consuming if none are present.
        for role_name, rel_key in (('producing', 'output'), ('consuming', 'input')):
            proxy_uuids: List[str] = []
            for rxn_stid in roles.get(rel_key, []):
                proxy_uuids.extend(reaction_stid_to_uuids.get(str(rxn_stid), []))
            if proxy_uuids:
                for puuid in dict.fromkeys(proxy_uuids):  # de-dup, keep order
                    rows.append({'entity_stable_id': str(entity_stid),
                                 'proxy_uuid': puuid,
                                 'proxy_role': role_name})
                break  # don't also emit consuming rows once producing matched

    out_df = pd.DataFrame(rows, columns=['entity_stable_id', 'proxy_uuid', 'proxy_role'])
    if not out_df.empty:
        out_df = out_df.sort_values(['entity_stable_id', 'proxy_uuid'])
    out_df.to_csv(output_file, index=False)
    n_entities = out_df['entity_stable_id'].nunique() if not out_df.empty else 0
    logger.info(
        f"Exported entity-reaction proxy mapping: {len(out_df)} rows "
        f"covering {n_entities} of {len(missing)} missing species"
    )


# ---------------------------------------------------------------------------
# Schema-backed provenance exports (see schema/logic_network.linkml.yaml)
# ---------------------------------------------------------------------------
def _uuid_to_stable_id_map(pathway_logic_network: pd.DataFrame,
                           uuid_mapping: Dict[str, str]) -> Dict[str, str]:
    """uuid -> node string id (stId or ``{parent}::variant::{members}``).

    ``uuid_mapping`` (reactome_id_to_uuid) can be stored either direction; detect
    it the same way :func:`export_uuid_to_reactome_mapping` does.
    """
    all_uuids: Set[str] = set()
    all_uuids.update(pathway_logic_network["source_id"].dropna().astype(str).unique())
    all_uuids.update(pathway_logic_network["target_id"].dropna().astype(str).unique())
    out: Dict[str, str] = {}
    if uuid_mapping:
        sample = next(iter(uuid_mapping.keys()))
        uuid_keyed = isinstance(sample, str) and sample.count("-") >= 4 and "::" not in sample
        if uuid_keyed:
            for u, s in uuid_mapping.items():
                if str(u) in all_uuids:
                    out[str(u)] = str(s)
        else:
            for s, u in uuid_mapping.items():
                if str(u) in all_uuids:
                    out[str(u)] = str(s)
    return out


_sets_chosen_cache: Dict[str, tuple] = {}


def _derive_sets_and_chosen(parent_complex: str, member_leaves: Set[str]) -> tuple:
    """(source_set_ids, chosen_member_ids) for a set_variant, from Neo4j."""
    if parent_complex in _sets_chosen_cache:
        sets, chooser = _sets_chosen_cache[parent_complex]
    else:
        from src.neo4j_connector import get_labels, get_complex_components, get_set_members
        sets = []
        chooser = {}
        try:
            comps = get_complex_components(parent_complex)
        except Exception:
            comps = {}
        for m in comps:
            try:
                labels = get_labels(m)
            except Exception:
                continue
            if any(lbl in labels for lbl in ("EntitySet", "DefinedSet", "CandidateSet")):
                sets.append(str(m))
                for sm in get_set_members(m):
                    chooser.setdefault(str(m), []).append(
                        (str(sm), frozenset(get_terminal_components(sm)))
                    )
        _sets_chosen_cache[parent_complex] = (sets, chooser)
    chosen: List[str] = []
    for _set_id, options in chooser.items():
        for sm_id, sm_leaves in options:
            if sm_leaves & member_leaves:
                chosen.append(sm_id)
    return sets, chosen


def export_nodes(pathway_logic_network: pd.DataFrame,
                 reaction_id_map: pd.DataFrame,
                 uuid_mapping: Dict[str, str],
                 output_file: str) -> None:
    """Write nodes.csv — one row per node (see schema class ``Node``).

    Provenance is derived post-hoc from the node id string, the edge topology,
    and Neo4j: node_kind, diagram_entity_id (the stId a diagram renders — parent
    complex for a variant), member_leaves, and the set decomposition.
    """
    from src.neo4j_connector import get_labels
    uuid_to_str = _uuid_to_stable_id_map(pathway_logic_network, uuid_mapping)
    vr_uids = set(reaction_id_map["uid"].astype(str))
    vr_to_reaction = dict(zip(reaction_id_map["uid"].astype(str),
                              reaction_id_map["reactome_id"].astype(str)))

    incoming_types: Dict[str, Set[str]] = {}
    has_outgoing: Set[str] = set()
    for _, e in pathway_logic_network.iterrows():
        s, t = e.get("source_id"), e.get("target_id")
        if pd.notna(s):
            has_outgoing.add(str(s))
        if pd.notna(t):
            incoming_types.setdefault(str(t), set()).add(str(e.get("edge_type")))

    all_uuids = set(uuid_to_str) | (vr_uids & (has_outgoing | set(incoming_types)))
    rows: List[Dict[str, Any]] = []
    for u in sorted(all_uuids):
        kind = "other"; diagram = None
        members: List[str] = []; sets: List[str] = []; chosen: List[str] = []
        if u in vr_uids and u not in uuid_to_str:
            kind, diagram = "reaction", vr_to_reaction.get(u)
        else:
            s = uuid_to_str.get(u)
            if s is None:
                pass
            elif "::variant::" in s:
                kind = "set_variant"
                diagram = s.split("::variant::")[0]
                members = [m for m in s.split("::variant::")[-1].split("_")
                           if m.startswith("R-")]
                sets, chosen = _derive_sets_and_chosen(diagram, set(members))
            else:
                diagram = s
                inc = incoming_types.get(u, set())
                if "dissociation" in inc and u not in has_outgoing:
                    kind, members = "dissociation_sink", [s]
                else:
                    try:
                        labels = get_labels(s)
                    except Exception:
                        labels = []
                    if "Complex" in labels:
                        kind = "simple_complex"
                        try:
                            members = sorted(get_terminal_components(s))
                        except Exception:
                            members = [s]
                    else:
                        kind, members = "simple_entity", [s]
        rows.append({
            "uuid": u,
            "node_kind": kind,
            "diagram_entity_id": diagram,
            "compartment": None,
            "member_leaves": "|".join(members),
            "source_sets": "|".join(sets),
            "chosen_members": "|".join(chosen),
        })
    cols = ["uuid", "node_kind", "diagram_entity_id", "compartment",
            "member_leaves", "source_sets", "chosen_members"]
    pd.DataFrame(rows, columns=cols).to_csv(output_file, index=False)
    logger.info(f"Exported {len(rows)} nodes with provenance: {output_file}")


def export_node_reaction_context(entity_uuid_registry: Dict[tuple, str],
                                 reaction_id_map: pd.DataFrame,
                                 catalyst_regulator_map: pd.DataFrame,
                                 output_file: str) -> None:
    """Write node_reaction_context.csv — (node, reaction, role) location rows."""
    vr_to_reaction = dict(zip(reaction_id_map["uid"].astype(str),
                              reaction_id_map["reactome_id"].astype(str)))
    seen: Set[tuple] = set()
    rows: List[Dict[str, Any]] = []

    for (eid, vr_uid, role), node_uuid in (entity_uuid_registry or {}).items():
        rid = vr_to_reaction.get(str(vr_uid))
        if rid is None or role not in ("input", "output"):
            continue
        key = (str(node_uuid), rid, role)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"context_node": str(node_uuid), "reaction_id": rid, "role": role})

    if catalyst_regulator_map is not None and not catalyst_regulator_map.empty:
        for _, r in catalyst_regulator_map.iterrows():
            cr_uuid = r.get("uuid"); rid = r.get("reaction_id"); et = str(r.get("edge_type"))
            if pd.isna(cr_uuid) or pd.isna(rid):
                continue
            role = "catalyst" if et == "catalyst" else "regulator"
            key = (str(cr_uuid), str(rid), role)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"context_node": str(cr_uuid), "reaction_id": str(rid), "role": role})

    cols = ["context_node", "reaction_id", "role"]
    pd.DataFrame(rows, columns=cols).to_csv(output_file, index=False)
    logger.info(f"Exported {len(rows)} node-reaction-context rows: {output_file}")
