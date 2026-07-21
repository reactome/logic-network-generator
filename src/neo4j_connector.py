import os
from typing import Any, Dict, List, Optional, Set, Union

import pandas as pd
from py2neo import Graph  # type: ignore

from src.argument_parser import logger

# Lazy module-level Neo4j connection. The Graph object is only created
# when first needed, so importing this module doesn't fail when Neo4j
# is unreachable (e.g. in CI or during pytest collection).
_graph: Optional[Graph] = None


def get_graph() -> Graph:
    global _graph
    if _graph is None:
        url = os.getenv("NEO4J_URL", "bolt://localhost:7687")
        user = os.getenv("NEO4J_USER", "neo4j")
        password = os.getenv("NEO4J_PASSWORD", "test")
        _graph = Graph(url, auth=(user, password))
    return _graph

# Module-level caches for bulk pre-fetched data
_labels_cache: Dict[str, List[str]] = {}
_components_cache: Dict[str, Dict[str, int]] = {}
_members_cache: Dict[str, Set[str]] = {}
_reference_entity_cache: Dict[str, Optional[str]] = {}
_reaction_io_cache: Dict[str, Dict[str, Set[str]]] = {}
_prefetch_done: bool = False


def clear_prefetch_cache() -> None:
    """Clear all pre-fetched caches. Call before processing a new pathway."""
    global _labels_cache, _components_cache, _members_cache
    global _reference_entity_cache, _reaction_io_cache, _prefetch_done
    _labels_cache.clear()
    _components_cache.clear()
    _members_cache.clear()
    _reference_entity_cache.clear()
    _reaction_io_cache.clear()
    _prefetch_done = False


def prefetch_entity_data(reaction_ids: List[str]) -> None:
    """Pre-fetch all entity data for a set of reactions in bulk.

    Replaces thousands of individual Neo4j queries with 5 bulk queries,
    dramatically improving performance for pathways with many entities.

    Args:
        reaction_ids: List of Reactome reaction stable IDs to pre-fetch data for
    """
    global _labels_cache, _components_cache, _members_cache
    global _reference_entity_cache, _reaction_io_cache, _prefetch_done

    clear_prefetch_cache()

    g = get_graph()
    rid_list = list(reaction_ids)

    # Query 1: Get all reaction inputs and outputs
    logger.info(f"Bulk pre-fetching data for {len(rid_list)} reactions...")
    io_results = g.run(
        """
        MATCH (r:ReactionLikeEvent)-[rel:input|output]->(e)
        WHERE r.stId IN $reaction_ids
        RETURN r.stId AS reaction_id, type(rel) AS rel_type, e.stId AS entity_id
        """,
        reaction_ids=rid_list,
    ).data()

    direct_entity_ids: Set[str] = set()
    for row in io_results:
        rid = row["reaction_id"]
        eid = row["entity_id"]
        rel = row["rel_type"]
        direct_entity_ids.add(eid)

        if rid not in _reaction_io_cache:
            _reaction_io_cache[rid] = {"input": set(), "output": set()}
        _reaction_io_cache[rid][rel].add(eid)

    logger.info(f"Found {len(direct_entity_ids)} direct input/output entities")

    if not direct_entity_ids:
        _prefetch_done = True
        return

    # Query 2: Discover all descendant entities and their labels.
    # The *0..10 depth bound covers every Reactome nesting we've seen
    # (real depth maxes out around 5); 10 is a generous safety margin
    # without being expensive. Bumping it has not been necessary.
    logger.info("Discovering all descendant entities...")
    desc_results = g.run(
        """
        MATCH (root)-[:hasComponent|hasCandidate|hasMember*0..10]->(entity)
        WHERE root.stId IN $entity_ids
        RETURN DISTINCT entity.stId AS entity_id, labels(entity) AS entity_labels
        """,
        entity_ids=list(direct_entity_ids),
    ).data()

    all_entity_ids: Set[str] = set()
    for row in desc_results:
        eid = row["entity_id"]
        all_entity_ids.add(eid)
        _labels_cache[eid] = row["entity_labels"]

    logger.info(f"Found {len(all_entity_ids)} total entities (including descendants)")
    all_ids_list = list(all_entity_ids)

    # Query 3: All hasComponent relationships (Complex → components) with stoichiometry
    logger.info("Bulk fetching component relationships...")
    comp_results = g.run(
        """
        MATCH (parent)-[rel:hasComponent]->(child)
        WHERE parent.stId IN $entity_ids
        RETURN parent.stId AS parent_id, child.stId AS child_id,
               rel.stoichiometry AS stoichiometry
        """,
        entity_ids=all_ids_list,
    ).data()
    for row in comp_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _components_cache:
            _components_cache[pid] = {}
        _components_cache[pid][cid] = row.get("stoichiometry") or 1
    logger.info(f"Cached {len(_components_cache)} complex -> component mappings")

    # Query 4: All hasCandidate|hasMember relationships (EntitySet → members)
    logger.info("Bulk fetching member relationships...")
    member_results = g.run(
        """
        MATCH (parent)-[:hasCandidate|hasMember]->(child)
        WHERE parent.stId IN $entity_ids
        RETURN parent.stId AS parent_id, child.stId AS child_id
        """,
        entity_ids=all_ids_list,
    ).data()
    for row in member_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _members_cache:
            _members_cache[pid] = set()
        _members_cache[pid].add(cid)
    logger.info(f"Cached {len(_members_cache)} set -> member mappings")

    # Query 5: All HGNC reference entity IDs
    logger.info("Bulk fetching reference entity IDs...")
    ref_results = g.run(
        """
        MATCH (rd:ReferenceDatabase)<-[:referenceDatabase]-(reg:ReferenceEntity)
              <-[:referenceGene]-(re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE rd.displayName = 'HGNC'
          AND pe.stId IN $entity_ids
        RETURN pe.stId AS entity_id, re.stId AS reference_id
        """,
        entity_ids=all_ids_list,
    ).data()
    for row in ref_results:
        _reference_entity_cache[row["entity_id"]] = row["reference_id"]
    logger.info(f"Cached {len(_reference_entity_cache)} reference entity mappings")

    _prefetch_done = True
    logger.info("Bulk pre-fetch complete")


def prefetch_entity_decomposition_data(entity_ids: List[str]) -> None:
    """Pre-fetch decomposition data (labels, components, members) for entity IDs.

    Unlike prefetch_entity_data which starts from reaction IDs and fetches
    inputs/outputs, this function starts from entity IDs directly and only
    fetches the data needed to recursively decompose them (labels, components,
    members). Used for catalyst/regulator entities that aren't covered by the
    main reaction-based prefetch.

    Args:
        entity_ids: List of Reactome entity stable IDs to pre-fetch decomposition data for
    """
    global _labels_cache, _components_cache, _members_cache

    # Filter out entities already in cache
    uncached = [eid for eid in entity_ids if eid not in _labels_cache]
    if not uncached:
        return

    g = get_graph()

    # Discover all descendant entities and their labels.
    # See prefetch_entity_data for why *0..10 is the depth bound.
    logger.info(f"Pre-fetching decomposition data for {len(uncached)} catalyst/regulator entities...")
    desc_results = g.run(
        """
        MATCH (root)-[:hasComponent|hasCandidate|hasMember*0..10]->(entity)
        WHERE root.stId IN $entity_ids
        RETURN DISTINCT entity.stId AS entity_id, labels(entity) AS entity_labels
        """,
        entity_ids=uncached,
    ).data()

    new_entity_ids: Set[str] = set()
    for row in desc_results:
        eid = row["entity_id"]
        if eid not in _labels_cache:
            new_entity_ids.add(eid)
            _labels_cache[eid] = row["entity_labels"]

    if not new_entity_ids:
        logger.info("No new entities to fetch decomposition data for")
        return

    new_ids_list = list(new_entity_ids)

    # hasComponent relationships (Complex → components) with stoichiometry
    comp_results = g.run(
        """
        MATCH (parent)-[rel:hasComponent]->(child)
        WHERE parent.stId IN $entity_ids
        RETURN parent.stId AS parent_id, child.stId AS child_id,
               rel.stoichiometry AS stoichiometry
        """,
        entity_ids=new_ids_list,
    ).data()
    for row in comp_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _components_cache:
            _components_cache[pid] = {}
        _components_cache[pid][cid] = row.get("stoichiometry") or 1

    # hasCandidate|hasMember relationships (EntitySet → members)
    member_results = g.run(
        """
        MATCH (parent)-[:hasCandidate|hasMember]->(child)
        WHERE parent.stId IN $entity_ids
        RETURN parent.stId AS parent_id, child.stId AS child_id
        """,
        entity_ids=new_ids_list,
    ).data()
    for row in member_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _members_cache:
            _members_cache[pid] = set()
        _members_cache[pid].add(cid)

    logger.info(
        f"Pre-fetched decomposition data: {len(new_entity_ids)} entities, "
        f"{len(comp_results)} component relations, {len(member_results)} member relations"
    )


def get_reaction_connections(pathway_id: str) -> pd.DataFrame:
    """Get reaction connections for a pathway from Neo4j.

    Args:
        pathway_id: Reactome pathway stable ID (e.g., "R-HSA-69620")

    Returns:
        DataFrame with preceding_reaction_id, following_reaction_id, and event_status columns

    Raises:
        ConnectionError: If Neo4j database is not accessible
        ValueError: If pathway_id is invalid or pathway not found
    """
    query: str = """
        MATCH (pathway:Pathway)-[:hasEvent*]->(r1:ReactionLikeEvent)
        WHERE pathway.stId = $pathway_id
        OPTIONAL MATCH (r1)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway)
        WHERE pathway.stId = $pathway_id
        RETURN DISTINCT r1.stId AS preceding_reaction_id,
               r2.stId AS following_reaction_id,
               CASE WHEN r2 IS NULL THEN 'No Preceding Event' ELSE 'Has Preceding Event' END AS event_status
        """

    try:
        result = get_graph().run(query, pathway_id=pathway_id).data()
        df: pd.DataFrame = pd.DataFrame(result)

        if df.empty:
            raise ValueError(
                f"No reactions found for pathway ID: {pathway_id}. "
                f"Verify the pathway exists in Reactome database and Neo4j is running."
            )

        logger.info(f"Found {len(df)} reaction connections for pathway {pathway_id}")
        return df

    except ValueError:
        raise
    except Exception as e:
        logger.error(f"Error querying Neo4j for pathway {pathway_id}", exc_info=True)
        raise ConnectionError(
            f"Failed to connect to Neo4j database at "
            f"{os.getenv('NEO4J_URL', 'bolt://localhost:7687')}. "
            f"Ensure Neo4j is running and accessible. Original error: {str(e)}"
        ) from e


def get_top_level_pathways() -> List[Dict[str, Any]]:
    """Get all top-level pathways for Homo sapiens from Reactome.

    Top-level pathways are those that are not contained within another pathway
    (i.e., no incoming hasEvent relationship from another pathway).

    Returns:
        List of dicts with 'stId' and 'name' keys for each top-level pathway

    Raises:
        ConnectionError: If Neo4j database is not accessible
    """
    query: str = """
        MATCH (p:TopLevelPathway)
        WHERE p.speciesName = 'Homo sapiens'
        RETURN p.stId AS stId, p.displayName AS name
        ORDER BY p.displayName
    """

    try:
        result = get_graph().run(query).data()
        logger.info(f"Found {len(result)} top-level pathways")
        return result
    except Exception as e:
        logger.error("Error in get_top_level_pathways", exc_info=True)
        raise ConnectionError(
            f"Failed to query top-level pathways from Neo4j at "
            f"{os.getenv('NEO4J_URL', 'bolt://localhost:7687')}. "
            f"Ensure Neo4j is running and accessible. Original error: {str(e)}"
        ) from e


def get_pathway_participating_entities(pathway_id: str) -> Set[str]:
    """Return every PhysicalEntity stId that participates in the pathway's reactions.

    Walks all ReactionLikeEvents under the pathway (via hasEvent) and collects
    the stIds of their inputs, outputs, catalysts, and regulators. This includes
    the *intact* Complex/EntitySet entities as Reactome curates them — which is
    exactly the set that may get decomposed into virtual variants in the logic
    network and therefore lose their original stId from the UUID mapping.

    Args:
        pathway_id: Reactome pathway stable ID (e.g., "R-HSA-69620")

    Returns:
        Set of PhysicalEntity stable IDs.

    Raises:
        ConnectionError: If Neo4j database is not accessible
    """
    # Explicit, indexed patterns per role. A blanket variable-length match
    # (e.g. ``-[*1..3]-``) over these relationship types is orders of magnitude
    # slower and can time out on large pathways.
    query: str = """
        MATCH (pathway:Pathway {stId: $pathway_id})-[:hasEvent*]->(r:ReactionLikeEvent)
        OPTIONAL MATCH (r)-[:input|output]->(io:PhysicalEntity)
        OPTIONAL MATCH (r)-[:catalystActivity]->(:CatalystActivity)
                         -[:physicalEntity]->(cat:PhysicalEntity)
        OPTIONAL MATCH (r)-[:regulatedBy]->(:Regulation)
                         -[:regulator]->(reg:PhysicalEntity)
        RETURN COLLECT(DISTINCT io.stId)
             + COLLECT(DISTINCT cat.stId)
             + COLLECT(DISTINCT reg.stId) AS stids
    """
    try:
        result = get_graph().run(query, pathway_id=pathway_id).data()
        if not result:
            return set()
        return {s for s in result[0]["stids"] if s}
    except Exception as e:
        logger.error(
            f"Error in get_pathway_participating_entities for {pathway_id}",
            exc_info=True,
        )
        raise ConnectionError(
            f"Failed to query participating entities from Neo4j at "
            f"{os.getenv('NEO4J_URL', 'bolt://localhost:7687')}. "
            f"Original error: {str(e)}"
        ) from e


def get_pathway_entity_reactions(
    pathway_id: str, entity_ids: List[str]
) -> Dict[str, Dict[str, List[str]]]:
    """For each entity, the pathway reactions that output (produce) or input it.

    Args:
        pathway_id: Reactome pathway stable ID (e.g., "R-HSA-195721")
        entity_ids: PhysicalEntity stIds to look up.

    Returns:
        ``{entity_stId: {"output": [reaction_stId, ...],
                          "input":  [reaction_stId, ...]}}``
        Entities with no participating reaction in the pathway are absent.

    Raises:
        ConnectionError: If Neo4j database is not accessible
    """
    if not entity_ids:
        return {}
    query: str = """
        MATCH (pathway:Pathway {stId: $pathway_id})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[rel:input|output]->(e:PhysicalEntity)
        WHERE e.stId IN $entity_ids
        RETURN e.stId AS entity, type(rel) AS rel,
               COLLECT(DISTINCT r.stId) AS reactions
    """
    try:
        result = get_graph().run(
            query, pathway_id=pathway_id, entity_ids=list(entity_ids)
        ).data()
        out: Dict[str, Dict[str, List[str]]] = {}
        for row in result:
            out.setdefault(row["entity"], {})[row["rel"]] = [
                rid for rid in row["reactions"] if rid
            ]
        return out
    except Exception as e:
        logger.error(
            f"Error in get_pathway_entity_reactions for {pathway_id}", exc_info=True
        )
        raise ConnectionError(
            f"Failed to query entity reactions from Neo4j at "
            f"{os.getenv('NEO4J_URL', 'bolt://localhost:7687')}. "
            f"Original error: {str(e)}"
        ) from e


def get_pathway_name(pathway_id: str) -> str:
    """Get the display name for a pathway by its stable ID.

    Args:
        pathway_id: Reactome pathway stable ID (e.g., "R-HSA-69620")

    Returns:
        The display name of the pathway

    Raises:
        ValueError: If pathway not found
        ConnectionError: If Neo4j database is not accessible
    """
    query: str = """
        MATCH (p:Pathway)
        WHERE p.stId = $pathway_id
        RETURN p.displayName AS name
    """

    try:
        result = get_graph().run(query, pathway_id=pathway_id).data()
        if not result:
            raise ValueError(f"Pathway with ID {pathway_id} not found")
        return result[0]["name"]
    except ValueError:
        raise
    except Exception as e:
        logger.error(f"Error in get_pathway_name for {pathway_id}", exc_info=True)
        raise ConnectionError(
            f"Failed to query pathway name from Neo4j at "
            f"{os.getenv('NEO4J_URL', 'bolt://localhost:7687')}. "
            f"Original error: {str(e)}"
        ) from e


def get_labels(entity_id: str) -> List[str]:
    if entity_id in _labels_cache:
        return _labels_cache[entity_id]

    try:
        result = get_graph().run(
            "MATCH (e) WHERE e.stId = $entity_id RETURN labels(e) AS labels",
            entity_id=entity_id,
        ).data()[0]["labels"]
        _labels_cache[entity_id] = result
        return result
    except Exception:
        logger.error("Error in get_labels", exc_info=True)
        raise


def get_complex_components(entity_id: str) -> Dict[str, int]:
    if entity_id in _components_cache:
        return _components_cache[entity_id]
    if _prefetch_done:
        return {}  # Not in bulk results means no components

    try:
        data = get_graph().run(
            """
            MATCH (entity)-[rel:hasComponent]->(component)
            WHERE entity.stId = $entity_id
            RETURN component.stId AS component_id, rel.stoichiometry AS stoichiometry
            """,
            entity_id=entity_id,
        ).data()
        result = {row["component_id"]: row.get("stoichiometry") or 1 for row in data}
        _components_cache[entity_id] = result
        return result
    except Exception:
        logger.error("Error in get_complex_components", exc_info=True)
        raise


def get_set_members(entity_id: str) -> Set[str]:
    if entity_id in _members_cache:
        return _members_cache[entity_id]
    if _prefetch_done:
        return set()  # Not in bulk results means no members

    try:
        data = get_graph().run(
            """
            MATCH (entity)-[:hasCandidate|hasMember]->(member)
            WHERE entity.stId = $entity_id
            RETURN collect(member.stId) AS member_ids
            """,
            entity_id=entity_id,
        ).data()
        result = set(data[0]["member_ids"])
        _members_cache[entity_id] = result
        return result
    except Exception:
        logger.error("Error in get_set_members", exc_info=True)
        raise


def get_reaction_input_output_ids(reaction_id: str, input_or_output: str) -> Set[str]:
    if reaction_id in _reaction_io_cache:
        return _reaction_io_cache[reaction_id].get(input_or_output, set())

    # Cypher syntax doesn't allow parameterizing relationship types, so the
    # input/output choice has to be embedded — but only after restricting to
    # the known vocabulary, so it can't be smuggled to anything else.
    if input_or_output not in {"input", "output"}:
        raise ValueError(f"input_or_output must be 'input' or 'output', got {input_or_output!r}")

    query: str = (
        f"MATCH (reaction)-[:{input_or_output}]-(io) "
        f"WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) "
        f"AND reaction.stId = $reaction_id "
        f"RETURN COLLECT(io.stId) AS io_ids"
    )

    try:
        return set(get_graph().run(query, reaction_id=reaction_id).data()[0]["io_ids"])
    except Exception:
        logger.error("Error in get_reaction_input_output_ids", exc_info=True)
        raise


def get_reaction_io_stoichiometry(reaction_id: str, input_or_output: str) -> Dict[str, int]:
    """Curated input/output stoichiometry for a reaction, keyed by entity stId.

    Reactome stores the coefficient on the input/output relationship (the same
    ``stoichiometry`` property that hasComponent uses). An entity linked with no
    explicit coefficient defaults to 1; where the same entity is linked by more
    than one relationship the coefficients are summed. Companion to
    :func:`get_reaction_input_output_ids`, which returns only the entity ids.
    """
    # Same relationship-type allowlist as get_reaction_input_output_ids: the
    # type can't be parameterized in Cypher, so restrict it before embedding.
    if input_or_output not in {"input", "output"}:
        raise ValueError(f"input_or_output must be 'input' or 'output', got {input_or_output!r}")

    query: str = (
        f"MATCH (reaction)-[rel:{input_or_output}]-(io) "
        f"WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) "
        f"AND reaction.stId = $reaction_id "
        f"RETURN io.stId AS io_id, rel.stoichiometry AS stoichiometry"
    )

    try:
        data = get_graph().run(query, reaction_id=reaction_id).data()
    except Exception:
        logger.error("Error in get_reaction_io_stoichiometry", exc_info=True)
        raise

    result: Dict[str, int] = {}
    for row in data:
        io_id = row["io_id"]
        if io_id is None:
            continue
        stoich_raw = row.get("stoichiometry")
        stoich = 1 if stoich_raw is None else int(stoich_raw)
        result[io_id] = result.get(io_id, 0) + stoich
    return result


def get_reference_entity_id(entity_id: str) -> Union[str, None]:
    if entity_id in _reference_entity_cache:
        return _reference_entity_cache[entity_id]
    if _prefetch_done:
        return None  # Not in bulk results means no HGNC reference

    try:
        data = get_graph().run(
            """
            MATCH (rd:ReferenceDatabase)<-[:referenceDatabase]-(reg:ReferenceEntity)
                  <-[:referenceGene]-(re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
            WHERE rd.displayName = 'HGNC'
              AND pe.stId = $entity_id
            RETURN re.stId AS id
            """,
            entity_id=entity_id,
        ).data()
        if len(data) == 0:
            _reference_entity_cache[entity_id] = None
            return None
        result = data[0]["id"]
        _reference_entity_cache[entity_id] = result
        return result
    except Exception:
        logger.error("Error in get_reference_entity_id", exc_info=True)
        raise


