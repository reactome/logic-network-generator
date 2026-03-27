from typing import Any, Dict, List, Optional, Set, Union

import pandas as pd
from py2neo import Graph  # type: ignore

from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))

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

    ids_str = ", ".join(f"'{rid}'" for rid in reaction_ids)

    # Query 1: Get all reaction inputs and outputs
    logger.info(f"Bulk pre-fetching data for {len(reaction_ids)} reactions...")
    query_io = f"""
        MATCH (r:ReactionLikeEvent)-[rel:input|output]->(e)
        WHERE r.stId IN [{ids_str}]
        RETURN r.stId as reaction_id, type(rel) as rel_type, e.stId as entity_id
    """
    io_results = graph.run(query_io).data()

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

    direct_ids_str = ", ".join(f"'{eid}'" for eid in direct_entity_ids)

    # Query 2: Discover all descendant entities and their labels
    # Follows hasComponent/hasCandidate/hasMember relationships up to 10 levels deep
    logger.info("Discovering all descendant entities...")
    query_descendants = f"""
        MATCH (root)-[:hasComponent|hasCandidate|hasMember*0..10]->(entity)
        WHERE root.stId IN [{direct_ids_str}]
        RETURN DISTINCT entity.stId as entity_id, labels(entity) as entity_labels
    """
    desc_results = graph.run(query_descendants).data()

    all_entity_ids: Set[str] = set()
    for row in desc_results:
        eid = row["entity_id"]
        all_entity_ids.add(eid)
        _labels_cache[eid] = row["entity_labels"]

    logger.info(f"Found {len(all_entity_ids)} total entities (including descendants)")

    all_ids_str = ", ".join(f"'{eid}'" for eid in all_entity_ids)

    # Query 3: All hasComponent relationships (Complex → components) with stoichiometry
    logger.info("Bulk fetching component relationships...")
    query_components = f"""
        MATCH (parent)-[rel:hasComponent]->(child)
        WHERE parent.stId IN [{all_ids_str}]
        RETURN parent.stId as parent_id, child.stId as child_id, rel.stoichiometry as stoichiometry
    """
    comp_results = graph.run(query_components).data()
    for row in comp_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _components_cache:
            _components_cache[pid] = {}
        _components_cache[pid][cid] = row.get("stoichiometry") or 1
    logger.info(f"Cached {len(_components_cache)} complex -> component mappings")

    # Query 4: All hasCandidate|hasMember relationships (EntitySet → members)
    logger.info("Bulk fetching member relationships...")
    query_members = f"""
        MATCH (parent)-[:hasCandidate|hasMember]->(child)
        WHERE parent.stId IN [{all_ids_str}]
        RETURN parent.stId as parent_id, child.stId as child_id
    """
    member_results = graph.run(query_members).data()
    for row in member_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _members_cache:
            _members_cache[pid] = set()
        _members_cache[pid].add(cid)
    logger.info(f"Cached {len(_members_cache)} set -> member mappings")

    # Query 5: All HGNC reference entity IDs
    logger.info("Bulk fetching reference entity IDs...")
    query_ref = f"""
        MATCH (rd:ReferenceDatabase)<-[:referenceDatabase]-(reg:ReferenceEntity)
              <-[:referenceGene]-(re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE rd.displayName = "HGNC"
          AND pe.stId IN [{all_ids_str}]
        RETURN pe.stId as entity_id, re.stId as reference_id
    """
    ref_results = graph.run(query_ref).data()
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

    ids_str = ", ".join(f"'{eid}'" for eid in uncached)

    # Discover all descendant entities and their labels
    logger.info(f"Pre-fetching decomposition data for {len(uncached)} catalyst/regulator entities...")
    query_descendants = f"""
        MATCH (root)-[:hasComponent|hasCandidate|hasMember*0..10]->(entity)
        WHERE root.stId IN [{ids_str}]
        RETURN DISTINCT entity.stId as entity_id, labels(entity) as entity_labels
    """
    desc_results = graph.run(query_descendants).data()

    new_entity_ids: Set[str] = set()
    for row in desc_results:
        eid = row["entity_id"]
        if eid not in _labels_cache:
            new_entity_ids.add(eid)
            _labels_cache[eid] = row["entity_labels"]

    if not new_entity_ids:
        logger.info("No new entities to fetch decomposition data for")
        return

    all_ids_str = ", ".join(f"'{eid}'" for eid in new_entity_ids)

    # hasComponent relationships (Complex → components) with stoichiometry
    query_components = f"""
        MATCH (parent)-[rel:hasComponent]->(child)
        WHERE parent.stId IN [{all_ids_str}]
        RETURN parent.stId as parent_id, child.stId as child_id, rel.stoichiometry as stoichiometry
    """
    comp_results = graph.run(query_components).data()
    for row in comp_results:
        pid = row["parent_id"]
        cid = row["child_id"]
        if pid not in _components_cache:
            _components_cache[pid] = {}
        _components_cache[pid][cid] = row.get("stoichiometry") or 1

    # hasCandidate|hasMember relationships (EntitySet → members)
    query_members = f"""
        MATCH (parent)-[:hasCandidate|hasMember]->(child)
        WHERE parent.stId IN [{all_ids_str}]
        RETURN parent.stId as parent_id, child.stId as child_id
    """
    member_results = graph.run(query_members).data()
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
        WHERE pathway.stId = '%s'
        OPTIONAL MATCH (r1)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway)
        WHERE pathway.stId = '%s'
        RETURN r1.stId AS preceding_reaction_id,
               r2.stId AS following_reaction_id,
               CASE WHEN r2 IS NULL THEN 'No Preceding Event' ELSE 'Has Preceding Event' END AS event_status
        """ % (pathway_id, pathway_id)

    try:
        result = graph.run(query).data()
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
            f"Failed to connect to Neo4j database at {uri}. "
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
        result = graph.run(query).data()
        logger.info(f"Found {len(result)} top-level pathways")
        return result
    except Exception as e:
        logger.error("Error in get_top_level_pathways", exc_info=True)
        raise ConnectionError(
            f"Failed to query top-level pathways from Neo4j at {uri}. "
            f"Ensure Neo4j is running and accessible. Original error: {str(e)}"
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
    query: str = f"""
        MATCH (p:Pathway)
        WHERE p.stId = '{pathway_id}'
        RETURN p.displayName AS name
    """

    try:
        result = graph.run(query).data()
        if not result:
            raise ValueError(f"Pathway with ID {pathway_id} not found")
        return result[0]["name"]
    except ValueError:
        raise
    except Exception as e:
        logger.error(f"Error in get_pathway_name for {pathway_id}", exc_info=True)
        raise ConnectionError(
            f"Failed to query pathway name from Neo4j at {uri}. "
            f"Original error: {str(e)}"
        ) from e


def get_labels(entity_id: str) -> List[str]:
    if entity_id in _labels_cache:
        return _labels_cache[entity_id]

    query_get_labels_template: str = """
       MATCH (e)
          WHERE e.stId = '%s'
       RETURN labels(e) AS labels
       """
    query: str = query_get_labels_template % entity_id

    try:
        result = graph.run(query).data()[0]["labels"]
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

    query_get_components_template: str = """
       MATCH (entity)-[rel:hasComponent]->(component)
           WHERE entity.stId = '%s'
       RETURN component.stId AS component_id, rel.stoichiometry AS stoichiometry
       """
    query: str = query_get_components_template % entity_id

    try:
        data = graph.run(query).data()
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

    query_get_members_template: str = """
        MATCH (entity)-[:hasCandidate|hasMember]->(member)
            WHERE entity.stId = '%s'
        RETURN collect(member.stId) as member_ids
        """
    query: str = query_get_members_template % entity_id

    try:
        result = set(graph.run(query).data()[0]["member_ids"])
        _members_cache[entity_id] = result
        return result
    except Exception:
        logger.error("Error in get_set_members", exc_info=True)
        raise


def get_reactions(pathway_id: str, taxon_id: str) -> List[str]:
    query_reaction_template: str = """
        MATCH (reaction)<-[:hasEvent*]-(pathway:Pathway)-[:species]->(species:Species)
             WHERE (reaction:Reaction OR reaction:ReactionLikeEvent)
                   AND pathway.stId='%s' AND species.taxId="%s"
        RETURN COLLECT(reaction.stId) AS reaction_ids
    """
    query: str = query_reaction_template % (pathway_id, taxon_id)

    try:
        return graph.run(query).data()[0]["reaction_ids"]
    except Exception:
        logger.error("Error in get_reactions", exc_info=True)
        raise


def get_reaction_input_output_ids(reaction_id: str, input_or_output: str) -> Set[str]:
    if reaction_id in _reaction_io_cache:
        return _reaction_io_cache[reaction_id].get(input_or_output, set())

    query_template: str = """
       MATCH (reaction)-[:%s]-(io)
           WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) AND reaction.stId='%s'
       RETURN COLLECT(io.stId) AS io_ids
    """
    relation_type: str = "input" if input_or_output == "input" else "output"
    query: str = query_template % (relation_type, reaction_id)

    try:
        return set(graph.run(query).data()[0]["io_ids"])
    except Exception:
        logger.error("Error in get_reaction_input_output_ids", exc_info=True)
        raise


def get_reference_entity_id(entity_id: str) -> Union[str, None]:
    if entity_id in _reference_entity_cache:
        return _reference_entity_cache[entity_id]
    if _prefetch_done:
        return None  # Not in bulk results means no HGNC reference

    query_template: str = """
        MATCH (reference_database:ReferenceDatabase)<-[:referenceDatabase]-(reference_entity_gene:ReferenceEntity)<-[:referenceGene]-(reference_entity:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE reference_database.displayName = "HGNC"
            AND pe.stId = '%s'
        RETURN reference_entity.stId as id
    """  # noqa
    query: str = query_template % entity_id

    try:
        data = graph.run(query).data()
        if len(data) == 0:
            _reference_entity_cache[entity_id] = None
            return None
        result = data[0]["id"]
        _reference_entity_cache[entity_id] = result
        return result
    except Exception:
        logger.error("Error in get_reference_entity_id", exc_info=True)
        raise


def contains_reference_gene_product_molecule_or_isoform(entity_id: str) -> bool:
    query_template = """
        MATCH (es:EntitySet)-[:hasCandidate|hasMember]->(pe:PhysicalEntity)
        WHERE es.stId = '%s'
            AND pe.referenceType IN ["ReferenceGeneProduct", "ReferenceIsoform", "ReferenceMolecule"]
        RETURN COUNT(pe) > 0 AS contains_reference
        """
    query = query_template % entity_id

    try:
        result = graph.run(query).data()
        return result[0]["contains_reference"]
    except Exception as e:
        logger.error(
            "Error in contains_reference_gene_product_molecule_or_isoform",
            exc_info=True,
        )
        raise e
