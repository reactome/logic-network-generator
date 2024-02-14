from typing import Any, Dict, List, Set

import pandas as pd
from py2neo import Graph

from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


def get_reaction_connections(pathway_id: str) -> pd.DataFrame:
    query: str = (
        """
       MATCH (pathway:Pathway)-[:hasEvent*]->(r1:ReactionLikeEvent)
           WHERE pathway.dbId = %s
       OPTIONAL MATCH (r1:ReactionLikeEvent)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway:Pathway)
           WHERE pathway.dbId = %s
       RETURN r1.dbId AS parent_reaction_id, r2.dbId AS child_reaction_id
    """
        % (pathway_id, pathway_id)
    )

    try:
        df: pd.DataFrame = pd.DataFrame(graph.run(query).data())
        df = df.astype({"parent_reaction_id": "Int64", "child_reaction_id": "Int64"})
        return df
    except Exception:
        logger.error("Error in get_reaction_connections", exc_info=True)
        raise


def get_all_pathways() -> List[Dict[str, Any]]:
    query: str = """
        MATCH (p:Pathway)
        WHERE p.speciesName='Homo sapiens'
        RETURN
            p.stId AS id,
            p.name[0] AS name
        LIMIT 10
        """

    try:
        return graph.run(query).data()
    except Exception:
        logger.error("Error in get_all_pathways", exc_info=True)
        raise


def get_labels(entity_id: int) -> List[str]:
    query_get_labels_template: str = """
       MATCH (e)
          WHERE e.dbId = %s
       RETURN labels(e) AS labels
       """
    query: str = query_get_labels_template % entity_id

    try:
        return graph.run(query).data()[0]["labels"]
    except Exception:
        logger.error("Error in get_labels", exc_info=True)
        raise


def get_complex_components(entity_id: int) -> Set[int]:
    query_get_components_template: str = """
       MATCH (entity)-[:hasComponent]->(component)
           WHERE entity.dbId = %s
       RETURN collect(component.dbId) AS component_ids
       """
    query: str = query_get_components_template % entity_id

    try:
        return set(graph.run(query).data()[0]["component_ids"])
    except Exception:
        logger.error("Error in get_complex_components", exc_info=True)
        raise


def get_set_members(entity_id: int) -> Set[int]:
    query_get_members_template: str = """
        MATCH (entity)-[:hasCandidate|hasMember]->(member)
            WHERE entity.dbId = %s
        RETURN collect(member.dbId) as member_ids
        """
    query: str = query_get_members_template % entity_id

    try:
        return set(graph.run(query).data()[0]["member_ids"])
    except Exception:
        logger.error("Error in get_set_members", exc_info=True)
        raise


def get_reactions(pathway_id: int, taxon_id: str) -> List[int]:
    query_reaction_template: str = """
        MATCH (reaction)<-[:hasEvent*]-(pathway:Pathway)-[:species]->(species:Species)
             WHERE (reaction:Reaction OR reaction:ReactionLikeEvent)
                   AND pathway.dbId=%s AND species.taxId="%s"
        RETURN COLLECT(reaction.dbId) AS reaction_ids
    """
    query: str = query_reaction_template % (pathway_id, taxon_id)

    try:
        return graph.run(query).data()[0]["reaction_ids"]
    except Exception:
        logger.error("Error in get_reactions", exc_info=True)
        raise


def get_reaction_input_output_ids(reaction_id: int, input_or_output: str) -> Set[int]:
    query_template: str = """
       MATCH (reaction)-[:%s]-(io)
           WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) AND reaction.dbId=%s
       RETURN COLLECT(io.dbId) AS io_ids
    """
    relation_type: str = "input" if input_or_output == "input" else "output"
    query: str = query_template % (relation_type, reaction_id)

    try:
        return set(graph.run(query).data()[0]["io_ids"])
    except Exception:
        logger.error("Error in get_reaction_input_output_ids", exc_info=True)
        raise


def do_entities_have_same_reference_entity(entity_id: int, entity_id2: int) -> bool:
    query: str = """
        MATCH (reference_database:ReferenceDatabase)<-[:referenceDatabase]-(reference_entity:ReferenceEntity)<-[:referenceGene]-(reference_entity_2:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE reference_database.displayName = "HGNC"
            AND (pe.dbId = $entity_id1 OR pe.dbId = $entity_id2)
        WITH reference_entity, count(DISTINCT pe) AS num_entities
        WHERE num_entities = 2
        RETURN count(reference_entity) > 0 AS result
        LIMIT 1
    """  # noqa

    try:
        df: pd.DataFrame = pd.DataFrame(
            graph.run(query, entity_id1=entity_id, entity_id2=entity_id2).data()
        )
        return df["result"].iloc[0]
    except Exception:
        logger.error("Error in do_entities_have_same_reference_entity", exc_info=True)
        raise
