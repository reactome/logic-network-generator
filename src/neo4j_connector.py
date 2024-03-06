from typing import Any, Dict, List, Set, Union

import pandas as pd
from py2neo import Graph  # type: ignore

from src.argument_parser import logger

uri: str = "bolt://localhost:7687"
graph: Graph = Graph(uri, auth=("neo4j", "test"))


def get_reaction_connections(pathway_id: str) -> pd.DataFrame:
    query: str = (
        """
        MATCH (pathway:Pathway)-[:hasEvent*]->(r1:ReactionLikeEvent)
        WHERE pathway.dbId = %s
        OPTIONAL MATCH (r1)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway)
        WHERE pathway.dbId = %s
        RETURN r1.dbId AS preceding_reaction_id,
               r2.dbId AS following_reaction_id,
               CASE WHEN r2 IS NULL THEN 'No Preceding Event' ELSE 'Has Preceding Event' END AS event_status
        """
        % (pathway_id, pathway_id)
    )

    try:
        df: pd.DataFrame = pd.DataFrame(graph.run(query).data())
        df["preceding_reaction_id"] = df["preceding_reaction_id"].astype("Int64")
        df["following_reaction_id"] = df["following_reaction_id"].astype("Int64")
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


def get_reference_entity_id(entity_id: int) -> Union[str, None]:
    query_template: str = """
        MATCH (reference_database:ReferenceDatabase)<-[:referenceDatabase]-(reference_entity_gene:ReferenceEntity)<-[:referenceGene]-(reference_entity:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE reference_database.displayName = "HGNC"
            AND pe.dbId = %s
        RETURN reference_entity.dbId as id
    """  # noqa
    query: str = query_template % entity_id

    try:
        data = graph.run(query).data()
        if len(data) == 0:
            return None
        return data[0]["id"]
    except Exception:
        logger.error("Error in get_reaction_input_output_ids", exc_info=True)
        raise


def contains_reference_gene_product_molecule_or_isoform(entity_id: int) -> bool:
    query_template = """
        MATCH (es:EntitySet)-[:hasCandidate|hasMember]->(pe:PhysicalEntity)
        WHERE es.dbId = %s
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
