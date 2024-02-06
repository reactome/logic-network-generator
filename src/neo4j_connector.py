from py2neo import Graph
import pandas as pd
from src.argument_parser import logger

uri = "bolt://localhost:7687"
graph = Graph(uri, auth=('neo4j', 'test'))


def get_reaction_connections(pathway_id):
    query = """
       MATCH (pathway:Pathway)-[:hasEvent*]->(r1:ReactionLikeEvent)
           WHERE pathway.dbId = %s
       OPTIONAL MATCH (r1:ReactionLikeEvent)<-[:precedingEvent]-(r2:ReactionLikeEvent)<-[:hasEvent*]-(pathway:Pathway)
           WHERE pathway.dbId = %s
       RETURN r1.dbId AS parent_reaction_id, r2.dbId AS child_reaction_id
    """ % (pathway_id, pathway_id)

    try:
        df = pd.DataFrame(graph.run(query).data())
        df = df.astype({'parent_reaction_id': 'Int64', 'child_reaction_id': 'Int64'})
        return df
    except Exception:
        logger.error("Error in get_reaction_connections", exc_info=True)
        raise


def get_all_pathways():
    query = """
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


def get_labels(entity_id):
    query_get_labels_template = """
       MATCH (e)
          WHERE e.dbId = %s
       RETURN labels(e) AS labels
       """
    query = query_get_labels_template % entity_id

    try:
        data = graph.run(query).data()
        return graph.run(query).data()[0]["labels"]
    except Exception:
        logger.error("Error in get_labels", exc_info=True)
        raise


def get_complex_components(entity_id):
    query_get_components_template = """
       MATCH (entity)-[:hasComponent]->(component)
           WHERE entity.dbId = %s
       RETURN collect(component.dbId) AS component_ids
       """
    query = query_get_components_template % entity_id

    try:
        return set(graph.run(query).data()[0]["component_ids"])
    except Exception:
        logger.error("Error in get_complex_components", exc_info=True)
        raise


def get_set_members(entity_id):
    query_get_members_template = """
        MATCH (entity)-[:hasCandidate|hasMember]->(member)
            WHERE entity.dbId = %s
        RETURN collect(member.dbId) as member_ids
        """
    query = query_get_members_template % entity_id

    try:
        return set(graph.run(query).data()[0]["member_ids"])
    except Exception:
        logger.error("Error in get_set_members", exc_info=True)
        raise


def get_reactions(pathway_id, taxon_id):
    query_reaction_template = """
        MATCH (reaction)<-[:hasEvent*]-(pathway:Pathway)-[:species]->(species:Species)
             WHERE (reaction:Reaction OR reaction:ReactionLikeEvent)
                   AND pathway.dbId=%s AND species.taxId="%s"
        RETURN COLLECT(reaction.dbId) AS reaction_ids
    """
    query = query_reaction_template % (pathway_id, taxon_id)

    try:
        return graph.run(query).data()[0]["reaction_ids"]
    except Exception:
        logger.error("Error in get_reactions", exc_info=True)
        raise


def get_reaction_input_output_ids(reaction_id, input_or_output):
    query_template = """
       MATCH (reaction)-[:%s]-(io)
           WHERE (reaction:Reaction OR reaction:ReactionLikeEvent) AND reaction.dbId=%s
       RETURN COLLECT(io.dbId) AS io_ids
    """
    relation_type = "input" if input_or_output == "input" else "output"
    query = query_template % (relation_type, reaction_id)

    try:
        return set(graph.run(query).data()[0]["io_ids"])
    except Exception:
        logger.error("Error in get_reaction_input_output_ids", exc_info=True)
        raise


def get_reference_entities(entity_id):
    query = """
        MATCH (e:Entity {entity_id: %s})-[:HAS_REFERENCE_ENTITY]->(ref:Entity)
        RETURN ref.entity_id AS reference_entity_id
    """ % entity_id

    try:
        df = pd.DataFrame(graph.run(query).data())
        df = df.astype({'reference_entity_id': 'str'})
        return df
    except Exception:
        logger.error("Error in get_reference_entities", exc_info=True)
        raise
