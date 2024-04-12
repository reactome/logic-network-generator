#!/usr/bin/python

from py2neo import Graph
import pandas as pd
import pprint
pp = pprint.PrettyPrinter(indent=4)

uri = "bolt://localhost:7687"
graph = Graph(uri, auth=('neo4j', 'test'))

query = """MATCH (d)
  WHERE d.dbId IS NOT NULL
  AND ("Event" IN labels(d) OR "PhysicalEntity" IN labels(d))
WITH d
OPTIONAL MATCH (d)--(species:Species)
WITH d, COLLECT(species.taxId) AS species_tax_ids
WITH d,
  CASE
    WHEN size(species_tax_ids) = 0 THEN TRUE
    WHEN "9606" IN species_tax_ids THEN
      CASE
        WHEN d.isChimeric IS NULL OR d.isChimeric = FALSE THEN TRUE
        ELSE FALSE
      END
    ELSE FALSE
  END AS is_human, species_tax_ids
WHERE is_human = TRUE
WITH d
OPTIONAL MATCH (d)-[:referenceEntity]->(reference_entity:ReferenceEntity)-[:referenceDatabase]->(reference_database:ReferenceDatabase)
RETURN
  d.dbId AS database_identifier,
  CASE
    WHEN "ReactionLikeEvent" IN labels(d) THEN "reaction-like-event"
    WHEN "Complex" IN labels(d) THEN "complex"
    WHEN "Drug" IN labels(d) THEN "drug"
    WHEN "EntitySet" IN labels(d) THEN "set"
    WHEN "Polymer" IN labels(d) THEN "polymer"
    WHEN "OtherEntity" IN labels(d) THEN "other-entity"
    WHEN "Pathway" IN labels(d) THEN "pathway"
    WHEN (reference_entity.databaseName = "UniProt") THEN "protein"
    WHEN (reference_entity.databaseName = "ENSEMBL") THEN
      CASE reference_entity.schemaClass
        WHEN "ReferenceGeneProduct" THEN "protein"
        WHEN "ReferenceRNASequence" THEN "rna"
        WHEN "ReferenceDNASequence" THEN "dna"
      END
    WHEN reference_entity.databaseName = "ChEBI" THEN "small-molecule"
    WHEN reference_entity.databaseName = "miRBase" THEN "miRNA"
    ELSE "N/A"
  END AS node_type,
  CASE
    WHEN d.displayName IS NOT NULL THEN d.displayName
    ELSE "N/A"
  END AS display_name,
  CASE
    WHEN reference_entity.name[0] IS NOT NULL THEN reference_entity.name[0]
    WHEN reference_entity.displayName IS NOT NULL THEN reference_entity.displayName
    ELSE "N/A"
  END AS reference_entity_name,
  CASE
    WHEN reference_entity.identifier IS NOT NULL THEN reference_entity.databaseName + ":" + reference_entity.identifier
    ELSE "N/A"
  END AS reference_entity_identifier,
  d.schemaClass AS instance_class"""

results = graph.run(query).data()
df = pd.DataFrame(results)

df.to_csv("db_id_to_name_mapping.tsv", sep="\t", index=False)
