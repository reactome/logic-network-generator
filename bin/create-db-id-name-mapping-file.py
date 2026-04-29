#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create database ID to name mapping file from Reactome Neo4j database.

This script extracts all human Event and PhysicalEntity nodes from the Reactome
database and creates a TSV mapping file containing:
- Database identifier (dbId)
- Node type (reaction-like-event, complex, protein, etc.)
- Display name
- Reference entity name
- Reference entity identifier
- Instance class

The mapping file is useful for converting Reactome database IDs to human-readable
names in downstream analysis.
"""

import argparse
import os
import sys
from typing import List, Dict, Any, Optional, Tuple

import pandas as pd
from py2neo import Graph
from py2neo.errors import ConnectionUnavailable

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.argument_parser import configure_logging, logger


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Create database ID to name mapping file from Reactome database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create mapping with default settings (no authentication)
  %(prog)s

  # Specify custom output file
  %(prog)s --output my_mapping.tsv

  # Use custom Neo4j connection
  %(prog)s --uri bolt://myserver:7687

  # Use authentication if required
  %(prog)s --username neo4j --password mypassword

  # Include all species (not just human)
  %(prog)s --all-species

  # Enable debug logging
  %(prog)s --debug
"""
    )

    parser.add_argument(
        "--output", "-o",
        default="db_id_to_name_mapping.tsv",
        help="Output TSV file path (default: db_id_to_name_mapping.tsv)"
    )

    parser.add_argument(
        "--uri",
        default="bolt://localhost:7687",
        help="Neo4j database URI (default: bolt://localhost:7687)"
    )

    parser.add_argument(
        "--username",
        default=None,
        help="Neo4j username (optional, only if authentication is enabled)"
    )

    parser.add_argument(
        "--password",
        default=None,
        help="Neo4j password (optional, only if authentication is enabled)"
    )

    parser.add_argument(
        "--all-species",
        action="store_true",
        help="Include all species (default: human only, taxId 9606)"
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )

    return parser.parse_args()


def build_query(all_species: bool = False) -> str:
    """Build the Cypher query for extracting database ID to name mappings.

    Args:
        all_species: If True, include all species; if False, only human (taxId 9606)

    Returns:
        Cypher query string
    """
    species_filter = ""
    if not all_species:
        species_filter = """
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
"""

    query = f"""MATCH (d)
  WHERE d.dbId IS NOT NULL
  AND ("Event" IN labels(d) OR "PhysicalEntity" IN labels(d))
{species_filter}
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

    return query


def fetch_mapping_data(
    graph: Graph,
    all_species: bool = False
) -> pd.DataFrame:
    """Fetch database ID to name mapping data from Neo4j.

    Args:
        graph: py2neo Graph instance connected to Neo4j
        all_species: If True, include all species; if False, only human

    Returns:
        DataFrame with mapping data

    Raises:
        ConnectionUnavailable: If Neo4j database is not accessible
        ValueError: If no data is returned from the query
    """
    logger.info("Building Cypher query...")
    query = build_query(all_species)

    logger.info("Executing query against Neo4j database...")
    logger.info("This may take several minutes for large databases...")

    try:
        results: List[Dict[str, Any]] = graph.run(query).data()
    except Exception as e:
        raise ConnectionUnavailable(
            f"Failed to execute query against Neo4j database. "
            f"Ensure Neo4j is running and accessible. Error: {str(e)}"
        ) from e

    if not results:
        raise ValueError(
            "Query returned no results. This may indicate:\n"
            "  1. The database is empty\n"
            "  2. No human entities exist (if using --all-species, check database content)\n"
            "  3. The database schema has changed"
        )

    logger.info(f"Retrieved {len(results)} entities from database")

    df = pd.DataFrame(results)

    # Validate DataFrame structure
    expected_columns = [
        "database_identifier",
        "node_type",
        "display_name",
        "reference_entity_name",
        "reference_entity_identifier",
        "instance_class"
    ]

    missing_columns = set(expected_columns) - set(df.columns)
    if missing_columns:
        raise ValueError(
            f"Query results missing expected columns: {missing_columns}"
        )

    return df


def save_mapping_file(df: pd.DataFrame, output_path: str) -> None:
    """Save mapping DataFrame to TSV file.

    Args:
        df: DataFrame to save
        output_path: Path to output TSV file

    Raises:
        IOError: If file cannot be written
    """
    logger.info(f"Writing mapping file to {output_path}...")

    try:
        df.to_csv(output_path, sep="\t", index=False)
    except IOError as e:
        raise IOError(
            f"Failed to write output file {output_path}. "
            f"Check permissions and disk space. Error: {str(e)}"
        ) from e

    logger.info(f"Successfully created mapping file: {output_path}")
    logger.info(f"File contains {len(df)} mappings")

    # Print statistics
    logger.info("\nMapping Statistics:")
    logger.info(f"  Total entities: {len(df)}")

    node_type_counts = df["node_type"].value_counts()
    logger.info("  Node types:")
    for node_type, count in node_type_counts.items():
        logger.info(f"    - {node_type}: {count}")


def main() -> None:
    """Main entry point for the script."""
    args = parse_arguments()
    configure_logging(args.debug, args.verbose)

    logger.info("="*70)
    logger.info("Database ID to Name Mapping Generator")
    logger.info("="*70)

    # Determine authentication
    auth: Optional[Tuple[str, str]] = None
    if args.username and args.password:
        auth = (args.username, args.password)
        logger.info(f"Using authentication (username: {args.username})")
    else:
        logger.info("Connecting without authentication")

    # Connect to Neo4j
    logger.info(f"Connecting to Neo4j at {args.uri}...")

    try:
        graph = Graph(args.uri, auth=auth)
        # Test connection
        graph.run("RETURN 1").data()
        logger.info("Successfully connected to Neo4j")
    except ConnectionUnavailable as e:
        logger.error(f"Failed to connect to Neo4j at {args.uri}")
        logger.error("Troubleshooting:")
        logger.error("  1. Ensure Neo4j is running: docker ps")
        logger.error("  2. Check Neo4j logs for errors")
        logger.error("  3. Verify connection details (URI)")
        if auth:
            logger.error("  4. Verify authentication credentials")
        logger.error(f"\nError: {str(e)}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error connecting to Neo4j: {str(e)}")
        sys.exit(1)

    # Fetch mapping data
    species_scope = "all species" if args.all_species else "human (taxId 9606)"
    logger.info(f"Fetching entities for {species_scope}...")

    try:
        df = fetch_mapping_data(graph, args.all_species)
    except ValueError as e:
        logger.error(f"Data validation error: {str(e)}")
        sys.exit(1)
    except ConnectionUnavailable as e:
        logger.error(f"Connection error: {str(e)}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error fetching data: {str(e)}")
        sys.exit(1)

    # Save mapping file
    try:
        save_mapping_file(df, args.output)
    except IOError as e:
        logger.error(f"File I/O error: {str(e)}")
        sys.exit(1)

    logger.info("\n" + "="*70)
    logger.info("Mapping file created successfully!")
    logger.info("="*70)


if __name__ == "__main__":
    main()
