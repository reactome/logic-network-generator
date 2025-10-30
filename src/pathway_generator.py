import os

import pandas as pd

from src.argument_parser import logger
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types
from src.logic_network_generator import create_pathway_logic_network
from src.neo4j_connector import get_reaction_connections
from src.reaction_generator import get_decomposed_uid_mapping


def generate_pathway_file(
    pathway_id: str, taxon_id: str, pathway_name: str, decompose: bool = False
) -> None:
    """Generate pathway logic network file with caching.

    Args:
        pathway_id: Reactome pathway database ID
        taxon_id: Taxonomy ID (currently unused)
        pathway_name: Human-readable pathway name
        decompose: Whether to decompose complexes/sets (default: False)

    Raises:
        ConnectionError: If Neo4j database is not accessible
        ValueError: If pathway data is invalid or pathway not found
        IOError: If cache files cannot be written
    """
    logger.info(f"Generating logic network for pathway {pathway_id}: {pathway_name}")

    # Define filenames for caching
    reaction_connections_file = f"reaction_connections_{pathway_id}.csv"
    decomposed_uid_mapping_file = f"decomposed_uid_mapping_{pathway_id}.csv"
    best_matches_file = f"best_matches_{pathway_id}.csv"

    try:
        # Load or fetch reaction connections
        if os.path.exists(reaction_connections_file):
            logger.info(f"Loading cached reaction connections from {reaction_connections_file}")
            reaction_connections = pd.read_csv(reaction_connections_file)
        else:
            logger.info(f"Fetching reaction connections from Neo4j for pathway {pathway_id}")
            reaction_connections = get_reaction_connections(pathway_id)
            try:
                reaction_connections.to_csv(reaction_connections_file, index=False)
                logger.info(f"Cached reaction connections to {reaction_connections_file}")
            except IOError as e:
                logger.warning(f"Could not cache reaction connections: {e}")
                # Continue without caching

        # Optional: Limit number of reactions for testing
        number_of_reaction_connections: int = -1
        if number_of_reaction_connections > 0:
            reaction_connections = reaction_connections.iloc[
                :number_of_reaction_connections
            ]

        # Load or generate decomposition and best matches
        if os.path.exists(decomposed_uid_mapping_file) and os.path.exists(best_matches_file):
            logger.info(f"Loading cached decomposition from {decomposed_uid_mapping_file}")
            decomposed_uid_mapping = pd.read_csv(
                decomposed_uid_mapping_file, dtype=decomposed_uid_mapping_column_types
            )
            best_matches = pd.read_csv(best_matches_file)
        else:
            logger.info("Decomposing complexes and entity sets...")
            [decomposed_uid_mapping, best_matches_list] = get_decomposed_uid_mapping(
                pathway_id, reaction_connections
            )
            best_matches = pd.DataFrame(
                best_matches_list, columns=["incomming", "outgoing"]
            )

            try:
                decomposed_uid_mapping.to_csv(decomposed_uid_mapping_file, index=False)
                best_matches.to_csv(best_matches_file, index=False)
                logger.info(f"Cached decomposition to {decomposed_uid_mapping_file}")
            except IOError as e:
                logger.warning(f"Could not cache decomposition results: {e}")
                # Continue without caching

        # Generate logic network
        logger.info("Creating pathway logic network...")
        pathway_logic_network = create_pathway_logic_network(
            decomposed_uid_mapping, reaction_connections, best_matches
        )

        # Save logic network
        output_file = f"pathway_logic_network_{pathway_id}.csv"
        try:
            pathway_logic_network.to_csv(output_file, index=False)
            logger.info(f"Successfully generated logic network: {output_file}")
            logger.info(f"Network contains {len(pathway_logic_network)} edges")
        except IOError as e:
            logger.error(f"Failed to write output file {output_file}: {e}")
            raise

    except (ConnectionError, ValueError) as e:
        logger.error(f"Failed to generate pathway {pathway_id}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error generating pathway {pathway_id}", exc_info=True)
        raise RuntimeError(f"Pathway generation failed: {str(e)}") from e
