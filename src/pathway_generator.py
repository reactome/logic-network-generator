import os
import re
from pathlib import Path

import pandas as pd

from src.argument_parser import logger
from src.decomposed_uid_mapping import decomposed_uid_mapping_column_types
from src.logic_network_generator import (
    create_pathway_logic_network,
    export_entity_reaction_proxy_mapping,
    export_node_reaction_context,
    export_nodes,
    export_uuid_to_reactome_mapping,
)
from src.neo4j_connector import get_reaction_connections
from src.reaction_generator import get_decomposed_uid_mapping


def sanitize_filename(name: str) -> str:
    """Sanitize a pathway name for use as a filename/directory name.

    Args:
        name: The pathway name to sanitize

    Returns:
        A sanitized version safe for filesystem use
    """
    # Replace spaces and special characters with underscores
    sanitized = re.sub(r'[^\w\-]', '_', name)
    # Replace multiple underscores with single
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    # Limit length to avoid filesystem issues
    if len(sanitized) > 100:
        sanitized = sanitized[:100]
    return sanitized


def generate_pathway_file(
    pathway_id: str,
    pathway_name: str,
    output_dir: str = "output",
) -> None:
    """Generate pathway logic network file with caching.

    Args:
        pathway_id: Reactome pathway database ID
        pathway_name: Human-readable pathway name
        output_dir: Base output directory (default: "output")

    Raises:
        ConnectionError: If Neo4j database is not accessible
        ValueError: If pathway data is invalid or pathway not found
        IOError: If cache files cannot be written

    Output files are organized as:
        {output_dir}/{pathway_name}_{pathway_id}/
            logic_network.csv           - Main logic network (what users need)
            stid_to_uuid_mapping.csv    - Stable ID to UUID mapping (what users need)
            cache/                      - Intermediate files
    """
    logger.info(f"Generating logic network for pathway {pathway_id}: {pathway_name}")

    # Create pathway-specific output directory
    base_output_dir = Path(output_dir)
    base_output_dir.mkdir(exist_ok=True)

    # Create pathway folder with sanitized name
    folder_name = f"{sanitize_filename(pathway_name)}_{pathway_id}" if pathway_name else f"pathway_{pathway_id}"
    pathway_output_dir = base_output_dir / folder_name
    pathway_output_dir.mkdir(exist_ok=True)

    # Create cache subdirectory for intermediate files
    cache_dir = pathway_output_dir / "cache"
    cache_dir.mkdir(exist_ok=True)

    # Define filenames for caching (in cache subdirectory)
    reaction_connections_file = cache_dir / "reaction_connections.csv"
    decomposed_uid_mapping_file = cache_dir / "decomposed_uid_mapping.csv"
    best_matches_file = cache_dir / "best_matches.csv"

    try:
        # Load or fetch reaction connections
        if os.path.exists(reaction_connections_file):
            logger.info(f"Loading cached reaction connections from {reaction_connections_file}")
            reaction_connections = pd.read_csv(reaction_connections_file, dtype=str)
        else:
            logger.info(f"Fetching reaction connections from Neo4j for pathway {pathway_id}")
            reaction_connections = get_reaction_connections(pathway_id)
            try:
                reaction_connections.to_csv(reaction_connections_file, index=False)
                logger.info(f"Cached reaction connections to {reaction_connections_file}")
            except IOError as e:
                logger.warning(f"Could not cache reaction connections: {e}")
                # Continue without caching

        # Load or generate decomposition and best matches
        if os.path.exists(decomposed_uid_mapping_file) and os.path.exists(best_matches_file):
            logger.info(f"Loading cached decomposition from {decomposed_uid_mapping_file}")
            decomposed_uid_mapping = pd.read_csv(
                decomposed_uid_mapping_file,
                dtype=decomposed_uid_mapping_column_types,  # type: ignore[arg-type]
            )
            best_matches = pd.read_csv(best_matches_file)
        else:
            logger.info("Decomposing complexes and entity sets...")
            [decomposed_uid_mapping, best_matches_list] = get_decomposed_uid_mapping(
                pathway_id, reaction_connections
            )
            best_matches = pd.DataFrame(
                best_matches_list,
                columns=["incomming", "outgoing", "reactome_id"],
            )

            try:
                decomposed_uid_mapping.to_csv(decomposed_uid_mapping_file, index=False)
                best_matches.to_csv(best_matches_file, index=False)
                logger.info(f"Cached decomposition to {decomposed_uid_mapping_file}")
            except IOError as e:
                logger.warning(f"Could not cache decomposition results: {e}")
                # Continue without caching

        # Augment connectivity with curator-drawn diagram flow (product->substrate
        # pairs the diagram links but precedingEvent may omit — esp. old pathways).
        # Only the connectivity/merge step sees this; the matching layer above
        # stays on pure precedingEvent. See reactome/logic-network-generator#39.
        from src.diagram_connectivity import augment_reaction_connections
        connectivity = augment_reaction_connections(pathway_id, reaction_connections)

        # Generate logic network
        logger.info("Creating pathway logic network...")
        result = create_pathway_logic_network(
            decomposed_uid_mapping, connectivity, best_matches
        )

        # Save logic network (main output file users need)
        output_file = pathway_output_dir / "logic_network.csv"
        try:
            result.logic_network.to_csv(output_file, index=False)
            logger.info(f"Successfully generated logic network: {output_file}")
            logger.info(f"Network contains {len(result.logic_network)} edges")
        except IOError as e:
            logger.error(f"Failed to write output file {output_file}: {e}")
            raise

        # Export UUID to Reactome stable ID mapping (main mapping file users need)
        uuid_to_reactome_file = pathway_output_dir / "stid_to_uuid_mapping.csv"
        try:
            export_uuid_to_reactome_mapping(
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                result.catalyst_regulator_map,
                str(uuid_to_reactome_file)
            )
            logger.info(f"Successfully exported stable ID to UUID mapping: {uuid_to_reactome_file}")
        except IOError as e:
            logger.error(f"Failed to write stable ID to UUID mapping file {uuid_to_reactome_file}: {e}")
            # Don't raise - this is supplementary

        # Export entity→reaction proxy mapping. Curated species (often Complexes)
        # that were expanded into virtual variants lose their own stId from the
        # UUID mapping; this file points each such species at the UUIDs of the
        # reactions that produce (or, failing that, consume) it, so consumers can
        # read reaction flux as a proxy for the species' state.
        proxy_mapping_file = pathway_output_dir / "entity_reaction_proxy_mapping.csv"
        try:
            export_entity_reaction_proxy_mapping(
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                pathway_id,
                str(proxy_mapping_file),
            )
            logger.info(f"Successfully exported entity-reaction proxy mapping: {proxy_mapping_file}")
        except Exception as e:
            logger.error(f"Failed to write entity-reaction proxy mapping file {proxy_mapping_file}: {e}")
            # Don't raise - this is supplementary

        # Schema-backed provenance files (schema/logic_network.linkml.yaml):
        # nodes.csv (node_kind, diagram_entity_id, member_leaves, set provenance)
        # and node_reaction_context.csv (node<->reaction location layer).
        try:
            export_nodes(
                result.logic_network,
                result.reaction_id_map,
                result.uuid_mapping,
                str(pathway_output_dir / "nodes.csv"),
            )
            export_node_reaction_context(
                result.entity_uuid_registry,
                result.reaction_id_map,
                result.catalyst_regulator_map,
                str(pathway_output_dir / "node_reaction_context.csv"),
            )
        except Exception as e:
            logger.error(f"Failed to write node provenance files: {e}", exc_info=True)
            # Don't raise - supplementary

        logger.info(f"Output directory: {pathway_output_dir}")

    except (ConnectionError, ValueError) as e:
        logger.error(f"Failed to generate pathway {pathway_id}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error generating pathway {pathway_id}", exc_info=True)
        raise RuntimeError(f"Pathway generation failed: {str(e)}") from e
