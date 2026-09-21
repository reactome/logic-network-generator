#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

# Determinism: the generator's output depends on set/dict iteration order, which
# Python randomizes per-process via hash seeding — so regenerating a pathway
# yields structurally different networks run-to-run (~5.8% of edges for TP53).
# Pin the hash seed so generation is reproducible. PYTHONHASHSEED must be set
# before the interpreter starts, so re-exec once if it isn't already fixed.
# Override with LNG_ALLOW_NONDETERMINISM=1. See reactome/logic-network-generator#42.
if os.environ.get("PYTHONHASHSEED") != "0" and os.environ.get("LNG_ALLOW_NONDETERMINISM") != "1":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

from typing import List, Tuple

import pandas as pd
from dotenv import dotenv_values, load_dotenv

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.argument_parser import configure_logging, logger, parse_args
from src.logic_network_generator import _reject_removed_env
from src.pathway_generator import generate_pathway_file
from src.neo4j_connector import get_top_level_pathways, get_pathway_name


def canonical_pathway_id(raw: str) -> str:
    """Accept 69620 or R-HSA-69620; return the stable id the queries need.

    The checked-in pathways.tsv carried bare numerics, and every one of them
    failed with "No reactions found ... Verify the pathway exists in Reactome
    database and Neo4j is running" -- a message that sends you to look at Neo4j
    when the id was simply the wrong shape. Normalising here makes that class of
    failure impossible rather than fixing one file.
    """
    pid = str(raw).strip()
    return f"R-HSA-{pid}" if pid.isdigit() else pid


def main() -> None:
    # Fail before ANY work if a removed flag is set. The per-pathway
    # guards inside the generator would each raise, but the loop below catches
    # every exception and continues, so without this check a stale flag would
    # produce "93 failed" buried in the log, leave the PREVIOUS run's
    # logic_network.csv files untouched on disk, and still exit 0 -- and the
    # benchmark that reads that directory would silently score the old catalog.
    _reject_removed_env()

    dotenv_path = os.path.join(os.path.dirname(__file__), "..", ".env")
    # load_dotenv populates os.environ (without clobbering already-set vars) so
    # NEO4J_URL/USER/PASSWORD in .env actually reach neo4j_connector.get_graph(),
    # which reads them via os.getenv. dotenv_values alone only builds a local
    # dict and would leave the connector on its hardcoded defaults.
    load_dotenv(dotenv_path)
    env_vars = dotenv_values(dotenv_path)
    args = parse_args()
    configure_logging(args.debug, args.verbose)

    output_dir = args.output_dir

    # Determine pathway source
    pathway_list_file = (
        args.pathway_list
        if args.pathway_list
        else env_vars.get("PATHWAY_LIST_FILE", None)
    )

    # Validate inputs
    if pathway_list_file:
        if not os.path.exists(pathway_list_file):
            logger.error(f"Pathway list file '{pathway_list_file}' does not exist.")
            sys.exit(1)
        elif not os.access(pathway_list_file, os.R_OK):
            logger.error(f"Pathway list file '{pathway_list_file}' is not readable.")
            sys.exit(1)
    elif not args.pathway_list and not args.pathway_id and not args.top_level_pathways:
        logger.error(
            "One of the following is required: '--pathway-id', '--pathway-list', '--top-level-pathways', or 'PATHWAY_LIST_FILE' environment variable."
        )
        sys.exit(1)

    pathway_list: List[Tuple[str, str]] = []

    if args.top_level_pathways:
        # Fetch all top-level pathways from the database
        logger.info("Fetching all top-level pathways from Reactome database...")
        try:
            top_level = get_top_level_pathways()
            pathway_list = [(p["stId"], p["name"]) for p in top_level]
            logger.info(f"Found {len(pathway_list)} top-level pathways")
        except Exception as e:
            logger.error(f"Error fetching top-level pathways: {e}")
            sys.exit(1)
    elif args.pathway_id:
        # Single pathway by ID - fetch name from database
        pathway_id = args.pathway_id
        try:
            pathway_name = get_pathway_name(pathway_id)
            logger.info(f"Found pathway: {pathway_name} (stId: {pathway_id})")
        except ValueError:
            logger.error(f"Pathway with ID {pathway_id} not found in database")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Error fetching pathway name: {e}")
            sys.exit(1)
        pathway_list = [(pathway_id, pathway_name)]
    elif pathway_list_file:
        try:
            pathways_df: pd.DataFrame = pd.read_csv(pathway_list_file, sep="\t")
            pathway_list = [(canonical_pathway_id(i), n)
                            for i, n in zip(pathways_df["id"].astype(str),
                                            pathways_df["pathway_name"])]
        except Exception as e:
            logger.error(f"Error reading pathway list file: {e}")
            sys.exit(1)

    logger.info(f"Processing {len(pathway_list)} pathway(s)")
    logger.info(f"Output directory: {output_dir}")

    successful = 0
    failed = 0

    for pathway_id, pathway_name in pathway_list:
        try:
            generate_pathway_file(pathway_id, pathway_name, output_dir)
            successful += 1
        except Exception as e:
            logger.error(f"Failed to process pathway {pathway_id} ({pathway_name}): {e}")
            failed += 1
            continue

    logger.info(f"Completed: {successful} successful, {failed} failed")
    if failed:
        # A partial catalog is not a catalog: downstream benchmarks glob the
        # output directory and would score whatever mix of new and stale
        # pathways happens to be there. Exit non-zero so a shell `&&` chain,
        # a Makefile or CI stops instead of proceeding.
        logger.error(f"{failed} pathway(s) failed; exiting non-zero")
        sys.exit(1)


if __name__ == "__main__":
    main()
