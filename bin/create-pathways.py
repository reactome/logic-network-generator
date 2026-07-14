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
from dotenv import dotenv_values

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.argument_parser import configure_logging, logger, parse_args
from src.pathway_generator import generate_pathway_file
from src.neo4j_connector import get_top_level_pathways, get_pathway_name


def main() -> None:
    dotenv_path = os.path.join(os.path.dirname(__file__), "..", ".env")
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
            return
        elif not os.access(pathway_list_file, os.R_OK):
            logger.error(f"Pathway list file '{pathway_list_file}' is not readable.")
            return
    elif not args.pathway_list and not args.pathway_id and not args.top_level_pathways:
        logger.error(
            "One of the following is required: '--pathway-id', '--pathway-list', '--top-level-pathways', or 'PATHWAY_LIST_FILE' environment variable."
        )
        return

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
            return
    elif args.pathway_id:
        # Single pathway by ID - fetch name from database
        pathway_id = args.pathway_id
        try:
            pathway_name = get_pathway_name(pathway_id)
            logger.info(f"Found pathway: {pathway_name} (stId: {pathway_id})")
        except ValueError:
            logger.error(f"Pathway with ID {pathway_id} not found in database")
            return
        except Exception as e:
            logger.error(f"Error fetching pathway name: {e}")
            return
        pathway_list = [(pathway_id, pathway_name)]
    elif pathway_list_file:
        try:
            pathways_df: pd.DataFrame = pd.read_csv(pathway_list_file, sep="\t")
            pathway_list = list(zip(pathways_df["id"].astype(str), pathways_df["pathway_name"]))
        except Exception as e:
            logger.error(f"Error reading pathway list file: {e}")
            return

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


if __name__ == "__main__":
    main()
