#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from typing import List, Tuple

import pandas as pd
from dotenv import dotenv_values

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.argument_parser import configure_logging, logger, parse_args
from src.pathway_generator import generate_pathway_file


def main() -> None:
    dotenv_path = os.path.join(os.path.dirname(__file__), "..", ".env")
    env_vars = dotenv_values(dotenv_path)
    args = parse_args()
    configure_logging(args.debug, args.verbose)

    pathway_list_file = (
        args.pathway_list
        if args.pathway_list
        else env_vars.get("PATHWAY_LIST_FILE", None)
    )
    if pathway_list_file:
        if not os.path.exists(pathway_list_file):
            logger.error(f"Pathway list file '{pathway_list_file}' does not exist.")
            return
        elif not os.access(pathway_list_file, os.R_OK):
            logger.error(f"Pathway list file '{pathway_list_file}' is not readable.")
            return
    elif not args.pathway_list and not args.pathway_id:
        logger.error(
            "Either '--pathway-list', '--pathway-id', or 'PATHWAY_LIST_FILE' environment variable is required."
        )
        return

    taxon_id = "9606"

    pathway_list: List[Tuple[str, str]] = []

    if args.pathway_id:
        pathway_list = [(args.pathway_id, "")]
    elif pathway_list_file:
        try:
            pathways_df: pd.DataFrame = pd.read_csv(pathway_list_file, sep="\t")
            pathway_list = list(zip(pathways_df["id"], pathways_df["pathway_name"]))
        except Exception as e:
            logger.error(f"Error reading pathway list file: {e}")
            return

    for pathway_id, pathway_name in pathway_list:
        generate_pathway_file(pathway_id, taxon_id, pathway_name)


if __name__ == "__main__":
    main()
