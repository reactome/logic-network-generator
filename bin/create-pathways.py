#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# flake8: noqa

import os
import sys
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.argument_parser import parse_args, configure_logging, logger
from src.pathway_generator import generate_pathway_file


def main():
    # parse command line arguments
    args = parse_args()

    # configure logging based on debug flag
    configure_logging(args.debug, args.verbose)

    taxon_id = "9606"

    if args.input_file:
        # Read pathways from the input file
        try:
            pathways_df = pd.read_csv(args.input_file, sep='\t')
            pathways = dict(zip(pathways_df['ID'], pathways_df['PathwayName']))
        except Exception as e:
            logger.error(f"Error reading input file: {e}")
            return
    else:
        logger.error("Input file (--input_file) is required.")
        return

    # create a .tsv file for pathways list
    pathways_list_df = pd.DataFrame(list(pathways.items()), columns=['ID', 'PathwayName'])
    pathways_list_df.to_csv(args.output, sep='\t', index=False)

    for pathway_id, pathway_name in pathways.items():
        generate_pathway_file(pathway_id, taxon_id, pathway_name, decompose=args.decompose)


if __name__ == "__main__":
    main()
