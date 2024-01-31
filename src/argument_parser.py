import argparse
import logging

def parse_args():
    parser = argparse.ArgumentParser(description='pathway_creation')
    parser.add_argument('--debug', action='store_true', help='Enable debugging')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    parser.add_argument('--input_file', type=str, help='Input file containing pathway information')

    return parser.parse_args()

# Configure the logging settings
def configure_logging(debug_flag, verbose_flag):
    if verbose_flag:
        log_level = logging.DEBUG
    elif debug_flag:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(filename='debug_log.txt', level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

