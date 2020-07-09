#########################################################################################
#
#
# parse config.yaml file containing pipeline's parameters
#
#
# ---------------------------------------------------------------------------------------
#
# Example: python yaml_parser.py -o <n_transfo>
# About the license: see the file LICENSE
###################################################################

import yaml
import argparse
import os
import sys

def get_parser():
    parser = argparse.ArgumentParser(
        description='output variable from yaml file',
        add_help=None,
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-o",
        required=True,
        help="parameters to output from config.yaml",
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-h',
        help='Help',
        nargs="*"
    )
    return parser


def main():
    """main function, reads and prints yaml config file parameters"""
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    # Read YAML file
    with open("config.yaml", 'r') as config_var:
        data_loaded = yaml.safe_load(config_var)
    # print yaml file parameters to be read by bash script
    print(data_loaded[arguments.o])
