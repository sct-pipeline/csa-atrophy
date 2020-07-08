
import yaml
import io
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
        help="variable to ouptut from config.yaml",
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-h',
        help='Help',
        nargs="*"
    )
    return parser

# get parser elements
parser = get_parser()
arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])


# Read YAML file
with open("config.yaml", 'r') as config_var:
    data_loaded = yaml.safe_load(config_var)


print(data_loaded[arguments.o])
