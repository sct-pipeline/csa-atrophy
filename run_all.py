# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Break down OpenMP jobs across sub-datasets
# example: python run_all.py -i config_sct_run_batch
#
#########################################################################################

import os
import argparse
import yaml


def get_parser(mandatory=None):
    """parser function"""
    parser = argparse.ArgumentParser(
        description="Break down OpenMP jobs across sub-datasets",
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )
    parser.add_argument(
        '-config',
        required=True,
        help='Path to config file, which contains parameters for the command sct_run_batch.',
    )
    parser.add_argument(
        '-o-shell',
        help='Path to the temporary batch script given to sbatch. By default basename is given a suffix _i with i the '
             'number of iterations over the batch script and extension .sh. Example: -o-shell job_csa_sublist will '
             'output job_csa_sublist_i.sh',
        default='job_csa_sublist'
    )
    return parser


def yaml_parser(config_file):
    """parse config_script.yml file containing pipeline's parameters"""
    with open(config_file, 'r') as config_var:
        config_param = yaml.safe_load(config_var)
    return config_param


# text for shell script
def bash_text(config_file, sublist, log_filename):
    bash_job = """#!/bin/sh
#SBATCH --account=def-jcohen
#SBATCH --time=0-08:00        # time (DD-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32    # number of OpenMP processes
#SBATCH --mem=128G

cd $SCRATCH
sct_run_batch -config {} -include-list {} -batch-log {}
""".format(config_file, str(sublist).replace("[", "").replace("]", "").replace("'", "").replace(",", ""), log_filename)
    return bash_job


def main():
    # Get parser arguments
    parser = get_parser()
    arguments = parser.parse_args()
    config_file = os.path.abspath(os.path.expanduser(arguments.config))
    config_param = yaml_parser(config_file)

    # get path for log file
    dir_output = config_param['path_output']
    path_output = os.path.abspath(os.path.expanduser(dir_output))
    # Get list of subjects in path data
    dir_data = config_param['path_data']
    path_data = os.path.abspath(os.path.expanduser(dir_data))
    list = os.listdir(path_data)
    list_subjects = [subject for subject in list if "sub" in subject]

    # Create X sublists of 32 subjects each
    n = 32
    sublists = [list_subjects[i:i + n] for i in range(0, len(list_subjects), n)]

    i = 0
    # Loop across the sublists
    for sublist in sublists:
        i = i + 1
        # Create temporary job shell script, default: job_csa_sublist_i.sh
        filename = os.path.abspath(os.path.expanduser(arguments.o_shell)) + str(i) + ".sh"
        log_filename = os.path.join(path_output, "log", "log_" + os.path.basename(filename).split(".")[1] + ".txt")
        # create shell script for sbatch
        with open(filename, 'w+') as temp_file:
            # bash content
            temp_file.write(bash_text(config_file, sublist, log_filename))
            temp_file.close()
        # Run it
        os.system('sbatch {}'.format(filename))


if __name__ == "__main__":
    main()
