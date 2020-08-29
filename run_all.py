# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Break down OpenMP jobs across sub-datasets
# example: python run_all.py -config config_sct_run_batch.yml
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
        '-job-template',
        help="""Path to sbatch config file containing sbatch options preceded of #SBATCH. Example: 
             #SBATCH --account=def-jcohen
             #SBATCH --time=0-08:00        # time (DD-HH:MM)
             #SBATCH --nodes=1
             #SBATCH --cpus-per-task=32    # number of OpenMP processes
             #SBATCH --mem=128G""",
    )
    return parser


def yaml_parser(config_file):
    """parse config_script.yml file containing pipeline's parameters"""
    with open(config_file, 'r') as config_var:
        config_param = yaml.safe_load(config_var)
    return config_param


# text for shell script
def bash_text(config_file, sublist, log_filename, job_template):
    bash_job = """#!/bin/sh
{}
cd $SCRATCH
sct_run_batch -config {} -include-list {} -batch-log {}
    """.format(job_template, config_file, str(sublist).replace("[", "").replace("]", "").replace("'", "").replace(",", ""), log_filename)
    return bash_job


def main():
    # Get parser arguments
    parser = get_parser()
    arguments = parser.parse_args()
    config_file = os.path.abspath(os.path.expanduser(arguments.config))

    # Check if sbatch config file was given with given with flag -job-template
    if arguments.job_template is not None:
        path_job_template = os.path.abspath(os.path.expanduser(arguments.job_template))
        job_template = open(path_job_template, 'r').read()
    else:
        job_template = """#SBATCH --account=def-jcohen
#SBATCH --time=0-08:00        # time (DD-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32    # number of OpenMP processes
#SBATCH --mem=128G
"""

    config_param = yaml_parser(config_file)
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
        filename = os.path.abspath(os.path.expanduser('tmp.job_csa_sublist')) + str(i) + ".sh"
        log_filename = os.path.join(os.path.dirname(filename), "log_" + os.path.basename(filename).split(".")[0] + ".txt")
        # create shell script for sbatch
        with open(filename, 'w+') as temp_file:
            # bash content
            temp_file.write(bash_text(config_file, sublist, log_filename, job_template))
            temp_file.close()
        # Run it
        os.system('sbatch {}'.format(filename))


if __name__ == "__main__":
    main()
