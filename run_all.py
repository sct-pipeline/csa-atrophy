# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Break down OpenMP jobs across sub-datasets
# example: python run_all.py
#
#########################################################################################

import os
import argparse

def get_parser(mandatory=None):
    """parser function"""
    parser = argparse.ArgumentParser(
        description="Compute statistics based on the csv files containing the CSA metrics:",
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )
    parser.add_argument(
        '-config',
        required=True,
        help='Path to config file, which contains parameters for the command sct_run_batch.',
    )
    return parser


# Get parser arguments
parser = get_parser()
arguments = parser.parse_args()
config_file = arguments.config
# Get list of subjects in path data
dir_data = 'data-multi-subject'
path_data = os.path.abspath(os.path.expanduser(dir_data))
list = os.listdir(path_data)
list_subjects = [subject for subject in list if "sub" in subject]

# Create X sublists of 32 subjects each
n = 32
sublists = [list_subjects[i:i + n] for i in range(0, len(list_subjects), n)]


# text for shell script
def bash_text(config_file, sublist):
    bash_job = """#!/bin/sh
#SBATCH --account=def-jcohen
#SBATCH --time=0-08:00        # time (DD-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32    # number of OpenMP processes
#SBATCH --mem=128G

cd $SCRATCH
sct_run_batch -config {} -include-list {}
""".format(config_file, str(sublist).replace("[", "").replace("]", "").replace("'", "").replace(",", ""))
    return bash_job


i = 0
# Loop across the sublists
for sublist in sublists:
    i = i + 1
    # Create temporary job shell script: tmp.job_csa_sublist_i.sh
    filename = os.path.abspath(os.path.expanduser('job_csa_sublist')) + str(i) + ".sh"
    # create shell script for sbatch
    with open(filename, 'w+') as temp_file:
        # bash content
        temp_file.write(bash_text(config_file, sublist))
        temp_file.close()
    # Run it
    os.system('sbatch {}'.format(filename))
    # remove tmp file
    # os.remove(filename_config)
    # os.remove(filename)
