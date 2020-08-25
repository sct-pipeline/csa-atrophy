# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Break down OpenMP jobs across sub-datasets
# example: python run_all.py
#
#########################################################################################

import os


# Get list of subjects in path data
dir_data = 'data-multi-subject'
path_data = os.path.abspath(os.path.expanduser(dir_data))
list = os.listdir(path_data)
list_subjects = [subject for subject in list if "sub" in subject]

# Create X sublists of 32 subjects each
n = 32
sublists = [list_subjects[i:i + n] for i in range(0, len(list_subjects), n)]


# text for shell script
def bash_text(config_file):
    bash_job = """#!/bin/sh
#SBATCH --account=def-jcohen
#SBATCH --time=0-08:00        # time (DD-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32    # number of OpenMP processes
#SBATCH --mem=128G

sct_run_batch -config """ + str(config_file) + """
"""
    return bash_job


# text for config_file
def config_text():
    config = """
# config file for sct_run_batch
path_data: data-multi-subject
path_output: csa_atrophy_results_t1
script: process_data.sh
script_args: /home/paul/Github/stats/config_script.yml
jobs: -1
batch_log: sct_run_batch_log.txt
continue_on_error: 1
subject_prefix: sub-
zip: true
include-list: """ + str(sublist) + """
"""
    return config


i = 0
# Loop across the sublists
for sublist in sublists:
    i = i + 1
    # Create temporary job shell script: tmp.job_csa_sublist_i.sh
    filename = os.path.abspath(os.path.expanduser('job_csa_sublist')) + str(i) + ".sh"
    filename_config = os.path.abspath(os.path.expanduser('config_sct_run_batch_')) + str(i) + ".yml"
    # create config file for sct_run_batch
    with open(filename_config, 'w+') as config_file:
        config_file.write(config_text())
        config_file.close()
    # create shell script for sbatch
    with open(filename, 'w+') as temp_file:
        # bash content
        temp_file.write(bash_text(filename_config))
        temp_file.close()
    # Run it
    command = 'sbatch ' + "tmp." + os.path.basename(filename)
    os.system(command)
    # remove tmp file
    os.remove(filename_config)
    os.remove(filename)
