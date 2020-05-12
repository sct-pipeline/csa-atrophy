#!/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Evaluate the robustness of automated CSA with global rescaled images
#
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
#
# About the license: see the file LICENSE.TXT
#########################################################################################
from __future__ import division


import sys, os
import argparse
import numpy as np
from numpy import genfromtxt
from scipy import stats
from spinalcordtoolbox.aggregate_slicewise import  save_as_csv # absence causes error ??
import sct_utils as sct
from spinalcordtoolbox.utils import Metavar, SmartFormatter




# Parser
#########################################################################################
def get_parser():
    parser = argparse.ArgumentParser(
        description='Compute statistics based on the csv file containing metrics:',
        add_help=None,
        formatter_class=SmartFormatter,
        prog=os.path.basename(__file__).strip(".py"))

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        '-i',
        required=True,
        metavar=Metavar.file,
        help=' metrics file csv',
        )
    mandatory.add_argument(
        '-r',
        required=True,
        help=' metrics file csv for rescaled',
        metavar=Metavar.file,
        )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-v',
        type=int,
        help='Verbose.',
        required=False,
        default=1,
        choices=(0, 1),)
    return parser

# Functions
#########################################################################################
# Computes mean of metric over all subjects
def mean_metric(metrics, key):
    mean_metric = np.mean(metrics[key])
    return mean_metric

# Computes percentage difference between theoric and measured CSA rescaled
def per_dif_metric(metrics, metrics_r, key, atrophy):
    perc_diff_metric = 100*(abs(np.subtract(metrics_r[key],(atrophy*metrics[key]))))/(atrophy*metrics[key])
    return perc_diff_metric

# Computes standard deviation between theoric and measured CSA rescaled
def std_metric(metrics):
    var = np.mean(abs(metrics - np.mean(metrics))**2)
    std_metric = np.sqrt(var)
    return std_metric

# computes pearson r between original and rescaled image
def pearsonr_metric(metrics, metrics_r, key):
    pearsonr_metric = stats.pearsonr(metrics_r[key], metrics[key])
    return pearsonr_metric

# Computes t test for significance of difference
def ttest_metric(metrics, metrics_r, key, atrophy):
    metric_ttest,metric_pval = stats.ttest_ind(metrics[key]*atrophy,metrics_r[key])
    return metric_ttest,metric_pval



# Main
#########################################################################################
def main():
    # get parser elements
    arguments = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    # fname
    fname = os.path.abspath(arguments.i)
    fname_r = os.path.abspath(arguments.r)

    # data extraction for numpy
    metrics_subj_o = genfromtxt(fname, delimiter=',', names=True, dtype=None, encoding='UTF-8')
    metrics_subj_r = genfromtxt(fname_r, delimiter=',', names=True, dtype=None, encoding='UTF-8')
    # metrics_perslice = genfromtxt(fname_perslice, delimiter=',', names=True, dtype=None, encoding='UTF-8')
    # metrics_perslice_r = genfromtxt(fname_perslice_r, delimiter=',', names=True, dtype=None, encoding='UTF-8')

    verbose = arguments.v
    sct.init_sct(log_level=verbose, update=True)  # Update log level

    # Computes mean of metric over all subjects
    mean_CSA = mean_metric(metrics_subj_o, 'MEANarea')
    mean_CSA_r = mean_metric(metrics_subj_r, 'MEANarea')
    sct.printv('the mean CSA on all original slices is ' + str(mean_CSA) + ' mm^2\n', 'info')
    sct.printv('the mean CSA fon all rescaled slices is ' + str(mean_CSA_r) + ' mm^2\n', 'info')

    # Computes standard deviation between theoric and measured CSA rescaled
    atrophy = (0.9) ** (2 / 3)
    CSA_subj_o_r = metrics_subj_o['MEANarea']*atrophy
    std_CSA = std_metric(CSA_subj_o_r)
    std_CSA_r = std_metric(metrics_subj_r['MEANarea'] )
    sct.printv('the std of CSA is ' + str(std_CSA) + ' mm^2\n', 'info')
    sct.printv('the std of CSA rescaled is ' + str(std_CSA_r) + ' mm^2\n', 'info')

    # computes pearson r between original and rescaled image
    pearsonr_CSA = pearsonr_metric(metrics_subj_o, metrics_subj_r, 'MEANarea')
    sct.printv('the pearson correlation coeff between original and rescaled image is ' + str(pearsonr_CSA) + '\n', 'info')

    # Computes percentage difference between theoric and measured CSA rescaled
    perc_diff_CSA = per_dif_metric(metrics_subj_o, metrics_subj_r, 'MEANarea', atrophy)
    sct.printv('Percentage difference from theoric atrophy is ' + str(perc_diff_CSA) + ' %\n', 'info')

    # Computes t test for significance of difference
    CSA_ttest,CSA_pval = ttest_metric(metrics_subj_o,metrics_subj_r, 'MEANarea', atrophy)
    sct.printv('p-value for the rescale test is ' + str(CSA_pval) + ' \n', 'info')

    verbose = arguments.v
    sct.init_sct(log_level=verbose, update=True)  # Update log level



# Run
#########################################################################################
if __name__ == "__main__":
    sct.init_sct()
    parser = get_parser()
    main()
