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

import pandas as pd
import numpy as np
import sys, os
import argparse
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from math import ceil


# Parser
#########################################################################################

def get_parser():
    parser = argparse.ArgumentParser(
        description='Compute statistics based on the csv file containing metrics:',
        add_help=None,
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py"))

    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-v',
        type=int,
        help='Verbose.',
        required=False,
        default=1,
        choices=(0, 1),)
    optional.add_argument(
        '-h',
        required=False,
        help=' help',)

    return parser

# Functions
#########################################################################################
    # data extraction for pandas
def get_data():
    files = []
    for file in os.listdir("results"):
        if file.endswith(".csv"):
            files.append(os.path.join(os.getcwd(),'results/',file))
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('r_')[1].split('.csv')[0]) for f in files ])
    metrics.to_csv("csa.csv")

def get_plot(atrophy, diff_arr):
    fig = plt.figure()
    y_pos = np.arange(len(atrophy))
    # plot
    plt.bar(y_pos, diff_arr, align='center', alpha=0.5)
    plt.xticks(y_pos, atrophy)
    plt.xlabel('rescaling factor')
    plt.title('error in function of rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    fig.savefig("err_plot.jpg")


def get_plot_sample(z, z_power, std):
    fig = plt.figure()
    # data for plotting
    i = np.arange(1.5, 8.0, 0.05)
    num_n = ((z+z_power)**2)*((2*std)**2)
    n = (num_n/((i)**2))
    # plot
    plt.plot(i, n)
    plt.ylabel('minimum number of participants')
    plt.xlabel('atrophy in mm^2')
    plt.title('minimum number of subjects to detect an atrophy ')
    plt.grid()
    fig.savefig("min_subj.jpg")



# Main
#########################################################################################
def main():
    # get parser elements
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])

    #read data
    get_data()
    data = pd.read_csv("csa.csv",decimal=".")

    #ground truth atrophy
    atrophy = sorted(set(data['rescale']))
    # Computes mean of metric over all subjects
    print("====================diff==========================\n")
    diff_arr = []
    for r in atrophy:
        gt_CSA = float(data.groupby('rescale')['MEAN(area)'].get_group(1).mean())
        r_CSA = data.groupby('rescale')['MEAN(area)'].get_group(r).mean()
        diff_perc = 100*(r_CSA - gt_CSA*(r**(2/3)))/(gt_CSA*(r**(2/3)))
        diff_arr.append(gt_CSA-r_CSA)
        print('the difference with ground truth for ',r,' rescaling is ',diff_perc,' %')

    # Computes standard deviation between theoric and measured CSA rescaled
    print("\n====================std==========================\n")
    std_arr = []
    for r in atrophy:
        std = data.groupby('rescale')['MEAN(area)'].get_group(r).std()
        print('CSA std on ',r,' rescaled image is ',float(std),' mm^2 ')

    # Computes t test for significance of difference
    print("\n====================ttest==========================\n")
    for r in atrophy:
        ttest,pvalue = stats.ttest_ind(data.groupby('rescale')['MEAN(area)'].get_group(r), data.groupby('rescale')['MEAN(area)'].get_group(1)*(r**(2/3)))
        print('p-value for ',r,' rescaled image is ',pvalue,' ')

    # calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    print("\n====================size==========================\n")
    # sample size with certainty 95% z(0.05)=1.645 power 0.8 zscore=1.282
    num_n = ((1.645+1.282)**2)*((2*std)**2)
    deno_n = (0.0178*80)**2
    n = ceil(num_n/deno_n)
    print('with 80% power, at 5% significance:')
    print('minimum sample size to detect annual mean MS atrophy (',deno_n,'mm^2): ',n )

    # plot grah if verbose is 2
    if arguments.v == 1:
        get_plot(atrophy, diff_arr)
        get_plot_sample(1.645, 1.282, std)
        print('\nfigures have been ploted in dataset')


# Run
#########################################################################################
if __name__ == "__main__":
#    sct.init_sct()
    parser = get_parser()
    main()
