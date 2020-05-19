#!/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Evaluate the robustness of automated CSA with global rescaled images
#
# example: python csa_rescale_stat.py -i <results>
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
        description='Compute statistics based on the csv file containing the metrics:',
        add_help=None,
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py"))

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        default='results',
        help='Input csv file path to results. (e.g. results)',
    )

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
#############################################################################
    # data extraction for pandas
def get_data(path_results):
    files = []
    for file in os.listdir(path_results):
        if file.endswith(".csv"):
            files.append(os.path.join(path_results,file))
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[2].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")

def get_plot(atrophy, diff_arr):
    fig = plt.figure()
    y_pos = np.arange(len(atrophy))
    # plot
    plt.bar(y_pos, np.absolute(diff_arr), align='center', alpha=0.5)
    plt.xticks(y_pos, atrophy)
    plt.xlabel('rescaling factor')
    plt.title('error in function of rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    fig.savefig("err_plot.jpg")


def get_plot_sample(z, z_power, std, mean_CSA):
    fig, ax = plt.subplots()
    # data for plotting
    n=[]
    for z_p in z_power:
        i = np.arange(1.5, 8.0, 0.05)
        i_perc = np.arange(1.5, 8.0, 0.05)
        num_n = ((z+z_p)**2)*((2*std)**2)
        n.append(num_n/((i)**2))
        # plot
    ax.plot(i, n[0], label=('80% power'))
    ax.plot(i, n[1], label=('90% power'))
    ax.set_ylabel('minimum number of participants')
    ax.set_xlabel('atrophy in mm^2')
    ax.set_title('minimum number of subjects to detect an atrophy ')
    ax.legend()
    ax.grid()
    # TODO: change this functions should take into account variable mean_CSA, default=80
    def forward(i):
        i2 = i/80*100
        return i2
    def inverse(i):
        return i/100*80
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel('atrophy in %')
    fig.savefig("min_subj.jpg", bbox_inches='tight')



# Main
########################################################################
def main():
    #read data
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

    # Computes standard deviation of subject mean CSA for each rescaling
    # TODO: normalization of CSA for intersubject studies
    print("\n====================std==========================\n")
    std_arr = []
    for r in atrophy:
        std = data.groupby('rescale')['MEAN(area)'].get_group(r).std()
        print('CSA std on ',r,' rescaled image is ',float(std),' mm^2 ')

    # Computes t test to measure the significance of the difference between rescaled CSA and original CSA * rescaling factor
    print("\n====================ttest==========================\n")
    for r in atrophy:
        ttest,pvalue = stats.ttest_ind(data.groupby('rescale')['MEAN(area)'].get_group(r), data.groupby('rescale')['MEAN(area)'].get_group(1)*(r**(2/3)))
        print('p-value for ',r,' rescaled image is ',pvalue,' ')

    # calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    print("\n====================size==========================\n")
    # sample size with certainty 95% z(0.05/2)=1.96, power 0.8 zscore=0.84, ratio patients/control 1:1
    # and with the assumption both samples have same std
    # (temp ref: the best option could be G*Power)
    num_n = ((1.96+0.84)**2)*((2*std)**2)
    deno_n = (0.1*80)**2
    n = ceil(num_n/deno_n)
    print('with 80% power, at 5% significance:')
    print('minimum sample size to detect mean 10% atrophy: ',n )

    # plot graph if verbose is 1
    if arguments.v == 1:
        get_plot(atrophy, diff_arr)
        get_plot_sample(1.96,(0.84, 1.282), std, 80)
        print('\nfigures have been ploted in dataset')


# Run
#########################################################################################
if __name__ == "__main__":
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    path_results = os.path.join(os.getcwd(),arguments.i)
    get_data(path_results)
    main()
