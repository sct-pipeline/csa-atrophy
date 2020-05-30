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
        description='Compute statistics based on the csv files containing the CSA metrics:',
        add_help=None,
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py"))

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        default='results',
        help='Input csv file path to results. (e.g. "results")',
    )

    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-v',
        help='Verbose, plotting figures',
        nargs="*"
    )

    optional.add_argument(
        '-h',
        help='Help',
        nargs="*"
    )

    return parser

# Functions
#############################################################################
    # data extraction for pandas
def get_data(path_results):
    files = []
    for file in os.listdir(path_results):
        if file.endswith(".csv"):
            files.append(os.path.join(path_results,file))
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[3].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")

#plot error in function of simulated atrophy
def get_plot(atrophy, diff_arr):
    fig = plt.figure()
    y_pos = np.arange(len(atrophy))
    # plot
    for i in range(len(diff_arr[1])):
        color = ['tab:red', 'tab:green', 'tab:blue']
        label = ('C2-C3','C2-C4','C2-C5')
        plt.bar(0.3*(i-1), np.absolute(diff_arr[0][i]), align='center', width=0.3, color=color[i], label=label[i])
        plt.bar(y_pos[1]+0.3*(i-1), np.absolute(diff_arr[1][i]), align='center', width=0.3, color=color[i])
        plt.bar(y_pos[2]+0.3*(i-1), np.absolute(diff_arr[2][i]), align='center', width=0.3, color=color[i])
        plt.bar(y_pos[3]+0.3*(i-1), np.absolute(diff_arr[3][i]), align='center', width=0.3, color=color[i])
        #TODO automatise the addition for more or less then 4 atrophies
        plt.bar(y_pos[4]+0.3*(i-1), np.absolute(diff_arr[4][i]), align='center', width=0.3, color=color[i])
        plt.bar(y_pos[5]+0.3*(i-1), np.absolute(diff_arr[5][i]), align='center', width=0.3, color=color[i])
    plt.legend()
    plt.xticks(y_pos, atrophy)
    plt.xlabel('rescaling factor')
    plt.title('error in function of rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    fig.savefig("err_plot.jpg")

#plot error in function of simulated atrophy
def get_plot2(df):
    fig = plt.figure()
    # create pandas group by
    r_CSA = df['CSA'].groupby(df['rescale']).mean()
    gt_CSA = (df.sort_values('rescale')['rescale'].unique()**2) * df.groupby('rescale')['CSA'].mean()[1]
    diff_CSA = r_CSA.sub(gt_CSA).abs()
    perc_diff_CSA = 100 * diff_CSA.div(gt_CSA)
    perc_diff_CSA.plot(kind='bar')
    plt.xlabel('rescaling factor')
    plt.title('error in function of rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    plt.savefig("csa_plot.jpg")

#plot error in function of simulated atrophy
def get_plot3(df):
    fig = plt.figure()
    df.boxplot(column=['CSA'], by='rescale')
    plt.ylabel('CSA in mm^2')
    plt.savefig("csa_plot2.jpg")

def get_plot4(df):
    fig = plt.figure()
    r_CSA = df['CSA'].groupby([df['rescale'],df['Filename']]).mean()
    gt_CSA = (df.groupby([df['rescale'],df['Filename']])['rescale'].unique()**2)* df.groupby([df['rescale'],df['Filename']])['CSA'].mean()[1]
    diff_CSA = r_CSA.sub(gt_CSA).abs()
    perc_diff_CSA = 100 * (diff_CSA).div(gt_CSA)
    (perc_diff_CSA.astype('float').groupby('rescale').mean()).plot(kind='bar')
    plt.xlabel('rescaling factor')
    plt.title('error in function of rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    plt.savefig("csa_persub_plot.jpg")

def get_plot5(df):
    fig = plt.figure()
    r_CSA = df['CSA'].groupby([df['rescale'],df['Filename']]).mean()
    print(r_CSA)
    gt_CSA = (df.groupby([df['rescale'],df['Filename']])['rescale'].unique()**2)* df.groupby([df['rescale'],df['Filename']])['CSA'].mean()[1]
    print(gt_CSA)
    diff_CSA = r_CSA.sub(gt_CSA)
    perc_diff_CSA = 100 * (diff_CSA).div(gt_CSA)
    perc_diff_CSA.to_frame(name='CSA').boxplot(column=['CSA'], by='rescale')
    plt.savefig("err_boxplot.jpg")


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
    # TODO: change: this functions should take into account variable mean_CSA, default=80
    def forward(i):
        i2 = i/80*100
        return i2
    def inverse(i):
        return i/100*80
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel('atrophy in %')
    fig.savefig("min_subj.jpg", bbox_inches='tight')



def diff(CSA_group, atrophy, Vert):
    print("====================diff==========================\n")
    # Computes mean of metric over all subjects
    diff_perc_arr = []
    gt_CSA = []
    r_CSA = []
    diff_arr1 = []
    diff_arr2 = []
    gt_arr =[]
    r_arr = []
    gt_arr2 =[]
    r_arr2 = []
    for r in atrophy:
        for n in Vert[1:len(Vert)]:
            # compute mean ground truth CSA for vertebras ranges
            gt_CSA.append(list((CSA_group.get_group((1,i)) for i in range(2, n+1))))
            gt_arr = np.concatenate(gt_CSA, axis=None).mean()*(r**(2))

            # compute mean rescaled CSA for vertebras ranges
            r_CSA.append(list(CSA_group.get_group((r,i)) for i in range(2,n+1)))
            r_arr = np.concatenate(r_CSA, axis=None).mean()

            # compute difference of mean CSA
            diff_perc_arr.append(100*abs(r_arr - gt_arr)/(gt_arr))
            diff_perc = np.array(diff_perc_arr).mean()
            diff_arr1.append(diff_perc)

            # print mean differnce
            print('the difference with ground truth for ',r,' rescaling on C2/C'+str(n)+' is ',round(diff_perc, 3),' %')
            diff_perc_arr = []
            gt_CSA = []
            r_CSA = []
        r_arr2.append(r_arr) # values of CSA for rescaled images
        gt_arr2.append(gt_arr)# values of CSA for original images *(r**(2/3))
        diff_arr2.append(diff_arr1)# values of difference of CSA
        r_arr = []
        gt_arr = []
        diff_arr1=[]
        print('\n')
    return diff_arr2

def std(CSA_group, atrophy, Vert):
    print("\n====================std==========================\n")
    # Computes standard deviation of subject mean CSA for each rescaling
    # TODO: normalization of CSA for intersubject studies
    std_arr = []
    for r in atrophy:
        for n in Vert[1:len(Vert)]:
            std = np.concatenate(list(CSA_group.get_group((r,i)) for i in range(2,n+1)), axis=None).std()
            print('CSA std on ',r,' rescaled image C2/C'+str(n)+' is ',round(std, 3),' mm^2 ')
            #std_arr.append(std)
        #std_arr.append(std_arr)
        print('\n')
        #return std_arr

def ttesst(CSA_group, atrophy, Vert):
    print("\n====================ttest==========================\n")
    # Computes t test to measure the significance of the difference between rescaled CSA and original CSA * rescaling factor
    for r in atrophy:
        for n in Vert[1:len(Vert)]:
            ttest,pvalue = stats.ttest_ind(np.array(list(CSA_group.get_group((1,i)) for i in (2,n))), np.array(list(CSA_group.get_group((1,i)) for i in (2,n)))*(r**2))
            print('p-value for ',r,' rescaled image C2/C'+str(n)+' is ',str(pvalue[0]),' ')
        print('\n')

def sample_size(std):
    print("\n====================size==========================\n")
    # calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    # sample size with certainty 95% z(0.05/2)=1.96, power 0.8 zscore=0.84, ratio patients/control 1:1
    # and with the assumption both samples have same std
    # (temp ref: the best option could be G*Power)
    num_n = ((1.96+0.84)**2)*((2*std[-1])**2)
    deno_n = (0.1*80)**2
    n = ceil(num_n/deno_n)
    print('with 80% power, at 5% significance, ratio 1:1 (patients/controls):')
    print('minimum sample size to detect mean 10% atrophy: ',n )

######################################################################################
def main():
    #read data
    data = pd.read_csv("csa.csv",decimal=".")
    data2 = {'Filename':data['Filename'],
             'VertLevel':data['VertLevel'],
             'CSA':data['MEAN(area)'],
             'rescale':data['rescale'],}
    df = pd.DataFrame(data2)

    # create pandas group by
    df['Filename']=list((os.path.basename(path).split('.')[0].split('_')[0]) for path in data['Filename'])
    CSA_group = df.groupby(['rescale','VertLevel'])['CSA']

    #ground truth atrophy
    atrophy = sorted(set(df['rescale'].values))
    # get vertebrae levels
    Vert = sorted(set(df['VertLevel'].values))
    std(CSA_group, atrophy, Vert)
    ttesst(CSA_group, atrophy, Vert)
    diff_err = diff(CSA_group, atrophy, Vert)

    # plot graph if verbose is present
    if arguments.v is not None:
        get_plot(atrophy, diff_err)
        get_plot2(df)
        get_plot4(df)
        get_plot5(df)
        get_plot3(df)
        # get_plot_sample(1.96,(0.84, 1.282), std, 80)
        print('\nfigures have been ploted in dataset')


# Run
#########################################################################################
if __name__ == "__main__":
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    if arguments.h is None:
        path_results = os.path.join(os.getcwd(),arguments.i)
        get_data(path_results)
        main()
    else:
        parser.print_help()
