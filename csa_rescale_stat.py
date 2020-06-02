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
    '''Fetch and concatenate all csv files into one
    :param path_results: path to folder containing csv files for statistics
    '''
    files = []
    for file in os.listdir(path_results):
        if file.endswith(".csv"):
            files.append(os.path.join(path_results,file))
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[3].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")

def get_plot(df_p1):
    '''plot percentage difference between measured simulated atrophy and ground truth atrophy
    for different vertebrae levels
    :param df_p1: dataframe for first plot
    '''
    fig= plt.figure()
    df_p1.groupby('rescale')[['perc_diff_C2_C3','perc_diff_C2_C4','perc_diff_C2_C5']].mean().plot(kind='bar')
    plt.xlabel('rescaling factor')
    plt.ylabel('error in %')
    plt.title('error in function of rescaling factor')
    plt.grid()
    plt.tight_layout()
    plt.savefig("err_plot.jpg")

#plot error in function of simulated atrophy
def get_plot2(df_p2):
    '''plot percentage difference between measured simulated atrophy and ground truth atrophy
    between C2 and C5
    :param df_p2: dataframe for second plot
    '''
    fig2 = plt.figure()
    df_p2.groupby('rescale').mean().perc_diff_C2_C5.plot(kind='bar')
    plt.xlabel('rescaling factor')
    plt.title('error in function of rescaling factor for vertebrae C2 to C5')
    plt.ylabel('error in %')
    plt.grid()
    plt.tight_layout()
    plt.savefig("csa_plot.jpg")


#plot error in function of simulated atrophy
def get_plot3(df_p3):
    '''plot CSA boxplot over different rescalings
    :param df_p3: dataframe for third plot
    '''
    fig3 = plt.figure()
    df_p3.boxplot(column=['CSA'], by='rescale')
    plt.ylabel('CSA in mm^2')
    plt.savefig("csa_boxplot.jpg")


def get_plot4(df_p4):
    '''plot percentage of error boxplot over different rescalings
    :param df_p4: dataframe for fourth plot
    '''
    fig4 = plt.figure()
    df_p4.boxplot(column=['perc_diff_C2_C5'], by='rescale')
    plt.ylabel('error in %')
    plt.savefig("err_boxplot.jpg")


def get_plot_sample(z, z_power, std, mean_CSA):
    '''plot minimum number of patients required to detect an atrophy of X
    :param z: z score for X % uncertainty
    :param z_power: z score for X % Power
    :param std: standard deviation for ground truth
    :param mean_CSA: mean value of CSA from which percentage atrophy is calculated (default = 80)
    '''
    fig6, ax = plt.subplots()
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


def std(CSA_group, atrophy, Vert):
    '''
    Computes standard deviation of subject mean CSA for each rescaling and different vertebrae levels
    '''
    print("\n====================std==========================\n")
    # TODO: normalization of CSA for intersubject studies
    std_arr = []
    for r in atrophy:
        for n in Vert[1:len(Vert)]:
            std = np.concatenate(list(CSA_group.get_group((r,i)) for i in range(2,n+1)), axis=None).std()
            print('CSA std on ',r,' rescaled image C2/C'+str(n)+' is ',round(std, 3),' mm^2 ')
        print('\n')


def ttesst(CSA_group, atrophy, Vert):
    '''
    Computes t test to measure the significance of the difference between rescaled CSA and original CSA * rescaling factor
    for each rescaling and different vertebrae levels
    '''
    print("\n====================ttest==========================\n")
    for r in atrophy:
        for n in Vert[1:len(Vert)]:
            ttest,pvalue = stats.ttest_ind(np.array(list(CSA_group.get_group((1,i)) for i in (2,n))), np.array(list(CSA_group.get_group((1,i)) for i in (2,n)))*(r**2))
            print('p-value for ',r,' rescaled image C2/C'+str(n)+' is ',str(pvalue[0]),' ')
        print('\n')

def sample_size(std):
    '''
    calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    sample size with certainty 95% z(0.05/2)=1.96, power 80% z(Power = 80%)=0.84, ratio patients/control 1:1
    and with the assumption both samples have same std
    '''
    print("\n====================size==========================\n")
    num_n = ((1.96+0.84)**2)*((2*std[-1])**2)
    deno_n = (0.1*80)**2
    n = ceil(num_n/deno_n)
    print('with 80% power, at 5% significance, ratio 1:1 (patients/controls):')
    print('minimum sample size to detect mean 10% atrophy: ',n )


def main():
    '''
    main function: gathers stats and calls plots
    '''
    #read data
    data = pd.read_csv("csa.csv",decimal=".")
    data2 = {'Filename':data['Filename'],
             'VertLevel':data['VertLevel'],
             'CSA':data['MEAN(area)'],
             'rescale':data['rescale'],}
    df = pd.DataFrame(data2)
    # create pandas group by
    df['Filename'] = list((os.path.basename(path).split('.')[0].split('_')[0]) for path in data['Filename'])
    df1 = df.copy()
    df2 = df.copy()
    df3 = df.copy()
    df4 = df.copy()

    df_gt2 = pd.DataFrame()


    # add NaN column to dataframe to insert GT values
    df['gt_CSA_C2_C5'] = np.nan
    df['gt_CSA_C2_C4'] = np.nan
    df['gt_CSA_C2_C3'] = np.nan

    # iterate throuh different vertebrae levels
    df_a = df.groupby(['rescale','Filename']).mean()
    n = []
    for i in range(5,2,-1):
        df_gt2 = pd.DataFrame()

        # get GT values
        if i==5:
            group_CSA_gt = df1.groupby('rescale').get_group(1).set_index('VertLevel').groupby('Filename').mean().CSA
        else:
            n.append(i+1)
            group_CSA_gt = df1.groupby('rescale').get_group(1).set_index('VertLevel').drop(index=n).groupby('Filename').mean().CSA

        # iterate through rescale groupby
        for name, group in df.groupby('rescale'):
            atrophy = group['rescale'].values
            # mean values from same subject  different vertebrae
            group2 = group.groupby('Filename').mean().reset_index()
            # to locate easily subjects put Filename as index
            group3 = group2.set_index('Filename')
            # iterate through dataframe subjects
            for subjectj in set(group2['Filename'].values):
                # if dataframe subject exist in GT (rescale = 1)
                if subjectj in group_CSA_gt.index.values:
                    group3.at[subjectj,'gt_CSA_C2_C'+str(i)] = group_CSA_gt.loc[subjectj] * (atrophy[0] ** 2)
            #print(group3['gt_CSA_C2_C'+str(i-1)])
            df_gt = df.groupby('rescale').get_group(atrophy[0]).groupby('Filename').mean()
            df_gt['gt_CSA_C2_C'+str(i)] = (group3['gt_CSA_C2_C'+str(i)].values)
            df_gt2 = pd.concat([df_gt2, df_gt])
        df_a['gt_CSA_C2_C'+str(i)] = df_gt2['gt_CSA_C2_C'+str(i)].values

    # add CSA vqlues for different vertebrae levels
    df_a['CSA_C2_C4'] = df2.set_index('VertLevel').drop(index=[5]).groupby(['rescale','Filename']).mean().values
    df_a['CSA_C2_C3'] = df3.set_index('VertLevel').drop(index=[4,5]).groupby(['rescale','Filename']).mean().values

    # difference between mean CSA and GT CSA for differennt vertebrae levels
    df_a['diff_C2_C5'] = df_a['CSA'].sub(df_a['gt_CSA_C2_C5']).abs()
    df_a['diff_C2_C4'] = df_a['CSA_C2_C4'].sub(df_a['gt_CSA_C2_C4']).abs()
    df_a['diff_C2_C3'] = df_a['CSA_C2_C3'].sub(df_a['gt_CSA_C2_C3']).abs()


    df_a['perc_diff_C2_C5'] = 100 * df_a['diff_C2_C5'].div(df_a['gt_CSA_C2_C5'])
    df_a['perc_diff_C2_C4'] = 100 * df_a['diff_C2_C4'].div(df_a['gt_CSA_C2_C4'])
    df_a['perc_diff_C2_C3'] = 100 * df_a['diff_C2_C3'].div(df_a['gt_CSA_C2_C3'])

    CSA_group = df4.groupby(['rescale','VertLevel'])['CSA']

    #ground truth atrophy
    atrophy = sorted(set(df4['rescale'].values))
    # get vertebrae levels
    Vert = sorted(set(df['VertLevel'].values))
    std(CSA_group, atrophy, Vert)
    ttesst(CSA_group, atrophy, Vert)

    # plot graph if verbose is present
    if arguments.v is not None:
        df_p1 = df_a.copy()
        get_plot(df_p1)
        df_p2 = df_a.copy()
        get_plot2(df_p2)
        df_p3 = df_a.copy()
        get_plot3(df_p3)
        df_p4 = df_a.copy()
        get_plot4(df_p4)
        #get_plot5(df5)

        get_plot_sample(1.96,(0.84, 1.282), std, 80)
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
