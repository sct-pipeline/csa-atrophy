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
    '''
    parser function
    '''
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
        '-l',
        help='Indicate vertebrae levels of interest \nexample: python csa_rescale_stat.py -i <results> -l 2 3 4 5 ',
        nargs="*",
    )

    optional.add_argument(
        '-h',
        help='Help',
        nargs="*"
    )

    return parser

# Functions
def get_data(path_results):
    '''Fetch and concatenate all csv files into one
    :param path_results: path to folder containing csv files for statistics
    '''
    files = []
    for file in os.listdir(path_results):
        if file.endswith(".csv"):
            files.append(os.path.join(path_results,file))
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[4].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")

def get_plot(df_p1, columns_to_plot):
    '''plot percentage difference between measured simulated atrophy and ground truth atrophy
    for different vertebrae levels
    :param df_p1: dataframe for first plot
    '''
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    print(df_p1)
    df_p1.groupby('Rescale')[columns_to_plot].mean().plot(kind='bar', ax=axes[0], grid=True)
    axes[0].set_title('mean error in function of rescaling factor');
    axes[0].set_ylabel('error in %')
    df_p1.groupby('Rescale')[columns_to_plot].std().plot(kind='bar', ax=axes[1], sharex=True, sharey=True, legend=False)
    axes[1].set_title('STD of error in function of rescaling factor');
    plt.xlabel('rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    plt.tight_layout()
    plt.savefig("err_plot.jpg")


def get_plot3(df_p3):
    '''plot CSA boxplot over different rescalings
    :param df_p3: dataframe for third plot
    '''
    fig3 = plt.figure()
    df_p3.boxplot(column=['CSA_original'], by='Rescale')
    plt.ylabel('CSA in mm^2')
    plt.savefig("csa_boxplot.jpg")


def get_plot4(df_p4, min_max_Vertlevels):
    '''plot percentage of error boxplot over different rescalings
    :param df_p4: dataframe for fourth plot
    '''
    fig4 = plt.figure()
    df_p4.boxplot(column=min_max_Vertlevels, by='Rescale')
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
        i = np.arange(1.5, 8.0, 0.05) # x_axis values ranging from 1.5 to 8.0 mm^2
        num_n = 2*((z+z_p)**2)*(std**2) # numerator of sample size equation
        n.append(num_n/((i)**2))
        # plot
    ax.plot(i, n[0], label=('80% power'))
    ax.plot(i, n[1], label=('90% power'))
    ax.set_ylabel('minimum number of participants per arm')
    ax.set_xlabel('atrophy in mm^2')
    ax.set_title('minimum number of subjects to detect an atrophy with 5% uncertainty\n std = '+str(round(std,2))+'mm², mean_CSA = '+str(mean_CSA)+'mm²')
    ax.legend()
    ax.grid()
    # create global variable for secax functions
    global mean_CSA_sample
    mean_CSA_sample = mean_CSA

    def forward(i):
        i2 = i/mean_CSA_sample*100
        return i2
    def inverse(i):
        return i/100*mean_CSA_sample
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel('atrophy in %')
    plt.savefig("min_subj.jpg", bbox_inches='tight')


def std(df_a, Vertlevels):
    '''
    Computes standard deviation of subject mean CSA for each rescaling and different vertebrae levels
    '''
    print("\n====================std==========================\n")
    # TODO: normalization of CSA for intersubject studies
    min_vert = min(list(Vertlevels))
    for i in Vertlevels[1:]:
        for name, group in df_a.groupby('Rescale'):
            gt_CSA = 'CSA_C'+str(min_vert)+'_C'+str(i)
            std = group[gt_CSA].std()
            atrophy = set(group.reset_index().Rescale)
            print('CSA std on '+str(atrophy)+'  rescaled image C'+str(min_vert)+'/C'+str(i)+' is ',round(std, 3),' mm^2 ')
        print('\n')


def sample_size(df_a, conf, power, mean_control, mean_patient):
    '''
    calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    sample size with certainty 95% z(0.05/2)=1.96, power 80% z(Power = 80%)=0.84, ratio patients/control 1:1
    and with the assumption both samples have same std
    :param df: dataframe containing csv files values
    :param conf: Confidence level (i.e. 0.60, 0.70, 0.8, 0.85, 0.90, 0.95)
    :param power: Power level (i.e. 0.60, 0.70, 0.8, 0.85, 0.90, 0.95)
    :param atrophy: expected atrophy to detect
    :param mean_control: mean value of control group
    :param mean_patient: mean value of patient group
    '''
    print("\n====================size==========================\n")
    data3 = {'Confidence_Level':[0.60, 0.70, 0.8, 0.85, 0.90, 0.95],
             'z_value':[0.842, 1.04, 1.28, 1.44, 1.64, 1.96],}

    dfs = pd.DataFrame(data3)
    dfs = dfs.set_index('Confidence_Level')
    std = df_a.groupby('Rescale').get_group(1)['CSA_original'].std()
    print('std: '+str(std))
    num_n = 2 * ((dfs.at[conf, 'z_value'] + dfs.at[power, 'z_value']) ** 2) * (std ** 2)
    deno_n = (abs(mean_control - mean_patient))**2
    sample = ceil(num_n/deno_n)
    print('with '+str(power*100)+'% power, at '+str(conf*100)+'% significance, ratio 1:1 (patients/controls):')
    print('minimum sample size to detect mean '+str(mean_control-mean_patient)+' mm² atrophy: '+ str(sample))


def add_to_dataframe(df, Vertlevels):
    '''dataframe column additions gt_CSA, diff_CSA, perc_diff_CSA for different vertbrae levels
    :param df: original dataframe
    :return df_a: modified dataframe with added gt_CSA, diff_CSA, perc_diff_CSA for different vertbrae levels
    '''
    # create dataframes
    df1 = df.copy()
    df2 = df.copy()
    df3 = df.copy()
    df_gt2 = pd.DataFrame()

    # iterate throuh different vertebrae levels
    # dataframe and variable for iteration
    df_a = df.groupby(['Rescale','Filename']).mean()
    n = []
    max_vert = max(list(Vertlevels))
    min_vert = min(list(Vertlevels))
    diff_vert = np.setdiff1d(list(set(df['VertLevel'].values)), list(Vertlevels))
    # iteration
    for i in range(max_vert,min_vert,-1):
        df_gt2 = pd.DataFrame()
        # get GT values
        if i==max_vert:
            group_CSA_gt = df1.groupby('Rescale').get_group(1).set_index('VertLevel').drop(index=diff_vert).groupby(['Filename']).mean().CSA_original
        else:
            n.append(i+1)
            group_CSA_gt = df1.groupby('Rescale').get_group(1).set_index('VertLevel').drop(index=n).groupby(['Filename']).mean().CSA_original
        # iterate through Rescale groupby
        for name, group in df.groupby('Rescale'):
            atrophy = group['Rescale'].values
            # mean values from same subject  different vertebrae
            group2 = group.groupby(['Filename']).mean().reset_index()
            # to locate easily subjects put Filename as index
            group3 = group2.set_index(['Filename'])
            # iterate through dataframe subjects
            for subjectj in set(group2['Filename'].values):
                # if dataframe subject exist in GT (Rescale = 1)
                if subjectj in group_CSA_gt.index.values:
                    group3.at[subjectj,'gt_CSA_C'+str(min_vert)+'_C'+str(i)] = group_CSA_gt.loc[subjectj] * (atrophy[0] ** 2)
            df_gt = df.groupby('Rescale').get_group(atrophy[0]).groupby('Filename').mean()
            df_gt['gt_CSA_C'+str(min_vert)+'_C'+str(i)] = (group3['gt_CSA_C'+str(min_vert)+'_C'+str(i)].values)
            df_gt2 = pd.concat([df_gt2, df_gt])
        # regroup to add row in new dataframe
        df_a['gt_CSA_C'+str(min_vert)+'_C'+str(i)] = df_gt2['gt_CSA_C'+str(min_vert)+'_C'+str(i)].values


    # add CSA values for vertebrae levels of interest
    m = []
    l = []
    max_vert2 = max(list(Vertlevels))
    min_vert2 = min(list(Vertlevels))
    diff_vert2 = np.setdiff1d(list(set(df['VertLevel'].values)), list(Vertlevels))
    # iterate throug diffent vertebrae levels
    for j in range(max_vert2, min_vert2, -1):
        df2 = df.copy()
        if j == max_vert2:
            df_a['CSA_C'+str(min_vert2)+'_C'+str(max_vert2)] = df2.set_index('VertLevel').drop(index=diff_vert).groupby(['Rescale','Filename']).mean().values
            df_a['diff_C'+str(min_vert2)+'_C'+str(j)] = df_a['CSA_C'+str(min_vert2)+'_C'+str(j)].sub(df_a['gt_CSA_C'+str(min_vert2)+'_C'+str(j)]).abs()
            df_a['perc_diff_C'+str(min_vert2)+'_C'+str(j)] = 100 * df_a['diff_C'+str(min_vert2)+'_C'+str(j)].div(df_a['gt_CSA_C'+str(min_vert2)+'_C'+str(j)])
        else:
            # iterate droping lower vertebrae
            m.append(j+1)
            df_a['CSA_C'+str(min_vert2)+'_C'+str(j)] = np.nan
            df_a['CSA_C'+str(min_vert2)+'_C'+str(j)] = df2.set_index('VertLevel').drop(index=m).groupby(['Rescale','Filename']).mean().values
            df_a['diff_C'+str(min_vert2)+'_C'+str(j)] = df_a['CSA_C'+str(min_vert2)+'_C'+str(j)].sub(df_a['gt_CSA_C'+str(min_vert2)+'_C'+str(j)]).abs()
            df_a['perc_diff_C'+str(min_vert2)+'_C'+str(j)] = 100 * df_a['diff_C'+str(min_vert2)+'_C'+str(j)].div(df_a['gt_CSA_C'+str(min_vert2)+'_C'+str(j)])

    return df_a


def main(Vertlevels):
    '''
    main function: gathers stats and calls plots
    '''
    #read data
    data = pd.read_csv("csa.csv",decimal=".")
    data2 = {'Filename':data['Filename'],
             'VertLevel':data['VertLevel'],
             'CSA_original':data['MEAN(area)'],
             'Rescale':data['rescale'],}
    df = pd.DataFrame(data2)
    pd.set_option('display.max_rows', None)



    # Use filename instead of path to file
    df['Filename'] = list(('_'.join(os.path.basename(path).split('.')[0].split('_')[0:3])) for path in data['Filename'])

    # dataframe column additions gt,diff,perc diff for different vertbrae levels
    max_vert = max(list(Vertlevels))
    min_vert = min(list(Vertlevels))
    if Vertlevels is None:
        Vertlevels = set(list(df[Vertlevels].values))
    df_a = add_to_dataframe(df, Vertlevels)

    # print mean CSA gt
    print("\n====================mean==========================\n")
    print(" mean CSA: " + str(df.groupby('Rescale').get_group(1)['CSA_original'].mean()))

    #compute sample size
    print(df_a)
    sample_size(df_a, 0.95,0.8, 7.77, 0)

    #ground truth atrophy
    df4 = df.copy()
    atrophies = sorted(set(df4['Rescale'].values))
    #display number of subjects in test
    print("\n====================number subjects==========================\n")
    for atrophy in atrophies:
        number_sub = df4.groupby('Filename')['CSA_original'].mean().count()
        print('For rescaling '+str(atrophy)+' number of subjects is ' +str(number_sub))

    # compute std for different vertebrae levels
    std(df_a, Vertlevels)
    std_v = df_a.groupby('Rescale').get_group(1)['CSA_original'].std()

    # plot graph if verbose is present
    if arguments.v is not None:
        df_p1 = df_a.copy()
        columns_to_plot = [i for i in df_p1.columns if 'perc_diff' in i]
        get_plot(df_p1, columns_to_plot)
        df_p3 = df_a.copy()
        get_plot3(df_p3)
        df_p4 = df_a.copy()
        min_max_Vert = ['perc_diff_C'+str(min_vert)+'_C'+str(max_vert)]
        get_plot4(df_p4, min_max_Vert)
        #get_plot5(df5)
        get_plot_sample(1.96,(0.84, 1.282), std_v, 80)
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
        main(list(map(int, arguments.l)))
    else:
        parser.print_help()
