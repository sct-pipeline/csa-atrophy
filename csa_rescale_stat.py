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
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[3].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")

def get_plot(df_p1):
    '''plot percentage difference between measured simulated atrophy and ground truth atrophy
    for different vertebrae levels
    :param df_p1: dataframe for first plot
    '''
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    df_p1.groupby('rescale')[['perc_diff_C2_C3','perc_diff_C2_C4','perc_diff_C2_C5','perc_diff_C3_C5']].mean().plot(kind='bar', ax=axes[0], grid=True)
    axes[0].set_title('mean error in function of rescaling factor');
    axes[0].set_ylabel('error in %')
    df_p1.groupby('rescale')[['perc_diff_C2_C3','perc_diff_C2_C4','perc_diff_C2_C5','perc_diff_C3_C5']].std().plot(kind='bar', ax=axes[1], sharex=True, sharey=True, legend=False)
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
    df_p3.boxplot(column=['CSA_C2_C5'], by='rescale')
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
        num_n = ((z+z_p)**2)*((2*5)**2)
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
    plt.savefig("min_subj.jpg", bbox_inches='tight')


def std(df_a):
    '''
    Computes standard deviation of subject mean CSA for each rescaling and different vertebrae levels
    '''
    print("\n====================std==========================\n")
    # TODO: normalization of CSA for intersubject studies
    for i in range(3,6):
        for name, group in df_a.groupby('rescale'):
            gt_CSA = 'CSA_C2_C'+str(i)
            std = group[gt_CSA].std()
            atrophy = set(group.reset_index().rescale)
            print('CSA std on '+str(atrophy)+'  rescaled image C2/C'+str(i)+' is ',round(std, 3),' mm^2 ')
        print('\n')
    for name, group in df_a.groupby('rescale'):
        gt_CSA = 'CSA_C3_C5'
        std = group[gt_CSA].std()
        atrophy = set(group.reset_index().rescale)
        print('CSA std on '+str(atrophy)+'  rescaled image C3/C5 is ',round(std, 3),' mm^2 ')


def cohen_d(df, small = None):
    '''compute cohen d to express sample effect
    :param df: original DataFrame
    :param small: correction of cohen d for small sample <50, not used by default
    '''
    print("\n====================cohen d==========================\n")
    group1_mean = df.groupby('rescale').get_group(1)['CSA_C2_C5'].mean()
    group1_std = df.groupby('rescale').get_group(1)['CSA_C2_C5'].std()
    group_r = df.groupby('rescale')['CSA_C2_C5']
    nb = df.groupby('rescale').get_group(1)['CSA_C2_C5'].count()
    atrophy = sorted(set(df['rescale'].values))
    if small is None:
        print('\n cohen d for samples >50')
        for r in atrophy:
            sd_pooled = np.sqrt((group1_std**2 + (group_r.get_group(r).std())**2)/2)
            d = (group1_mean - group_r.get_group(r).mean())/sd_pooled
            print('cohen d: '+str(d)+'  for atrophy simultation  '+str(r))
    else:
        print('\n cohen d for small samples <50')
        for r in atrophy:
            sd_pooled = np.sqrt((group1_std**2 + (group_r.get_group(r).std())**2)/2)
            d = ((group1_mean - group_r.get_group(r).mean())/sd_pooled)*((nb-3)/(nb-2.25))*np.sqrt((nb-2)/nb)
            print('cohen d: '+str(round(d, 3))+'  for atrophy simultation  '+str(r))


def sample_size(df, conf, power, mean_control, mean_patient):
    '''
    calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    sample size with certainty 95% z(0.05/2)=1.96, power 80% z(Power = 80%)=0.84, ratio patients/control 1:1
    and with the assumption both samples have same std
    :param df: dataframe containing csv files values
    :param conf: Confdence level in % (i.e. 0.60, 0.70, 0.8, 0.85, 0.90, 0.95)
    :param power: Power level in % (i.e. 0.60, 0.70, 0.8, 0.85, 0.90, 0.95)
    :param atrophy: expected atrophy to detect
    :param mean_control: mean value of control group
    :param mean_patient: mean value of patient group
    '''
    print("\n====================size==========================\n")
    data3 = {'Confidence_Level':[0.60, 0.70, 0.8, 0.85, 0.90, 0.95],
             'z_value':[0.842, 1.04, 1.28, 1.44, 1.64, 1.96],}

    dfs = pd.DataFrame(data3)
    dfs = dfs.set_index('Confidence_Level')
    num_n = 2 * ((dfs.at[conf, 'z_value'] + dfs.at[power, 'z_value']) ** 2) * ((df.groupby('rescale').get_group(1)['CSA_C2_C5'].std()) ** 2)
    # mean CSA is approximatively 80 so %atrophy * 80 is the difference between study means
    deno_n = (abs(mean_control - mean_patient))**2
    sample = ceil(num_n/deno_n)
    print('with 80% power, at 5% significance, ratio 1:1 (patients/controls):')
    print('minimum sample size to detect mean 10% atrophy: '+ str(sample))


def dataframe_add(df):
    '''dataframe column additions gt_CSA, diff_CSA, perc_diff_CSA for different vertbrae levels
    :param df: original dataframe
    :return df_a: modified dataframe with added gt_CSA, diff_CSA, perc_diff_CSA for different vertbrae levels
    '''
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
    for i in range(5,1,-1):
        df_gt2 = pd.DataFrame()
        # get GT values
        if i==5:
            group_CSA_gt = df1.groupby('rescale').get_group(1).set_index('VertLevel').groupby('Filename').mean().CSA_C2_C5
        elif i==2:
            group_CSA_gt = df1.groupby('rescale').get_group(1).set_index('VertLevel').drop(index=i).groupby('Filename').mean().CSA_C2_C5
        else:
            n.append(i+1)
            group_CSA_gt = df1.groupby('rescale').get_group(1).set_index('VertLevel').drop(index=n).groupby('Filename').mean().CSA_C2_C5
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
    df_a['CSA_C2_C3'] = df2.set_index('VertLevel').drop(index=[4,5]).groupby(['rescale','Filename']).mean().values
    df_a['CSA_C3_C5'] = df2.set_index('VertLevel').drop(index=[2]).groupby(['rescale','Filename']).mean().values

    # difference between mean CSA and GT CSA for differennt vertebrae levels
    df_a['diff_C2_C5'] = df_a['CSA_C2_C5'].sub(df_a['gt_CSA_C2_C5']).abs()
    df_a['diff_C2_C4'] = df_a['CSA_C2_C4'].sub(df_a['gt_CSA_C2_C4']).abs()
    df_a['diff_C2_C3'] = df_a['CSA_C2_C3'].sub(df_a['gt_CSA_C2_C3']).abs()
    df_a['diff_C3_C5'] = df_a['CSA_C3_C5'].sub(df_a['gt_CSA_C2_C2']).abs()


    df_a['perc_diff_C2_C5'] = 100 * df_a['diff_C2_C5'].div(df_a['gt_CSA_C2_C5'])
    df_a['perc_diff_C2_C4'] = 100 * df_a['diff_C2_C4'].div(df_a['gt_CSA_C2_C4'])
    df_a['perc_diff_C2_C3'] = 100 * df_a['diff_C2_C3'].div(df_a['gt_CSA_C2_C3'])
    df_a['perc_diff_C3_C5'] = 100 * df_a['diff_C3_C5'].div(df_a['gt_CSA_C2_C2'])
    return df_a



def main():
    '''
    main function: gathers stats and calls plots
    '''
    #read data
    data = pd.read_csv("csa.csv",decimal=".")
    data2 = {'Filename':data['Filename'],
             'VertLevel':data['VertLevel'],
             'CSA_C2_C5':data['MEAN(area)'],
             'rescale':data['rescale'],}
    df = pd.DataFrame(data2)
    pd.set_option('display.max_rows', None)

    # Use filename instead of path to file
    df['Filename'] = list((os.path.basename(path).split('.')[0].split('_')[0]) for path in data['Filename'])

    # dataframe column additions gt,diff,perc diff for different vertbrae levels
    df_a = dataframe_add(df)

    # print mean CSA gt
    print("\n====================mean==========================\n")
    print(" mean CSA: " + str(df.groupby('rescale').get_group(1)['CSA_C2_C5'].mean()))

    # compute cohen d
    df5 = df.copy()
    cohen_d(df5)

    #compute sample size
    sample_size(df, 0.95,0.8, 80, 78)


    #ground truth atrophy
    df4 = df.copy()
    atrophy = sorted(set(df4['rescale'].values))
    #display number of subjects in test
    print("\n====================number subjects==========================\n")
    for atr in atrophy:
        number = df4.groupby('rescale').get_group(atr).groupby('Filename')['CSA_C2_C5'].mean().count()
        print('number of subjects for '+str(atr)+' rescaling is ' +str(number))

    # compute std for different vertebrae levels
    std(df_a)

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
