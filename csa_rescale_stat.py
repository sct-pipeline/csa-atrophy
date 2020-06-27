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
        description='Compute statistics based on the csv files containing the csa metrics:',
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
def concatenate_csv_files(path_results):
    '''Fetch and concatenate all csv files into one
    :param path_results: path to folder containing csv files for statistics
    '''
    files = []
    for file in os.listdir(path_results):
        if file.endswith(".csv"):
            files.append(os.path.join(path_results,file))
    metrics = pd.concat([pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[4].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")

def plot_perc_err(df, columns_to_plot):
    '''plot percentage difference between measured simulated atrophy and ground truth atrophy
    for different vertebrae levels
    :param df: dataframe for first plot
    :param columns_to_plot: perc_diff dataframe columns that will be plot
    '''
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    df.groupby('Rescale')[columns_to_plot].mean().plot(kind='bar', ax=axes[0], grid=True)
    axes[0].set_title('mean error in function of rescaling factor');
    axes[0].set_ylabel('error in %')
    df.groupby('Rescale')[columns_to_plot].std().plot(kind='bar', ax=axes[1], sharex=True, sharey=True, legend=False)
    axes[1].set_title('STD of error in function of rescaling factor');
    plt.xlabel('rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    plt.tight_layout()
    plt.savefig("err_plot.png")


def boxplot_csa(df):
    '''plot csa boxplot over different rescalings
    :param df: dataframe for third plot
    '''
    fig3 = plt.figure()
    df.boxplot(column=['csa_original'], by='Rescale')
    plt.ylabel('csa in mm^2')
    plt.savefig("csa_boxplot.png")


def boxplot_perc_err(df, min_max_vertlevels):
    '''plot percentage of error boxplot over different rescalings
    :param df: dataframe for fourth plot
    :param min_max_vertlevels: uses dataframe columns with most distant vertebrae levels. Example: perc_diff_c2_c5
    '''
    fig4 = plt.figure()
    df.boxplot(column=min_max_vertlevels, by='Rescale')
    plt.ylabel('error in %')
    plt.savefig("err_boxplot.png")


def plot_sample_size(z_conf, z_power, std, mean_csa):
    '''plot minimum number of patients required to detect an atrophy of X
    :param z_conf: z score for X % uncertainty. Example: z_conf=1.96
    :param z_power: z score for X % Power. Example: z_power=(0.84, 1.282)
    :param std: standard deviation around mean csa of control subjects,
    csa STD for atrophied subject and control subject are considered equivalent
    :param mean_csa: mean value of csa from which percentage atrophy is calculated. Example: 80
    '''
    fig6, ax = plt.subplots()
    # data for plotting
    n=[]
    for z_p in z_power:
        i = np.arange(1.5, 8.0, 0.05) # x_axis values ranging from 1.5 to 8.0 mm^2
        num_n = 2*((z_conf+z_p)**2)*(std**2) # numerator of sample size equation
        n.append(num_n/((i)**2))
        # plot
    ax.plot(i, n[0], label=('80% power'))
    ax.plot(i, n[1], label=('90% power'))
    ax.set_ylabel('number of participants per arm of study \n(patients or controls) with ratio 1:1')
    ax.set_xlabel('atrophy in mm^2')
    ax.set_title('minimum number of participants to detect an atrophy with 5% uncertainty\n std = '+str(round(std,2))+'mm², mean_csa = '+str(mean_csa)+'mm²')
    ax.legend()
    ax.grid()
    # create global variable for secax functions
    global mean_csa_sample
    mean_csa_sample = mean_csa

    def forward(i):
        i2 = i/mean_csa_sample*100
        return i2
    def inverse(i):
        return i/100*mean_csa_sample
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel('atrophy in %')
    plt.savefig("min_subj.png", bbox_inches='tight')


def std(df_a, vertlevels):
    '''
    Computes standard deviation of subject mean csa for each rescaling and different vertebrae levels
    :param df_a: dataframe grouped by subject containing information from add_to_dataframe
    :param vertlevels: vertebrae levels of interest list, by default list contains all vertebrae levels present in csv files
    '''
    print("\n====================std==========================\n")
    # TODO: possible normalization of csa for intersubject studies
    min_vert = min(list(vertlevels))
    for i in vertlevels[1:]:
        for name, group in df_a.groupby('Rescale'):
            gt_csa = 'csa_c'+str(min_vert)+'_c'+str(i)
            std = group[gt_csa].std()
            atrophy = set(group.reset_index().Rescale)
            print('csa std on '+str(atrophy)+'  rescaled image c'+str(min_vert)+'/c'+str(i)+' is ',round(std, 3),' mm^2 ')
        print('\n')


def sample_size(df_a, conf, power, mean_control=None, mean_patient=None, atrophy=None):
    '''
    calculate the minimum number of patients required to detect an atrophy of X (i.e. power analysis)
    sample size with certainty 95% z(0.05/2)=1.96, power 80% z(Power = 80%)=0.84, ratio patients/control 1:1
    and with the assumption both samples have same std. Example sample_size(df_a, 0.95, 0.8, mean_control=None, mean_patient=None, atrophy=7.7))
    :param df_a: dataframe grouped by subject containing information from add_to_dataframe
    :param conf: Confidence level. Example 0.8
    :param power: Power level. Example 0.9
    :param mean_control: mean csa value of control group
    :param mean_patient: mean csa value of patient group
    :param atrophy: expected atrophy to detect in mm^2. Example atrophy=7.7
    '''
    print("\n====================size==========================\n")
    z_score_dict = {'confidence_Level':[0.60, 0.70, 0.8, 0.85, 0.90, 0.95],
             'z_value':[0.842, 1.04, 1.28, 1.44, 1.64, 1.96],}

    df_sample = pd.DataFrame(z_score_dict)
    df_sample = df_sample.set_index('confidence_Level')
    std = df_a.groupby('Rescale').get_group(1)['csa_original'].std()
    print('std: '+str(std))
    num_n = 2 * ((df_sample.at[conf, 'z_value'] + df_sample.at[power, 'z_value']) ** 2) * (std ** 2)
    if atrophy:
        atrophy = atrophy
    elif mean_control is not None & mean_patient is not None:
        atrophy = abs(mean_control - mean_patient)
    else:
        print('input error: either input mean_control and mean_patient or atrophy in mm^2')
    deno_n = (abs(atrophy))**2
    sample = ceil(num_n/deno_n)
    print('with '+str(power*100)+'% power, at '+str(conf*100)+'% significance, ratio 1:1 (patients/controls):')
    print('minimum sample size to detect mean '+str(atrophy)+' mm² atrophy: '+ str(sample))


def add_to_dataframe(df, vertlevels):
    '''dataframe column additions gt_csa, diff_csa, perc_diff_csa for different vertbrae levels
    :param df: original dataframe with csv file data
    :param vertlevels: vertebrae levels of interest list, by default list contains all vertebrae
    levels present in csv files
    :return df_a: modified dataframe with added gt_csa, diff_csa, perc_diff_csa for different vertebrae levels
    '''
    # create dataframes
    df1 = df.copy()
    df2 = df.copy()
    df_gt2 = pd.DataFrame()

    # variables for iteration
    df_a = df.groupby(['Rescale','Filename']).mean()
    n = []
    max_vert = max(vertlevels)
    min_vert = min(vertlevels)
    diff_vert = np.setdiff1d(list(set(df['VertLevel'].values)), list(vertlevels))
    # iterate across different vertebrae levels and add ground truth values to each subject in dataframe
    for i in range(max_vert,min_vert,-1):
        df_gt2 = pd.DataFrame()
        # get GT values
        if i==max_vert:
            group_csa_gt = df1.groupby('Rescale').get_group(1).set_index('VertLevel').drop(index=diff_vert).groupby(['Filename']).mean().csa_original
        else:
            n.append(i+1)
            group_csa_gt = df1.groupby('Rescale').get_group(1).set_index('VertLevel').drop(index=n).groupby(['Filename']).mean().csa_original
        # iterate across Rescale groupby
        for name, group in df.groupby('Rescale'):
            atrophy = group['Rescale'].values
            # mean values from same subject  different vertebrae
            group2 = group.groupby(['Filename']).mean().reset_index()
            # to locate easily subjects put Filename as index
            group3 = group2.set_index(['Filename'])
            # iterate across dataframe subjects
            for subject_j in set(group2['Filename'].values):
                # if dataframe subject exist in GT (Rescale = 1)
                if subject_j in group_csa_gt.index.values:
                    group3.at[subject_j,'gt_csa_c'+str(min_vert)+'_c'+str(i)] = group_csa_gt.loc[subject_j] * (atrophy[0] ** 2)
            df_gt = df.groupby('Rescale').get_group(atrophy[0]).groupby('Filename').mean()
            df_gt['gt_csa_c'+str(min_vert)+'_c'+str(i)] = (group3['gt_csa_c'+str(min_vert)+'_c'+str(i)].values)
            df_gt2 = pd.concat([df_gt2, df_gt])
        # regroup to add row in new dataframe
        df_a['gt_csa_c'+str(min_vert)+'_c'+str(i)] = df_gt2['gt_csa_c'+str(min_vert)+'_c'+str(i)].values

    # add csa, diff and perc_diff values for vertebrae levels of interest for each subject
    m = []
    l = []
    max_vert2 = max(list(vertlevels))
    min_vert2 = min(list(vertlevels))
    diff_vert2 = np.setdiff1d(list(set(df['VertLevel'].values)), list(vertlevels))
    # iterate across different vertebrae levels
    for j in range(max_vert2, min_vert2, -1):
        if j == max_vert2:
            df_a['csa_c'+str(min_vert2)+'_c'+str(max_vert2)] = df2.set_index('VertLevel').drop(index=diff_vert).groupby(['Rescale','Filename']).mean().values
            df_a['diff_c'+str(min_vert2)+'_c'+str(j)] = df_a['csa_c'+str(min_vert2)+'_c'+str(j)].sub(df_a['gt_csa_c'+str(min_vert2)+'_c'+str(j)]).abs()
            df_a['perc_diff_c'+str(min_vert2)+'_c'+str(j)] = 100 * df_a['diff_c'+str(min_vert2)+'_c'+str(j)].div(df_a['gt_csa_c'+str(min_vert2)+'_c'+str(j)])
        else:
            m.append(j+1)
            df_a['csa_c'+str(min_vert2)+'_c'+str(j)] = np.nan
            df_a['csa_c'+str(min_vert2)+'_c'+str(j)] = df2.set_index('VertLevel').drop(index=m).groupby(['Rescale','Filename']).mean().values
            df_a['diff_c'+str(min_vert2)+'_c'+str(j)] = df_a['csa_c'+str(min_vert2)+'_c'+str(j)].sub(df_a['gt_csa_c'+str(min_vert2)+'_c'+str(j)]).abs()
            df_a['perc_diff_c'+str(min_vert2)+'_c'+str(j)] = 100 * df_a['diff_c'+str(min_vert2)+'_c'+str(j)].div(df_a['gt_csa_c'+str(min_vert2)+'_c'+str(j)])
    return df_a


def main(vertlevels_input):
    '''
    main function: gathers stats and calls plots
    :param vertlevels_input: vertebrae levels of interest, arguments of flag -l
    '''
    #read data
    data = pd.read_csv("csa.csv",decimal=".")
    data2 = {'Filename':data['Filename'],
             'VertLevel':data['VertLevel'],
             'csa_original':data['MEAN(area)'],
             'Rescale':data['rescale'],}
    df = pd.DataFrame(data2)
    pd.set_option('display.max_rows', None)

    # Change dataframe['Filenamee'] to basename,
    # Remove rescale suffix
    df['Filename'] = list((os.path.basename(path).split('_r')[0]+'_'+os.path.basename(path).split('_')[3].split('.nii.gz')[0]) for path in data['Filename'])

    # verify if vertlevels were given in input by user
    if vertlevels_input is None:
        vertlevels = list(set(df['VertLevel'].values))
        print(vertlevels)
    elif vertlevels_input is not None:
        vertlevels = list(map(int, vertlevels_input))
        if all(elem in set(list(df['VertLevel'].values)) for elem in vertlevels):
            pass
        else:
            print('error: Input vertebrae levels ',vertlevels,' do not exist in csv files')
            exit()
    # dataframe column additions gt,diff, perc diff for different vertbrae levels
    df_a = add_to_dataframe(df, vertlevels)

    # print mean csa without rescaling
    print("\n====================mean==========================\n")
    mean_csa = df.groupby('Rescale').get_group(1)['csa_original'].mean()
    print(" mean csa: " + str(mean_csa))

    #compute sample size
    # conf = confidence level
    # power = power level
    sample_size(df_a, conf=0.95, power=0.8, mean_control=None, mean_patient=None, atrophy=7.7)

    #ground truth atrophy
    atrophies = sorted(set(df['Rescale'].values))
    # display number of subjects in test (transformed subjects are considered different)
    print("\n====================number subjects==========================\n")
    for atrophy in atrophies:
        number_sub = df.groupby('Filename')['csa_original'].mean().count()
        print('For rescaling '+str(atrophy)+' number of subjects is ' +str(number_sub))

    # compute std for different vertebrae levels
    std(df_a, vertlevels)
    std_v = df_a.groupby('Rescale').get_group(1)['csa_original'].std()

    # plot graph if verbose is present
    if arguments.v is not None:
        columns_to_plot = [i for i in df_a.columns if 'perc_diff' in i]
        plot_perc_err(df_a, columns_to_plot)
        boxplot_csa(df_a)
        max_vert = max(vertlevels)
        min_vert = min(vertlevels)
        min_max_Vert = ['perc_diff_c'+str(min_vert)+'_c'+str(max_vert)]
        boxplot_perc_err(df_a, min_max_Vert)
        # z = z_score for confidence level,
        # z_power = z_score for confidence level,
        # std_v = STD of subjects without rescaling CSA values
        # mean_csa =  mean CSA value of subjects without resclaing
        plot_sample_size(z_conf=1.96, z_power=(0.84, 1.282), std=std_v, mean_csa=mean_csa)
        print('\nfigures have been ploted in dataset')


# Run
#########################################################################################
if __name__ == "__main__":
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    if arguments.h is None:
        path_results = os.path.join(os.getcwd(),arguments.i)
        concatenate_csv_files(path_results)
        main(arguments.l)
    else:
        parser.print_help()
