from __future__ import division

# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Evaluate the robustness of automated CSA with global rescaled images
#
# example: python csa_rescale_stat.py -i <results>
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
#
# About the license: see the file LICENSE
#########################################################################################

import pandas as pd
import numpy as np
import os
import argparse
from scipy import stats
import matplotlib.pyplot as plt
from math import ceil
import yaml


# Parser
#########################################################################################

def get_parser():
    """parser function"""
    parser = argparse.ArgumentParser(
        description="Compute statistics based on the csv files containing the CSA metrics:",
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        default='csa_atrophy_results',
        help='Path to folder that contains output csv files (e.g. "csa_atrophy_results/results")',
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-v',
        help='Verbose, plotting figures',
        nargs="*"
    )
    optional.add_argument(
        '-l',
        help='Vertebrae levels on which to compute the statistics. \nExample: -l 2 3 4 5',
        nargs="*",
    )
    optional.add_argument(
        '-o',
        help='Path where figures will be saved. By default, they will be saved in the current directory.',
        default=""
    )
    optional.add_argument(
        '-config',
        help='Path to config file, which contains parameters for the statistics and figures.',
    )
    return parser


# Functions
def concatenate_csv_files(path_results):
    """Fetch and concatenate data from all csv files in results/csa_data to compute statistics with pandas
    :param path_results: path to folder containing csv files for statistics
    """
    files = []
    for file in os.listdir(path_results):
        if ".csv" in file and "csa" in file:
            files.append(os.path.join(path_results, file))
    metrics = pd.concat(
        [pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[4].split('.csv')[0]) for f in files])
    metrics.to_csv("csa.csv")


def yaml_parser(config_file):
    """parse config_script.yml file containing pipeline's parameters"""
    with open(config_file, 'r') as config_var:
        config_param = yaml.safe_load(config_var)
    return config_param


def plot_perc_err(df, columns_to_plot, path_output):
    """plot percentage difference between simulated atrophy and ground truth atrophy
    for different vertebrae levels
    :param df: dataframe for first plot
    :param columns_to_plot: perc_diff dataframe columns for plotting
    :param path_output: directory in which plot is saved
    """
    df = df.reset_index().set_index('Rescale')
    df['subject'] = list(tf.split('_T2w')[0] for tf in df['Filename'])
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    df.groupby('Rescale')[columns_to_plot].mean().plot(kind='bar', ax=axes[0], grid=True)
    axes[0].set_title('mean error function of rescaling factor')
    axes[0].set_ylabel('error in %')
    df.groupby(['Rescale', 'subject']).mean().groupby('Rescale').std()[columns_to_plot].plot(kind='bar', ax=axes[1],
                                                                                             sharex=True, sharey=True,
                                                                                             legend=False)
    axes[1].set_title('STD of error function of rescaling factor')
    plt.xlabel('rescaling factor')
    plt.ylabel('error in %')
    plt.grid()
    plt.tight_layout()
    output_file = path_output + "/err_plot.png"
    plt.savefig(output_file)


def boxplot_csa(df, path_output):
    """boxplot CSA for different rescaling values
    :param df: dataframe for second plot
    :param path_output: directory in which plot is saved
    """
    fig2 = plt.figure()
    df.boxplot(column=['csa_original'], by='Rescale')
    plt.ylabel('csa in mm^2')
    output_file = path_output + "/csa_boxplot.png"
    plt.savefig(output_file)


def boxplot_perc_err(df, min_max_vertlevels, path_output):
    """boxplot percentage of error over different rescaling values
    :param df: dataframe for third plot
    :param min_max_vertlevels: uses dataframe column with most distant vertebrae levels. Example: perc_diff_c2_c5
    :param path_output: directory in which plot is saved
    """
    fig3 = plt.figure()
    df.boxplot(column=min_max_vertlevels, by='Rescale')
    plt.ylabel('error in %')
    output_file = path_output + "/err_boxplot.png"
    plt.savefig(output_file)


def plot_sample_size(z_conf, z_power, std, mean_csa, path_output):
    """plot minimum number of patients required to detect an atrophy of a given value
    :param z_conf: z score for X % uncertainty. Example: z_conf=1.96
    :param z_power: z score for X % Power. Example: z_power=(0.84, 1.282)
    :param std: STD around mean CSA of control subjects (without rescaling),
    CSA STD for atrophied subjects and control subjects are considered equal
    :param mean_csa: mean value of CSA to compute the atrophy percentage. Example: 80
    :param path_output: directory in which plot is saved
    """
    fig_sample, ax = plt.subplots()
    # data for plotting
    n = []
    for z_p in z_power:
        atrophy = np.arange(1.5, 8.0, 0.05)  # x_axis values ranging from 1.5 to 8.0 mm^2
        num_n = 2 * ((z_conf + z_p) ** 2) * (std ** 2)  # numerator of sample size equation
        n.append(num_n / ((atrophy) ** 2))
    # plot
    ax.plot(atrophy, n[0], label='80% power')
    ax.plot(atrophy, n[1], label='90% power')
    ax.set_ylabel('number of participants per group of study \n(patients or controls) with ratio 1:1')
    ax.set_xlabel('atrophy in mm^2')
    # create global variable for secax (second axis) functions
    global mean_csa_sample
    mean_csa_sample = mean_csa

    ax.set_title('minimum number of participants to detect an atrophy with 5% uncertainty\n std = ' + str(
        round(std, 2)) + 'mm², mean_csa = ' + str(mean_csa_sample) + 'mm²')
    ax.legend()
    ax.grid()
    def forward(atrophy):
        return atrophy / mean_csa_sample * 100

    def inverse(atrophy):
        return atrophy / 100 * mean_csa_sample

    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel('atrophy in %')
    output_file = path_output + "/min_subj.png"
    plt.savefig(output_file, bbox_inches='tight')


def std(df, vertlevels):
    """
    Compute STD of inter-subject mean CSA for each rescaling across different vertebrae levels
    :param df_a: dataframe grouped by subject containing information from add_to_dataframe
    :param vertlevels: vertebrae levels of interest list, by default list contains all vertebrae levels
    present in csv files
    """
    # TODO: possible normalization of csa for inter-subject studies
    df = df.reset_index().set_index('Rescale')
    df['subject'] = list(tf.split('_T')[0] for tf in df['Filename'])
    min_vert = min(list(vertlevels))
    for i in vertlevels[1:]:
        for name, group in df.groupby('Rescale'):
            csa_col = 'csa_c' + str(min_vert) + '_c' + str(i)
            std = group.groupby('subject')[csa_col].mean().std()
            cov = stats.variation(group.groupby('subject')[csa_col].mean()) * 100
            atrophy = set(group.reset_index().Rescale)
            print('csa std on ', atrophy, '  rescaled image c' + str(min_vert) + '/c' + str(i) + ' is ',
                  round(std, 3), ' mm^2 and cov is ', round(cov, 3), '%')
            if group.index[0] == 1:
                std_v = group.groupby('subject')['csa_original'].mean().std()
        print('\n')
    return std_v


def std_suject(df, vertlevels):
    """
    Compute mean STD of intra-subject (with transformation) CSA for each rescaling across different vertebrae levels
    :param df_a: dataframe grouped by subject containing information from add_to_dataframe
    :param vertlevels: vertebrae levels of interest list, by default list contains all vertebrae levels
    present in csv files
    """
    df = df.reset_index().set_index('Rescale')
    df['subject'] = list(tf.split('_T')[0] for tf in df['Filename'])
    min_vert = min(list(vertlevels))
    for i in vertlevels[1:]:
        for name, group in df.groupby('Rescale'):
            csa_col = 'csa_c' + str(min_vert) + '_c' + str(i)
            mean_cov = (group.groupby('subject')[csa_col].std() / group.groupby('subject')[csa_col].mean()).mean() * 100
            mean_std = group.groupby('subject')[csa_col].agg([np.mean, np.std]).mean()
            atrophy = set(group.reset_index().Rescale)
            print('mean csa std with ', atrophy, 'rescaling on c' + str(min_vert) + '/c' + str(i) + ' vertebrae is ',
                  round(mean_std['mean'], 3), ' mm^2, std', round(mean_std['std'], 3), ' cov is ', round(mean_cov, 3),
                  '%')
        print('\n')


def sample_size(df_a, atrophy, conf, power, mean_control=None, mean_patient=None):
    """
    Calculate the minimum number of patients required to detect an atrophy of a given value (i.e. power analysis),
    ratio patients/control 1:1 and with the assumption that both samples have the same STD.
    ref: Ard, M Colin, and Steven D Edland. “Power calculations for clinical trials in Alzheimer's disease.”
    Journal of Alzheimer's disease : JAD vol. 26 Suppl 3,Suppl 3 (2011): 369-77. doi:10.3233/JAD-2011-0062
    Example: sample_size(df_a, 0.95, 0.8, mean_control=None, mean_patient=None, atrophy=7.7)
    :param df_a: dataframe grouped by subject containing information from add_to_dataframe
    :param conf: Confidence level. Example 0.8
    :param power: Power level. Example 0.9
    :param mean_control: mean CSA value of control group
    :param mean_patient: mean csa value of patient group
    :param atrophy: expected atrophy in mm^2. Example atrophy=7.7
    """
    z_score_dict = {'confidence_Level': [0.60, 0.70, 0.8, 0.85, 0.90, 0.95],
                    'z_value': [0.842, 1.04, 1.28, 1.44, 1.64, 1.96], }

    df_sample = pd.DataFrame(z_score_dict)
    df_sample = df_sample.set_index('confidence_Level')
    std_sample = df_a.groupby('Rescale').get_group(1)['csa_original'].std()
    print('std: ' + str(std_sample))
    num_n = 2 * ((df_sample.at[conf, 'z_value'] + df_sample.at[power, 'z_value']) ** 2) * (std_sample ** 2)
    if atrophy:
        atrophy = atrophy
    elif mean_control is not None & mean_patient is not None:
        atrophy = abs(mean_control - mean_patient)
    else:
        print('input error: either input mean_control and mean_patient or atrophy in mm^2')
    deno_n = (abs(atrophy)) ** 2
    sample = ceil(num_n / deno_n)
    print('with ' + str(power * 100) + '% power, at ' + str(
        conf * 100) + '% significance, ratio 1:1 (patients/controls):')
    print('minimum sample size to detect mean ' + str(atrophy) + ' mm² atrophy: ' + str(sample))


def add_to_dataframe(df, vertlevels):
    """ dataframe column additions gt_csa, diff_csa, perc_diff_csa for different vertebrae levels
    :param df: original dataframe with csv file data
    :param vertlevels: vertebrae levels of interest list, by default list contains all vertebrae
    levels present in csv files
    :return df_a: modified dataframe with added gt_csa, diff_csa, perc_diff_csa for different vertebrae levels
    """
    # create dataframes
    df1 = df.copy()
    df2 = df.copy()

    # variables for iteration
    df_a = df.groupby(['Rescale', 'Filename']).mean()
    n = []
    max_vert = max(vertlevels)
    min_vert = min(vertlevels)
    diff_vert = np.setdiff1d(list(set(df['VertLevel'].values)), list(vertlevels))
    # iterate across different vertebrae levels and add ground truth values to each subject in dataframe
    for i in range(max_vert, min_vert, -1):
        # df_gt2 is used as an intermediary dataframe to store values from df_gt
        df_gt2 = pd.DataFrame()
        # get GT values
        if i == max_vert:
            group_csa_gt = df1.groupby('Rescale').get_group(1).set_index('VertLevel').drop(index=diff_vert).groupby(
                ['Filename']).mean().csa_original
        else:
            n.append(i + 1)
            group_csa_gt = df1.groupby('Rescale').get_group(1).set_index('VertLevel').drop(index=n).groupby(
                ['Filename']).mean().csa_original
        # iterate across Rescale groupby
        for name, group in df.groupby('Rescale'):
            atrophy = group['Rescale'].values[0]
            # mean CSA values across vertebrae for each transformaion and subject
            group2 = group.groupby(['Filename']).mean()
            # iterate across dataframe subjects
            for subject_j in set(group2.index.values):
                # if dataframe subject exist in GT (without rescaling)
                if subject_j in group_csa_gt.index.values:
                    group2.at[subject_j, 'gt_csa_c' + str(min_vert) + '_c' + str(i)] = group_csa_gt.loc[subject_j] * (atrophy ** 2)
            # df_gt is used to store values from group2 if present in GT
            df_gt = df.groupby('Rescale').get_group(atrophy).groupby('Filename').mean()
            df_gt['gt_csa_c' + str(min_vert) + '_c' + str(i)] = (
                group2['gt_csa_c' + str(min_vert) + '_c' + str(i)].values)
            df_gt2 = pd.concat([df_gt2, df_gt])
        # add column ground truth to final dataframe
        df_a['gt_csa_c' + str(min_vert) + '_c' + str(i)] = df_gt2['gt_csa_c' + str(min_vert) + '_c' + str(i)].values

    # add columns csa, diff and perc_diff values for vertebrae levels of interest for each subject
    m = []
    max_vert2 = max(list(vertlevels))
    min_vert2 = min(list(vertlevels))
    # iterate across different vertebrae levels
    for j in range(max_vert2, min_vert2, -1):
        if j == max_vert2:
            # add column csa
            df_a['csa_c' + str(min_vert2) + '_c' + str(max_vert2)] = df2.set_index('VertLevel').drop(
                index=diff_vert).groupby(['Rescale', 'Filename']).mean().values
            # add column diff
            df_a['diff_c' + str(min_vert2) + '_c' + str(j)] = df_a['csa_c' + str(min_vert2) + '_c' + str(j)].sub(
                df_a['gt_csa_c' + str(min_vert2) + '_c' + str(j)]).abs()
            # add column perc_diff
            df_a['perc_diff_c' + str(min_vert2) + '_c' + str(j)] = 100 * df_a[
                'diff_c' + str(min_vert2) + '_c' + str(j)].div(df_a['gt_csa_c' + str(min_vert2) + '_c' + str(j)])
        else:
            m.append(j + 1)
            # If a subject or transformation is present in GT but not rescaled fill value with nan
            df_a['csa_c' + str(min_vert2) + '_c' + str(j)] = np.nan
            # add column csa
            df_a['csa_c' + str(min_vert2) + '_c' + str(j)] = df2.set_index('VertLevel').drop(index=m).groupby(
                ['Rescale', 'Filename']).mean().values
            # add column diff
            df_a['diff_c' + str(min_vert2) + '_c' + str(j)] = df_a['csa_c' + str(min_vert2) + '_c' + str(j)].sub(
                df_a['gt_csa_c' + str(min_vert2) + '_c' + str(j)]).abs()
            # add column perc_diff
            df_a['perc_diff_c' + str(min_vert2) + '_c' + str(j)] = 100 * df_a[
                'diff_c' + str(min_vert2) + '_c' + str(j)].div(df_a['gt_csa_c' + str(min_vert2) + '_c' + str(j)])
    return df_a


def main():
    """
    main function, gather stats and call plots
    :param vertlevels_input: vertebrae levels of interest, arguments of flag -l
    """
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    path_results = os.path.join(os.getcwd(), arguments.i)
    concatenate_csv_files(path_results)
    vertlevels_input = arguments.l
    path_output = arguments.o
    # read data
    data = pd.read_csv("csa.csv", decimal=".")
    data2 = {'Filename': data['Filename'],
             'VertLevel': data['VertLevel'],
             'csa_original': data['MEAN(area)'],
             'Rescale': data['rescale'], }
    df = pd.DataFrame(data2)
    pd.set_option('display.max_rows', None)

    # fetch parameters from config.yaml file
    config_param = yaml_parser(arguments.config)

    # Change dataframe['Filename'] to basename and remove rescale suffix
    df['Filename'] = list(
        (os.path.basename(path).split('_r')[0] + '_' + os.path.basename(path).split('_')[3].split('.nii.gz')[0]) for
        path in data['Filename'])

    # verify if vertlevels of interest were given in input by user
    if vertlevels_input is None:
        vertlevels = list(set(df['VertLevel'].values))
    elif vertlevels_input is not None:
        vertlevels = list(map(int, vertlevels_input))
        if all(elem in set(list(df['VertLevel'].values)) for elem in vertlevels):
            pass
        else:
            print('error: Input vertebrae levels ', vertlevels, ' do not exist in csv files')
            exit()
    # dataframe column additions gt, diff, perc_diff for different vertebrae levels
    df_a = add_to_dataframe(df, vertlevels)

    # print mean CSA without rescaling
    print("\n====================mean==========================\n")
    mean_csa = df.groupby('Rescale').get_group(1)['csa_original'].mean()
    print(" mean csa: " + str(mean_csa))

    # compute sample size
    print("\n==================sample_size======================\n")
    # configuration parameters can be modified in config.yaml file
    atrophy = config_param['stats']['sample_size']['atrophy_sample']
    # conf = confidence level
    conf = config_param['stats']['sample_size']['conf']
    # power = power level
    power = config_param['stats']['sample_size']['power']
    sample_size(df_a, atrophy, conf, power, mean_control=None, mean_patient=None)

    # ground truth atrophy
    atrophies = sorted(set(df['Rescale'].values))
    # display number of subjects in test (multiple transformations of the same subjects are considered different)
    print("\n=================number subjects=======================\n")
    df['subject'] = list(tf.split('_T')[0] for tf in df['Filename'])
    for atrophy in atrophies:
        number_sub = df.groupby('subject')['csa_original'].mean().count()
        number_tf = df.groupby(['Filename', 'subject'])['csa_original'].mean().count()
        print('For rescaling ' + str(atrophy) + ' number of subjects is ' + str(number_sub) +
        ' and number of transformations per subject ' + str(number_tf))


    # compute STD for different vertebrae levels
    print("\n=======================std============================\n")
    std_v = std(df_a, vertlevels)
    print("\n===================std_subject=========================\n")
    std_suject(df_a, vertlevels)
    print("\n===================original_csa========================\n")
    print(df_a['csa_original'])

    # plot graph if verbose is present
    if arguments.v is not None:
        columns_to_plot = [i for i in df_a.columns if 'perc_diff' in i]
        plot_perc_err(df_a, columns_to_plot, path_output)
        boxplot_csa(df_a, path_output)
        max_vert = max(vertlevels)
        min_vert = min(vertlevels)
        min_max_vert = ['perc_diff_c' + str(min_vert) + '_c' + str(max_vert)]
        boxplot_perc_err(df_a, min_max_vert, path_output)
        # z_conf = z_score for confidence level,
        # z_power = z_score for power level,
        # std_v = STD of subjects without rescaling CSA values
        # mean_csa =  mean CSA value of subjects without rescaling
        plot_sample_size(z_conf=1.96, z_power=(0.84, 1.282), std=std_v, mean_csa=mean_csa, path_output=path_output)
        print('\nfigures have been plotted in dataset')


# Run
#########################################################################################
if __name__ == "__main__":
    main()
