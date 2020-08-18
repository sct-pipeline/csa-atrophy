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
from ruamel.yaml import YAML


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
    mandatory.add_argument(
        '-config',
        required=True,
        help='Path to config file, which contains parameters for the statistics and figures.',
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-fig',
        help='Generate figures',
        action='store_true'
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
    return parser


# Functions
def concatenate_csv_files(path_results):
    """Fetch and concatenate data from all csv files in results/csa_data to compute statistics with pandas
    :param path_results: path to folder containing csv files for statistics
    """
    files = []
    for file in os.listdir(path_results):
        if ".csv" in file and "csa" and "sub" in file:
            files.append(os.path.join(path_results, file))
    if not files:
        raise FileExistsError("Folder {} does not contain any results csv file.".format(path_results))
    metrics = pd.concat(
        [pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[4].split('.csv')[0]) for f in files])
    # output csv file in PATH_RESULTS
    metrics.to_csv(os.path.join(path_results, r'csa_all.csv'))


def yaml_parser(config_file):
    """parse config_script.yml file containing pipeline's parameters"""
    with open(config_file, 'r') as config_var:
        yaml = YAML(typ='safe')
        config_param = yaml.load(config_var)
    return config_param


def plot_perc_err(df, columns_to_plot, df_results, path_output):
    """plot percentage difference between simulated atrophy and ground truth atrophy
    for different vertebrae levels
    :param df: dataframe for first plot
    :param columns_to_plot: perc_diff dataframe columns for plotting
    :param path_output: directory in which plot is saved
    """
    df = df.reset_index().set_index('rescale_in_percent')
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    df.groupby('rescale_in_percent')[columns_to_plot].mean().plot(kind='bar', ax=axes[0], grid=True)
    axes[0].set_title('mean error function of area rescaling')
    axes[0].set_ylabel('error in %')
    df_results['intra_sub_cov'].plot(kind='bar', ax=axes[1], sharex=True, sharey=True, legend=False)
    axes[1].set_title('COV of intra-subject CSA in function of area rescaling')
    plt.xlabel('area rescaling in %')
    plt.ylabel('COV in %')
    plt.grid()
    plt.tight_layout()
    output_file = path_output + "/err_plot.png"
    plt.savefig(output_file)


def boxplot_csa(df, path_output):
    """boxplot CSA for different rescaling values
    :param df: dataframe with csv files data
    :param path_output: directory in which plot is saved
    """
    fig2 = plt.figure()
    df.boxplot(column=['MEAN(area)'], by='rescale_in_percent', showmeans=True, meanline=True)
    plt.title('Boxplot of CSA in function of area rescaling')
    plt.suptitle("")
    plt.ylabel('csa in mm^2')
    plt.xlabel('area rescaling in %')
    output_file = path_output + "/boxplot_csa.png"
    plt.savefig(output_file)


def boxplot_atrophy(df, column_ratio, path_output):
    """boxplot error for different rescaling values
    :param df: dataframe with csv files data
    :param min_max_vertlevels: uses dataframe column with most distant vertebrae levels. Example: perc_diff_c2_c5
    :param path_output: directory in which plot is saved
    """
    fig3 = plt.figure()
    df.boxplot(column=column_ratio, by='rescale_in_percent', showmeans=True, meanline=True)
    plt.title('boxplot of measured atrophy in function of area rescaling')
    plt.ylabel('measured atrophy in %')
    plt.xlabel('area rescaling in %')
    plt.suptitle("")
    output_file = path_output + "/boxplot_atrophy.png"
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


def inter_subject_mean(df, df_results):
    """
    Compute inter-subject CSA average for each rescaling
    :param df: dataframe with csv files data
    :param df_results: results dataframe
    """
    df = df.reset_index()
    for name, group in df.groupby('rescale_in_percent'):
        inter_sub_mean = group.groupby('subject')['MEAN(area)'].mean().mean()
        atrophy = set(group.reset_index().rescale_in_percent)
        print("For rescaling: {} MEAN inter-subject CSA is {} mm^2".format(atrophy, round(inter_sub_mean, 3)))
        # fill results dataframe
        df_results.loc[atrophy, 'inter_sub_mean'] = inter_sub_mean
    print('\n')


def inter_subject_std(df, df_results):
    """
    Compute inter-subject CSA STD for each rescaling
    :param df: dataframe with csv files data
    :param df_results: results dataframe
    """
    # TODO: possible normalization of csa for inter-subject studies
    df = df.reset_index()
    for name, group in df.groupby('rescale_in_percent'):
        inter_sub_std = group.groupby('subject')['MEAN(area)'].mean().std()
        inter_sub_cov = stats.variation(group.groupby('subject')['MEAN(area)'].mean()) * 100
        atrophy = set(group.reset_index().rescale_in_percent)
        print("For rescaling: {} STD of inter-subject CSA is {} mm^2 and COV is {} %".format(atrophy, round(inter_sub_std, 3), round(inter_sub_cov, 3)))
        # fill results dataframe
        df_results.loc[atrophy, 'inter_sub_std'] = inter_sub_std
        df_results.loc[atrophy, 'inter_sub_cov'] = inter_sub_cov
    print('\n')


def intra_subject_std(df, df_results):
    """
    Compute intra-subject CSA STD for each rescaling
    :param df: dataframe with csv files data
    :param df_results: results dataframe
    """
    df = df.reset_index()
    for name, group in df.groupby('rescale_in_percent'):
        intra_sub_std = group.groupby('subject')['MEAN(area)'].std().mean()
        intra_sub_cov = (group.groupby('subject')['MEAN(area)'].std() / group.groupby('subject')['MEAN(area)'].mean()).mean() * 100
        atrophy = set(group.reset_index().rescale_in_percent)
        print("For rescaling: {} STD of intra-subject CSA is {} mm^2 and COV is {} %".format(atrophy, round(intra_sub_std, 3), round(intra_sub_cov, 3)))
        # fill results dataframe
        df_results.loc[atrophy, 'intra_sub_std'] = intra_sub_std
        df_results.loc[atrophy, 'intra_sub_cov'] = intra_sub_cov
    print('\n')



def sample_size(df, atrophy, conf, power, mean_control=None, mean_patient=None):
    """
    Calculate the minimum number of patients required to detect an atrophy of a given value (i.e. power analysis),
    ratio patients/control 1:1 and with the assumption that both samples have the same STD.
    ref: Suresh and Chandrashekara 2012. “Sample size estimation and power analysis for clinical research studies”
    doi: 10.4103/0974-1208.97779
    Example: sample_size(df_a, 0.95, 0.8, mean_control=None, mean_patient=None, atrophy=7.7)
    :param df: dataframe with csv files data
    :param atrophy: expected atrophy in mm^2. Example atrophy=7.7
    :param conf: Confidence level. Example 0.8
    :param power: Power level. Example 0.9
    :param mean_control: mean CSA value of control group
    :param mean_patient: mean csa value of patient group
    """
    z_score_dict = {'confidence_Level': [0.60, 0.70, 0.8, 0.85, 0.90, 0.95],
                    'z_value': [0.842, 1.04, 1.28, 1.44, 1.64, 1.96], }

    df_z = pd.DataFrame(z_score_dict)
    df_z = df_z.set_index('confidence_Level')
    std_sample = df.groupby('rescale').get_group(1).groupby('subject')['MEAN(area)'].mean().std()
    # Check if user input atrophy percentage or mean_control and mean_patient
    if atrophy:
        atrophy = atrophy
    elif mean_control is not None & mean_patient is not None:
        atrophy = abs(mean_control - mean_patient)
    else:
        print('input error: either input mean_control and mean_patient or atrophy in mm^2')
    # sample size calculation
    num_n = 2 * ((df_z.at[conf, 'z_value'] + df_z.at[power, 'z_value']) ** 2) * (std_sample ** 2)
    deno_n = (abs(atrophy)) ** 2
    sample = ceil(num_n / deno_n)
    print("Minimum sample size to detect an atrophy of {} mm² is: {}".format(atrophy, sample))
    print("With parameters: \n - STD: {} \n - power: {} %\n - significance: {} %\n - ratio 1:1 (patients/controls)".format(round(std_sample, 3), power*100, conf*100))


def add_gt_to_dataframe(df):
    """  Add theoretic CSA values (rX^2 * MEAN(area)) to dataframe
    :param df: dataframe with csv files data
    :return df: modified dataframe with added theoretic_csa
    """
    # iterate across different vertebrae levels and add ground truth values to each subject in dataframe
    # get CSA values without rescale
    csa_without_rescale = df.groupby('rescale').get_group(1).groupby(['basename']).mean()['MEAN(area)']
    df['theoretic_csa_c' + str(min_vert) + '_c' + str(max_vert)] = np.nan
    # iterate across rescaling coefficients
    for name, group in df.groupby('rescale'):
        # get group rescale value
        group = group.reset_index().set_index('basename')
        atrophy = group['rescale'].values[0]
        # iterate across dataframe subjects
        for subject in group.index.values:
            # if dataframe subject exist in csa_without_rescale register theoretic csa value in th_csa_cX_cY
            if subject in csa_without_rescale.index.values:
                group.loc[subject, 'theoretic_csa_c' + str(min_vert) + '_c' + str(max_vert)] = csa_without_rescale.loc[subject] * (atrophy ** 2)
        df.loc[atrophy, 'theoretic_csa_c' + str(min_vert) + '_c' + str(max_vert)] = group['theoretic_csa_c' + str(min_vert) + '_c' + str(max_vert)].values
    return df


def add_stat_to_dataframe(df):
    """ add columns diff_csa, perc_diff_csa to dataframe
        :param df: original dataframe with csv file data
        :return df_a: modified dataframe with added gt_csa, diff_csa, perc_diff_csa for different vertebrae levels
        """
    # add column diff
    df['diff_c' + str(min_vert) + '_c' + str(max_vert)] = df['MEAN(area)'].sub(
        df['theoretic_csa_c' + str(min_vert) + '_c' + str(max_vert)]).abs()
    # add column perc_diff
    df['perc_diff_c' + str(min_vert) + '_c' + str(max_vert)] = 100 * df[
        'diff_c' + str(min_vert) + '_c' + str(max_vert)].div(df['theoretic_csa_c' + str(min_vert) + '_c' + str(max_vert)])
    # add column ratio
    for name, group in df.reset_index().groupby('rescale'):
        atrophy = group['rescale'].values[0]
        df.loc[atrophy, 'csa_without_rescale'] = df.groupby('rescale').get_group(1)['MEAN(area)'].values
    df['ratio_c' + str(min_vert) + '_c' + str(max_vert)] = 100 - 100 * df['MEAN(area)'].div(df['csa_without_rescale'])
    return df


def main():
    """
    main function, gather stats and call plots
    """
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    path_results = os.path.join(os.getcwd(), arguments.i)

    # aggregate all csv results files
    concatenate_csv_files(path_results)
    vertlevels_input = arguments.l
    path_output = arguments.o

    # read data
    data = pd.read_csv(os.path.join(path_results, r'csa_all.csv'), decimal=".")

    # create a dataframe from the csv files
    df_vert = pd.DataFrame(data)
    pd.set_option('display.max_rows', None)

    # fetch parameters from config.yaml file
    config_param = yaml_parser(arguments.config)

    # add column 'basename' to dataframe
    df_vert['basename'] = list(
        (os.path.basename(path).split('_r')[0] + '_' + os.path.basename(path).split('_')[3].split('.nii.gz')[0]) for
        path in df_vert['Filename'])

    # verify if vertlevels of interest were given in input by user
    if vertlevels_input is None:
        vertlevels = list(set(df_vert['VertLevel'].values))
    elif vertlevels_input:
        vertlevels = list(map(int, vertlevels_input))
        if not all(elem in set(list(df_vert['VertLevel'].values)) for elem in vertlevels):
            raise ValueError("\nInput vertebral levels '{}' do not exist in csv files".format(vertlevels))
    # register vertebrae levels of interest (Default: all vertebrae levels in csv files)
    # diff_vert = np.setdiff1d(list(set(df_vert['VertLevel'].values)), list(vertlevels))
    print("Stats are averaged across vertebral levels: {}".format(vertlevels))

    # Create new dataframe with only selected vertebral levels
    df = df_vert[df_vert['VertLevel'].isin(vertlevels)]

    # Drop column VertLevel (no more used)
    df = df.drop('VertLevel', 1)

    # Average values across levels, for each subject
    df.groupby(['rescale', 'basename']).mean()

    # Add column "transfo" to DF

    # TODO update code below

    df_sub = pd.DataFrame
    df_sub['']

    # df = df_vert.set_index('VertLevel').drop(index=diff_vert).groupby(['rescale', 'basename']).mean()
    # df = add_gt_to_dataframe(df, max_vert, min_vert)

    # Dataframe column additions diff, perc_diff for different vertebrae levels
    df = add_stat_to_dataframe(df, max_vert, min_vert)
    # add column 'subject' to dataframe
    df['subject'] = list(tf.split('_T')[0] for tf in df.reset_index()['basename'])
    df['rescale_in_percent'] = 100 - (df.reset_index()['rescale'] ** 2).values * 100
    df['rescale_in_percent'] = df['rescale_in_percent'].round(2)

    # Create results dataframe
    rescaling = sorted(set(df.reset_index()['rescale_in_percent'].values))
    data_results = {'rescaling':rescaling}
    df_results = pd.DataFrame(data_results).set_index('rescaling')

    # Print mean CSA for each rescaling
    print("\n====================mean==========================\n")
    inter_subject_mean(df, df_results)

    # compute sample size
    print("==================sample_size======================\n")
    # configuration parameters can be modified in config.yaml file
    atrophy = config_param['stats']['sample_size']['atrophy_sample']
    # conf = confidence level
    conf = config_param['stats']['sample_size']['conf']
    # power = power level
    power = config_param['stats']['sample_size']['power']
    sample_size(df, atrophy, conf, power, mean_control=None, mean_patient=None)

    # display number of subjects in test (multiple transformations of the same subjects are considered different)
    print("\n=================number subjects=======================\n")
    df['subject'] = list(tf.split('_T')[0] for tf in df.reset_index()['basename'])
    for r in rescaling:
        number_sub = df.groupby('subject')['MEAN(area)'].mean().count()
        number_tf = df.groupby(['rescale', 'subject'])['MEAN(area)'].count().values[0]
        print('For rescaling ' + str(r) + ' number of subjects is ' + str(number_sub) +
              ' and number of transformations per subject ' + str(number_tf))
        # fill results dataframe
        df_results.loc[r, 'number_sub'] = number_sub
        df_results.loc[r, 'number_tf'] = number_tf

    # compute STD for different vertebrae levels
    print("\n=======================inter_subject_std============================\n")
    inter_sub_std = inter_subject_std(df, df_results)
    print("\n===================intra_subject_std=========================\n")
    intra_subject_std(df, df_results)
    # register results dataframe in a csv file
    df_results.to_csv(os.path.join(path_results, r'stats_results.csv'))
    # round dataframe
    df = df.round(2)

    # plot graph if verbose is present
    if arguments.fig:
        if not os.path.isdir(arguments.o):
            os.makedirs(arguments.o)

        # plot percentage difference between simulated atrophy and ground truth atrophy
        column_perc_diff = [i for i in df.columns if 'perc_diff' in i]
        plot_perc_err(df, column_perc_diff, df_results, path_output)

        # boxplot CSA across different rescaling values
        boxplot_csa(df, path_output)

        # boxplot of atrophy across different rescaling values
        column_ratio = [j for j in df.columns if 'ratio' in j]
        boxplot_atrophy(df, column_ratio, path_output)

        # plot minimum number of patients required to detect an atrophy of a given value
        # z_conf = z_score for confidence level,
        # z_power = z_score for power level,
        # std = STD of subjects without rescaling CSA values
        # mean_csa =  mean CSA value of subjects without rescaling
        plot_sample_size(z_conf=1.96, z_power=(0.84, 1.282), std=df_results.loc[0, 'inter_sub_std'], mean_csa=df_results.loc[0, 'inter_sub_mean'], path_output=path_output)
        print('\nfigures have been plotted in dataset')


# Run
#########################################################################################
if __name__ == "__main__":
    main()
