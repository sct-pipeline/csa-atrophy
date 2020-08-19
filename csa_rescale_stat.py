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


def sample_size(df, config_param):
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
    sample_size = []
    # configuration parameters can be modified in config.yaml file
    # conf = confidence level
    conf = config_param['stats']['sample_size']['conf']
    # power = power level
    power = config_param['stats']['sample_size']['power']
    z_score_dict = {'confidence_Level': [0.60, 0.70, 0.8, 0.85, 0.90, 0.95],
                    'z_value': [0.842, 1.04, 1.28, 1.44, 1.64, 1.96], }

    df_z = pd.DataFrame(z_score_dict)
    df_z = df_z.set_index('confidence_Level')
    for name, group in df.groupby('rescale'):
        std = group['std_inter'].values[0]
        mean_patient = group['mean_inter'].values[0]
        mean_control = df.groupby('rescale').get_group(1)['mean_inter'].values[0]
        atrophy = mean_control - mean_patient
        if atrophy != 0:
            num_n = 2 * ((df_z.at[conf, 'z_value'] + df_z.at[power, 'z_value']) ** 2) * (std ** 2)
            deno_n = (abs(atrophy)) ** 2
            sample_size.append(ceil(num_n / deno_n))
        else:
            print('hello')
            sample_size.append('inf')
    return sample_size
        #print("Minimum sample size to detect an atrophy of {} mm² is: {}".format(atrophy, sample))
        #print("With parameters: \n - STD: {} \n - power: {} %\n - significance: {} %\n - ratio 1:1 (patients/controls)".format(round(std_sample, 3), power*100, conf*100))


def add_gt_to_dataframe(df):
    """  Add theoretic CSA values (rX^2 * MEAN(area)) to dataframe
    :param df: dataframe with csv files data
    :return df: modified dataframe with added theoretic_csa
    """
    # iterate across different vertebrae levels and add ground truth values to each subject in dataframe
    # get CSA values without rescale
    df = df.set_index('rescale')
    csa_without_rescale = df.groupby('rescale').get_group(1)
    csa_without_rescale = csa_without_rescale.set_index('subject')
    df['theoretic_csa'] = np.nan
    # iterate across rescaling coefficients
    for rescale, group in df.groupby('rescale'):
        # get group rescale value
        group = group.reset_index().set_index('subject')
        # iterate across dataframe subjects
        for subject in group.index.values:
            # if dataframe subject exist in csa_without_rescale register theoretic csa value in th_csa_cX_cY
            if subject in csa_without_rescale.index.values:
                group.loc[subject, 'theoretic_csa'] = csa_without_rescale.loc[subject, 'mean'] * (rescale ** 2)
        df.loc[rescale, 'theoretic_csa'] = group['theoretic_csa'].values
        df.loc[rescale, 'csa_without_rescale'] = csa_without_rescale['mean'].values
    df = df.reset_index()
    return df

def main():
    """
    main function, gather stats and call plots
    """
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    path_results = os.path.join(os.getcwd(), arguments.i)
    vertlevels_input = arguments.l
    path_output = arguments.o

    # aggregate all csv results files
    concatenate_csv_files(path_results)

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
    df = df.groupby(['rescale', 'basename']).mean()
    df = df.reset_index()
    df['subject'] = list(tf.split('_T')[0] for tf in df['basename'])
    df['transfo'] = list(tf.split('w_')[1] for tf in df['basename'])
    df = df.drop('basename', 1)


    # Ceate dataframe for computing stats per subject: df_sub
    print("\n====================subject_dataframe==========================\n")
    df_sub = pd.DataFrame()
    # add necessary columns to df_sub dataframe
    df_sub['rescale'] = df.groupby(['rescale','subject']).mean().reset_index()['rescale']
    df_sub['rescale_area'] = 100 * (1 - df.groupby(['rescale', 'subject']).mean().reset_index()['rescale'] ** 2)
    df_sub['subject'] = df.groupby(['rescale', 'subject']).mean().reset_index()['subject']
    df_sub['num_tf'] = df.groupby(['rescale', 'subject'])['transfo'].count().values
    # add stats to per subject datframe
    df_sub['mean'] = df.groupby(['rescale', 'subject']).mean()['MEAN(area)'].values
    df_sub['std'] = df.groupby(['rescale', 'subject']).std()['MEAN(area)'].values
    df_sub['cov'] = 100 * df_sub['std'].div(df_sub['mean'])
    df_sub = add_gt_to_dataframe(df_sub)
    df_sub['rescale_estimated'] = 100 * df_sub['mean'].div(df_sub['csa_without_rescale'])
    df_sub['error'] = (df_sub['mean'] - df_sub['theoretic_csa']).abs()
    df_sub['perc_error'] = 100 * (df_sub['mean'] - df_sub['theoretic_csa']).abs().div(df_sub['mean'])
    print(round(df_sub, 2))


    # Ceate dataframe for computing stats across subject: df_rescale
    print("\n====================rescaling_dataframe==========================\n")
    df_rescale = pd.DataFrame()
    df_rescale['rescale'] = df_sub.groupby(['rescale']).mean().reset_index()['rescale']
    df_rescale['rescale_area'] = 100 * (1 - df_sub.groupby('rescale').mean().reset_index()['rescale'] ** 2)
    df_rescale['num_sub'] = df_sub.groupby('rescale')['mean'].count().values
    df_rescale['mean_inter'] = df_sub.groupby('rescale').mean()['mean'].values
    df_rescale['std_intra'] = df_sub.groupby('rescale').mean()['std'].values
    df_rescale['cov_intra'] = df_sub.groupby('rescale').mean()['cov'].values
    df_rescale['std_inter'] = df_sub.groupby('rescale').std()['mean'].values
    df_rescale['mean_rescale_estimated'] = df_sub.groupby('rescale').mean()['rescale_estimated'].values
    df_rescale['mean_rescale_estimated'] = df_sub.groupby('rescale').std()['rescale_estimated'].values
    df_rescale['mean_error'] = df_sub.groupby('rescale').mean()['error'].values
    df_rescale['std_error'] = df_sub.groupby('rescale').std()['error'].values
    df_rescale['sample_size'] = sample_size(df_rescale, config_param)
    print(round(df_rescale,2))


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
