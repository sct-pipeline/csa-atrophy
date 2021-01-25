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

# TODO: add plot of STD across transfo

import pandas as pd
import numpy as np
import os
import argparse
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
        help='Path to config file, which contains parameters for the statistics and figures. Example: config_script.yml',
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
        path = os.path.join(path_results, file)
        if ".csv" in file and "csa" and "sub" in file:
            files.append(path)
    if not files:
        raise FileExistsError("Folder {} does not contain any results csv file.".format(path_results))
    #metrics = pd.concat(
       #[pd.read_csv(f).assign(rescale=os.path.basename(f).split('_')[4].split('.csv')[0]) for f in files])
    print("Concatenate csv files. This will take a few seconds...")
    metrics = pd.concat(pd.read_csv(f) for f in files)
    # output csv file in PATH_RESULTS
    metrics.to_csv(os.path.join(path_results, r'csa_all.csv'))


def yaml_parser(config_file):
    """parse config_script.yml file containing pipeline's parameters"""
    with open(config_file, 'r') as config_var:
        yaml = YAML(typ='safe')
        config_param = yaml.load(config_var)
    return config_param


def plot_perc_err(df, path_output):
    """plot percentage difference between measured rescaling and rescaling
    :param df: dataframe for computing stats across subject: df_rescale
    :param path_output: directory in which plot is saved
    """
    df = df.reset_index().set_index('rescale_area')
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    df.groupby('rescale_area')['mean_perc_error'].mean().plot(kind='bar', ax=axes[0], grid=True)
    axes[0].set_title('mean error function of area rescaling')
    axes[0].set_ylabel('error in %')
    df['std_perc_error'].plot(kind='bar', ax=axes[1], sharex=True, sharey=True, legend=False)
    axes[1].set_title('STD of error in function of area rescaling')
    plt.xlabel('area rescaling in %')
    plt.ylabel('STD in %')
    plt.grid()
    plt.tight_layout()
    output_file = os.path.join(path_output, "fig_err.png")
    plt.savefig(output_file)
    print("--> Created figure: {}".format(output_file))


def boxplot_csa(df, path_output):
    """boxplot CSA for different rescaling values
    :param df: dataframe with csv files data: df_sub
    :param path_output: directory in which plot is saved
    """
    # TODO: round xlabel
    # TODO: find a way to display ylabel title with superscript
    fig2 = plt.figure()
    # round rescale area
    df['rescale_area'] = round(df['rescale_area'], ndigits=0).astype(int)
    df.boxplot(column=['mean'], by='rescale_area', showmeans=True, meanline=True)
    plt.title('Boxplot of CSA in function of area rescaling')
    plt.suptitle("")
    plt.ylabel('CSA in mm^2')
    plt.xlabel('CSA scaling')
    output_file = os.path.join(path_output, "fig_boxplot_csa.png")
    plt.savefig(output_file)
    print("--> Created figure: {}".format(output_file))


def boxplot_atrophy(df, path_output):
    """boxplot error for different rescaling values
    :param df: dataframe for computing stats per subject: df_sub
    :param path_output: directory in which plot is saved
    """
    fig = plt.figure()
    # convert to percentage
    df['rescale_estimated'] = df['rescale_estimated'] * 100
    # round rescale area
    df['rescale_area'] = round(df['rescale_area'], ndigits=0).astype(int)
    df.boxplot(column='rescale_estimated', by='rescale_area', positions=sorted(set(df['rescale_area'].values)), showmeans=True, meanline=True, figsize=(10,8))
    min_rescale = min(df['rescale_area'].values)
    max_rescale = max(df['rescale_area'].values)
    plt.plot([min_rescale, max_rescale], [min_rescale, max_rescale], ls="--", c=".3")
    plt.title('Boxplot of measured atrophy in function of CSA scaling')
    plt.ylabel('Atrophied CSA in % of the un-rescaled CSA')
    plt.xlabel('CSA scaling')
    plt.axis('scaled')
    output_file = os.path.join(path_output, "fig_boxplot_atrophy.png")
    plt.savefig(output_file)
    print("--> Created figure: {}".format(output_file))


def plot_sample_size(z_conf, z_power, std_arr, mean_csa, path_output):
    """plot minimum number of patients required to detect an atrophy of a given value
    :param z_conf: z score for X % uncertainty. Example: z_conf=1.96
    :param z_power: z score for X % Power. Example: z_power=(0.84, 1.282)
    :param std_arr: STD around mean CSA of control subjects (without rescaling),
    CSA STD for atrophied subjects and control subjects are considered equal
    :param mean_csa: mean value of CSA to compute the atrophy percentage. Example: 80
    :param path_output: directory in which plot is saved
    """
    fig_sample, ax = plt.subplots()
    # data for plotting
    n_t1 = []
    n_t2 =[]
    for z_p in z_power:
        # x_axis values ranging from 1.5 to 8.0 mm^2
        atrophy = np.arange(1.5, 8.0, 0.05)
        # numerator of sample size equation T1w
        num_n_t1 = 2 * ((z_conf + z_p) ** 2) * (std_arr[0] ** 2)
        n_t1.append(num_n_t1 / (atrophy ** 2))
        # numerator of sample size equation T2w
        num_n_t2 = 2 * ((z_conf + z_p) ** 2) * (std_arr[1] ** 2)
        n_t2.append(num_n_t2 / (atrophy ** 2))
    # plot
    ax.plot(atrophy / mean_csa[0] * 100, n_t1[0], 'b.', markevery=3, markersize=10, label='t1 80% power')
    ax.plot(atrophy / mean_csa[0] * 100, n_t1[1], 'r.', markevery=3, markersize=10, label='t1 90% power')
    ax.plot(atrophy / mean_csa[1] * 100, n_t2[0], 'c', label='t2 80% power')
    ax.plot(atrophy / mean_csa[1] * 100, n_t2[1], 'm', label='t2 90% power')
    ax.set_ylabel('number of participants per group of study \n(patients or controls) with ratio 1:1')
    ax.set_xlabel('atrophy in %')
    plt.suptitle('minimum number of participants to detect an atrophy with 5% uncertainty', fontsize=16, y=1.05)
    ax.set_title('T1w (1.0 mm iso): SD = {} mm², mean CSA= {} mm² \nT2w (0.8 mm iso): SD = {} mm², mean CSA= {} mm²'.format(
            str(
                round(std_arr[0], 2)), str(mean_csa[0]), str(round(std_arr[1], 2)) , str(mean_csa[1])))
    ax.legend()
    ax.grid()
    output_file = os.path.join(path_output, "fig_min_subj.png")
    plt.savefig(output_file, bbox_inches='tight')
    print("--> Created figure: {}".format(output_file))


def error_function_of_csa(df, path_output):
    """Scatter plot of CSA in function of error %
    :param df: dataframe for computing stats per subject: df_sub
    :param path_output: directory in which plot is saved
    """
    fig, ax = plt.subplots(figsize=(7, 7))
    df['Normalized CSA in mm2'] = df['mean'].div(df['rescale'])
    df['error %'] = df['perc_error']
    # compute linear regression
    z = np.polyfit(x=df.loc[:, 'error %'], y=df.loc[:, 'Normalized CSA in mm2'], deg=1)
    p = np.poly1d(z)
    # plot
    df.plot.scatter(x='error %', y='Normalized CSA in mm2', c='rescale', colormap='viridis')
    min_err = min(df['error %'].values)
    max_err = max(df['error %'].values)
    plt.plot([min_err, max_err], [min_err*z[0]+z[1], max_err*z[0]+z[1]], ls="--", c=".3")
    plt.title('Normalized CSA in function of error %,\n linear regression: {}'.format(p))
    # plt.title("CSA in function of % error")
    ax.legend(loc='upper right')
    plt.grid()
    output_file = os.path.join(path_output, "fig_err_in_function_of_csa.png")
    plt.savefig(output_file, bbox_inches='tight')
    print("--> Created figure: {}".format(output_file))


def sample_size(df, config_param):
    """
    Calculate the minimum number of patients required to detect an atrophy of a given value (i.e. power analysis),
    ratio patients/control 1:1 and with the assumption that both samples have the same STD.
    ref: Suresh and Chandrashekara 2012. “Sample size estimation and power analysis for clinical research studies”
    doi: 10.4103/0974-1208.97779
    :param df: dataframe for computing stats across subject: df_rescale
    :param config_param: configuration parameters can be modified in config.yaml file. Example conf = 0.95
    :return sample_size: sample size for each rescaling
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
            sample_size.append('inf')
    return sample_size

def error_function_of_intra_cov(df, path_output):
    """Scatter plot of intra-subject COV in function of error %
    :param df: dataframe for computing stats per subject: df_sub
    :param path_output: directory in which plot is saved
    """
    fig, ax = plt.subplots(figsize=(7, 7))
    # remove rescale=1 because error=0
    df = df[df['rescale'] != 1]
    # compute linear regression
    z = np.polyfit(x=df.loc[:, 'perc_error'], y=df.loc[:, 'cov'], deg=1)
    p = np.poly1d(z)
    # plot
    df.plot.scatter(x='perc_error', y='cov', c='rescale', colormap='viridis')
    min_err = min(df['perc_error'].values)
    max_err = max(df['perc_error'].values)
    plt.plot([min_err, max_err], [min_err*z[0]+z[1], max_err*z[0]+z[1]], ls="--", c=".3")
    ax.set_xlabel('perc_error')
    ax.set_ylabel('cov')
    plt.title('COV in function of % error,\n linear regression: {}'.format(p))
    plt.title("COV in function of % error")
    ax.legend(loc='upper right')
    plt.grid()
    output_file = os.path.join(path_output, "fig_err_in_function_of_cov.png")
    plt.savefig(output_file, bbox_inches='tight')
    print("--> Created figure: {}".format(output_file))


def error_function_of_intra_cov_outlier(df, path_output):
    """Scatter plot of intra-subject COV in function of error % to identify outliers with either high error %
     or high intra-subject COV
    :param df: dataframe for computing stats per subject: df_sub
    :param path_output: directory in which plot is saved
    """
    fig, ax = plt.subplots(figsize=(7, 7))
    # identified outliers either high error % or high intra-subject COV of t1 images
    outliers_t1_all = ['sub-brnoUhb01', 'sub-brnoUhb02', 'sub-brnoUhb03', 'sub-brnoUhb06', 'sub-brnoUhb07', 'sub-brnoUhb08',
                       'sub-barcelona04', 'sub-barcelona06', 'sub-beijingPrisma03', 'sub-milan03', 'sub-oxfordFmrib01',
                       'sub-tokyo750w03', 'sub-pavia04']
    # identified t1 outliers remaining after subjects removed due to missing vertebrae levels missing CSA
    outliers_t1 = ['sub-brnoUhb01', 'sub-brnoUhb08', 'sub-milan03', 'sub-oxfordFmrib01', 'sub-cmrrb05',
                   'sub-tokyo750w03', 'sub-pavia04']
    # identified outliers either high error % or high intra-subject COV of t2 images
    outliers_t2_all = ['sub-tokyo750w02', 'sub-tokyo750w04', 'sub-tokyo750w06', 'sub-beijingVerio02',
                       'sub-beijingVerio03', 'sub-ubc02', 'sub-vuiisIngenia05']
    # identified t2 outliers remaining after subjects removed due to missing vertebrae levels missing CSA
    outliers_t2 = ['sub-tokyo750w02', 'sub-tokyo750w04', 'sub-tokyo750w06', 'sub-beijingVerio02', 'sub-beijingVerio03',
                   'sub-ubc02', 'sub-vuiisIngenia05']
    # remove rescale=1 because error=0
    df = df[df['rescale'] != 1]
    ax.scatter(df['perc_error'], df['cov'], color='tab:blue', label='others')
    # plot scatter outliers of t1w images
    for outlier_t1 in outliers_t1:
        df_t1 = df.groupby(['subject']).get_group(outlier_t1)
        ax.scatter(df_t1['perc_error'], df_t1['cov'], color='tab:red', label=outlier_t1)
        df_t1 = []
    # plot scatter outliers of t2w images
    for outlier_t2 in outliers_t2:
        df_t2 = df.groupby(['subject']).get_group(outlier_t2)
        ax.scatter(df_t2['perc_error'], df_t2['cov'], color='tab:olive', label=outlier_t2)
        df_t2 = []
    # plot
    ax.set_xlabel('perc_error')
    ax.set_ylabel('cov')
    plt.title("COV in function of error %")
    ax.legend(loc='upper right')
    plt.grid()
    # save image
    output_file = os.path.join(path_output, "fig_err_in_function_of_cov_outlier.png")
    plt.savefig(output_file, bbox_inches='tight')
    print("--> Created figure: {}".format(output_file))


def add_columns_df_sub(df):
    """  Add columns theoretic CSA values (rX^2 * MEAN(area)) and CSA without rescaling to dataframe
    :param df: dataframe for computing stats per subject: df_sub
    :return df: modified dataframe with added theoretic_csa and csa_without_rescale
    """
    # get CSA values without rescale
    df = df.set_index('rescale')
    csa_without_rescale = df.groupby('rescale').get_group(1)
    csa_without_rescale = csa_without_rescale.set_index('subject')
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
    path_results = os.path.abspath(os.path.expanduser(arguments.i))
    vertlevels_input = arguments.l
    path_output = os.path.abspath(arguments.o)

    # aggregate all csv results files
    concatenate_csv_files(path_results)

    # read data
    data = pd.read_csv(os.path.join(path_results, r'csa_all.csv'), decimal=".")

    # create a dataframe from the csv files
    df_vert = pd.DataFrame(data)
    pd.set_option('display.max_rows', None)

    # identify rows with missing values
    print("Remove rows with missing values...")
    lines_to_drop = df_vert[df_vert['MEAN(area)'] == 'None'].index
    df_vert['subject'] = list(sub.split('data_processed/')[1].split('/anat')[0] for sub in df_vert['Filename'])

    # remove rows with missing values
    df_vert = df_vert.drop(df_vert.index[lines_to_drop])
    df_vert['MEAN(area)'] = pd.to_numeric(df_vert['MEAN(area)'])
    print("  Rows removed: {}".format(lines_to_drop))

    # fetch parameters from config.yaml file
    config_param = yaml_parser(arguments.config)

    # add useful columns to dataframe
    df_vert['basename'] = list(os.path.basename(path).split('.nii.gz')[0] for path in df_vert['Filename'])
    df_vert['rescale'] = list(float(b.split('RPI_r_r')[1].split('_')[0]) for b in df_vert['basename'])
    df_vert['slices'] = list(int(slices.split(':')[1]) - int(slices.split(':')[0]) + 1 for slices in df_vert['Slice (I->S)'])

    # verify if vertlevels of interest were given in input by user
    if vertlevels_input is None:
        vertlevels = list(set(df_vert['VertLevel'].values))
    elif vertlevels_input:
        vertlevels = list(map(int, vertlevels_input))
        if not all(elem in set(list(df_vert['VertLevel'].values)) for elem in vertlevels):
            raise ValueError("\nInput vertebral levels '{}' do not exist in csv files".format(vertlevels))
    # register vertebrae levels of interest (Default: all vertebrae levels in csv files)
    print("Stats are averaged across vertebral levels: {}".format(vertlevels))

    # Create new dataframe with only selected vertebral levels
    df = df_vert[df_vert['VertLevel'].isin(vertlevels)]
    # Drop column VertLevel (no more used)
    df = df.drop('VertLevel', 1)
    # Average values across levels, for each subject
    df = df.groupby(['rescale', 'basename']).mean()
    # Reset index because the groupby assigned rescale and basename as the new indexes. We want to re-generate a
    # number-based index
    df = df.reset_index()
    df['subject'] = list(tf.split('_T')[0] for tf in df['basename'])
    df['transfo'] = list(tf.split('_t')[1].split('_seg')[0] for tf in df['basename'])
    # Sum number of slices across selected vertebrae
    df['num_slices'] = df_vert.groupby(['rescale', 'basename'])['slices'].sum().values
    df = df.drop('basename', 1)

    # Create dataframe for computing stats per subject: df_sub
    print("\n==================== subject_dataframe ==========================\n")
    df_sub = pd.DataFrame()
    # add necessary columns to df_sub dataframe
    df_sub['rescale'] = df.groupby(['rescale', 'subject']).mean().reset_index()['rescale']
    df_sub['rescale_area'] = 100 * (df.groupby(['rescale', 'subject']).mean().reset_index()['rescale'] ** 2)
    df_sub['subject'] = df.groupby(['rescale', 'subject']).mean().reset_index()['subject']
    df_sub['num_tf'] = df.groupby(['rescale', 'subject'])['transfo'].count().values
    df_sub['num_slices'] = df.groupby(['rescale', 'subject'])['num_slices'].mean().values
    # add stats to per subject dataframe
    df_sub['mean'] = df.groupby(['rescale', 'subject']).mean()['MEAN(area)'].values
    df_sub['std'] = df.groupby(['rescale', 'subject']).std()['MEAN(area)'].values
    df_sub['cov'] = df_sub['std'].div(df_sub['mean'])
    df_sub = add_columns_df_sub(df_sub)
    df_sub['rescale_estimated'] = df_sub['mean'].div(df_sub['csa_without_rescale'])
    df_sub['error'] = (df_sub['mean'] - df_sub['theoretic_csa']).abs()
    df_sub['perc_error'] = 100 * (df_sub['mean'] - df_sub['theoretic_csa']).abs().div(df_sub['theoretic_csa'])
    print(df_sub)
    # save dataframe in a csv file
    df_sub.to_csv(os.path.join(path_output, r'csa_sub.csv'))

    # Create dataframe for computing stats across subject: df_rescale
    print("\n==================== rescaling_dataframe ==========================\n")
    df_rescale = pd.DataFrame()
    df_rescale['rescale'] = df_sub.groupby(['rescale']).mean().reset_index()['rescale']
    df_rescale['rescale_area'] = df_sub.groupby('rescale_area').mean().reset_index()['rescale_area']
    df_rescale['mean_slices'] = df_sub.groupby(['rescale']).mean()['num_slices'].values
    df_rescale['std_slices'] = df_sub.groupby(['rescale']).std()['num_slices'].values
    df_rescale['num_sub'] = df_sub.groupby('rescale')['mean'].count().values
    df_rescale['mean_inter'] = df_sub.groupby('rescale').mean()['mean'].values
    df_rescale['std_intra'] = df_sub.groupby('rescale').mean()['std'].values
    df_rescale['cov_intra'] = df_sub.groupby('rescale').mean()['cov'].values
    df_rescale['std_inter'] = df_sub.groupby('rescale').std()['mean'].values
    df_rescale['cov_inter'] = df_sub.groupby('rescale').std()['mean'].div(
        df_sub.groupby('rescale').mean()['mean']).values
    df_rescale['mean_rescale_estimated'] = df_sub.groupby('rescale').mean()['rescale_estimated'].values
    df_rescale['std_rescale_estimated'] = df_sub.groupby('rescale').std()['rescale_estimated'].values
    df_rescale['mean_perc_error'] = df_sub.groupby('rescale').mean()['perc_error'].values
    df_rescale['mean_error'] = df_sub.groupby('rescale').mean()['error'].values
    df_rescale['std_perc_error'] = df_sub.groupby('rescale').std()['perc_error'].values
    df_rescale['sample_size'] = sample_size(df_rescale, config_param)
    print(df_rescale)
    # save dataframe in a csv file
    df_rescale.to_csv(os.path.join(path_output, r'csa_rescale.csv'))

    # plot graph if verbose is present
    if arguments.fig:
        if path_output:
            os.makedirs(path_output, exist_ok=True)

        # plot percentage difference between simulated atrophy and ground truth atrophy
        plot_perc_err(df_rescale, path_output)

        # boxplot CSA across different rescaling values
        boxplot_csa(df_sub, path_output)

        # boxplot of atrophy across different rescaling values
        boxplot_atrophy(df_sub, path_output)

        # plot minimum number of patients required to detect an atrophy of a given value
        # z_score for confidence level,
        z_score_confidence = config_param['fig']['sample_size']['conf']
        # z_score for power level,
        z_score_power = config_param['fig']['sample_size']['power']
        # std = STD of subjects without rescaling CSA values
        # mean_csa =  mean CSA value of subjects without rescaling
        plot_sample_size(z_conf=z_score_confidence , z_power=z_score_power , std_arr=[7.56 , 8.29] ,
                         mean_csa=[69.70 , 76.12] , path_output=path_output)
        # scatter plot of COV in function of error %
        error_function_of_intra_cov(df_sub, path_output=path_output)
        # scatter plot of COV in function of error % to identify outliers
        error_function_of_intra_cov_outlier(df_sub, path_output=path_output)
        # scatter plot of CSA in function of error %
        error_function_of_csa(df_sub, path_output=path_output)

if __name__ == "__main__":
    main()
