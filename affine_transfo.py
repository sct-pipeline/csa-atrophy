#!/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# implement random transformations to mimic subject repositioning, for each rescaling,
#
# generate randomized transformations per subject, and freeze them in the repos,
# so that results can be reproduced by using the frozen params in subsequent runs of the pipeline,
# Example:
# shifting ±5 voxels in each direction,
# rotation ±10° in each direction,
#
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
# SCT source: spinalcordtoolbox/testing/create_test_data.py
#
# example python affine_transfo.py -i <sub-amu01 sub-amu02>
# About the license: see the file LICENSE
###################################################################

# TODO: add feature to re-use transformation if csv file exists
# TODO: clarify what i_dir is. v
# TODO: save image as FLOAT32 v
# TODO: add flag interp to apply transfo to seg-labeled data --> in that case, interp=0 (nearestneighbor)
# TODO (less priority): check padding (seems unecessary)

import glob, os, sys
import math
from numpy.random import rand
import numpy as np
import argparse
import csv
import pandas as pd
import yaml

import nibabel as nib

from scipy.ndimage import affine_transform


def get_parser():
    parser = argparse.ArgumentParser(
        description='apply random rotation and translation with values following a gaussian distribution:',
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        help="Input nifti image to apply the transformation on.",
    )
    mandatory.add_argument(
        "-config",
        required=True,
        help="Path to yml config file that includes parameters for the transformation.",
    )
    mandatory.add_argument(
        '-transfo',
        help="Path to csv file that contains the transformation. If the transformation for the specific subject already "
             "exists, it will read the transformation instead of creating a new random one.",
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-o',
        help="Suffix for output file name.\nExample: '-i MYFILE.nii -o _t' would output MYFILE_t.nii",
        default='_t',
    )
    optional.add_argument(
        '-interpolation',
        help="Transformation spline interpolation order. Order must be in the range 0-5",
        default=5,
        type=int
    )
    return parser


def yaml_parser(config_file):
    """parse config.yaml file containing pipeline's parameters"""
    print(config_file)
    with open(config_file, 'r') as config_var:
        config_param = yaml.safe_load(config_var)
    return config_param


def random_values(df, subject_name, config_param):
    """generate gaussian distribution random values to simulate subject repositioning,
    transformation values are kept in a pandas dataframe that is afterwards converted to a csv file
     :return angle_IS: angle of rotation around Inferior/Superior axis
     :return angle_PA: angle of rotation around Posterior/Anterior axis
     :return angle_LR: angle of rotation around Left/Right axis
     :return shift_LR: value of shift along Left/Right axis
     :return shift_PA: value of shift along Posterior/Anterior axis
     :return shift_IS: value of shift along Inferior/Superior axis
     """
    values = 2 * rand(6) - 1
    # transformations bounds can be modified in config_script.yml
    # random angle values are within ±angle_bound° around each axis,
    angle_bound = config_param['transfo']['bounds']['angle_bound']
    # random shift values are within ±shift_bound voxels in each direction,
    shift_bound = config_param['transfo']['bounds']['shift_bound']

    angle_IS = angle_bound * values[0]
    shift_IS = shift_bound * values[1]

    angle_PA = angle_bound * values[2]
    shift_PA = shift_bound * values[3]

    angle_LR = angle_bound * values[4]
    shift_LR = shift_bound * values[5]

    transfo_dict = {
        'subjects': subject_name,
        'angle_IS': angle_IS,
        'angle_PA': angle_PA,
        'angle_LR': angle_LR,
        'shift_LR': shift_LR,
        'shift_PA': shift_PA,
        'shift_IS': shift_IS
    }
    df = df.append(transfo_dict, ignore_index=True)
    return df, angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS


def get_image(img, angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS):
    """fetch nibabel image and calculate minimum padding necessary for transformed image,
    minimum padding is computed to contain maximum rotation and shift.
     :param img: nibabel image
     :param angle_IS: angle of rotation around Inferior/Superior axis
     :param angle_PA: angle of rotation around Posterior/Anterior axis
     :param angle_LR: angle of rotation around Left/Right axis
     :param shift_LR: value of shift along Left/Right axis
     :param shift_PA: value of shift along Posterior/Anterior axis
     :param shift_IS: value of shift along Inferior/Superior axis
     :return data: image data with a padding
     :return min_pad: number of voxels added on each side of the image
     """
    # upload and pad image to avoid edge overflow during transformation
    data = img.get_fdata()
    data = data.astype(np.int32)
    max_shift = np.max(np.abs((shift_LR, shift_PA, shift_IS)))
    max_axes = np.max(data.shape)
    max_angle = np.deg2rad(np.max(np.abs((angle_IS, angle_PA, angle_LR))))
    # estimate upper limit of largest frame increase due to transformation
    increase_angle = (max_axes / 2) * (np.sqrt(2) * np.sin(((np.pi / 4) - max_angle)) - 1)
    increase_axes = np.sqrt(2 * increase_angle ** 2)
    min_pad = math.ceil(increase_axes + np.sqrt(3 * (max_shift ** 2)))
    # padding choices: {‘constant’, ‘edge’, ‘symmetric’, ‘reflect’, ‘wrap’}
    data = np.pad(data, min_pad, 'constant')
    print('min padding:', min_pad)
    print('data shape (image with padding):', data.shape)
    return data, min_pad


def transfo(angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS, data, interpolation):
    """apply rotation and translation on image
     :param angle_IS: angle of rotation around Inferior/Superior axis
     :param angle_PA: angle of rotation around Posterior/Anterior axis
     :param angle_LR: angle of rotation around Left/Right axis
     :param shift_LR: value of shift along Left/Right axis
     :param shift_PA: value of shift along Posterior/Anterior axis
     :param shift_IS: value of shift along Inferior/Superior axis
     :param data: padded image data
     :return data: image data with a padding
     :return data_rot: return image data after random transformation
     """
    # print angles and shifts
    print('angles of rotation IS:', angle_IS, ' PA:', angle_PA, ' LR:', angle_LR)
    print('number of pixel shift LR:', shift_LR, ' PA:', shift_PA, ' IS:', shift_IS)

    # find center of data
    c_in = 0.5 * np.array(data.shape)

    # rotation matrix around IS
    cos_theta = np.cos(np.deg2rad(-angle_IS))
    sin_theta = np.sin(np.deg2rad(-angle_IS))
    rotation_affine_IS = np.array([[cos_theta, -sin_theta, 0],
                                   [sin_theta, cos_theta, 0],
                                   [0, 0, 1]])
    affine_arr_rotIS = rotation_affine_IS.dot(np.eye(3))

    # rotation matrix around PA
    cos_fi = np.cos(np.deg2rad(-angle_PA))
    sin_fi = np.sin(np.deg2rad(-angle_PA))
    rotation_affine_PA = np.array([[cos_fi, 0, sin_fi],
                                   [0, 1, 0],
                                   [-sin_fi, 0, cos_fi]])
    affine_arr_rotIS_rotPA = rotation_affine_PA.dot(affine_arr_rotIS)

    # rotation matrix around LR
    cos_gamma = np.cos(np.deg2rad(-angle_LR))
    sin_gamma = np.sin(np.deg2rad(-angle_LR))
    rotation_affine_LR = np.array([[1, 0, 0],
                                   [0, cos_gamma, -sin_gamma],
                                   [0, sin_gamma, cos_gamma]])
    # affine array for rotation around IS, AP and RL
    affine_arr_rotIS_rotPA_rotLR = rotation_affine_LR.dot(affine_arr_rotIS_rotPA)

    print('rotation matrix: \n', affine_arr_rotIS_rotPA_rotLR)

    # offset to shift the center of the old grid to the center of the new grid + random shift
    shift = c_in.dot(affine_arr_rotIS_rotPA_rotLR) - c_in - np.array([shift_LR, shift_PA, shift_IS])
    # resampling data
    # TODO; check if order=3 is much faster
    data_shift_rot = affine_transform(data, affine_arr_rotIS_rotPA_rotLR, offset=shift, order=interpolation)

    return data_shift_rot


def main():
    """Main function, transform and save image"""
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    suffix = arguments.o
    interpolation = arguments.interpolation
    # fetch parameters from config_script.yml file
    path_config_file = arguments.config
    config_param = yaml_parser(path_config_file)

    # Images of selected subject chosen by user in command line instructions, are copied and transformed
    # if a csv file containing transformation values does not yet exist it will be created,
    # otherwise command line input csv file is used and a new subject is added in a new row of a pandas dataframe
    if os.path.isfile(arguments.transfo):
        df = pd.read_csv(arguments.transfo, delimiter=',')
        print(df)
    else:
        print('error', arguments.transfo, ' is not present in current directory')
        print('creating new file named ', arguments.transfo)
        transfo_column_name = ['subjects', 'angle_IS', 'angle_PA', 'angle_LR', 'shift_LR', 'shift_PA', 'shift_IS']
        df = pd.DataFrame(columns=transfo_column_name)

    # transformations are applied on selected subject
    fname_path = arguments.i
    fname = os.path.basename(fname_path)
    if fname_path:
        # create new path to save data
        path_tf = fname_path.split('.nii.gz')[0] + str(suffix) + '.nii.gz'
        subject = os.path.basename(path_tf).split('.nii.gz')[0]
        # if path to transformed file already exists overwrite
        if os.path.isfile(path_tf):
            os.remove(path_tf)
        # load image
        img = nib.load(fname_path)
        print('\n----------affine transformation subject: ' + fname + '------------')
        # check if the subject is not already in dataframe, otherwise use dataframe subject values
        if subject not in df['subjects'].values:
            df, angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS = random_values(df, subject, config_param)
        else:
            print(df.set_index('subjects').loc[subject].values)
            angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS = df.set_index('subjects').loc[
                subject].values

        # nibabel data follows the RAS+ convention (Right, Anterior, Superior in the ascending direction)
        data, min_pad = get_image(img, angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS)
        data_shift_rot = transfo(angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS, data, interpolation)
        # load data back to nifti format
        img_t = nib.Nifti1Image(data_shift_rot, img.affine)
        print('new image shape: ', img_t.shape)
        print('new image path: ' + path_tf)
        img_t.to_filename(path_tf)
        # raise output error if the subject does not exist
    else:
        print('error: ' + fname_path + ' is not a valid subject')
    df.set_index('subjects').to_csv(arguments.transfo)


if __name__ == "__main__":
    main()
