#!/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# implement random transformations to mimic subject repositioning, for each rescaling,
#
# generate randomized transfo per subject, and freeze them in the repos,
# so that we can reproduce the results by inputting the frozen params in subsequent runs of the pipeline,
# shifting (±5 voxels in each direction),
# rotation (±10° in each direction),
#
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
# SCT source: spinalcordtoolbox/testing/create_test_data.py
#
# examle python affine_transfo.py -i <sub-amu01 sub-amu02>
# About the license: see the file LICENSE
###################################################################

import glob, os, sys
import math
from numpy.random import randn
import numpy as np
import argparse

import nibabel as nib

from scipy.ndimage import affine_transform


def get_parser():
    parser = argparse.ArgumentParser(
        description='apply random rotation and translation with values following a gaussian distritbution:',
        add_help=None,
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py"))

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        help="Input T2w images to be transformed specify 'all' to transform entire dataset",
        type=str,
        nargs="*"
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-h',
        help='Help',
        nargs="*"
    )
    optional.add_argument(
        '-o',
        help="Suffix for output file name.\nexample: '-i MYFILE.nii -o _t' would output MYFILE_t.nii",
        default='_t',
    )
    optional.add_argument(
        '-p',
        help='Path to subject image directory, default: ./data',
        default= 'data',
    )
    return parser


def random_values():
    """generate gaussian distribution random values to simulate subject repositioning
     :return angle_IS: angle of rotation around Superior/Inferior axis
     :return angle_PA: angle of rotation around Anterior/Superior axis
     :return angle_LR: angle of rotation around Right/Left axis
     :return shift_LR: value of shift along Left/Right axis
     :return shift_PA: value of shift along Anterior/Superior axis
     :return shift_IS: value of shift along Inferior/Superior axis
     """
    values = randn(6)
    # for 95% of subjects repositioning (2*std away)
    std_angle = 5 # rotation (±10° in each direction),
    std_shift = 2.5 # shifting (±5 voxels in each direction)
    angle_IS = std_angle*values[0]
    shift_IS = std_shift*values[1]

    angle_PA = std_angle*values[2]
    shift_PA = std_shift*values[3]

    angle_LR = std_angle*values[4]
    shift_LR = std_shift*values[5]
    return angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS


def get_image(img, angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS):
     """fetch nibabel image and calculate minimum padding necessary for rotation boundaries
     :param img: nibabel image
     :param angle_IS: angle of rotation around Superior/Inferior axis
     :param angle_PA: angle of rotation around Anterior/Superior axis
     :param angle_LR: angle of rotation around Right/Left axis
     :param shift_LR: value of shift along Left/Right axis
     :param shift_PA: value of shift along Anterior/Superior axis
     :param shift_IS: value of shift along Inferior/Superior axis
     :return data: image data with a padding
     :return min_pad: number of voxels added on each side of the image
     """
     # upload and pad image
     data = img.get_fdata()
     # pad image to avoid edge effect during rotation
     max_shift = np.max(np.abs((shift_LR, shift_PA, shift_IS)))
     max_axes = np.max(data.shape)
     max_angle = np.deg2rad(np.max(np.abs((angle_IS, angle_PA, angle_LR))))
     # estimate upper limit of largest increase with rotation
     increase_angle = (max_axes/2)*(np.sqrt(2)*np.sin(((np.pi/4)-max_angle))-1)
     increase_axes = np.sqrt(2*increase_angle**2)
     min_pad = math.ceil(increase_axes+np.sqrt(3*(max_shift**2)))
     # padding choices: {‘constant’, ‘edge’, ‘symmetric’, ‘reflect’, ‘wrap’}
     data = np.pad(data, min_pad, 'constant')
     print('min padding:',min_pad)
     print('data shape (image with padding):',data.shape)
     return data, min_pad


def transfo(angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS, data):
     """apply rotation and translation on image
     :param angle_IS: angle of rotation around Superior/Inferior axis
     :param angle_PA: angle of rotation around Anterior/Superior axis
     :param angle_LR: angle of rotation around Right/Left axis
     :param shift_LR: value of shift along Left/Right axis
     :param shift_PA: value of shift along Anterior/Superior axis
     :param shift_IS: value of shift along Inferior/Superior axis
     :param data: padded image data
     :return data: image data with a padding
     :return data_rot: returns image data after applied random rotation and translation
     """
     # print angles and shifts
     print('angles for rotation IS:',angle_IS,' PA:', angle_PA, ' LR:',angle_LR)
     print('number of pixel shift LR:',shift_LR,' PA:',shift_PA,' IS:', shift_IS)

     c_in=0.5*np.array(data.shape) # find center of data

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

     #rotation matrix around LR
     cos_gamma = np.cos(np.deg2rad(-angle_LR))
     sin_gamma = np.sin(np.deg2rad(-angle_LR))
     rotation_affine_LR = np.array([[1, 0, 0],
                                    [0, cos_gamma, -sin_gamma],
                                    [0, sin_gamma, cos_gamma]])
     affine_arr_rotIS_rotPA_rotLR = rotation_affine_LR.dot(affine_arr_rotIS_rotPA) # affine array for rotation auround IS, AP and RL

     print('rotation matrix: \n', affine_arr_rotIS_rotPA_rotLR)

     # offset to shift the center of the old grid to the center of the new new grid + random shift
     shift = c_in.dot(affine_arr_rotIS_rotPA_rotLR)-c_in - np.array([shift_LR, shift_PA, shift_IS])
     # resampling data
     data_shift_rot = affine_transform(data, affine_arr_rotIS_rotPA_rotLR, offset=shift, order=5)

     return data_shift_rot


def main(fname, suffix):
    """Main function, crop and save image"""
    name = os.path.basename(fname).split(fname)[0]
    path = os.path.join(os.getcwd(), fname) # get file path
    print(path.split('.nii.gz')[0])
    path_tf = os.path.join(path, path.split('.nii.gz')[0]+str(suffix)+'.nii.gz')# create new path to save data
    print(path_tf)
    if os.path.isfile(path_tf):
        os.remove(path_tf)
    img = nib.load(fname) # load image
    print('\n----------affine transformation subject: '+name+'------------')
    angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS = random_values()
    # nibabel data follows the RAS+ (Right, Anterior, Superior are in the ascending direction) convention,
    data, min_pad = get_image(img, angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS)
    data_shift_rot = transfo(angle_IS, angle_PA, angle_LR, shift_LR, shift_PA, shift_IS, data)
    # load data back to nifti format
    img_t = nib.Nifti1Image(data_shift_rot, img.affine)
    print('new image shape: ',img_t.shape)
    print('new image path: '+path_tf)
    img_t.to_filename(path_tf)


if __name__ == "__main__":
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    if arguments.h is None:
        for subject in arguments.i:
            if os.path.isdir(arguments.p+'/'+subject): # verify presence of subject in dataset
                path = glob.glob(arguments.p+'/'+str(subject)+'/anat/*T2w.nii.gz') # find subject in dataset
                for fnames in path:
                    main(fnames, arguments.o)
            elif subject == 'all':
                for fnames in glob.glob(arguments.p+'/*/*/*T2w.nii.gz'): # apply transformation on all subjects of dataset
                    main(fnames, arguments.o)
            else:
                print('error: '+subject+' is not a valide subject')
    else:
        parser.print_help()
