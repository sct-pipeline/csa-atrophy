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

from skimage.transform import rotate

import nibabel as nib

from  scipy.ndimage import shift

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
        '-p',
        help='path to subject image directory',
        default='data',
    )
    return parser


def random_values():
    """generate gaussian distribution random values to simulate subject repositioning
     :return angle_IS: angle of rotation around Superior/Inferior axis
     :return angle_AP: angle of rotation around Anterior/Superior axis
     :return angle_RL: angle of rotation around Right/Left axis
     :return shift_RL: value of shift along Left/Right axis
     :return shift_PA: value of shift along Anterior/Superior axis
     :return shift_IS: value of shift along Inferior/Superior axis
     """
    values = randn(6)
    np.set_printoptions(precision=3, suppress=True)
    # for 95% of subjects repositioning (2*std away)
    std_angle = 5 # rotation (±10° in each direction),
    std_shift = 2.5 # shifting (±5 voxels in each direction)
    angle_IS = std_angle*values[0]
    shift_IS = std_shift*values[1]

    angle_AP = std_angle*values[2]
    shift_PA = std_shift*values[3]

    angle_RL = std_angle*values[4]
    shift_LR = std_shift*values[5]
    return angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS


def get_image(img, angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS):
     """fetch nibabel image and calculate minimum padding necessary for rotation boundaries
     :param img: nibabel image
     :param angle_IS: angle of rotation around Superior/Inferior axis
     :param angle_AP: angle of rotation around Anterior/Superior axis
     :param angle_RL: angle of rotation around Right/Left axis
     :param shift_RL: value of shift along Left/Right axis
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
     max_angle = np.deg2rad(np.max(np.abs((angle_IS, angle_AP, angle_RL))))
     # estimate upper limit of largest increase with rotation
     increase_angle = (max_axes/2)*(np.sqrt(2)*np.sin(((np.pi/4)-max_angle))-1)
     increase_axes = np.sqrt(2*increase_angle**2)
     min_pad = math.ceil(increase_axes+np.sqrt(3*(max_shift**2)))
     # padding choices: {‘constant’, ‘edge’, ‘symmetric’, ‘reflect’, ‘wrap’}
     data = np.pad(data, min_pad, 'constant')
     print('min padding:',min_pad)
     print('data shape (image with padding):',data.shape)
     return data, min_pad


def transfo(angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS, data):
     """apply rotation and translation on image
     :param angle_IS: angle of rotation around Superior/Inferior axis
     :param angle_AP: angle of rotation around Anterior/Superior axis
     :param angle_RL: angle of rotation around Right/Left axis
     :param shift_RL: value of shift along Left/Right axis
     :param shift_PA: value of shift along Anterior/Superior axis
     :param shift_IS: value of shift along Inferior/Superior axis
     :param data: padded image data     
     :return data: image data with a padding
     :return data_rot: returns image data after applied random rotation and translation
     """
     # print angles and shifts
     print('angles for rotation IS:',angle_IS,' RL:',angle_RL,' AP:', angle_AP)
     print('number of pixel shift LR:',shift_LR,' PA:',shift_PA,' IS:', shift_IS)

      # shift in RAS dimensions
     data_shiftRA = shift(data, np.array([shift_LR, shift_PA, shift_IS]))

      # ROTATION AROUND IS AXIS
     # rotate (in deg in counter-clockwise direction.), and re-grid using linear interpolation
     data_rotIS = rotate(data_shiftRA, angle_IS, resize=False, center=None, order=1, mode='constant', cval=0, clip=False, preserve_range=False)

      # ROTATION AROUND RL AXIS
     # Swap x-z axes (to make a rotation within y-z plane, as rotate will apply rotation on the first 2 dims)
     data_rotIS_swap = data_rotIS.swapaxes(0, 2)

      # rotate (in deg), and re-grid using linear interpolation
     data_rotIS_swap_rotRL = rotate(data_rotIS_swap, angle_RL, resize=False, center=None, order=1, mode='constant',cval=0, clip=False, preserve_range=False)
     # swap back
     data_rotIS_rotRL = data_rotIS_swap_rotRL.swapaxes(0, 2)

      # ROTATION AROUND AP AXIS
     # Swap y-z axes (to make a rotation within x-z plane)
     data_rotIS_rotRL_swap = data_rotIS_rotRL.swapaxes(1, 2)
     # rotate (in deg), and re-grid using linear interpolation
     data_rotIS_rotRL_swap_rotAP = rotate(data_rotIS_rotRL_swap, angle_AP, resize=False, center=None, order=1,
                                          mode='constant', cval=0, clip=False, preserve_range=False)

      # swap back
     data_rot = data_rotIS_rotRL_swap_rotAP.swapaxes(1, 2)
     return data_rot


def main(fname):
    """Main function, crop and save image"""
    # iterate transformation for each subject
    #for fname == fnames:
    name = os.path.basename(fname).split(fname)[0]
    path = os.path.join(os.getcwd(), fname) # get file path
    img = nib.load(fname) # load image
    print('\n----------affine transformation subject: '+name+'------------')
    angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS = random_values()
    # nibabel data follows the RAS+ (Right, Anterior, Superior are in the ascending direction) convention,
    data, min_pad = get_image(img, angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS)
    data_rot = transfo(angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS, data)
    # Crop image (to remove padding)
    data_crop = data_rot[min_pad:min_pad+img.shape[0], min_pad:img.shape[1]+min_pad, min_pad:img.shape[2]+min_pad]
    # load data back to nifti format
    img_t = nib.Nifti1Image(data_crop, img.affine)
    print('new image shape: ',img_t.shape)
    # create new path to save data
    new_path = path.split('.nii.gz')[0] + '-t.nii.gz'
    print('new image path: '+new_path)
    img_t.to_filename(new_path)


if __name__ == "__main__":
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    if arguments.h is None:
        for subject in arguments.i:
            if os.path.isdir(arguments.p+'/'+subject):
                path = glob.glob(arguments.p+'/'+str(subject)+'/anat/*T2w.nii.gz')
                for fnames in path:
                    main(fnames)
            elif subject == 'all':
                for fnames in glob.glob(arguments.p+'/*/*/*T2w.nii.gz'):
                    main(fnames)
            else:
                print('error: '+subject+' is not a valide subject')
    else:
        parser.print_help()
