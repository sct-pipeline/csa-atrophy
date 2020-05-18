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
#
# About the license: see the file LICENSE.TXT
###################################################################
import glob, os
import math
from numpy.random import randn
import numpy as np
import matplotlib.pyplot as plt

from skimage.transform import EuclideanTransform
from skimage.transform import rotate
from skimage.transform import warp

import nibabel as nib

from  scipy.ndimage import shift

##############################VALUES###############################
# generate Gaussian values
def random_values():
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
##############################IMAGE###############################
def get_image(img, angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS):
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

##############################ROTATION###############################
def transfo(angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS, data):
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

##############################MAIN###############################
def main():
    i=1
    # iterate transformation for each subject,
    for fname in glob.glob('data/*/*/*T2w.nii.gz'):
        fname = os.path.join(os.getcwd(), fname) # get file path
        img = nib.load(fname) # load image
        name = os.path.basename(fname)
        print('----------affine transformation subject '+name+'------------')
        angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS = random_values()
        # nibabel data follows the RAS+ (Right, Anterior, Superior are in the ascending direction) convention,
        data, min_pad = get_image(img, angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS)
        data_rot = transfo(angle_IS, angle_AP, angle_RL, shift_LR, shift_PA, shift_IS, data)
        # Crop image (to remove padding)
        data_crop = data_rot[min_pad:min_pad+img.shape[0], min_pad:img.shape[1]+min_pad, min_pad:img.shape[2]+min_pad]
        # load data back to nifti format
        img_t = nib.Nifti1Image(data_crop, img.affine)
        print('new image shape',img_t.shape)
        nib.save(img_t, fname)
        i += 1

##############################RUN###############################
if __name__ == "__main__":
    main()
