#!/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Rescale image by changing subject image affine ndarray
# 
#
# example: python csa_rescale_stat.py -i <results>
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
#
# About the license: see the file LICENSE.TXT
#########################################################################################
import glob, os, sys
from numpy.random import randn
import numpy as np
import argparse

import nibabel as nib

def get_parser():
    parser = argparse.ArgumentParser(
        description='apply isotropic rescaling to image:',
        add_help=None,
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py"))

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        help='path to T2w mri image',
    )
    mandatory.add_argument(
        "-r",
        required=True,
        help='rescaling coefficient',
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-h',
        help='Help',
        nargs="*"
    )
    return parser


#MAIN
############################################################
def main(fname, coef_r):
    print(coef_r)
    # iterate transformation for each subject,
    img = nib.load(fname) # load image
    data = img.get_fdata()# load data
    img.affine[:3,:3] = img.affine[:3,:3]*float(coef_r) # apply rescale
    img_t = nib.Nifti1Image(data, img.affine) # changge affine for data

    # save rescaled image
    fname_out = fname.split('.nii.gz')[0] + '_r'+str(coef_r)+'.nii.gz'
    nib.save(img_t, fname_out)

#RUN
############################################################
if __name__ == "__main__":
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    fname = arguments.i
    coef_r = arguments.r
    main(fname, coef_r)
