#!/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Rescale image by changing subject image affine ndarray,
# rescaled image is saved with suffix _r + rescaling factor. Example: <input file name>_r0.90.nii.gz
#
# Example: python affine_rescale -i <data/sub-amu01/anat/sub-amu01_T2w.nii.gz> -r 0.9
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
#
# About the license: see the file LICENSE
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
        help='path to T2w MRI data',
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
def main():
    """Main function, isotropic rescale of input images according to coef_r"""

    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args(args=None if sys.argv[0:] else ['--help'])
    # fname is the name of input image
    fname = arguments.i
    # coef_r is image rescaling coefficient
    coef_r = arguments.r
    name_coef = arguments.r
    if coef_r == "gt":
        coef_r = 1

    # load image
    img = nib.load(fname)
    data = img.get_fdata()# load data
    img.affine[:3,:3] = img.affine[:3,:3]*float(coef_r) # apply rescale
    img_t = nib.Nifti1Image(data, img.affine) # change affine for data

    # save rescaled image
    fname_out = fname.split('.nii.gz')[0] + '_r'+str(name_coef)+'.nii.gz'
    nib.save(img_t, fname_out)

#RUN
############################################################
if __name__ == "__main__":
    main()
