from __future__ import division

# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# detect Ponto-Medullary Junction location and pass parameters to process_data.sh
#
# example: python get_pmj.py -i sub-amu01_T1w_pmj.nii.gz -o nx
# ---------------------------------------------------------------------------------------
# Authors: Paul Bautin
#
# About the license: see the file LICENSE
#########################################################################################

import os
import argparse
import nibabel as nib
import numpy as np

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
        help='Path to image with Ponto-Medullary Junction label (e.g. "csa_atrophy_results/data_processed/sub-amu01/anat/sub-amu01_T1w_pmj.nii.gz")',
    )
    mandatory.add_argument(
        "-o",
        required=True,
        help='Return desired image parameter(e.g. image shape: nx)',
    )
    return parser

# check for QC, if problem, do it manually, push to the dataset, and run a 2nd pass
# run a small python script that will fetch the `z_pmj` coordinate of the label
# note: you can also run: sct_label_utils -i IMAGE_pmj.nii.gz -display
# fetch n_x dim
# Crop image


def main():
    """
    main function
    """
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    img_cropped = nib.load(arguments.i)
    data_cropped = img_cropped.get_fdata()
    # first 3 values of shape respectively refer to space x, y and z
    nx, ny, nz = img_cropped.shape
    # Ponto-Medullary Junction location
    x_pmj, y_pmj, z_pmj = np.where(data_cropped)
    # dictionnary to pass elements to process_data.sh
    dict_param = {
                  'nx' : nx,
                  'ny' : ny,
                  'nz' : nz,
                  'x_pmj': x_pmj[0],
                  'y_pmj': y_pmj[0],
                  'z_pmj': z_pmj[0],
    }
    print(dict_param[arguments.o])


# Run
#########################################################################################
if __name__ == "__main__":
    main()
