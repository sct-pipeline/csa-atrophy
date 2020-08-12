#!/bin/bash
#
# Process data. This script should be run within the dataset folder.
#
# Usage:
#   ./process_data.sh <SUBJECT>
#
# Example:
#   ./process_data.sh sub-03
#
# Author: Julien Cohen-Adad, Paul Bautin
###################################################

# Uncomment for full verbose
set -v

# Immediately exit if error
set -e

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Global variables
SUBJECT=$1
# The following global variables are retrieved from config_script.yml file
n_transfo=$(yaml_parser -o n_transfo -i config_script.yml)
rescaling=$(yaml_parser -o rescaling -i config_script.yml)
R_COEFS=$(echo $rescaling | tr '[]' ' ' | tr ',' ' ' | tr "'" ' ')
contrast=$(yaml_parser -o contrast -i config_script.yml)
if [ $contrast == "t2" ]; then
  contrast_str="T2w"
fi
if [ $contrast == "t1" ]; then
  contrast_str="T1w"
fi


# FUNCTIONS
# ==============================================================================
# Check if manual label already exists. If it does, copy it locally. If it does
# not, perform automatic labeling.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  # TODO: this function should ONLY be applied for the scale=1 stage.
  local contrast="$3"
  local contrast_str="$4"
  # Update global variable with segmentation file name
  FILELABEL="${file}_labels"
  FILELABELMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELABEL}-manual.nii.gz"
  if [ -e "$FILELABELMANUAL" ]; then
    echo "manual labeled file was found: $FILELABELMANUAL"
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -discfile ${file}_labels-manual.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}

    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz
  else
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # Create labels in the Spinal Cord
    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}

# Check if manual segmentation already exists. If it does, copy it locally. If
# it does not, perform seg.
segment_if_does_not_exist(){
  local file="$1"
  local contrast="$2"
  # Update global variable with segmentation file name
  FILESEG="${file}_seg"
  FILESEGMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILESEG}-manual.nii.gz"
  if [ -e $FILESEGMANUAL ]; then
    echo "Found! Using manual segmentation."
    rsync -avzh $FILESEGMANUAL ${FILESEG}.nii.gz
    sct_qc -i ${file}.nii.gz -s ${FILESEG}.nii.gz -p sct_deepseg_sc -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    # Segment spinal cord
    sct_deepseg_sc -i ${file}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}


# SCRIPT STARTS HERE
# ==============================================================================
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to results folder, where most of the outputs will be located
#TODO: replace PATH_RESULTS by $PATH_DATA_PROCESSED
cd $PATH_DATA_PROCESSED
# Copy source images
cp -r $PATH_DATA/${SUBJECT} .
cd $SUBJECT
# we don't need dwi data, so let's remove it
rm -r dwi

# Image analysis
#=============================================================================
# seg cord on anat
# TODO
# label cord on anat
# TODO

# iterate across rescaling
for r_coef in ${R_COEFS[@]}; do
  # if data already exists, remove it
  if [ -d "anat_r${r_coef}" ]; then
    echo "anat_r${r_coef} already exists: removing folder..."
    rm -r "anat_r${r_coef}"
  fi
  if [ -f "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv" ]; then
    echo "csa_perlevel_${SUBJECT}_${r_coef}.csv already exists: remove it..."
    rm "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv"
  fi

  cp -r anat anat_r$r_coef
  cd anat_r${r_coef}

  # Rescale header of nifti file
  affine_rescale -i ${SUBJECT}_${contrast_str}.nii.gz -r ${r_coef}

  # create list of array to iterate on (e.g.: seq_transfo = 1 2 3 4 5 if n_transfo=5)
  seq_transfo=$(seq ${n_transfo})
  for i_transfo in ${seq_transfo[@]}; do
    # Image random transformation (rotation, translation). By default transformation values are taken from
    # "transfo_values.csv" file if it already exists.
    # We keep a transfo_values.csv file, so that after first pass of the pipeline and QC, if segmentations
    # need to be manually-corrected, we want the transformations to be the same for the 2nd pass of the pipeline.
    affine_transfo -i ${SUBJECT}_${contrast_str}_r${r_coef}.nii.gz -i_dir $PATH_RESULTS -o _t${i_transfo} -o_file "$PATH_RESULTS"/transfo_values.csv
    file_c=${SUBJECT}_${contrast_str}_r${r_coef}_t${i_transfo}
    # Segment spinal cord (only if it does not exist)
    segment_if_does_not_exist ${file_c} ${contrast}
    # name segmented file
    file_c_seg=${FILESEG}
    # Rescale and apply transformation on the reference label under anat/
    # TODO

    # TODO: remove gt
    if [ $r_coef == "gt" ] && [ ${i_transfo} == 1 ];then
      label_if_does_not_exist $file_c $file_c_seg $contrast $contrast_str
      cp ${file_c_seg}_labeled.nii.gz ${PATH_RESULTS}/${SUBJECT}
      mv ${PATH_RESULTS}/${SUBJECT}/${file_c_seg}_labeled.nii.gz ${PATH_RESULTS}/${SUBJECT}/${file_c_seg%%_rgt*}.nii.gz
      path_label=${PATH_RESULTS}/${SUBJECT}/${file_c_seg%%_rgt*}
    else
      echo "---------------------"$path_label
      affine_rescale -i ${path_label}.nii.gz -r ${r_coef}
      affine_transfo -i ${path_label}_r${r_coef}.nii.gz -i_dir $PATH_RESULTS -o _t${i_transfo} -o_file "$PATH_DATA_PROCESSED"/transfo_values.csv
      mv ${path_label}_r${r_coef}_t${i_transfo}.nii.gz ${path_label}_r${r_coef}_t${i_transfo}_seg_labeled.nii.gz
      cp ${path_label}_r${r_coef}_t${i_transfo}_seg_labeled.nii.gz ${PATH_RESULTS}/${SUBJECT}/anat_r${r_coef}/${file_c_seg}_labeled.nii.gz
    fi

    # Compute average CSA between C2 and C5 levels (append across subjects)
    sct_process_segmentation -i $file_c_seg.nii.gz -vert 2:5 -perlevel 1 -vertfile ${file_c_seg}_labeled.nii.gz -o $PATH_DATA_PROCESSED/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv -qc ${PATH_QC}
    # add files to check
    FILES_TO_CHECK+=(
    "$PATH_RESULTS/csa_data/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv"
    "$PATH_RESULTS/${SUBJECT}/anat_r${r_coef}/${file_c_seg}.nii.gz"
    )
  done
  cd ../
  cp -r $PATH_DATA/${SUBJECT}/anat .
done



# Verify presence of output files and write log file if error
# =============================================================================
for file in ${FILES_TO_CHECK[@]}; do
  if [ ! -e $file ]; then
    echo "$file does not exist" >> $PATH_LOG/error.log
  fi
done
