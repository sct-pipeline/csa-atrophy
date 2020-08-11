#!/bin/bash
#
# Process data. This script should be run within the dataset folder.
#
# Usage:
#   ./process_data.sh <SUBJECT> <FILEPARAM>
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
# The following global variables are retrieved from config.yaml file
# set n_transfo to the desired number of transformed images of same subject for segmentation,
# n_transfo also represents the number of iterations of the transformation, segmentation and labeling process
n_transfo=$(yaml_parser -o n_transfo)
# define rescaling coefficients (always keep value 1 for reference)
rescaling=$(yaml_parser -o rescaling)
R_COEFS=$(echo $rescaling | tr '[]' ' ' | tr ',' ' ' | tr "'" ' ')
contrast=$(yaml_parser -o contrast)
if [ $contrast == "t2" ]; then
  contrast_str="T2w"
fi
if [ $contrast == "t1" ]; then
  contrast_str="T1w"
fi
path_output=$(yaml_parser -o path_output)


# FUNCTIONS
# ==============================================================================
# Check if manual label already exists. If it does, copy it locally. If it does
# not, perform automatic labeling.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  local scale="$3"

  # Update global variable with segmentation file name
  FILELABEL="${file}_labels"
  FILELABELMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELABEL}-manual.nii.gz"
  if [ -e "$FILELABELMANUAL" ]; then
    echo "manual labeled file was found: $FILELABELMANUAL"
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -discfile ${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${file}_labels-manual.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz
  else
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -scale-dist ${scale} -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # Create labels in the Spinal Cord
    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz
  fi
}

# Check if manual segmentation already exists. If it does, copy it locally. If
# it does not, perform seg.
segment_if_does_not_exist(){
  local file="$1"
  local contrast="$2"
  # Update global variable with segmentation file name
  FILESEG="${file}_seg"
  if [ -e "${PATH_SEGMANUAL}/${FILESEG}-manual.nii.gz" ]; then
    echo "Found manual segmentation: ${PATH_SEGMANUAL}/${FILESEG}-manual.nii.gz"
    rsync -avzh "${PATH_SEGMANUAL}/${FILESEG}-manual.nii.gz" ${FILESEG}.nii.gz
    sct_qc -i ${file}.nii.gz -s ${FILESEG}.nii.gz -p sct_deepseg_sc -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    # Segment spinal cord
    sct_deepseg_sc -i ${file}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}

# SCRIPT STARTS HERE
# ==============================================================================
# Go to results folder, where most of the outputs will be located
cd $PATH_RESULTS
# Copy ###source images
cp -r $PATH_DATA/${SUBJECT} $PATH_RESULTS
cd $SUBJECT
rm -r dwi

# T2w resampling
#=============================================================================
# iterate rescaling and transformation on subject
for r_coef in ${R_COEFS[@]}; do
    if [ -d "anat_r${r_coef}" ]; then
   rm -r "anat_r${r_coef}"
   echo "anat_r${r_coef} already exists: creating folder"
 fi
 if [ -f "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv" ]; then
   rm "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv"
   echo "csa_perlevel_${SUBJECT}_${r_coef}.csv already exists: overwriting current csv file"
 fi
  # rename anat to explicit rescaling coefficient
  mv anat anat_r$r_coef
  cd anat_r${r_coef}

  seq_transfo=$(seq ${n_transfo})
  for i_transfo in ${seq_transfo[@]}; do
    # Image homothetic rescaling
    affine_rescale -i ${SUBJECT}_${contrast_str}.nii.gz -r ${r_coef}
    # Image random transformation (rotation, translation). By default transformation values are taken from
    # "transfo_values.csv" file if it already exists
    affine_transfo -i ${SUBJECT}_${contrast_str}_r${r_coef}.nii.gz -i_dir $path_output -o _t${i_transfo} -o_file "$PATH_DATA_PROCESSED"/transfo_values.csv

    file_t2=${SUBJECT}_${contrast_str}_r${r_coef}_t${i_transfo}
    # Segment spinal cord (only if it does not exist)
    segment_if_does_not_exist ${file_t2} ${contrast}
    # name segmented file
    file_t2_seg=${FILESEG}

    # Create labels in the cord, function uses by default labels file in directory seg_manual
    label_if_does_not_exist $file_t2 $file_t2_seg $R_COEFS $contrast
    file_label=$FILELABEL
    # Compute average CSA between C2 and C5 levels (append across subjects)
    sct_process_segmentation -i $file_t2_seg.nii.gz -vert 2:5 -perlevel 1 -vertfile ${file_t2_seg}_labeled.nii.gz -o $PATH_DATA_PROCESSED/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv -qc ${PATH_QC}
    # add files to check
    FILES_TO_CHECK+=(
    "$PATH_RESULTS/csa_data/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv"
    "$PATH_RESULTS/${SUBJECT}/anat_r${r_coef}/${file_t2_seg}.nii.gz"
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
