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
config_script=$3
# The following global variables are retrieved from config_script.yml file
n_transfo=$(yaml_parser -o n_transfo -i $config_script)
rescaling=$(yaml_parser -o rescaling -i $config_script)
R_COEFS=$(echo $rescaling | tr '[]' ' ' | tr ',' ' ' | tr "'" ' ')
contrast=$(yaml_parser -o contrast -i $config_script)
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
#=======================================================================
# Segment spinal cord (only if it does not exist) in dir anat
cd anat
echo ${SUBJECT}_${contrast_str}
echo $contrast
file_c=${SUBJECT}_${contrast_str}
segment_if_does_not_exist $file_c ${contrast}
# name segmented file
file_c_seg=${FILESEG}
# Label spinal cord (only if it does not exist) in dir anat
label_if_does_not_exist $file_c $file_c_seg $contrast $contrast_str
file_label=${file_c_seg}_labeled
cd ../

# iterate across rescaling
for r_coef in ${R_COEFS[@]}; do
  mkdir anat_r$r_coef
  cd anat_r${r_coef}

  # Rescale header of nifti file
  # TODO: pass variable to point to -config yml file
  # rescale nifti native image
  affine_rescale -i ../anat/${file_c}.nii.gz -r ${r_coef} -o ${file_c}_r${r_coef}.nii.gz
  file_c_r=${file_c}_r${r_coef}
  #rescale nifti segmented and labled image
  affine_rescale -i ../anat/${file_label}.nii.gz -r ${r_coef} -o ${SUBJECT}_${contrast_str}_r${r_coef}_seg_labeled.nii.gz
  file_label_c_r=${SUBJECT}_${contrast_str}_r${r_coef}_seg_labeled

  # create list of array to iterate on (e.g.: seq_transfo = 1 2 3 4 5 if n_transfo=5)
  seq_transfo=$(seq ${n_transfo})
  for i_transfo in ${seq_transfo[@]}; do
    # Image random transformation (rotation, translation). By default transformation values are taken from
    # "transfo_values.csv" file if it already exists.
    # We keep a transfo_values.csv file, so that after first pass of the pipeline and QC, if segmentations
    # need to be manually-corrected, we want the transformations to be the same for the 2nd pass of the pipeline.
    affine_transfo -i ${file_c_r}.nii.gz -transfo $PATH_RESULTS/transfo_values.csv -config ../../../../$config_script -o _t${i_transfo}
    file_c_r_t=${SUBJECT}_${contrast_str}_r${r_coef}_t${i_transfo}
    affine_transfo -i ${file_label_c_r}.nii.gz -transfo $PATH_RESULTS/transfo_values.csv -config ../../../../$config_script -o _t${i_transfo} -interpolation 0
    file_label_c_r_t=${file_label_c_r}_t${i_transfo}
    # Segment spinal cord (only if it does not exist)
    segment_if_does_not_exist ${file_c_r_t} ${contrast}
    # name segmented file
    file_c_r_t_seg=${FILESEG}
    # Compute average CSA between C2 and C5 levels (append across subjects)
    sct_process_segmentation -i $file_c_r_t_seg.nii.gz -vert 2:5 -perlevel 1 -vertfile $file_label_c_r_t.nii.gz -o $PATH_RESULTS/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv -qc ${PATH_QC}
    # add files to check
    FILES_TO_CHECK+=(
    "$PATH_RESULTS/csa_data/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv"
    "$PATH_RESULTS/${SUBJECT}/anat_r${r_coef}/${file_c_seg}.nii.gz"
    )
  done
  cd ../
done



# Verify presence of output files and write log file if error
# =============================================================================
for file in ${FILES_TO_CHECK[@]}; do
  if [ ! -e $file ]; then
    echo "$file does not exist" >> $PATH_LOG/error.log
  fi
done
