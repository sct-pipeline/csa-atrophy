#!/bin/bash
#
# Process data. This script should be run within the subject's folder.
#
# Usage:
#   ./process_data.sh <SUBJECT> <FILEPARAM>
#
# Example:
#   ./process_data.sh sub-03
#
# Author: Julien Cohen-Adad, Paul Bautin

# The following global variables are retrieved from parameters.sh but could be
# overwritten here:


# Uncomment for full verbose
set -v

# Immediately exit if error
set -e

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Retrieve input params
SUBJECT=$1
#FILEPARAM=$2


# FUNCTIONS
# ==============================================================================
# Check if manual label already exists. If it does, copy it locally. If it does
# not, perform labeling.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  local scale="$3"
  # Update global variable with segmentation file name
  FILELABEL="${file}_labels"
  if [ -e "${PATH_SEGMANUAL}/${file}_labels-manual.nii.gz" ]; then
    echo "Found manual label: ${PATH_SEGMANUAL}/${file}_labels-manual.nii.gz"
    rsync -avzh "${PATH_SEGMANUAL}/${file}_labels-manual.nii.gz" ${FILELABEL}.nii.gz
  else
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c t2 -scale-dist ${scale} -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # Create labels in the cord
    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz

    # If automatic labeling did not work, you can initialize with manual identification of C2-C3 disc:
    # sct_label_utils -i ${file}.nii.gz -create-viewer 3 -o label_vert.nii.gz -msg "Click at the posterior tip of inter-vertebral disc"
    # sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c t2 -scale-dist ${scale} -initlabel label_vert.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz

    # If manual labeling did not work, you can initialize with manual identification discs:
    #sct_label_utils -i ${file}.nii.gz -create-viewer 2,3,4,5,6 -o label_vert.nii.gz -msg "Click at the posterior tip of inter-vertebral disc"
    #sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c t2 -discfile label_vert.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz
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
    sct_deepseg_sc -i ${file}.nii.gz -c $contrast -qc ${PATH_QC} -qc-subject ${SUBJECT}
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
# define resampling coefficients (always keep value 1 for reference)
R_COEFS=(0.85 0.90 0.95 0.97 0.99 1)
# iterate resample on subject
for r_coef in ${R_COEFS[@]}; do
  if [ -d "anat_r${r_coef}" ]; then
    rm -r "anat_r${r_coef}"
    echo "anat_r${r_coef} allready exists: creating folder"
  fi
  if [ -f "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv" ]; then
    rm "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv"
    echo "csa_perlevel_${SUBJECT}_${r_coef}.csv already exists: overwriting current csv file"
  fi
  # rename anat to explicit resampling coefficient
  mv anat anat_r$r_coef
  cd anat_r${r_coef}
  for file_tf in ${SUBJECT}_T2w*.nii.gz; do
    #subject_tf = basename "file_tf"
    # Image homothetic rescaling
    echo ${file_tf%%.*}
    python3 ../../../affine_rescale.py -i ${file_tf%%.*}.nii.gz -r ${r_coef}
    # sct_resample -i ${SUBJECT}_T2w.nii.gz -o ${SUBJECT}_T2w_r${r_coef}.nii.gz -f ${r_coef}x${r_coef}x${r_coef}
    file_t2=${file_tf%%.*}_r${r_coef}
    # Segment spinal cord (only if it does not exist)
    segment_if_does_not_exist $file_t2 "t2"
    # name segmented file
    file_t2_seg=$FILESEG
    # Create labels in the cord at C3 and C5 cervical vertebral levels (only if it does not exist)
    label_if_does_not_exist $file_t2 $file_t2_seg ${r_coef}
    file_label=$FILELABEL
    # Compute average CSA between C2 and C5 levels (append across subjects)
    # sct_process_segmentation -i $file_t2_seg.nii.gz -vert 1:3 -vertfile ${file_t2_seg}_labeled.nii.gz -o $PATH_RESULTS/csa_${SUBJECT}_${r_coef}.csv -qc ${PATH_QC}
    sct_process_segmentation -i $file_t2_seg.nii.gz -vert 2:5 -perlevel 1 -vertfile ${file_t2_seg}_labeled.nii.gz -o $PATH_RESULTS/CSA_perlevel_${file_tf%%.*}_${r_coef}.csv -qc ${PATH_QC}
  done
  cd ../
  cp -r $PATH_DATA/${SUBJECT}/anat .
  # add files to check
  FILES_TO_CHECK+=(
  "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv"
  "$PATH_RESULTS/${SUBJECT}/anat_r${r_coef}/${file_t2_seg}.nii.gz"
  )

done


# Verify presence of output files and write log file if error
# =============================================================================
for file in ${FILES_TO_CHECK[@]}; do
  if [ ! -e $file ]; then
    echo "$file does not exist" >> $PATH_LOG/error.log
  fi
done
