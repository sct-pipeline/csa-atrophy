#!/bin/bash
#
# Process data. This script should be run within the subject's folder.
#
# Usage:
#   ./process_data.sh <SUBJECT> <FILEPARAM>
#
# Example:
#   ./process_data.sh sub-03 parameters.sh
#
# Author: Julien Cohen-Adad

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
  # Update global variable with segmentation file name
  FILELABEL="${file}_labels"
  if [ -e "${PATH_SEGMANUAL}/${file}_labels-manual.nii.gz" ]; then
    echo "Found manual label: ${PATH_SEGMANUAL}/${file}_labels-manual.nii.gz"
    rsync -avzh "${PATH_SEGMANUAL}/${file}_labels-manual.nii.gz" ${FILELABEL}.nii.gz
  else
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c t2 -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # Create labels in the cord at C3 and C5 mid-vertebral levels
    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 3,5 -o ${FILELABEL}.nii.gz
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


# T2w resampling
# =============================================================================

# define resamplling coefficients
R_COEFS=(0.85 0.9 0.95 1)
# iterate resample on subject
for r_coef in ${R_COEFS[@]}; do
  # rename anat to explicit resampling coefficient
  mv anat anat_r$r_coef
  cd anat_r${r_coef}
  # Image homothetic rescaling
  sct_resample -i ${SUBJECT}_T2w.nii.gz -f ${r_coef}x${r_coef}x${r_coef}
  # rename file to explicit resampling
  mv ${SUBJECT}_T2w_r.nii.gz ${SUBJECT}_T2w_r${r_coef}.nii.gz
  file_t2=${SUBJECT}_T2w_r${r_coef}
  # Segment spinal cord (only if it does not exist)
  segment_if_does_not_exist $file_t2 "t2"
  # name segmented file
  file_t2_seg=$FILESEG
  # Create a close mask around the spinal cord for more accurate registration
  # sct_create_mask -i $file_t2_r.nii.gz -p centerline,$file_t2_seg_r.nii.gz -size 35mm -f cylinder -o mask_t2_r.nii.gz
  # Create labels in the cord at C1 and C3 upper cervical vertebral levels (only if it does not exist)
  label_if_does_not_exist $file_t2 $file_t2_seg
  file_label=$FILELABEL
  # Register to template
  sct_register_to_template -i $file_t2.nii.gz -s $file_t2_seg.nii.gz -l $file_label.nii.gz -c t2 -param step=1,type=seg,algo=centermassrot:step=2,type=im,algo=syn,iter=5,slicewise=1,metric=CC,smooth=0 -qc $PATH_QC
  # Warp template
  # Note: we don't need the white matter atlas at this point, therefore use flag "-a 0"
  sct_warp_template -d $file_t2.nii.gz -w warp_template2anat.nii.gz -a 0 -ofolder label_T2w -qc ${PATH_QC}
  # Compute average CSA between C1 and C3 levels (append across subjects)
  sct_process_segmentation -i $file_t2_seg.nii.gz -vert 1:3 -vertfile label_T2w/template/PAM50_levels.nii.gz -o $PATH_RESULTS/csa_r_${r_coef}.csv -append 1 -qc ${PATH_QC}
  # sct_process_segmentation -i $file_t2_seg_r.nii.gz -vert 1:3 -perslice 1 -vertfile label_T2w/template/PAM50_levels.nii.gz -o $PATH_RESULTS/CSA_perslice_r.csv -append 1 -qc ${PATH_QC}
  cd ../
  cp -r $PATH_DATA/${SUBJECT}/anat .
done


# Verify presence of output files and write log file if error
# =============================================================================
FILES_TO_CHECK=(
  "$file_t2_seg.nii.gz"
)
for file in ${FILES_TO_CHECK[@]}; do
  if [ ! -e $file ]; then
    echo "$SUBJECT/$file does not exist" >> $PATH_LOG/error.log
  fi
done
