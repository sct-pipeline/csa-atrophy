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

# TODO: crop data to save computation time

# Uncomment for full verbose
set -x

# Immediately exit if error
set -e

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# get starting time:
start=`date +%s`

# Global variables
SUBJECT=$1
config_script=$2
echo "SUBJECT: $SUBJECT"
echo "config_script: $config_script"

# The following global variables are retrieved from config_script.yml file
n_transfo=$(yaml_parser -o n_transfo -i $config_script)
rescaling=$(yaml_parser -o rescaling -i $config_script)
R_COEFS=$(echo $rescaling | tr '[]' ' ' | tr ',' ' ' | tr "'" ' ')
contrast=$(yaml_parser -o contrast -i $config_script)
# TODO: enable to input a list of contrast and loop across contrasts
transfo_file=$(yaml_parser -o transfo_file -i $config_script)
echo "transfo_file: $transfo_file"


# FUNCTIONS
# ==============================================================================
crop_image(){
  local file=$1
  local file_seg=$2
  local contrast=$3
  FILE_PMJ=${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${file}_pmj
  # Verify if a manually detected pmj is present
  if [ -e "${FILE_PMJ}-manual.nii.gz" ]; then
    echo "Found manual PMJ detection; file: ${FILE_PMJ}-manual.nii.gz"
    # parameters to crop image
    local nx=$(get_pmj -i ${FILE_PMJ}-manual.nii.gz -o nx)
    local z_pmj=$(get_pmj -i ${FILE_PMJ}-manual.nii.gz -o z_pmj)
    local x_min=$(get_pmj -i ${FILE_PMJ}-manual.nii.gz -o x_min)
    local x_max=$(get_pmj -i ${FILE_PMJ}-manual.nii.gz -o x_max)
    # crop original image and segmented image
    sct_crop_image -i ${file}.nii.gz -xmin ${x_min} -xmax ${x_max} -zmax ${z_pmj}
    sct_crop_image -i ${file_seg}.nii.gz -xmin ${x_min} -xmax ${x_max} -zmax ${z_pmj}
  else
    # Detect ponto-medullary junction automatically
    sct_detect_pmj -i ${file}.nii.gz -s ${file_seg}.nii.gz -c $contrast -qc ${PATH_QC}
    # parameters to crop image
    local nx=$(get_pmj -i ${file}_pmj.nii.gz -o nx)
    local z_pmj=$(get_pmj -i ${file}_pmj.nii.gz -o z_pmj)
    local x_min=$(get_pmj -i ${file}_pmj.nii.gz -o x_min)
    local x_max=$(get_pmj -i ${file}_pmj.nii.gz -o x_max)
    # crop original image and segmented image
    sct_crop_image -i ${file}.nii.gz -xmin ${x_min} -xmax ${x_max} -zmax ${z_pmj}
    sct_crop_image -i ${file_seg}.nii.gz -xmin ${x_min} -xmax ${x_max} -zmax ${z_pmj}
  fi
}


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
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -discfile ${FILELABELMANUAL} -qc ${PATH_QC} -qc-subject ${SUBJECT}

    sct_label_utils -i ${file_seg}_labeled.nii.gz -vert-body 0 -o ${FILELABEL}.nii.gz
  else
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
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
# Copy config files to output results folder
mkdir ${PATH_RESULTS}/$SUBJECT/
cp -u $config_script ${PATH_RESULTS}/$SUBJECT/
# Go to results folder, where most of the outputs will be located
cd $PATH_DATA_PROCESSED
# Copy source images
cp -r $PATH_DATA/${SUBJECT} .
cd $SUBJECT
# we don't need dwi data, so let's remove it
rm -r dwi


# Image analysis in folder anat
#=======================================================================
cd anat
# Reorient to RPI and resample file
if [ $contrast == "t2" ]; then
  contrast_str="T2w"
  file_c=${SUBJECT}_${contrast_str}
  # Reorient to RPI
  sct_image -i ${file_c}.nii.gz -setorient RPI -o ${file_c}_RPI.nii.gz
  file_c_o=${file_c}_RPI
  # Resample to 0.8mm iso
  sct_resample -i ${file_c}_RPI.nii.gz -mm 0.8x0.8x0.8 -o ${file_c_o}r.nii.gz
  file_c_or=${file_c_o}r
fi
if [ $contrast == "t1" ]; then
  contrast_str="T1w"
  file_c=${SUBJECT}_${contrast_str}
  # Reorient to RPI
  sct_image -i ${file_c}.nii.gz -setorient RPI -o ${file_c}_RPI.nii.gz
  file_c_o=${file_c}_RPI
  # Resample to 1mm iso
  sct_resample -i ${file_c}_RPI.nii.gz -mm 1x1x1 -o ${file_c_o}r.nii.gz
  file_c_or=${file_c_o}r
fi

# Segment spinal cord (only if it does not exist) in dir anat
segment_if_does_not_exist $file_c_or ${contrast}
# name segmented file
file_c_or_seg=${FILESEG}
# crop image
crop_image ${file_c_or} ${file_c_or_seg} ${contrast}
file_c_or_crop=${file_c_or}_crop
file_c_or_seg_crop=${file_c_or_seg}_crop
# Label spinal cord (only if it does not exist) in dir anat
label_if_does_not_exist $file_c_or_crop $file_c_or_seg_crop $contrast $contrast_str
file_c_or_crop_label=${file_c_or_seg_crop}_labeled
cd ../


# CSA measure iterated across rescaling in folder anat_r{r_coef}
# ====================================================================
for r_coef in ${R_COEFS[@]}; do
  # If directory exists (e.g. 2nd pass after QC and manual correction), we should remove it
  if [ -d "anat_r${r_coef}" ]; then
    rm -r "anat_r${r_coef}"
    echo "anat_r${r_coef} already exists: removing folder"
  fi
  if [ -f "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv" ]; then
    rm "$PATH_RESULTS/csa_perlevel_${SUBJECT}_${r_coef}.csv"
    echo "csa_perlevel_${SUBJECT}_${r_coef}.csv already exists: removing file"
  fi
  mkdir anat_r$r_coef
  cd anat_r${r_coef}

  # Rescale header of native nifti file
  file_c_or_crop_r=${file_c_or_crop}_r${r_coef}
  affine_rescale -i ../anat/${file_c_or_crop}.nii.gz -r ${r_coef} -o ${file_c_or_crop_r}.nii.gz
  # rescale labeled segmentation
  file_c_or_crop_label_r=${file_c_or_crop_r}_seg_labeled
  affine_rescale -i ../anat/${file_c_or_crop_label}.nii.gz -r ${r_coef} -o ${file_c_or_crop_label_r}.nii.gz

  # create list of array to iterate on (e.g.: seq_transfo = 1 2 3 4 5 if n_transfo=5)
  seq_transfo=$(seq ${n_transfo})
  # Iterate across transformations
  for i_transfo in ${seq_transfo[@]}; do
    # Image random transformation (rotation, translation). By default transformation values are taken from
    # "transfo_values.csv" file if it already exists.
    # We keep a transfo_values.csv file, so that after first pass of the pipeline and QC, if segmentations
    # need to be manually-corrected, we want the transformations to be the same for the 2nd pass of the pipeline.
    affine_transfo -i ${file_c_or_crop_r}.nii.gz -transfo ${PATH_RESULTS}/$transfo_file -config ${PATH_RESULTS}/$SUBJECT/$config_script -o _t${i_transfo}
    file_c_or_crop_r_t=${file_c_or_crop_r}_t${i_transfo}
    # transform the labeled segmentation with same transfo values
    affine_transfo -i ${file_c_or_crop_label_r}.nii.gz -transfo ${PATH_RESULTS}/$transfo_file -config ${PATH_RESULTS}/$SUBJECT/$config_script -o _t${i_transfo} -interpolation 0
    file_c_or_crop_label_r_t=${file_c_or_crop_label_r}_t${i_transfo}
    # Segment spinal cord (only if it does not exist)
    segment_if_does_not_exist ${file_c_or_crop_r_t} ${contrast}
    # name segmented file
    file_c_or_crop_r_t_seg=${FILESEG}
    # Compute average CSA between C2 and C5 levels (append across subjects)
    sct_process_segmentation -i $file_c_or_crop_r_t_seg.nii.gz -vert 2:5 -perlevel 1 -vertfile $file_c_or_crop_label_r_t.nii.gz -o $PATH_RESULTS/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv -qc ${PATH_QC}
    # add files to check
    FILES_TO_CHECK+=(
    "$PATH_RESULTS/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv"
    "$PATH_DATA_PROCESSED/${SUBJECT}/anat_r${r_coef}/${file_c_or_crop_r_t}.nii.gz"
    "$PATH_DATA_PROCESSED/${SUBJECT}/anat_r${r_coef}/${file_c_or_crop_r_t_seg}.nii.gz"
    "$PATH_DATA_PROCESSED/${SUBJECT}/anat_r${r_coef}/${file_c_or_crop_label_r_t}.nii.gz"
    )
  done
  cd ../
done

# Check the presence of output files and write log file if error
# ==============================================================
for file in ${FILES_TO_CHECK[@]}; do
  if [ ! -e $file ]; then
    echo "$file does not exist" >> $PATH_LOG/error.log
  fi
done

# Display useful info for the log
# ===============================
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"