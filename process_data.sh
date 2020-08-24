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

# TODO: simplify variable names: no need to call file_or_c_blablabla. Instead: update file

# Uncomment for full verbose
set -x

# Immediately exit if error
set -e

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# get starting time:
start_all=`date +%s`
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

# Check if manual label already exists. If it does, copy it locally. If it does
# not, perform automatic labeling.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  local contrast="$3"
  local contrast_str="$4"
  # Update global variable with segmentation file name
  FILELABEL="${file}_labels-disc"
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
  local qc=$3
  # Update global variable with segmentation file name
  FILESEG="${file}_seg"
  FILESEGMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILESEG}-manual.nii.gz"
  if [ -e $FILESEGMANUAL ]; then
    echo "Found! Using manual segmentation."
    rsync -avzh $FILESEGMANUAL ${FILESEG}.nii.gz
    sct_qc -i ${file}.nii.gz -s ${FILESEG}.nii.gz -p sct_deepseg_sc $qc
  else
    # Segment spinal cord
    sct_deepseg_sc -i ${file}.nii.gz -c ${contrast} $qc
  fi
}


# SCRIPT STARTS HERE
# ==============================================================================
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short
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
  interp="0.8x0.8x0.8"
elif [ $contrast == "t1" ]; then
  contrast_str="T1w"
  interp="1x1x1"
fi
file=${SUBJECT}_${contrast_str}
# Reorient to RPI
sct_image -i ${file}.nii.gz -setorient RPI -o ${file}_RPI.nii.gz
file=${file}_RPI
# Resample isotropically
sct_resample -i ${file}.nii.gz -mm $interp -o ${file}_r.nii.gz
file=${file}_r
end=`date +%s`
runtime=$((end-start))
echo "+++++++++++ TIME: Duration of reorienting and resampling:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"

# Segment spinal cord (only if it does not exist) in dir anat
start=`date +%s`
segment_if_does_not_exist $file ${contrast} "-qc ${PATH_QC} -qc-subject ${SUBJECT}"
# name segmented file
file_seg=${FILESEG}
end=`date +%s`
runtime=$((end-start))
echo "+++++++++++ TIME: Duration of segmentation on non cropped image:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"

# Label spinal cord (only if it does not exist) in dir anat
start=`date +%s`
label_if_does_not_exist $file $file_seg $contrast $contrast_str
file_label=${file_seg}_labeled
end=`date +%s`
runtime=$((end-start))
echo "+++++++++++ TIME: Duration of labelling:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"

# dilate segmentation (for cropping)
sct_maths -i ${file_seg}.nii.gz -dilate 15 -shape cube -o ${file_seg}_dil.nii.gz
# crop image
sct_crop_image -i ${file}.nii.gz -m ${file_seg}_dil.nii.gz
file=${file}_crop
# crop segmentation
sct_crop_image -i ${file_seg}.nii.gz -m ${file_seg}_dil.nii.gz
file_seg=${file_seg}_crop
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
  file_r=${file}_r${r_coef}
  affine_rescale -i ../anat/${file}.nii.gz -r ${r_coef} -o ${file_r}.nii.gz
  # rescale labeled segmentation
  file_label_r=${file_r}_seg_labeled
  affine_rescale -i ../anat/${file_label}.nii.gz -r ${r_coef} -o ${file_label_r}.nii.gz

  # create list of array to iterate on (e.g.: seq_transfo = 1 2 3 4 5 if n_transfo=5)
  seq_transfo=$(seq ${n_transfo})
  # Iterate across transformations
  for i_transfo in ${seq_transfo[@]}; do
    # Image random transformation (rotation, translation). By default transformation values are taken from
    # "transfo_values.csv" file if it already exists.
    # We keep a transfo_values.csv file, so that after first pass of the pipeline and QC, if segmentations
    # need to be manually-corrected, we want the transformations to be the same for the 2nd pass of the pipeline.
    start=`date +%s`
    affine_transfo -i ${file_r}.nii.gz -transfo ${PATH_RESULTS}/transfo_${file_r} -config $config_script -o _t${i_transfo}
    file_r_t=${file_r}_t${i_transfo}
    end=`date +%s`
    runtime=$((end-start))
    echo "+++++++++++ TIME: Duration of of image transfo t${i_transfo}:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
    # transform the labeled segmentation with same transfo values

    start=`date +%s`
    affine_transfo -i ${file_label_r}.nii.gz -transfo ${PATH_RESULTS}/transfo_${file_r} -config $config_script -o _t${i_transfo} -interpolation 0
    file_label_r_t=${file_label_r}_t${i_transfo}
    end=`date +%s`
    runtime=$((end-start))
    echo "+++++++++++ TIME: Duration of labelling transfo t${i_transfo}:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
    # Segment spinal cord (only if it does not exist)
    start=`date +%s`
    segment_if_does_not_exist ${file_r_t} ${contrast}
    end=`date +%s`
    runtime=$((end-start))
    echo "+++++++++++ TIME: Duration of segmentation t${i_transfo}:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
    # name segmented file
    file_r_t_seg=${FILESEG}
    # Compute average CSA between C2 and C5 levels (append across subjects)
    start=`date +%s`
    sct_process_segmentation -i $file_r_t_seg.nii.gz -vert 2:5 -perlevel 1 -vertfile $file_label_r_t.nii.gz -o $PATH_RESULTS/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv
    end=`date +%s`
    runtime=$((end-start))
    echo "+++++++++++ TIME: Duration of sct_process_segmentation t${i_transfo}:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
    # add files to check
    FILES_TO_CHECK+=(
    "$PATH_RESULTS/csa_perlevel_${SUBJECT}_t${i_transfo}_${r_coef}.csv"
    "$PATH_DATA_PROCESSED/${SUBJECT}/anat_r${r_coef}/${file_r_t}.nii.gz"
    "$PATH_DATA_PROCESSED/${SUBJECT}/anat_r${r_coef}/${file_r_t_seg}.nii.gz"
    "$PATH_DATA_PROCESSED/${SUBJECT}/anat_r${r_coef}/${file_label_r_t}.nii.gz"
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
runtime=$((end-start_all))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"