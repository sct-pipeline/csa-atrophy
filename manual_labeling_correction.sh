#!/bin/bash
#########################################################################################
#
# If you spot issues (wrong labeling), add their filename in the config.yaml file.
# Then manually create labels in the cord on the posterior tip of inter-vertebral discs from C2 to C5
#
# The bash script outputs all effectuated manual labellings to 'results_corrected/seg_manual'.
# It is now possible to re-run the whole process, pointing to the manual corrections.
#
# ---------------------------------------------------------------------------------------
#
# Example: ./manual_labeling_corrections.sh
# About the license: see the file LICENSE
###################################################################

# Local folder to output the manual labels
PATH_DATA=$(yaml_parser -o path_data -i config_sct_run_batch.yml)
PATH_ORIGINAL_RESULTS=$(yaml_parser -o PATH_ORIGINAL_RESULTS -i config_script.yml)
PATH_ORIGINAL_CSV=$(yaml_parser -o PATH_ORIGINAL_CSV -i config_script.yml)
PATH_DESTINATION_RESULTS=$(yaml_parser -o PATH_DESTINATION_RESULTS -i config_script.yml)
PATH_SEGMANUAL=$(yaml_parser -o PATH_SEGMANUAL -i config_script.yml)
mkdir -p ${PATH_SEGMANUAL}
mkdir -p ${PATH_DESTINATION_RESULTS}
cp ${PATH_ORIGINAL_CSV}/transfo_values.csv ${PATH_DESTINATION_RESULTS}

FILES_SEGMANUAL=$(yaml_parser -o FILES_SEG -i config_correction.yml)
FILES_S=$(echo $FILES_SEGMANUAL | tr '[]' ' ' | tr "'" ' ' | tr ',' ' ')
# Loop across files
if [ "$FILES_S" == "None" ]
then
  echo "no files found for manual segmentation"
else
  for file in ${FILES_S[@]}; do
    subject=${file%%_*}
    rescaling_split=${file%%_t*}
    rescaling=${rescaling_split#*_r}
    mkdir -p data_to_correct/$subject/anat
    cp $PATH_ORIGINAL_RESULTS/$subject/anat_r${rescaling}/$file data_to_correct/$subject/anat
  done
fi

# List of subjects to create manual label
FILES_LABELMANUAL=$(yaml_parser -o FILES_LABEL -i config_correction.yml)
FILES_L=$(echo $FILES_LABELMANUAL | tr '[]' ' ' | tr "'" ' ' | tr ',' ' ')
# Loop across files
if [ "$FILES_L" == "None" ]
then
  echo "no files found for labeling correction"
else
  for file in ${FILES_L[@]}; do
    subject=${file%%_*}
    rescaling_split=${file%%_t*}
    rescaling=${rescaling_split#*_r}
    mkdir -p data_to_correct/$subject/anat
    cp $PATH_ORIGINAL_RESULTS/$subject/anat_r${rescaling}/$file data_to_correct/$subject/anat
  done
fi


FILES_GMMANUAL=$(yaml_parser -o FILES_GMSEG -i config_correction.yml)
FILES_G=$(echo $FILES_GMMANUAL | tr '[]' ' ' | tr "'" ' ' | tr ',' ' ')
# Loop across files
if [ "$FILES_G" == "None" ]
then
  echo "no files found for grey matter manual segmentation"
else
  for file in ${FILES_G[@]}; do
    subject=${file%%_*}
    rescaling_split=${file%%_t*}
    rescaling=${rescaling_split#*_r}
    mkdir -p data_to_correct/$subject/anat
    cp $PATH_ORIGINAL_RESULTS/$subject/anat_r${rescaling}/$file data_to_correct/$subject/anat
  done
fi

sg_manual_correction -config config_correction.yml -path-in data_to_correct -path-out $PATH_DATA
