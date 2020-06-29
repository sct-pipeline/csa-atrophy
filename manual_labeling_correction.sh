#!/bin/bash
# Local folder to output the manual labels
PATH_ORIGINAL_RESULTS="csa_atrophy_results/results"
PATH_DESTINATION_RESULTS="csa_atrophy_results_corrected/results"
PATH_SEGMANUAL="csa_atrophy_results_corrected/seg_manual"
mkdir -p ${PATH_SEGMANUAL}
mkdir -p ${PATH_DESTINATION_RESULTS}
cp ${PATH_ORIGINAL_RESULTS}/transfo_values.csv ${PATH_DESTINATION_RESULTS}
# List of subjects to create manual labels
FILES=(
  "sub-amu01_T2w_r0.90_t1.nii.gz"
)
# Loop across files
for file in ${FILES[@]}; do
  subject=${file%%_*}
  rescaling_split=${file%%_t*}
  rescaling=${rescaling_split#*_r}
  sct_label_utils -i $PATH_ORIGINAL_RESULTS/$subject/anat_r${rescaling}/$file -create-viewer 2,3,4,5,6 -o ${PATH_SEGMANUAL}/${file%%.nii.gz}_labels-manual.nii.gz -msg "Click at the posterior tip of inter-vertebral disc"
done
