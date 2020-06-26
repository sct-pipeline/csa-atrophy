#!/bin/bash
# Local folder to output the manual labels
PATH_SEGMANUAL="results_corrected/seg_manual"
mkdir -p ${PATH_SEGMANUAL}
# List of subjects to create manual labels
FILES=(
  "sub-amu01_T2w_r0.85_t1.nii.gz"
)
# Loop across subjects
for file in ${FILES[@]}; do
  subject=${file%%_*}
  rescaling_split=${file%%_t*}
  rescaling=${rescaling_split#*_r}
  sct_label_utils -i results/$subject/anat_r${rescaling}/$file -create-viewer 2,3,4,5 -o ${PATH_SEGMANUAL}/${file%%.nii.gz}_labels-manual.nii.gz -msg "Click at the posterior tip of inter-vertebral disc"
done
