# csa-atrophy
Evaluate the sensitivity of atrophy detection with SCT

# Data
~~~
Dataset/
└── csa_rescale_stat.py
└── process_data.sh
└── csa_fetch_dataset.sh
└── data
    └── sub-subj01
        └── anat
            └── sub-subj01_T2w.nii.gz
            └── sub-subj01_T2w.json
     └── sub-subj02
         └── anat
            └── sub-subj02_T2w.nii.gz
            └── sub-subj02_T2w.json
 └── results

~~~

# How to run
Download (or git clone) this repository:
~~~
git clone https://github.com/sct-pipeline/csa-atrophy.git
cd csa-atrophy
~~~
Installation:
csa-atrophy requires specific python packages for computing statistics and processing images. If not already present on the computer's python environment such packages will automatically be installed by running pip command:
~~~
pip install -e .
~~~
Fetch dataset:
Suggested testing dataset must be downloaded from "Spine Generic Public Database". To download latest version of the whole multi-subject dataset run commands:
  ~~~
  curl -o spinegeneric.zip -L https://github.com/spine-generic/data-multi-subject/archive/r20200801.zip
  unzip spinegeneric.zip
  ~~~
Run the following script within the Dataset folder to extract CSA. This script can be run on desired subjects using flag -include and in parallel processing using flag -jobs.
~~~
sct_run_batch -path-data data-multi-subject-master -path-output csa_atrophy_results process_data.sh
~~~
To output statistics, run in Dataset
~~~
python csa_rescale_stat.py -i csa_atrophy_results/results/csa_data -o csa_atrophy_results -v
~~~

# Quality Control

After running the analysis, check your Quality Control (QC) report by opening the file qc/index.html. Use the “Search” feature of the QC report to quickly jump to segmentations or labeling results. If you spot issues (wrong labeling), add their filename in the variable array "FILES" of the 'manual_labeling_correction.sh' script. Then manually create labels in the cord on the posterior tip of inter-vertebral discs from C2 to C5 with command:
~~~
./manual_labeling_correction.sh
~~~
The bash script outputs all effectuated manual labelings to 'results_corrected/seg_manual'.
It is now possible to re-run the whole process, pointing to the manual corrections. With the command below labeling will use the manual corrections present in 'seg_manual', otherwise labeling will be done automatically.
~~~
sct_run_batch -path-data data-multi-subject-master -path-output csa_atrophy_results_corrected -path-segmanual seg_manual process_data.sh
~~~
