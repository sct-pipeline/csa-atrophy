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
# Installation
csa-atrophy requires specific python packages for computing statistics and processing images. If not already present on the computer's python environment such packages will automatically be installed by running pip command:
~~~
pip install -e .
~~~
# How to run
Download (or git clone) this repository:
~~~
git clone https://github.com/sct-pipeline/csa-atrophy.git
cd csa-atrophy
~~~
Fetch dataset (2 choices):
  - To fetch/sync dataset run command: (this file should be edited according to your needs)
  ~~~
  ./csa_fetch_dataset.sh
  ~~~
  - To fetch specific version of openneuro repository (version: 1.0.5, Files: 2501, Size: 7.9GB, Subjects: 248) follow instructions to set openeuro CLI: https://www.npmjs.com/package/openneuro-cli
  ~~~
  openneuro download --snapshot 1.0.5 ds001919 data
  ~~~
Run the following script within the Dataset folder to extract CSA. This script can be run on desired subjects using flag -include and in parallel processing using flag -jobs.
~~~
sct_run_batch -config_sct_run_batch.yml
~~~
To output statistics, run in Dataset
~~~
python csa_rescale_stat.py -i csa_atrophy_results/results/csa_data -o csa_atrophy_results -v
~~~

# Quality Control

After running the analysis, check your Quality Control (QC) report by opening the file qc/index.html. Use the “Search” feature of the QC report to quickly jump to segmentations or labeling results. If you spot issues (wrong labeling), add their filename in the 'config_correction.yml' file (see https://spine-generic.readthedocs.io/en/latest/analysis_pipeline.html for further indications). Then manually create labels in the cord on the posterior tip of inter-vertebral discs from C2 to C5 with command:
~~~
manual_correction -config config_correction.yml -path-in csa_atrophy_results/results -path-out data
~~~
The bash script outputs all effectuated manual corrections to 'data/derivatives/labels'.
It is now possible to re-run the whole process, pointing to the manual corrections. With the command below labeling will use the manual corrections present in 'data/derivatives/labels', otherwise labeling will be done automatically.
~~~
sct_run_batch -path-data data -path-output csa_atrophy_results_corrected -path-segmanual data/derivatives/labels -script process_data.sh
~~~
