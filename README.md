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
csa-atrophy requires specific python packages for computing statistics and processing images. If not already present on the computers python environment such packages will automatically be installed by running pip command:
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
Run the following script within the Dataset folder to extract CSA. This script can be ran on desired subjects using flag -include and in parallel processing using flag -jobs.
~~~
sct_run_batch -path-data data process_data.sh
~~~
To output statistics, run in Dataset
~~~
python csa_rescale_stat.py -i results -v
~~~

