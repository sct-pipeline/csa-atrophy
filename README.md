# csa-atrophy
Evaluate the sensitivity of atrophy detection with SCT

# Data
~~~
Dataset/
└── participants.json
└── participants.tsv
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
        └── CSA_rescale_stat.py
~~~
# How to run
Download (or git clone) this repository:
~~~
git clone https://github.com/sct-pipeline/csa-atrophy.git
cd csa-atrophy
~~~
to fetch dataset run command: (this file should be edited according to your needs and chmod permission should be given)
~~~
./csa-fetch-dataset
~~~
run the script within the Dataset folder (using sct venv if needed)
~~~
sct_run_batch -path-data data process_data.sh -jobs 2
~~~
to output statistics, run in Dataset
~~~
python CSA_rescale_stat.py -i results/CSA.csv -r results/CSA_r.csv
~~~
