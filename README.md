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
run the script within the Dataset folder
~~~
sct_run_batch -path-data data process_data.sh -jobs 2
~~~
