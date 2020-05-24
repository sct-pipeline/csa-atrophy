# csa-atrophy
Evaluate the sensitivity of atrophy detection with SCT

# Data
~~~
Dataset/
└── csa_rescale_stat.py
└── process_data.sh
└── csa-fetch-dataset.sh
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

~~~
# How to run
Download (or git clone) this repository:
~~~
git clone https://github.com/sct-pipeline/csa-atrophy.git
cd csa-atrophy
~~~
To fetch/sync dataset run command: (this file should be edited according to your needs)
~~~
./csa-fetch-dataset.sh
~~~
Run the script within the Dataset folder (using sct venv if needed)
~~~
sct_run_batch -path-data data process_data.sh -jobs 2
~~~
To output statistics, run in Dataset (if needed run requirements file)
~~~
pip install -r requirements.txt
python csa_rescale_stat.py -i results
~~~
