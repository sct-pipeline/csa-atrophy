# csa-atrophy

Evaluate the sensitivity of atrophy detection with SCT

# How to run

This code has been tested using Python 3.7.

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
curl -o spinegeneric_r20200801.zip -L https://github.com/spine-generic/data-multi-subject/archive/r20200801.zip
unzip spinegeneric_r20200801.zip
~~~
To run csa-atrophy pipeline on a sub-dataset of spine-generic there are 2 options: 
1. Create a new folder with a sub-dataset (make sure to update config file path_data);
2. Use flag -include in `sct_run_batch` to select one subject at a time.

Run the following script within the Dataset folder to extract CSA. This script can be run on desired subjects using flag -include and in parallel processing using flag -jobs.
~~~
sct_run_batch -config config_sct_run_batch.yml
~~~
To output statistics, run in Dataset
~~~
python csa_rescale_stat.py -i csa_atrophy_results/data_processed -o csa_atrophy_results -config config_script.yml -v
~~~

# Quality Control

After running the analysis, check your Quality Control (QC) report by opening the file qc/index.html. Use the “Search” feature of the QC report to quickly jump to segmentations or labeling results. If you spot issues (wrong labeling), add their filename in the variable array "FILES_SEGMANUAL" of the 'config.yaml' file. Then manually create labels in the cord on the posterior tip of inter-vertebral discs from C2 to C5 with command:
~~~
manual_correction -config config_correction.yml -path-in csa_atrophy_results/results -path-out spinegeneric_r20200801
~~~
The bash script outputs all effectuated manual labelings to 'results_corrected/seg_manual'.
It is now possible to re-run the whole process, pointing to the manual corrections. With the command below labeling will use the manual corrections present in 'seg_manual', otherwise labeling will be done automatically.
~~~
sct_run_batch -path-data spinegeneric_r20200801 -path-output csa_atrophy_results_corrected -path-segmanual seg_manual -script process_data.sh
~~~
