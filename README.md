# csa-atrophy

Evaluate the sensitivity of atrophy detection with SCT. The algorithm works as follows:
- Consider subject I --> sI
- Applies a rescaling on the native image (e.g. 1, 0.95, 0.8) --> rX
- Applies random affine transformation --> tY
- Segment the cord
- Compute CSA --> CSA(sI, rX, tY)

After everything is done, compute stats:
- intra-subject STD: STD[CSA(sI, rX, :)] --> STD_intrasub(sI, rX)
- inter-subject STD: STD_intrasub(:, rX) --> STD_intersub(rX)
- intra-subject error MEAN: MEAN[CSA(sI, rX, :)] - rX^2*MEAN[CSA(sI, 1, :)] (for rX included in list of rX with rX=1 excluded)
  --> MEAN(sI, rX)
- inter-subject error STD: STD[MEAN(:, rX)]
- TODO: sample size

Plot results:
- STD_intersub
- TODO: sample size

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

Edit the file `config_sct_run_batch.yml` according to your setup. Notable flags include:
- `path_data`: If you downloaded the spine-generic data at another location, make sure to update the path;
- `include_list`: If you only want to run the script in a few subjects, list them here. Example:
  `include_list: ['sub-unf04', 'sub-unf05']`

See `sct_run_batch -h` to look at the available options.

Run the analysis:
~~~
sct_run_batch -config config_sct_run_batch.yml
~~~

:note: **desired subjects using flag -include and in parallel processing using flag -jobs.**

To output statistics, run in Dataset
~~~
csa_rescale_stat -i csa_atrophy_results/data_processed -o csa_atrophy_results -config config_script.yml -v
~~~

# Quality Control

After running the analysis, check your Quality Control (QC) report by opening the file qc/index.html. Use the “Search” feature of the QC report to quickly jump to segmentations or labeling results. If you spot issues (wrong labeling), add their filename in the 'config_correction.yml' file (see https://spine-generic.readthedocs.io/en/latest/analysis_pipeline.html for further indications). Then manually create labels in the cord on the posterior tip of inter-vertebral discs from C2 to C5 with command:
~~~
manual_correction -config config_correction.yml -path-in csa_atrophy_results/data_processed -path-out spinegeneric_r20200801
~~~
The bash script outputs all effectuated manual labelings to the derivatives directory in 'spinegeneric_r20200801'.
It is now possible to re-run the whole process, pointing to the manual corrections. With the command below labeling will use the manual corrections present in 'seg_manual', otherwise labeling will be done automatically.
~~~
sct_run_batch -path-data spinegeneric_r20200801 -path-output csa_atrophy_results_corrected -script process_data.sh
~~~
