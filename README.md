![csa-atrophy](https://github.com/sct-pipeline/csa-atrophy/blob/master/csa_atrophy_scheme3.png)

[![Documentation Status](https://readthedocs.org/projects/sphinx/badge/?version=master)](https://csa-atrophy.readthedocs.io/en/latest/)
   
# csa-atrophy

The csa-atrophy framework aims to evaluate the robustness and the sensitivity of an automated analysis pipeline for detecting SC atrophy.

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

Download the [Spine Generic Multi-Subject dataset](https://github.com/spine-generic/data-multi-subject#download). 

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
csa_rescale_stat -i csa_atrophy_results/results -o csa_atrophy_results -config config_script.yml -fig
~~~



# Statistics

Plot results:
- STD_intersub
- Mean and STD inter-subject error percentage in function of rescaling
- sample size: minimum number of patients to detect an atrophy of X with Y% power and Z% uncertainty
- CSA values boxplot in function of rescaling
- Error values boxplot in function of rescaling

