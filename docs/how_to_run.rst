# csa-atrophy

The csa-atrophy framework aims to evaluate the robustness and the sensitivity of an automated analysis pipeline for detecting SC atrophy.

# How to run

This code has been tested using Python 3.7.

Download (or git clone) this repository:

.. code-block:: python

    git clone https://github.com/sct-pipeline/csa-atrophy.git
    cd csa-atrophy


Installation:
csa-atrophy requires specific python packages for computing statistics and processing images. If not already present on the computer's python environment such packages will automatically be installed by running pip command:

.. code-block:: python

    pip install -e .

Download the results file from [Spine Generic Multi-Subject dataset](https://github.com/spine-generic/data-multi-subject/releases/tag/r20201130).

Edit the file `config_sct_run_batch.yml` according to your setup. Notable flags include:

*  `path_data`: If you downloaded the spine-generic data at another location, make sure to update the path;
* `include_list`: If you only want to run the script in a few subjects, list them here. Example:
  `include_list: ['sub-unf04', 'sub-unf05']`

See `sct_run_batch -h` to look at the available options.

Run the analysis:

.. code-block:: python

    sct_run_batch -config config_sct_run_batch.yml


note: desired subjects using flag -include and in parallel processing using flag -jobs.

To output statistics, run in Dataset

.. code-block:: python

    csa_rescale_stat -i csa_atrophy_results/results -o csa_atrophy_results -config config_script.yml -fig

