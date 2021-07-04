Quality control
================

After running the analysis, check your Quality Control (QC) report by opening the file qc/index.html. Use the
“Search” feature of the QC report to quickly jump to segmentations or labeling results. If you spot issues
(wrong labeling), add their filenames in the 'config_correction.yml' file
(see https://spine-generic.rtfd.io/en/latest/analysis-pipeline.html for further indications). Then, manually create
labels in the cord at the level of inter-vertebral discs C1-C2, C2-C3, ..., C4-C5 with the command:

.. code-block:: python

    manual_correction -config config_correction.yml -path-in csa_atrophy_results/data_processed -path-out PATH_DATA

The bash script outputs all manual labelings to the derivatives directory in the dataset path defined in `path_data`.
It is now possible to re-run the whole process. With the command below labeling will use the manual corrections that
are present in the derivatives/ folder of the dataset, otherwise labeling will be done automatically.

.. code-block:: python

    sct_run_batch -config config_sct_run_batch.yml
