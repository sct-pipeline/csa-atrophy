CSA-atrophy
============
CSA-atrophy evaluates the robustness and the sensitivity of an automated analysis pipeline for detecting SC atrophy. Notably, the proposed framework utilizes image scaling and applies a random rigid transformation to mimic subject repositioning (scan-rescan). This enables the quantification of the accuracy and precision of the estimated CSA across various degrees of simulated atrophy. As presented in section statistics, statistics from these experiments such as power analyses and minimum sample sizes are derived.

.. image:: csa_atrophy_scheme3.png


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Statistics

   statistics/introduction.rst
   statistics/intra_subject.rst
   statistics/inter_subject.rst
   statistics/sample_size.rst


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Quality Control

    quality_control.rst