.. csa-atrophy documentation master file, created by
   sphinx-quickstart on Fri Mar 12 16:44:50 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to csa-atrophy's documentation!
=========================================
csa-atrophy evaluates the robustness and the sensitivity of an automated analysis pipeline for detecting SC atrophy. Notably, the proposed framework utilizes image scaling and applies a random rigid transformation to mimic subject repositioning (scan-rescan). This enables the quantification of the accuracy and precision of the estimated CSA across various degrees of simulated atrophy. As presented below, statistics from these experiments such as power analyses and minimum sample sizes are derived.

.. image:: csa_atrophy_scheme3.png

.. toctree::
   :maxdepth: 2
   :caption: Statistics:
   Intro.rst
   Intra_subject.rst
   Inter_subject.rst
   Sample_size.rst
