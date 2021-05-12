.. csa-atrophy documentation master file, created by
   sphinx-quickstart on Fri Mar 12 16:44:50 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to csa-atrophy's documentation!
=========================================

STATISTICS
============
Statistics documentation to evaluate the sensitivity of atrophy detection with SCT. Each image CSA is indexed as follows: :math:`CSA(sI, rX, tY)` where

- :math:`(sI)` corresponds to subject :math:`I`
- :math:`(rX)` corresponds to the applied rescaling :math:`X` on the native image (e.g. 1, 0.95, 0.8)
- :math:`(tY)` corresponds to the applied random affine transformation Y on the native image

.. toctree::
   About <README.md>
   subject_dataframe
   rescale_dataframe
   sample_size
   :maxdepth: 2
   :caption: Contents:
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
