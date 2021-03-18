STATISTICS
============
A module about statistics.

Subject dataframe
==================

After everything is done, compute stats: Per-subject stat: Panda dataframe ``df_sub``:

CSA estimation
"""""""""""""""

intra-subject MEAN: MEAN[CSA(sI, rX, :)] --> MEAN_intra(sI, rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 10
    
Intra-subject SD
""""""""""""""""""

intra-subject STD: STD[CSA(sI, rX, :)] --> STD_intra(sI, rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 11
    
    
Intra-subject COV
""""""""""""""""""

intra-subject COV: STD[CSA(sI, rX, :)] / MEAN[CSA(sI, rX, :)] --> COV_intra(sI, rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 12
    
Rescale estimation
""""""""""""""""""""

rescale_estimated_subject MEAN: MEAN[CSA(sI, rX, :) / CSA(sI, 1, :)] --> MEAN_rescale_estimated_subject(sI, rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 14 
    
Error
""""""""""""""""""

intra-subject error MEAN: MEAN[CSA(sI, rX, :)] - (rX^2 * MEAN[CSA(sI, 1, :)]) --> MEAN_error_intra(sI, rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 15
    
Rescale dataframe
=====================

Across-subject stats: Panda dataframe ``df_rescale``
 
Mean SD of CSA across Monte-Carlo samples
""""""""""""""""""""""""""""""""""""""""""""

intra-subject STD: MEAN[STD_intra(:, rX)] --> STD_intra(rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 9
   
Mean COV of CSA across Monte-Carlo samples
""""""""""""""""""""""""""""""""""""""""""""

intra-subject COV: MEAN[COV_intra_sub(:, rX)] --> COV_intra(rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 10
   
Mean COV of CSA across subjects
""""""""""""""""""""""""""""""""""""""""""""

inter-subject STD: STD[MEAN_intra(:, rX)] --> STD_inter(rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 11
   
Mean rescale estimated (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

rescale_estimated (across subjects) MEAN: MEAN[MEAN_rescale_estimated_subject(:, rX)] --> MEAN_rescale_estimated(rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 12
   
SD rescale estimated (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

rescale_estimated (across subjects) STD: STD[MEAN_rescale_estimated_subject(:, rX)] --> STD_rescale_estimated(rX)

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 13
   
Mean error in percentage (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

error in percentage (across subjects) MEAN: MEAN[MEAN_error_intra(:, rX)]

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 14
   
SD od error in percentage (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

error in percentage (across subjects) STD: STD[MEAN_error_intra(:, rX)]

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 15
   

Power analysis
=====================
   
Between group minimum sample size
"""""""""""""""""""""""""""""""""""

the minimum sample size (number of subjects per study arm) necessary to detect an atrophy between groups was computed based on a two-sample (unpaired) bilateral t-test using the following formula (Wang and Ji 2020; Wittes 2002):

:math:`n_{unpaired} = (z_{α/2} + z_{β})^2(SD + SD)^2 / diff_{group}^2`.

Where n_{unpaired} is the minimum sample size required to differentiate between groups with a given power (z_{β} corresponds to the power z score, e.g. 80% power gives β=0.2 and z_{β}= -0.84) and level of significance (z_{α/2} corresponds to the significance level z score, e.g. 5% level of significance gives 𝛂=0.05 and z_{α/2}=-1.96), SD is the inter-subject standard deviation of the mean CSA (which was calculated by taking the mean CSA across Monte Carlo samples). diff_{group} group is the difference of the mean CSA between the groups.

.. autofunction:: csa_rescale_stat.sample_size  

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 5-7
   
   
