STATISTICS
============
A module about statistics to evaluate the sensitivity of atrophy detection with SCT. The algorithm works as follows:

- Consider subject :math:`I \to sI`
- Applies a rescaling on the native image (e.g. 1, 0.95, 0.8) :math:`\to rX`
- Applies random affine transformation :math:`\to tY`
- Segment the cord
- Compute CSA :math:`\to CSA(sI, rX, tY)`

Subject dataframe
==================

After everything is done, compute stats: Per-subject stat: Panda dataframe ``df_sub``:

CSA estimation
"""""""""""""""

intra-subject MEAN: :math:`MEAN[CSA(sI, rX, :)] \to MEAN_intra(sI, rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 10
    
Intra-subject SD
""""""""""""""""""

intra-subject STD: :math:`STD[CSA(sI, rX, :)] \to STD_intra(sI, rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 11
    
    
Intra-subject COV
""""""""""""""""""

intra-subject COV: :math:`\frac{STD[CSA(sI, rX, :)]}{MEAN[CSA(sI, rX, :)]} \to COV_intra(sI, rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 12
    
Rescale estimation
""""""""""""""""""""

rescale_estimated_subject MEAN: :math:`MEAN \left[\frac{CSA(sI, rX, :)}{CSA(sI, 1, :)}\right] \to MEAN_{rescale\_estimated\_subject}(sI, rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 14 
    
Error
""""""""""""""""""

intra-subject error MEAN: :math:`MEAN[CSA(sI, rX, :)] - rX^2  \times MEAN[CSA(sI, 1, :)] \to MEAN_{error\_intra}(sI, rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 15
    
Rescale dataframe
=====================

Across-subject stats: Panda dataframe ``df_rescale``
 
Mean SD of CSA across Monte-Carlo samples
""""""""""""""""""""""""""""""""""""""""""""

intra-subject STD: :math:`MEAN[STD_intra(:, rX)] \to STD_intra(rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 9
   
Mean COV of CSA across Monte-Carlo samples
""""""""""""""""""""""""""""""""""""""""""""

intra-subject COV: :math:`MEAN[COV_{intra\_sub}(:, rX)] \to COV_intra(rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 10
   
Mean COV of CSA across subjects
""""""""""""""""""""""""""""""""""""""""""""

inter-subject STD: :math:`STD[MEAN_intra(:, rX)] \to STD_inter(rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 11
   
Mean rescale estimated (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

rescale_estimated (across subjects) MEAN: :math:`MEAN[MEAN_{rescale\_estimated\_subject}(:, rX)] \to MEAN_{rescale\_estimated}(rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 12
   
SD of rescale estimated (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

rescale_estimated (across subjects) STD: :math:`STD[MEAN_{rescale\_estimated\_subject}(:, rX)] \to STD_{rescale\_estimated}(rX)`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 13
   
Mean error in percentage (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

error in percentage (across subjects) MEAN: :math:`MEAN[MEAN_{error\_intra}(:, rX)]`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 14
   
SD of error in percentage (across subjects)
""""""""""""""""""""""""""""""""""""""""""""

error in percentage (across subjects) STD: :math:`STD[MEAN_{error\_intra}(:, rX)]`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 15
   

Power analysis
=====================
   
Between-group minimum sample size
"""""""""""""""""""""""""""""""""""

The minimum sample size (number of subjects per study arm) necessary to detect an atrophy between groups was computed based on a two-sample (unpaired) bilateral t-test using the following formula (Wang and Ji 2020; Wittes 2002):

:math:`n_{unpaired} = \frac{(z_{α/2} + z_{β})^2(\sigma_{100}+\sigma_{rX})^2}{\Delta_{sub} ^2}`

Where :math:`n_{unpaired}` is the minimum sample size required to differentiate between groups with a given power (:math:`z_{β}` corresponds to the power z score, e.g. 80% power gives β=0.2 and :math:`z_{β}`= -0.84) and level of significance (:math:`z_{α/2}` corresponds to the significance level z score, e.g. 5% level of significance gives 𝛂=0.05 and :math:`z_{α/2}`=-1.96), SD is the inter-subject standard deviation of the mean CSA (which was calculated by taking the mean CSA across Monte Carlo samples). :math:`diff_{group}` group is the difference of the mean CSA between the groups.

.. autofunction:: csa_rescale_stat.sample_size  

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 5-7

Within-subject minimum sample size
"""""""""""""""""""""""""""""""""""

the minimum sample size necessary to detect an atrophy in a within-subject (repeated-measures) study was computed based on a two-sample bilateral paired t-test using the following formula (Altmann et al. 2009):

:math:`n_{within-sub} = \frac{(z_{α/2} + z_{β})^2(\sigma_{diff})^2}{\Delta_{sub} ^2}`
   
Where :math:`\sigma_{diff}` is the standard deviation between longitudinal CSA measures across  subjects and :math:`\Delta_{sub}` is the mean of the difference between longitudinal CSA measures.

.. autofunction:: csa_rescale_stat.sample_size

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 9-10