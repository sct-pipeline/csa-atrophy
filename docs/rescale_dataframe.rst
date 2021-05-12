STATISTICS
============
Statistics documentation to evaluate the sensitivity of atrophy detection with SCT. Each image computed CSA is indexed as follows: :math:`CSA(sI, rX, tY)` where

- :math:`(sI)` corresponds to subject :math:`I`
- :math:`(rX)` corresponds to the applied rescaling :math:`X` on the native image (e.g. 1, 0.95, 0.8)
- :math:`(tY)` corresponds to the applied random affine transformation Y on the native image

Rescale dataframe
=====================

Across-subject statistics are regrouped in the Panda dataframe ``df_rescale``
 
Mean intra-subject SD
""""""""""""""""""""""""""""""""""""""""""""

Per scaling mean intra-subject SD: :math:`MEAN[\sigma_{(rX,sI,:)}(:, rX)] \to \overline{\sigma}_{(:,rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 9
   
Mean intra-subject COV
""""""""""""""""""""""""""""""""""""""""""""

Per scaling mean intra-subject COV: :math:`MEAN[COV_{(rX,sI,:)}(:, rX)] \to \overline{COV}_{(:,rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 10
   
Inter-subject SD
""""""""""""""""""""""""""""""""""""""""""""

Per scaling SD of the intra-subject CSA: :math:`STD[\overline{CSA}_{(sI, rX, :)}(:, rX)] \to \sigma_{(:,rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 11
   
Mean rescale estimated (RE)
""""""""""""""""""""""""""""""""""""""""""""

Per scaling mean rescale_estimated across subjects: :math:`MEAN[\overline{RE}_{(sI, rX, :)}(:, rX)] \to \overline{RE}_{(:, rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 12
   
SD of rescale estimated
""""""""""""""""""""""""""""""""""""""""""""

Per scaling SD of the rescale_estimated across subjects: :math:`STD[\overline{RE}_{(sI, rX, :)}(:, rX)] \to (\sigma_{RE_{(:, rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 13
   
Mean error
""""""""""""""""""""""""""""""""""""""""""""

Per scaling mean error on intra-subject CSA estimation: :math:`MEAN[Error_{(sI, rX, :)}(:, rX)] \to \overline{Error}_{(:, rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 14
   
SD of error
""""""""""""""""""""""""""""""""""""""""""""

Per scaling SD of error on intra-subject CSA estimation: :math:`STD[Error_{(sI, rX, :)}(:, rX)] \to \sigma_{Error_{(:, rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 15
   

Power analysis
=====================
   
Between-group minimum sample size
"""""""""""""""""""""""""""""""""""

The minimum sample size (number of subjects per study arm) necessary to detect an atrophy between groups was computed based on a two-sample (unpaired) bilateral t-test using the following formula (Wang and Ji 2020; Wittes 2002):

:math:`n_{unpaired} = \frac{(z_{α/2} + z_{β})^2(\sigma_{inter\_sI\_r1}+\sigma_{inter\_sI\_rX})^2}{\Delta_{sub} ^2}`

Where :math:`n_{unpaired}` is the minimum sample size required to differentiate between groups with a given power (:math:`z_{β}` corresponds to the power z score, e.g. 80% power gives β=0.2 and :math:`z_{β}`= -0.84) and level of significance (:math:`z_{α/2}` corresponds to the significance level z score, e.g. 5% level of significance gives 𝛂=0.05 and :math:`z_{α/2}`=-1.96), SD is the inter-subject standard deviation of the mean CSA (which was calculated by taking the mean CSA across Monte Carlo samples). :math:`diff_{group}` group is the difference of the mean CSA between the groups.

.. autofunction:: csa_rescale_stat.sample_size  

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 5-7

Within-subject minimum sample size
"""""""""""""""""""""""""""""""""""

the minimum sample size necessary to detect an atrophy in a within-subject (repeated-measures) study was computed based on a two-sample bilateral paired t-test using the following formula (Altmann et al. 2009):

:math:`n_{within\_sub} = \frac{(z_{α/2} + z_{β})^2(\sigma_{diff})^2}{\Delta_{sub} ^2}`
   
Where :math:`\sigma_{diff}` is the standard deviation between longitudinal CSA measures across  subjects and :math:`\Delta_{sub}` is the mean of the difference between longitudinal CSA measures.

.. autofunction:: csa_rescale_stat.sample_size

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 9-10


RESULTS
============

