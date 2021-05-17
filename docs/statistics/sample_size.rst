Power analysis
=====================
   
Between-group minimum sample size
"""""""""""""""""""""""""""""""""""

The minimum sample size, number of subjects per group (study arm), necessary to detect an atrophy between groups was computed based on a two-sample (unpaired) bilateral t-test using the following formula (Wang and Ji 2020; Wittes 2002):

:math:`n_{unpaired} = \frac{(z_{α/2} + z_{β})^2(\sigma_{(:,r1)}+\sigma_{(:,rX)})^2}{\Delta_{sub} ^2}`

Where :math:`n_{unpaired}` is the minimum sample size required to differentiate between groups with a given power (:math:`z_{β}` corresponds to the power z score, e.g. 80% power gives β=0.2 and :math:`z_{β}`= -0.84) and level of significance (:math:`z_{α/2}` corresponds to the significance level z score, e.g. 5% level of significance gives 𝛂=0.05 and :math:`z_{α/2}`=-1.96), SD is the inter-subject standard deviation of the mean CSA (which was calculated by taking the mean CSA across Monte Carlo samples). :math:`diff_{group}` group is the difference of the mean CSA between the groups.


.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 5-7

Within-subject minimum sample size
"""""""""""""""""""""""""""""""""""

the minimum sample size necessary to detect an atrophy in a within-subject (repeated-measures) study was computed based on a two-sample bilateral paired t-test using the following formula (Altmann et al. 2009):

:math:`n_{within\_sub} = \frac{(z_{α/2} + z_{β})^2(\sigma_{diff})^2}{\Delta_{sub} ^2}`
   
Where :math:`\sigma_{diff}` is the standard deviation between longitudinal CSA measures across  subjects and :math:`\Delta_{sub}` is the mean of the difference between longitudinal CSA measures.

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 352-363
   :emphasize-lines: 9-10


