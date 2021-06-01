Intra-subject
=============

Intra-subject statistics. These statistics are gathered per rescaling and per subject in the Panda dataframe ``df_sub``:

Intra-subject CSA (CSA estimation)
""""""""""""""""""""""""""""""""""

CSA averaged across transformations.

:math:`\mu_t \{{CSA_{sI, rX}}}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 10
    
Intra-subject SD
""""""""""""""""

SD of CSA across transformations.

:math:`\sigma_t \{CSA_{sI,rX}}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 11
    
    
Intra-subject COV
""""""""""""""""""

COV of CSA across transformations.

:math:`COV_t \{CSA_{sI,rX}}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 12
    
Rescale estimation (RE)
"""""""""""""""""""""""

ratio of the atrophied CSA divided by the un-rescaled CSA averaged across transformations (gives an estimation of the applied scaling).

:math:`\mu_t \left\{ \frac{CSA_{sI, rX}}{CSA_{sI, rX}} \right\}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 14 
    
Error
"""""

mean absolute error on CSA estimation averaged across transformations.

:math:`\mu_t \{{CSA_{sI, rX}}} - \mu_t\{{CSA_{sI, r1}}} \dot (rX)^2`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 15
