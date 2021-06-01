Inter-subject
=============

Inter-subject statistics. These statistics are gathered per scaling in the Panda dataframe ``df_rescale``
 
Mean intra-subject SD
"""""""""""""""""""""

Intra-subject SD averaged across subjects.

:math:`\mu_s \{ \sigma_t \{ CSA_{rX} \} \}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 9
   
Mean intra-subject COV
""""""""""""""""""""""

Intra-subject COV averaged across subjects.

:math:`\mu_s \{ COV_t \{ CSA_{rX} \} \}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 10
   
Inter-subject SD
""""""""""""""""

SD of intra-subject CSA across subjects.

:math:`\sigma_s \{ \mu_t \{ CSA_{rX} \} \}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 11
   
Mean rescale estimated (RE)
"""""""""""""""""""""""""""

rescale_estimated averaged across subjects.

:math:`\mu_s \left \{ \mu_t \left\{ \frac{CSA_{rX}}{CSA_{r1}} \right\}\right\}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 12
   
SD of rescale estimated
"""""""""""""""""""""""

SD of rescale_estimated across subjects.

:math:`\sigma_s \left\{\mu_t \left\{ \frac{CSA_{rX}}{CSA_{r1}} \right\}\right\}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 13
   
Mean error
""""""""""

error on the intra-subject CSA estimation averaged across subjects.

:math:`\mu_s \{ \mu_t \{ CSA_{rX} \} - \mu_t \{ CSA_{r1} \cdot (rX)^2 \} \}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 14
   
SD of error
"""""""""""

SD of error on intra-subject CSA estimation across subjects.

:math:`\sigma_s \{ \mu_t \{ CSA_{rX} \} - \mu_t \{ CSA_{r1} \cdot (rX)^2 \} \}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 15
