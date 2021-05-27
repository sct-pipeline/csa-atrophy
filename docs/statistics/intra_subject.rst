Intra-subject
=============

Intra-subject statistics. These statistics are gathered per rescaling and per subject in the Panda dataframe ``df_sub``:

Intra-subject CSA (CSA estimation)
""""""""""""""""""""""""""""""""""

CSA averaged across transformations.

:math:`MEAN[CSA(sI, rX, :)] \to \overline{CSA}_{(sI,rX,:)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 10
    
Intra-subject SD
""""""""""""""""

SD of CSA across transformations.

:math:`STD[CSA(sI, rX, :)] \to \sigma_{(sI,rX,:)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 11
    
    
Intra-subject COV
""""""""""""""""""

COV of CSA across transformations.

:math:`\frac{STD[CSA(sI, rX, :)]}{MEAN[CSA(sI, rX, :)]} \to COV_{(sI,rX,:)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 12
    
Rescale estimation (RE)
"""""""""""""""""""""""

ratio of the atrophied CSA divided by the un-rescaled CSA averaged across transformations (gives an estimation of the applied scaling).

:math:`MEAN \left[\frac{CSA(sI, rX, :)}{CSA(sI, r1, :)}\right] \to \overline{RE}_{(sI,rX,:)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 14 
    
Error
"""""

mean absolute error on CSA estimation averaged across transformations.

:math:`MEAN[CSA(sI, rX, :)] - rX^2  \times MEAN[CSA(sI, r1, :)] \to Error_{(sI,rX,:)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 15
