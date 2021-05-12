Subject dataframe
==================

Per-subject statistics are regrouped in the Panda dataframe ``df_sub``:

Intra-subject CSA (CSA estimation)
""""""""""""""""""""""""""""""""""""

Per rescaling and per subject, computation of the mean CSA across transformations: :math:`MEAN[CSA(sI, rX, :)] \to \overline{CSA}_{(sI,rX,:)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 10
    
Intra-subject SD
""""""""""""""""""

Per rescaling and per subject, SD computation of CSA across transformations: :math:`STD[CSA(sI, rX, :)] \to \sigma_{(sI,rX,:)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 11
    
    
Intra-subject COV
""""""""""""""""""

Per rescaling and per subject, COV computation of CSA across transformations: :math:`\frac{STD[CSA(sI, rX, :)]}{MEAN[CSA(sI, rX, :)]} \to COV_{(sI,rX,:)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 12
    
Rescale estimation (RE)
"""""""""""""""""""""""""

Per rescaling and per subject, mean ratio computation of the atrophied CSA divided by the un-rescaled CSA across transformations: :math:`MEAN \left[\frac{CSA(sI, rX, :)}{CSA(sI, r1, :)}\right] \to \overline{RE}_{(sI,rX,:)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 14 
    
Error
""""""""""""""""""

Per rescaling and per subject, computation of the mean absolute error on CSA estimation across transformations: :math:`MEAN[CSA(sI, rX, :)] - rX^2  \times MEAN[CSA(sI, r1, :)] \to Error_{(sI,rX,:)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 401-419
   :emphasize-lines: 15
