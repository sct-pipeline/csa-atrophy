Rescale dataframe
=====================

Across-subject statistics are regrouped in the Panda dataframe ``df_rescale``
 
Mean intra-subject SD
""""""""""""""""""""""""""""""""""""""""""""

Per scaling mean intra-subject SD: :math:`MEAN[\sigma_{(sI, rX, :)}(:, rX)] \to \overline{\sigma}_{(:,rX)}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 9
   
Mean intra-subject COV
""""""""""""""""""""""""""""""""""""""""""""

Per scaling mean intra-subject COV: :math:`MEAN[COV_{(sI, rX, :)}(:, rX)] \to \overline{COV}_{(:,rX)}`

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

Per scaling SD of the rescale_estimated across subjects: :math:`STD[\overline{RE}_{(sI, rX, :)}(:, rX)] \to \sigma_{RE_{(:, rX)}}`

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

Per scaling SD of error on intra-subject CSA estimation: :math:`STD[Error_{(sI, rX, :)}(:, rX)] \to \sigma_{Error_{(:, rX)}}`

.. literalinclude:: ../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 15
