Inter-subject
=============

Inter-subject statistics. These statistics are gathered per scaling in the Panda dataframe ``df_rescale``
 
Mean intra-subject SD
"""""""""""""""""""""

Intra-subject SD averaged across subjects.
:math:`MEAN[\sigma_{(sI, rX, :)}(:, rX)] \to \overline{\sigma}_{(:,rX)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 9
   
Mean intra-subject COV
""""""""""""""""""""""

Intra-subject COV averaged across subjects.
:math:`MEAN[COV_{(sI, rX, :)}(:, rX)] \to \overline{COV}_{(:,rX)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 10
   
Inter-subject SD
""""""""""""""""

SD of intra-subject CSA across subjects.
:math:`STD[\overline{CSA}_{(sI, rX, :)}(:, rX)] \to \sigma_{(:,rX)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 11
   
Mean rescale estimated (RE)
"""""""""""""""""""""""""""

rescale_estimated averaged across subjects.
:math:`MEAN[\overline{RE}_{(sI, rX, :)}(:, rX)] \to \overline{RE}_{(:, rX)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 12
   
SD of rescale estimated
"""""""""""""""""""""""

SD of rescale_estimated across subjects.
:math:`STD[\overline{RE}_{(sI, rX, :)}(:, rX)] \to \sigma_{RE_{(:, rX)}}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 13
   
Mean error
""""""""""

error on the intra-subject CSA estimation averaged across subjects.
:math:`MEAN[Error_{(sI, rX, :)}(:, rX)] \to \overline{Error}_{(:, rX)}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 14
   
SD of error
"""""""""""

SD of error on intra-subject CSA estimation across subjects.
:math:`STD[Error_{(sI, rX, :)}(:, rX)] \to \sigma_{Error_{(:, rX)}}`

.. literalinclude:: ../../csa_rescale_stat.py
   :language: python
   :lines: 422-439
   :emphasize-lines: 15
