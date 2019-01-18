.. _hicNormalize:

hicNormalize
============

.. argparse::
   :ref: hicexplorer.hicNormalize.parse_arguments
   :prog: hicNormalize


Background
^^^^^^^^^^

To be able to compare different Hi-C interaction matrices the matrices need to be normalized to a equal level of read coverage or
value ranges. hicNormalize accomplish this by offering two modes: 0-1 range normalization or a read count normalization.

Usage example
^^^^^^^^^^^^^

Normalize to 0-1 range
""""""""""""""""""""""

norm_range mode of hicNormalize normalizes all reads to the 0 to 1 range i.e. the maximum value of the interaction matrix
becomes 1 and the minimum value 0.

.. code:: bash

    $ hicNormalize -m matrix.cool --normalize norm_range -o matrix_0_1_range.cool

Normalize to smallest read count
""""""""""""""""""""""""""""""""


All matrices are normalized in the way the total read count of each matrix is equal to the read count
of the matrix with the smallest read count of all input matrices.

Example
-------

- matrix.cool with a read count of 10000
- matrix2.cool with a read count of 12010
- matrix3.cool with a read count of 11000

In this example each entry in matrix2.cool and matrix3.cool are normalized with a factor of 12010 / 10000 respective with 11000 / 10000.

.. code:: bash

    $ hicNormalize -m matrix.cool matrix2.cool matrix3.cool --normalize smallest 
      -o matrix_normalized.cool matrix2_normalized.cool matrix3_normalized.cool