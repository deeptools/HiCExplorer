.. _hicNormalize:

hicNormalize
============

.. argparse::
   :ref: hicexplorer.hicNormalize.parse_arguments
   :prog: hicNormalize


Background
^^^^^^^^^^

To be able to compare different Hi-C interaction matrices the matrices need to be normalized to a equal level of read coverage or
value ranges. This tool helps to acomplish this. 
Each Hi-C interaction matrix has a different read coverage 

Usage example
^^^^^^^^^^^^^

Running 0-1 range
""""""""""""""""""""""""""

hicNormalize -m matrix.cool 