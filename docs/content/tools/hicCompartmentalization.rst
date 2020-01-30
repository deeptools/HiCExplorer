.. _hicCompartmentalization:

hicCompartmentalization
=======================

.. argparse::
   :ref: hicexplorer.hicCompartmentalization.parse_arguments
   :prog: hicCompartmentalization


PCA to compute the global compartmentalization signal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To our knowledge this method has been first introduced by Wibke Schwarzer
et al. 2017 (Nature. 2017 Nov 2; 551(7678): 51–56). In this method, a
global (genome-wide) strength for compartmentalization is computed as
(AA + BB) / (AB + BA)
after rearranging the bins of obs/exp based on their corresponding pc1 values.
For this purpose, first pc1 values are reordered incrementally, then the same
order of bins is used to rearrange the bins in obs/exp matrix.

.. code:: bash

    $ _hicCompartmentalization --obsexp_matrices obsExpMatrix.h5 --pca pc1.bedgraph
      -o global_signal.png
