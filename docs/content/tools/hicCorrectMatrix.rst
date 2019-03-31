.. _hicCorrectMatrix:

hicCorrectMatrix
================

.. contents::
    :local:

.. argparse::
   :ref: hicexplorer.hicCorrectMatrix.parse_arguments
   :prog: hicCorrectMatrix


With HiCExplorer version 3.0 we offer an additional Hi-C interaction matrix correction algorithm: Knight-Ruiz.


.. code:: bash

    $ hicCorrectMatrix correct --matrix matrix.cool --correctionMethod KR --chromosomes chrUextra chr3LHet --outFileName corrected_KR.cool


The iterative correction can be used via:

.. code:: bash

    $ hicCorrectMatrix correct --matrix matrix.cool --correctionMethod ICE --chromosomes chrUextra chr3LHet --iterNum 500  --outFileName corrected_ICE.cool --filterThreshold -1.5 5.0