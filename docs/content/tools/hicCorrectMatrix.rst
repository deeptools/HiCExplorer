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


HiCExplorer version 3.1 changes the way data is transfered from Python to C++ for the KR correction algorithm. With these changes 
the following runtime and peak memory usage on Rao 2014 GM12878 primary + replicate data is possible:

- KR on 25kb: 165 GB, 1:08 h 
- ICE on 25kb: 224 GB, 3:10 h 
- KR on 10kb: 228 GB, 1:42 h
- ICE on 10kb: 323 GB, 4:51 h

- KR on 1kb: 454 GB, 16:50 h
- ICE on 1kb: >600 GB, > 2.5 d (we interrupted the computation and strongly recommend to use KR on this resolution)

For HiCExplorer versions <= 3.0 KR performs as follows:

- KR on 25kb: 159 GB, 57:11 min
- KR on 10kb: >980 GB, -- (out of memory on 1TB node, we do not have access to a node with more memory on our cluster)