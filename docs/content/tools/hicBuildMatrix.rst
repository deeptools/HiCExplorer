.. _hicBuildMatrix:

hicBuildMatrix
==============

.. argparse::
   :ref: hicexplorer.hicBuildMatrix.parse_arguments
   :prog: hicBuildMatrix


Please note that the file type extension for the output matrix (``--outFileName``) must be given! This can be **.h5**, **.cool** or the specializations of cool 
**.mcool**, if the path is given! For **.scool** files please create one **.cool** file per cell and merge it together with scHiCExplorer's scHicMergeToSCool.

Building multicooler matrices
------------------------------

``hicBuildMatrix`` supports building multicooler matrices which are for example needed for visualization with `HiGlass <https://higlass.io/>`__.
To do so, use as out file format either .cool or .mcool and define the desired resolutions as `--binSize`.
``hicBuildMatrix`` builds the interaction matrix for the highest resolution and merges the bins for the lower resolutions.
The lower resolutions need to be an integer multiplicative of the highest resolution.

.. code:: bash

    $ hicBuildMatrix -s forward.bam reverse.bam -o multi_resolution.cool 
      --binSize 10000 20000 50000 100000 --QCfolder QC

Introducing with version 3.5 we support multiple restriction and dangling end sequences, and multiple restriction cut site files. 
Hi-C protocols that use multiple restriction cut enzymes benefit from this and get now an improved QC report.
Version 3.5 adds also the support for a chromosome size file which can help to get interaction matrices with a predefined size. Capture Hi-C or 
single-cell Hi-C data, where it is not guaranteed that reads from all areas of the chromosome are present benefit from this latest improvement.