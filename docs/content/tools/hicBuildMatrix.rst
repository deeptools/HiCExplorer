.. _hicBuildMatrix:

hicBuildMatrix
==============

.. argparse::
   :ref: hicexplorer.hicBuildMatrix.parse_arguments
   :prog: hicBuildMatrix


Multiple interaction matrices
-----------------------------

hicBuildMatrix supports the build of multicooler matrices, to do so, use as outfile format either .cool or .mcool and define
for `--binSize` the desired resolutions. hicBuildMatrix builds the interaction matrix for the highest resolution and
merges the bins to achieve the lower resolutions. The lower resolutions need to be an integer multiplicative of the highest
resolution.

.. code:: bash

    $ hicBuildMatrix -s forward.bam reverse.bam -o multi_resolution.cool 
      --binSize 10000 20000 50000 100000 --QCfolder QC