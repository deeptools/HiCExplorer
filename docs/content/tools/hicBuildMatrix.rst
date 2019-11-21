.. _hicBuildMatrix:

hicBuildMatrix
==============

.. argparse::
   :ref: hicexplorer.hicBuildMatrix.parse_arguments
   :prog: hicBuildMatrix


Building multicooler matrices
------------------------------

``hicBuildMatrix`` supports building multicooler matrices which are for example needed for visualization with HiGlass https://higlass.io/
To do so, use as outfile format either .cool or .mcool and define the desired resolutions as `--binSize`.
``hicBuildMatrix`` builds the interaction matrix for the highest resolution and merges the bins for the lower resolutions.
The lower resolutions need to be an integer multiplicative of the highest resolution.

.. code:: bash

    $ hicBuildMatrix -s forward.bam reverse.bam -o multi_resolution.cool 
      --binSize 10000 20000 50000 100000 --QCfolder QC
