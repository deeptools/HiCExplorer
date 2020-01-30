.. _hicAdjustMatrix:

hicAdjustMatrix
===============

.. argparse::
   :ref: hicexplorer.hicAdjustMatrix.parse_arguments
   :prog: hicAdjustMatrix

hicAdjustMatrix can mask, remove or keep defined regions from a BED file or given chromosomes. This can be useful to create a smaller Hi-C matrix of a single chromosome for e.g. testing purposes or to remove repetitive chromosome ends before calculating A/B compartments using :doc:`hicPCA`.

Example usages
--------------

.. code:: bash

    $ hicAdjustMatrix -m matrix.cool --action keep --chromosomes chr1 -o matrix_chr1.cool

.. code:: bash

    $ hicAdjustMatrix -m matrix.cool --action mask --regions mask_regions.bed -o matrix_masked.cool

mask_regions.bed

.. code::

    chr1    10  30
    chr1    50  300
