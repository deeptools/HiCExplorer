.. _hicQC:

hicQC
=====

.. argparse::
   :ref: hicexplorer.hicPrepareQCreport.parse_arguments
   :prog: hicQC

Details
^^^^^^^

hicQC can be used to generate summary tables and plots of the QC measures from several samples generated using :doc:`hicBuildMatrix`. Additionally, an HTML output will be generated.

Examples
^^^^^^^^

:doc:`hicBuildMatrix` generates a QC.log file per processed Hi-C sample in a folder specified in the ``--QCfolder`` argument. The QC measures from several samples can be merged in summary tables and plots. An example usage is:

.. code-block:: bash

    $ hicQC --logfiles ./sample_1/QC.log ./sample_2/QC.log /sample_3/QC.log \
    --labels "Sample 1" "Sample 2" "Sample3" \
    -o QC_plots

.. image:: ../../images/pairs_sequenced.png

Here we can see how many reads were sequenced per sample (pairs considered), how many reads were mappable, unique and of high quality and how many reads passed all quality controls and are thus useful for further analysis (pairs used). All quality controls used for read filtering are explained below.

.. image:: ../../images/unmappable_and_non_unique.png

The next figure contains the fraction of reads with respect to the total number of reads that did not map, that have a low quality score or that didn't map uniquely to the genome.
In our example we can see that sample 3 has the highest fraction of pairs used. We explain the differences between the 3 samples below.

.. image:: ../../images/pairs_discarded.png

This figure contains the fraction of read pairs (with respect to mappable and unique reads) that were discarded when building the Hi-C matrix.

**Dangling ends**
These are reads that start with the restriction site and constitute reads that were digested but no ligated.
Sample 1 in our example has a high fraction of dangling ends (and thus less pairs used). Reasons for this can be inefficient ligation or removal of danging ends during sample processing.
**Duplicated pairs**
These are reads that have the same sequence coming from PCR amplification.
Sample 2 in our example was amplified too much and thus has a very high fraction of duplicated pairs.
**Same fragment**
These are read mates, facing inward, separated by up to 800 bp that do not have a restriction enzyme in between. These read pairs are not valid Hi-C pairs.
**Self circle**
Self circles are defined as pairs within 25kb with 'outward' read orientation.
**Self ligation**
These are read pairs with a restriction site in between that are within 800 bp.

.. image:: ../../images/distance.png

This figure contains the fraction of read pairs (with respect to mappable reads) that compose inter chromosomal, short range (< 20kb) or long range contacts.
Inter chromosomal reads of a wild-type sample are expected to be low. A high fraction of inter chromosomal reads is an indicator of low sample quality or for example cell cycle changes.

.. image:: ../../images/read_orientation.png

The last figure shows the fractions of inward, outward, left or right read pairs (with respect to mappable reads). Deviations from an equal distribution indicates problems during sample preparation.
