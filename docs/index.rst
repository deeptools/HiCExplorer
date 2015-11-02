HiCExplorer
===========

Set of programs to process, analyze and visualize Hi-C data
-----------------------------------------------------------

HiCexplorer addresses the common tasks of Hi-C analysis from processing to visualization.

.. image:: ./images/hicexplorer.png


Usage of HiCExplorer
--------------------

Reads mapping
^^^^^^^^^^^^^

Mates have to be mapped individually to avoid mapper specific heuristics designed
for standard paired-end libraries.

We have used the HiCExplorer sucessfuly with `bowtie2` and `hisat2`. And important parameter is
`--reorder` which tells bowtie2 or hisat2 to ouput the *sam* files in the **exact** same order as
in the *.fastq* files.


.. code-block:: bash

   # map the reads, each mate individually using
   # for example bowtie2
   $ bowtie2 -p 16 --local --reorder -x index_path \
       -U mate_R1.fastq.gz ) 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam

   $ bowtie2 -p 16 --local --reorder -x index_path \
       -U mate_R2.fastq.gz ) 2>>mate_R2.log | samtools view -Shb - > mate_R2.bam


Creation of a Hi-C matrix
^^^^^^^^^^^^^^^^^^^^^^^^^

Once the reads have been mapped the Hi-C matrix can be built. For this, the minimal
extra information required is the `binSize` used for the matrix. Is it best
to enter a low number like 10.000 because lower
resolution matrices can be easily constructed using `hicMergeMatrixBins`. Matrices
at restriction fragment resolution can be created by providing a file
containing the restriction sites, this file can be created with the tool
`findRestSite` that is part of HiCExplorer.

.. code-block:: bash

   # build matrix from independently mated read pairs
   $ hicBuildMatrix --samFiles mate_R1.bam mate_R2.bam \
                    --binSize 10000 \
                    --restrictionSequence GATC \
                    --outBam hic.bam \
                    -o hic_matrix.npz > build_matrix.log


`hicBuildMatrix` creates two files, a bam file containing only the valid Hi-C read pairs and a matrix containing the
Hi-C contacts at the given resolution. The bam file is useful to check the quality of the
Hi-C library on the genome browser. A good Hi-C library should contain piles of reads near
the restriction fragment sites. The log output contains the number of valid pairs, duplicated pairs and
other useful information. Usually, only 25% of the reads are valid and used to build the Hi-C matrix mostly
because of the reads that are on repetitive regions that need to be discarded.


Correction of Hi-C matrix
^^^^^^^^^^^^^^^^^^^^^^^^^

The Hi-C matrix has to be corrected to remove GC, open chromatin biases and, most importantly,
to normalize the number of restriction sites per bin. Because a fraction of bins from repetitive regions
contain few contacts it is necessary to filter those regions first. Also, in mammalian genomes some regions
enriched by reads should be discarded. To aid in the filtering of regions `hicCorrectMatrix` generates a
diagnostic plot as follows:

.. code-block:: bash

   $ hicCorrectMatrix diagnostic_plot -m hic_matrix.npz -o hic_corrected.npz


The plot should look like this:

.. figure:: ./images/diagnostic_plot.png
    :scale: 70 %
    :align: center

    Histogram of the number of counts per bin.



For the upper threshold is only important to remove very high outliers and thus a value of 5 could be used.
For the lower threshold it is recommended to use a value between -2 and -1. What it not desired is to
try to correct low count bins which could result simply in an amplification of noise. For the upper threshold
is not so concerning because those bins will be scaled down.

Once the thresholds have been decided, the matrix can be corrected

.. code-block:: bash

   # correct Hi-C matrix
   $ hicCorrectMatrix -m hic_matrix.npz --filterThreshold -1.5 5 -o hic_corrected.npz


Visualization of results
^^^^^^^^^^^^^^^^^^^^^^^^

There are two ways to see the resulting matrix, one using `hicPlotMatrix` and the
other is using `hicPlotTADs`. The first one allows the visualization over large regions
while the second one is preferred to see specific parts together with other information,
for example genes or bigwig tracks.

Because of the large differences in counts found int he matrix, it is better to
plot the counts using the `--log1p` option.

.. code-block:: bash

   $ hicPlotMatrix -m hic_corrected.npz -o hic_plot.png --region 1:20000000-80000000 --log1p



.. figure:: ./images/corrected_matrix_example.png
    :scale: 90 %
    :align: center

    Corrected Hi-C counts in log scale.


The following is the list of tools available in HiCExplorer


=============================== ===========================================================================================
tool                            description
=============================== ===========================================================================================
:ref:`findRestSites`              Identifies the genomic locations of restriction sites
:ref:`hicBuildMatrix`            Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads
:ref:`hicCorrectMatrix`          Uses iterative correction to remove biases from a Hi-C matrix
:ref:`hicFindEnrichedContacts`   Identifies enriched Hi-C contacts
:ref:`hicCorrelate`              Computes and visualises the correlation of Hi-C matrices
:ref:`hicFindTADs`               Identifies Topologically Associating Domains (TADs)
:ref:`hicMergeMatrixBins`        Merges consecutives bins on a Hi-C matrix to reduce resolution
:ref:`hicPlotMatrix`             Plots a Hi-C matrix as a heatmap
:ref:`hicPlotTADs`               Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)
:ref:`hicSumMatrices`            Adds Hi-C matrices of the same size
:ref:`hicExport`                 Export matrix to text formats
=============================== ===========================================================================================


Installation
------------


Assuming that python is already installed, the easiest way to get the latest HiCExplorer
version is to use the ``pip`` command. If the  ``scipy``, ``numpy`` and ``matplotlib`` dependencies are not yet installed
they can be easily installed using `anaconda <http://continuum.io/downloads>`_.

.. code-block:: bash

   $ pip install git+ssh://git@github.com/maxplanck-ie/HiCExplorer.git


Otherwise, HiCExplorer can be downloaded and installed:

.. code-block:: bash

   $ git clone https://github.com/maxplanck-ie/HiCExplorer.git
   $ cd HiCExplorer
   $ python setup.py install

By default, the script will install the python library and executable
codes globally, which means you need to be root or administrator of
the machine to complete the installation. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script.


Mailing list
------------

If you have questions, requests, or bugs to report, please email the
`deepTools mailing list <https://groups.google.com/forum/#!forum/deeptools>`_



This tool suite is developed by the `Bioinformatics Unit <http://www.ie-freiburg.mpg.de/bioinformaticsfac>`_
at the `Max Planck Institute for Immunobiology and Epigenetics <http://www.ie-freiburg.mpg.de/>`_, Freiburg.


Contents:
---------

.. toctree::
   :maxdepth: 2

   content/list-of-tools
   content/mES-HiC_analysis


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

