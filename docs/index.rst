.. HiCExplorer documentation master file, created by
   sphinx-quickstart on Wed Sep 23 13:37:43 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


HiCExplorer
===========

Set of programs to process, analyze and visualize Hi-C data
-----------------------------------------------------------

HiCexplorer addresses the common tasks of Hi-C analysis from processing to visualization.

![gallery](https://raw.githubusercontent.com/maxplanck-ie/HiCExplorer/master/examples/images/hicexplorer.png?token=AEu_1VmdSzz0lipVV1DMKuYgYcIjUb4qks5U6zbwwA%3D%3D)

Usage of HiCExplorer
--------------------

Examples of usage

.. code-block:: bash

   # map the reads, each mate individually using for example 
   # bowtie2
   $ bowtie2 -p 16 --local --reorder -x indes_path \
       -U mate_R1.fastq.gz ) 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam
   
   # build matrix from idependently mated read pairs
   $ hicBuildMatrix --samFiles mate1.sam mate2.sam \
                    --binSize 10000 \
                    --restrictionSequence GATC \
                    --outBam hic.bam \
                    -o hic_matrix.npz

   # this creates two files, a bam file containing
   # only the valid Hi-C read pairs
   # and a matrix containing the
   # Hi-C contacts at the given resolution.

   # correct Hi-C matrix
   $ hicCorrectMatrix -m hic_matrix.npz -o hic_corrected.npz

   # visualize the corrected matrix
   ## This needs different options now?? --tracks
   $ hicPlotMatrix -m hic_corrected.npz -o hic_plot.png



The following is the list of tools available in HiCExplorer


========================= ===========================================================================================
tool                      description
========================= ===========================================================================================
hicBuildMatrix            Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads
hicCorrectMatrix          Uses iterative correction to remove biases from a Hi-C matrix
hicFindEnrichedContacts   Identifies enriched Hi-C contacts
hicCorrelate              Computes and visualises the correlation of Hi-C matrices
hicFindTADs               Identifies Topologically Associating Domains (TADs)
hicMergeMatrixBins        Merges consecutives bins on a Hi-C matrix to reduce resolution
hicPlotMatrix             Plots a Hi-C matrix as a heatmap
hicPlotTADs               Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)
hicSumMatrices            Adds Hi-C matrices of the same size
========================= ===========================================================================================


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



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

