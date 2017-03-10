HiCExplorer
===========

Set of programs to process, normalize, analyze and visualize Hi-C data
----------------------------------------------------------------------

HiCexplorer addresses the common tasks of Hi-C analysis from processing to visualization.

.. image:: ./images/hicex2.png


The following is the list of tools available in HiCExplorer
-----------------------------------------------------------


=============================== ===========================================================================================
tool                            description
=============================== ===========================================================================================
:ref:`findRestSites`             Identifies the genomic locations of restriction sites
:ref:`hicBuildMatrix`            Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads
:ref:`hicCorrectMatrix`          Uses iterative correction to remove biases from a Hi-C matrix
:ref:`hicFindEnrichedContacts`   Identifies enriched Hi-C contacts
:ref:`hicCorrelate`              Computes and visualises the correlation of Hi-C matrices
:ref:`hicFindTADs`               Identifies Topologically Associating Domains (TADs)
:ref:`hicMergeMatrixBins`        Merges consecutives bins on a Hi-C matrix to reduce resolution
:ref:`hicPlotDistVsCounts`       Plot the decay in interaction frequency with distance
:ref:`hicPlotMatrix`             Plots a Hi-C matrix as a heatmap
:ref:`hicPlotTADs`               Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)
:ref:`hicSumMatrices`            Adds Hi-C matrices of the same size
:ref:`hicPlotDistVsCounts`       Plots distance vs. Hi-C counts of corrected data
:ref:`hicExport`                 Export matrix to text formats
=============================== ===========================================================================================


Getting Help
------------

* For general questions, please use Biostars with Tag `hicexplorer` : `Biostars <https://www.biostars.org/t/hicexplorer/>`_
* For specific questions and feature requests, use the `deepTools mailing list <https://groups.google.com/forum/#!forum/deeptools>`_
* For suggesting changes/enhancements and to report bugs, please create an issue on `our GitHub repository <https://github.com/maxplanck-ie/HiCExplorer>`_


Contents:
---------

.. toctree::
   :maxdepth: 2

   content/installation
   content/list-of-tools
   content/example_usage


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
---------

Please cite HiCExplorer as follows:

Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening, Kin Chung Lam, Bianca Habermann, Asifa Akhtar, Thomas Manke
**"High-resolution TADs reveal DNA sequences underlying genome organization in flies". bioRxiv 115063; doi: https://doi.org/10.1101/115063 **

.. image:: images/logo_mpi-ie.jpg

This tool suite is developed by the `Bioinformatics Unit <http://www.ie-freiburg.mpg.de/bioinformaticsfac>`_
at the `Max Planck Institute for Immunobiology and Epigenetics <http://www.ie-freiburg.mpg.de/>`_, Freiburg.
