HiCExplorer
===========

Set of programs to process, normalize, analyze and visualize Hi-C data
----------------------------------------------------------------------

HiCexplorer addresses the common tasks of Hi-C analysis from processing to visualization.

.. image:: ./images/hicex2.png


The following is the list of tools available in HiCExplorer
-----------------------------------------------------------


=============================== ==========================================================================================================================================================
tool                            description
=============================== ==========================================================================================================================================================
:ref:`findRestSite`              Identifies the genomic locations of restriction sites
:ref:`hicBuildMatrix`            Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads
:ref:`hicQC`                     Plots QC measures from the output of hicBuildMatrix
:ref:`hicCorrectMatrix`          Uses iterative correction to remove biases from a Hi-C matrix
:ref:`hicFindEnrichedContacts`   Identifies enriched Hi-C contacts
:ref:`hicCorrelate`              Computes and visualises the correlation of Hi-C matrices
:ref:`hicFindTADs`               Identifies Topologically Associating Domains (TADs)
:ref:`hicPCA`                    Computes for A / B compartments the eigenvectors 
:ref:`hicTransform`              Computes a obs_exp matrix like Lieberman-Aiden (2009), a pearson correlation matrix and or a covariance matrix. These matrices can be used for plotting.
:ref:`hicMergeMatrixBins`        Merges consecutive bins on a Hi-C matrix to reduce resolution
:ref:`hicMergeTADbins`           Uses a BED file of domains or TAD boundaries to merge the bin counts of a Hi-C matrix.
:ref:`hicPlotDistVsCounts`       Plot the decay in interaction frequency with distance
:ref:`hicPlotMatrix`             Plots a Hi-C matrix as a heatmap
:ref:`hicPlotTADs`               Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)
:ref:`hicPlotViewpoint`          A plot with the interactions around a reference point or region.  
:ref:`hicSumMatrices`            Adds Hi-C matrices of the same size
:ref:`hicPlotDistVsCounts`       Plots distance vs. Hi-C counts of corrected data
:ref:`hicExport`                 Export matrix to text formats
:ref:`hicInfo`                   Shows information about a Hi-C matrix file (no. of bins, bin length, sum, max, min, etc)
:ref:`hicCompareMatrices`        Computes difference or ratio between two matrices
:ref:`hicLog2Ratio`              Computes the log2 ratio between two matrices.     
=============================== ==========================================================================================================================================================


Getting Help
------------

* For general questions, please use Biostars with Tag `hicexplorer` : `Biostars <https://www.biostars.org/t/hicexplorer/>`_
* For specific questions and feature requests, use the `deepTools mailing list <https://groups.google.com/forum/#!forum/deeptools>`_
* For suggesting changes/enhancements and to report bugs, please create an issue on `our GitHub repository <https://github.com/deeptools/HiCExplorer>`_


Contents:
---------

.. toctree::
   :maxdepth: 2

   content/installation
   content/list-of-tools
   content/example_usage
   content/News

Citation
---------

Please cite HiCExplorer as follows:

Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening,Kin Chung Lam, Bianca Habermann, Asifa Akhtar, Thomas Manke.
**"High-resolution TADs reveal DNA sequences underlying genome organization in flies". Nature Communications**, Volume 9, Article number: 189 (2018), doi: https://doi.org/10.1038/s41467-017-02525-w

.. image:: images/logo_mpi-ie.jpg

This tool suite is developed by the `Bioinformatics Unit <http://www.ie-freiburg.mpg.de/bioinformaticsfac>`_
at the `Max Planck Institute for Immunobiology and Epigenetics <http://www.ie-freiburg.mpg.de/>`_, Freiburg.
