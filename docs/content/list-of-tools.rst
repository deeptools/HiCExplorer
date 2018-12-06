HiCExplorer tools
=================

.. contents::
    :local:


+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
| tool                           | type             | input files                       | main output file(s)                         | application                                                                       |
+================================+==================+===================================+=============================================+===================================================================================+
|:ref:`findRestSite`             | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicBuildMatrix`           | preprocessing    | 2 BAM/SAM files                   | hicMatrix object                            | Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads    |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCorrectMatrix`         | preprocessing    | hicMatrix object                  | normalized hicMatrix object                 | Uses iterative correction to remove biases from a Hi-C matrix                     |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicMergeMatrixBins`       | preprocessing    | hicMatrix object                  | hicMatrix object                            | Merges consecutives bins on a Hi-C matrix to reduce resolution                    |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicSumMatrices`           | preprocessing    | 2 or more hicMatrix objects       | hicMatrix object                            | Adds Hi-C matrices of the same size                                               |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicNormalize`             | preprocessing    | multiple Hi-C matrices            | multiple Hi-C matrices                      | Normalize data to 0 to 1 range or to smallest total read count                    |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicFindEnrichedContacts`  | analysis         | hicMatrix object                  | hicMatrix object                            | Identifies enriched Hi-C contacts                                                 |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCorrelate`             | analysis         | 2 or more hicMatrix objects       | a heatmap/scatterplot                       | Computes and visualises the correlation of Hi-C matrices                          |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicFindTADs`              | analysis         | hicMatrix object                  | bedGraph file (TAD score), a boundaries.bed | Identifies Topologically Associating Domains (TADs)                               |
|                                |                  |                                   | file, a domains.bed file (TADs)             |                                                                                   |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotMatrix`            | visualization    | hicMatrix object                  | a heatmap of Hi-C contacts                  | Plots a Hi-C matrix as a heatmap                                                  |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotTADs`              | visualization    | hicMatrix object, a config file   | Hi-C contacts on a given region, along with | Plots TADs as a track that can be combined with other tracks                      |
|                                |                  |                                   | other provided signal (bigWig) or regions   | (genes, signal, interactions)                                                     |
|                                |                  |                                   | (bed) file                                  |                                                                                   |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotDistVsCounts`      | visualization    | hicMatrix object                  | log log plot of Hi-C contacts per distance  | Quality control                                                                   |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicConvertFormat`         | data integration | one/multiple Hi-C file formats    | Hi-C matrices/outputs in several formats    | Convert matrix to different formats                                               |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicAdjustMatrix`          | data integration | one Hi-C file formats             | Hi-C matrix                                 | Removes, masks or keeps specified regions of a matrix                             |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicInfo`                  | information      | one or more hicMatrix objects     | Screen info                                 | Prints information about  matrices, like size, maximum, minimux, bin size, etc.   |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPCA`                   | analysis         | one Hi-C matrix                   | bedgraph or bigwig file(s) for each         | Computes for A / B compartments the eigenvectors                                  |
|                                |                  |                                   | eigenvector                                 |                                                                                   |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicTransform`             | analysis         | one Hi-C matrix                   | Hi-C matrix                                 | Computes a obs_exp matrix like Lieberman-Aiden (2009), a pearson correlation      |
|                                |                  |                                   |                                             | matrix and or a covariance matrix. These matrices can be used for plotting.       |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotViewpoint`         | visualization    | one Hi-C matrix                   | A viewpoint plot                            | A plot with the interactions around a reference point or region.                  |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicQC`                    | information      | log files from hicBuildMatrix     | A quality control report                    | Quality control of the created contact matrix.                                    |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCompareMatrices`       | analysis         | two Hi-C matrices                 | one Hi-C matrix                             | Applies diff, ratio or log2ratio on matrices to compare them.                     |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicAverageRegions`        | analysis         | multiple Hi-C matrices            | one npz object                              | Averages the given locations. Visualization with hicPlotAverageRegions            |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotAverageRegions`    | visualization    | one npz file                      | one image                                   | Visualization of hicAverageRegions.                                               |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicMergeTADbins`          | preprocessing    | one Hi-C matrix, one BED file     | one Hi-C matrix                             | Uses a BED file of domains or TAD boundaries to merge the                         |
|                                |                  |                                   |                                             | bin counts of a Hi-C matrix.                                                      |
+--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+


General principles
^^^^^^^^^^^^^^^^^^

A typical HiCExplorer command could look like this:

.. code:: bash

 $ hicPlotMatrix -m myHiCmatrix.h5 \
 -o myHiCmatrix.pdf \
 --clearMaskedBins \
 --region chrX:10,000,000-15,000,000 \
 --vMin -4 --vMax 4 \


You can always see all available command-line options via --help:

.. code:: bash

 $ hicPlotMatrix --help

- Output format of plots should be indicated by the file ending, e.g. ``MyPlot.pdf`` will return a pdf file, ``MyPlot.png`` a png-file.
- Most of the tools that produce plots can also output the underlying data - this can be useful in cases where you don't like the HiCExplorer visualization, as you can then use the data matrices produced by deepTools with your favorite plotting tool, such as R.
- The vast majority of command line options are also available in Galaxy (in a few cases with minor changes to their naming).


Tools for Hi-C data pre-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`findRestSite`
"""""""""""""""""""
:ref:`hicBuildMatrix`
"""""""""""""""""""""
:ref:`hicSumMatrices`
"""""""""""""""""""""
:ref:`hicMergeMatrixBins`
"""""""""""""""""""""""""
:ref:`hicCorrectMatrix`
"""""""""""""""""""""""
:ref:`hicNormalize`
"""""""""""""""""""

Tools for Hi-C QC
^^^^^^^^^^^^^^^^^

:ref:`hicQC`
""""""""""""
:ref:`hicCorrelate`
"""""""""""""""""""
:ref:`hicPlotDistVsCounts`
""""""""""""""""""""""""""
:ref:`hicInfo`
""""""""""""""

Tools for Hi-C data analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`hicCompareMatrices`
"""""""""""""""""""""""""
:ref:`hicFindEnrichedContacts`
""""""""""""""""""""""""""""""
:ref:`hicPCA`
"""""""""""""
:ref:`hicTransform`
"""""""""""""""""""
:ref:`hicAverageRegions`
""""""""""""""""""""""""

Tools for TADs processing
^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`hicFindTADs`
""""""""""""""""""
:ref:`hicMergeTADbins`
""""""""""""""""""""""

Tools for Hi-C and TADs visualization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`hicPlotMatrix`
""""""""""""""""""""
:ref:`hicPlotTADs`
""""""""""""""""""
:ref:`hicPlotViewpoint`
"""""""""""""""""""""""
:ref:`hicAggregateContacts`
"""""""""""""""""""""""""""
:ref:`hicPlotAverageRegions`
""""""""""""""""""""""""""""


Hi-C contact matrix handling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`hicConvertFormat`
"""""""""""""""""""""""
:ref:`hicAdjustMatrix`
""""""""""""""""""""""
