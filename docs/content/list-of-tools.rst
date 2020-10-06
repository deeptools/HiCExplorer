HiCExplorer tools
=================

Tools for Hi-C data pre-processing
----------------------------------
.. toctree::
    :maxdepth: 1

    tools/hicFindRestSite
    tools/hicBuildMatrix
    tools/hicSumMatrices
    tools/hicMergeMatrixBins
    tools/hicCorrectMatrix
    tools/hicNormalize

Tools for Hi-C QC
-----------------
.. toctree::
    :maxdepth: 1  

    tools/hicQuickQC
    tools/hicQC
    tools/hicCorrelate
    tools/hicPlotDistVsCounts
    tools/hicInfo

Tools for Hi-C data analysis
----------------------------
.. toctree::
    :maxdepth: 1  

    tools/hicCompareMatrices
    tools/hicPCA
    tools/hicTransform
    tools/hicAverageRegions
    tools/hicDetectLoops
    tools/hicValidateLocations
    tools/hicMergeLoops
    tools/hicHyperoptDetectLoops
    tools/hicHyperoptDetectLoopsHiCCUPS
    tools/hicCompartmentalization
    tools/hicPlotSVL

Tools for TADs processing
-------------------------
.. toctree::
    :maxdepth: 1  

    tools/hicFindTADs
    tools/hicMergeDomains
    tools/hicDifferentialTAD
    tools/hicMergeTADbins

Tools for Hi-C and TADs visualization
-------------------------------------
.. toctree::
    :maxdepth: 1  

    tools/hicPlotMatrix
    tools/hicPlotTADs
    tools/hicPlotViewpoint
    tools/hicAggregateContacts
    tools/hicPlotAverageRegions

Hi-C contact matrix handling
----------------------------
.. toctree::
    :maxdepth: 1  

    tools/hicConvertFormat
    tools/hicAdjustMatrix

Capture Hi-C analysis
---------------------
.. toctree::
    :maxdepth: 1  

    tools/chicQualityControl
    tools/chicViewpointBackgroundModel
    tools/chicViewpoint
    tools/chicSignificantInteractions
    tools/chicAggregateStatistic
    tools/chicDifferentialTest
    tools/chicPlotViewpoint

.. contents::
    :local:


+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
| tool                                 | type             | input files                       | main output file(s)                         | application                                                                       |
+======================================+==================+===================================+=============================================+===================================================================================+
|:ref:`hicFindRestSite`                | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicBuildMatrix`                 | preprocessing    | 2 BAM/SAM files                   | hicMatrix object                            | Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads    |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCorrectMatrix`               | preprocessing    | hicMatrix object                  | normalized hicMatrix object                 | Uses iterative correction or Knight-Ruiz to remove biases from a Hi-C matrix      |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicMergeMatrixBins`             | preprocessing    | hicMatrix object                  | hicMatrix object                            | Merges consecutive bins on a Hi-C matrix to reduce resolution                     |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicSumMatrices`                 | preprocessing    | 2 or more hicMatrix objects       | hicMatrix object                            | Adds Hi-C matrices of the same size                                               |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicNormalize`                   | preprocessing    | multiple Hi-C matrices            | multiple Hi-C matrices                      | Normalize data to 0 to 1 range or to smallest total read count                    |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCorrelate`                   | analysis         | 2 or more hicMatrix objects       | a heatmap/scatter plot                      | Computes and visualizes the correlation of Hi-C matrices                          |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicFindTADs`                    | analysis         | hicMatrix object                  | bedGraph file (TAD score), a boundaries.bed | Identifies Topologically Associating Domains (TADs)                               |
|                                      |                  |                                   | file, a domains.bed file (TADs)             |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicMergeDomains`                | analysis         | multiple TAD domain files         | tad domain file with merged tad locations   | Merges detect TADs locations of different resolutions; hierarchical relation      |
|                                      |                  |                                   | multiple files with plotted TAD relations   | between TADs as multiple plots                                                    |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicDifferentialTAD`             | analysis         | two Hi-C matrices                 | two diff_tad files: accepted H0 and         | Identifies differential Topologically Associating Domains (TADs) between          |
|                                      |                  | one TAD domain file               | rejected H0. Similar to BED6                | two Hi-C matrices                                                                 |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotMatrix`                  | visualization    | hicMatrix object                  | a heatmap of Hi-C contacts                  | Plots a Hi-C matrix as a heatmap                                                  |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotTADs`                    | visualization    | hicMatrix object, a config file   | Hi-C contacts on a given region, along with | Plots TADs as a track that can be combined with other tracks                      |
|                                      |                  |                                   | other provided signal (bigWig) or regions   | (genes, signal, interactions)                                                     |
|                                      |                  |                                   | (bed) file                                  |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotDistVsCounts`            | visualization    | hicMatrix object                  | log log plot of Hi-C contacts per distance  | Quality control                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicConvertFormat`               | data integration | one/multiple Hi-C file formats    | Hi-C matrices/outputs in several formats    | Convert matrix to different formats                                               |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicAdjustMatrix`                | data integration | one Hi-C file formats             | Hi-C matrix                                 | Removes, masks or keeps specified regions of a matrix                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicInfo`                        | information      | one or more hicMatrix objects     | Screen info                                 | Prints information about  matrices, like size, maximum, minimum, bin size, etc.   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPCA`                         | analysis         | one Hi-C matrix                   | bedgraph or bigwig file(s) for each         | Computes for A / B compartments the eigenvectors                                  |
|                                      |                  |                                   | eigenvector                                 |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicTransform`                   | analysis         | one Hi-C matrix                   | Hi-C matrix                                 | Computes a obs_exp matrix like Lieberman-Aiden (2009), a pearson correlation      |
|                                      |                  |                                   |                                             | matrix and or a covariance matrix. These matrices can be used for plotting.       |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotViewpoint`               | visualization    | one Hi-C matrix                   | A viewpoint plot                            | A plot with the interactions around a reference point or region.                  |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicQC`                          | information      | log files from hicBuildMatrix     | A quality control report                    | Quality control of the created contact matrix.                                    |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicQuickQC`                     | information      | 2 BAM/SAM files                   | An estimated quality control report         | Estimated quality report of the Hi-C data.                                        |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCompareMatrices`             | analysis         | two Hi-C matrices                 | one Hi-C matrix                             | Applies diff, ratio or log2ratio on matrices to compare them.                     |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicAverageRegions`              | analysis         | multiple Hi-C matrices            | one npz object                              | Averages the given locations. Visualization with hicPlotAverageRegions            |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicDetectLoops`                 | analysis         | one Hi-C matrices                 | bedgraph file with loop locations           | Detects enriched regions. Visualization with hicPlotmatrix and --loop parameter.  |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicValidateLocations`           | analysis         | one loop, one protein peak file   | bedgraph file with matched loop locations,  | Matches loop locations with protein peak positions                                |
|                                      |                  |                                   | one file with loop / protein statistics     |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicMergeLoops`                  | analysis         | multiple loop files               | bedgraph file with merged loop locations    | Merges detect loop locations of different resolutions                             |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicHyperoptDetectLoops`         | analysis         | one Hi-C matrix, protein peaks    | best parameter setting                      | Search for best parameter setting for hicDetectLoops                              |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicHyperoptDetectLoopsHiCCUPS`  | analysis         | one Hi-C matrix, protein peaks    | best parameter setting                      | Search for best parameter setting for Juicer's HiCCUPS                            |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicCompartmentalization`        | visualization    | one Hi-C interaction matrix       | one image                                   | The global compartmentalization signal.                                           |
|                                      |                  | one PCA bedgraph file             | polarization plot                           |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotAverageRegions`          | visualization    | one npz file                      | one image                                   | Visualization of hicAverageRegions.                                               |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicPlotSVL`                     | analysis         | one / multiple Hi-C matrices      | one image, p-values file, raw data file     | Computes short/long range contacts; a box plot, a p-value and raw data file       |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`hicMergeTADbins`                | preprocessing    | one Hi-C matrix, one BED file     | one Hi-C matrix                             | Uses a BED file of domains or TAD boundaries to merge the                         |
|                                      |                  |                                   |                                             | bin counts of a Hi-C matrix.                                                      |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicQualityControl`             | preprocessing    | Hi-C matrices                     | two plots                                   | Checks for sparsity of viewpoints and removes them if too sparse.                 |
|                                      |                  | reference point BED file          | accepted reference point BED file           |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicViewpointBackgroundModel`   | preprocessing    | Hi-C matrices                     | background model file                       | Creates a background model for all given samples and reference points.            |
|                                      |                  | reference point BED file          |                                             |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicViewpoint`                  | preprocessing    | Hi-C matrices                     | viewpoint file(s)                           | Creates per sample per viewpoint one viewpoint file.                              |
|                                      |                  | background model file             |                                             |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicSignificantInteractions`    | preprocessing    | viewpoint file(s)                 | significant interaction file(s)             | Detects significant interactions per viewpoint based on the background and        |
|                                      | analysis         | background model file             | target file(s)                              | neighborhood merging via x-fold and loose p-values.                               |
|--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicAggregateStatistic`         | preprocessing    | viewpoint files(s)                | aggregated file(s) for differential test    | Aggregates for one viewpoint of two samples via a target file the locations to    |
|                                      |                  | target file (s)                   |                                             | test for differential interactions.                                               |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicDifferentialTest`           | analysis         | aggregated file(s) of two samples | H0_accepted-, H0_rejected-files             | Tests with chi2 or fisher for differential interactions of two samples.           |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
|:ref:`chicPlotViewpoint`              | visualization    | viewpoint file(s)                 | one image per viewpoint                     | Visualization of a viewpoint.                                                     |
|                                      |                  | differential expression file(s)   |                                             |                                                                                   |
|                                      |                  | significant interactions file(s)  |                                             |                                                                                   |
+--------------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+

 
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
