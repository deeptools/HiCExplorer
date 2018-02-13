HiCExplorer tools
=================


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
   |:ref:`hicExport`                | data integration | multiple Hi-C file formats        | Hi-C matrices/outputs in several formats    | Export matrix to different formats                                                |
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
   |:ref:`hicLog2Ratio`             | analysis         | two Hi-C matrices                 | one Hi-C matrix                             | Computes the log2 ratio between two matrices.                                     |
   +--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
   |:ref:`hicMergeTADbins`          | preprocessing    | one Hi-C matrix, one BED file     | one Hi-C matrix                             | Uses a BED file of domains or TAD boundaries to merge the                         |
   |                                |                  |                                   |                                             | bin counts of a Hi-C matrix.                                                      |
   +--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
  
.. toctree::
    :glob:
    
    tools/*
