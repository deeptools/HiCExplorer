HiCExplorer tools
=================

   .. contents::
       :local:

   +--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
   | tool                           | type             | input files                       | main output file(s)                         | application                                                                       |
   +================================+==================+===================================+=============================================+===================================================================================+
   |:ref:`findRestSites`            | preprocessing    | 1 genome FASTA file               | bed file with restriction site coordinates  | Identifies the genomic locations of restriction sites                             |
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
   |:ref:`hicExport`                | data integration | multiple hiC file formats         | HiC matices/outputs in several formats      | Export matrix to different formats                                                |
   +--------------------------------+------------------+-----------------------------------+---------------------------------------------+-----------------------------------------------------------------------------------+
