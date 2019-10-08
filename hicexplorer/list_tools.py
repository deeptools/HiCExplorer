#!/usr/bin/env python
# -*- coding: utf-8 -*-
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import sys
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
HiCExplorer addresses the common tasks of Hi-C analysis from processing to visualization.
Each tool should be called by its own name as in the following example:

 $ hicPlotMatrix -m hic_matrix.h5 -o plot.pdf

If you find HiCExplorer useful for your research please cite as:

Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening, Kin Chung Lam,
Bianca Habermann, Asifa Akhtar, Thomas Manke.
"High-resolution TADs reveal DNA sequences underlying genome organization in flies".
Nature Communications, Volume 9, Article number: 189 (2018), doi: https://doi.org/10.1038/s41467-017-02525-w

Joachim Wolff, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach,
Thomas Manke, Rolf Backofen, Fidel Ramírez, Björn A Grüning.
"Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and visualization",
Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W11–W16, doi: https://doi.org/10.1093/nar/gky504

The following is the list of tools:

findRestSite                 Identifies the genomic locations of restriction sites
hicBuildMatrix               Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads
hicQuickQC                   Estimates the quality of Hi-C dataset
hicQC                        Plots QC measures from the output of hicBuildMatrix
hicCorrectMatrix             Uses iterative correction to remove biases from a Hi-C matrix
hicDetectLoops               Identifies enriched Hi-C contacts
hicCorrelate                 Computes and visualizes the correlation of Hi-C matrices
hicFindTADs                  Identifies Topologically Associating Domains (TADs)
hicPCA                       Computes for A / B compartments the eigenvectors
hicTransform                 Computes a obs_exp matrix like Lieberman-Aiden (2009), a pearson correlation matrix and or a covariance matrix. These matrices can be used for plotting.
hicMergeMatrixBins           Merges consecutive bins on a Hi-C matrix to reduce resolution
hicMergeTADbins              Uses a BED file of domains or TAD boundaries to merge the bin counts of a Hi-C matrix.
hicPlotDistVsCounts          Plot the decay in interaction frequency with distance
hicPlotMatrix                Plots a Hi-C matrix as a heatmap
hicPlotTADs                  Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)
hicPlotViewpoint             A plot with the interactions around a reference point or region.
hicAggregateContacts         A tool that allows plotting of aggregated Hi-C sub-matrices of a specified list of positions.
hicSumMatrices               Adds Hi-C matrices of the same size
hicPlotDistVsCounts          Plots distance vs. Hi-C counts of corrected data
hicInfo                      Shows information about a Hi-C matrix file (no. of bins, bin length, sum, max, min, etc)
hicCompareMatrices           Computes difference or ratio between two matrices
hicAverageRegions            Computes the average of multiple given regions, usually TAD regions
hicPlotAverageRegions        visualization of hicAverageRegions
hicNormalize                 Normalizes the given matrices to 0-1 range or the smallest read coverage
hicConvertFormat             Converts between different Hi-C interaction matrices
hicAdjustMatrix              Keeps, removes or masks regions in a Hi-C matrix
hicValidateLocations         Compare the loops with known peak protein locations
hicMergeLoops                Merges loops of different resolutions
hicCompartmentsPolarization  Compute the global compartmentalization signal
chicQualityControl           Quality control for cHi-C data
chicViewpointBackgroundModel Background model computation for cHi-C analysis
chicViewpoint                Computation of all viewpoints based on background model for cHi-C analysis
chicSignificantInteractions  Detection of significant interactions per viewpoint based on background model
chicAggregateStatistic       Compiling of target regions for two samples as input for differential analysis
chicDifferentialTest         Differential analysis of interactions of two samples
chicPlotViewpoint            Plotting of viewpoint with background model and highlighting of significant and differential regions
hicPlotSVL                   Computing short vs long range contacts and plotting the results

For more information visit: http://hicexplorer.readthedocs.org
""")

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    return args


def main(args=None):
    if args is None and len(sys.argv) == 1:
        args = ["--help"]
    process_args(args)
