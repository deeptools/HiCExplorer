#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
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

   findRestSites            Identifies the genomic locations of restriction sites
   hicBuildMatrix           Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads
   hicQC                    Plots QC measures from the output of hicBuildMatrix
   hicCorrectMatrix         Uses iterative correction to remove biases from a Hi-C matrix
   hicFindEnrichedContacts  Identifies enriched Hi-C contacts
   hicCorrelate             Computes and visualizes the correlation of Hi-C matrices
   hicFindTADs	            Identifies Topologically Associating Domains (TADs)
   hicMergeMatrixBins	    Merges consecutive bins on a Hi-C matrix to reduce resolution
   hicPCA                   Computes the principal components (eigenvectors) for A/B compartment analysis
   hicTransform             Computes obs_exp (like Lieberman-Aiden), pearson and covariance matrix for A/B compartment analysis
   hicPlotDistVsCounts	    Plot the decay in interaction frequency with distance
   hicPlotMatrix	        Plots a Hi-C matrix as a heatmap and can add a pca track to it
   hicPlotTADs	            Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)
   hicPlotViewpoint         Plots the number of interactions around a reference point
   hicSumMatrices	        Adds Hi-C matrices of the same size
   hicPlotDistVsCounts	    Plots distance vs. Hi-C counts of corrected data
   hicInfo                  Shows information about a Hi-C matrix (no. of bins, bin length, sum, max, min, etc)
   hicCompareMatrices       Computes difference or ratio between two matrices
   hicAdjustMatrix          Remove, keep or mask specific regions in a matrix
   hicAverageRegions        Sum up defined reference points and a range up/downstream of it. Useful for TAD analysis
   hicPlotAverageRegions    Plot tool for hicAverageRegions
   hicNormalize             Normalizes matrices to 0 1 range or to lowest read coverage
   hicConvertFormat         Converter for different interaction matrix file formats.

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
