"""hicCompartmentalization: a thin layer over the C++ executable hicCompartmentalization.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import numpy as np
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import logging
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import convertNansToZeros
matplotlib.use('Agg')
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicCompartmentalization, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicCompartmentalization')


def get_indices(obs_exp, row):
    indices = np.arange(obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[0],
                        obs_exp.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)[1] + 1)  # end is inclusive
    return indices


def count_interactions(obs_exp, pc1, quantiles_number, offset):
    "Counts the total interaction on obs_exp matrix per quantile and "
    "normalizes it by the number of bins per quantile."
    chromosomes = pc1["chr"].unique()

    interaction_sum = np.zeros((quantiles_number, quantiles_number))
    number_of_bins = np.zeros((quantiles_number, quantiles_number))
    if offset:
        for dist in offset:
            assert (dist >= 0)
            indices = np.arange(0, obs_exp.matrix.shape[0] - dist)
            obs_exp.matrix[indices, indices + dist] = np.nan
            obs_exp.matrix[indices + dist, indices] = np.nan
    for chrom in chromosomes:  # It is only handeling cis contacts
        # pc1_chr = pc1.loc[pc1["chr"] == chrom].reset_index(drop=True)
        # chr_range = obs_exp.getChrBinRange(chrom)
        # chr_submatrix = obs_exp.matrix[chr_range[0]:chr_range[1],
        #                                chr_range[0]:chr_range[1]]
        for qi in range(0, quantiles_number):
            row_indices = pc1.loc[pc1["quantile"] == qi]["bin_id"]
            if row_indices.empty:
                continue
            row_indices = [i for list in row_indices for i in list]
            for qj in range(0, quantiles_number):
                col_indices = pc1.loc[pc1["quantile"] == qj]["bin_id"]
                if col_indices.empty:
                    continue
                col_indices = [i for list in col_indices for i in list]
                submatrix = obs_exp.matrix[np.ix_(row_indices, col_indices)]
                submatrix = submatrix.todense()
                submatrix = submatrix[~np.isnan(submatrix)]  # remove nans
                submatrix = submatrix[~np.isinf(submatrix)]
                interaction_sum[qi, qj] += np.sum(submatrix)
                interaction_sum[qj, qi] += np.sum(submatrix)
                number_of_bins[qi, qj] += submatrix.shape[1]
                number_of_bins[qj, qi] += submatrix.shape[1]
    return interaction_sum / number_of_bins


def within_vs_between_compartments(normalised_sum_per_quantile,
                                   quantiles_number):
    """
    This function computes the interaction between two compartments and whithin
    each compartment. It adds one quantile at the time to compute these values
    on a larger block of the matrix (normalised_sum_per_quantile)until the
    maximum number of quantiles reached.
    """
    within_to_between = []
    for q in range(1, quantiles_number):
        within_comps = normalised_sum_per_quantile[0:q, 0:q].sum() + \
            normalised_sum_per_quantile[quantiles_number - q:quantiles_number,
                                        quantiles_number - q:quantiles_number].sum()

        between_comps = normalised_sum_per_quantile[0:q,
                                                    quantiles_number - q:quantiles_number].sum() + \
            normalised_sum_per_quantile[
            quantiles_number - q:quantiles_number,
            0:q].sum()

        within_to_between.append(within_comps / between_comps)
    return within_to_between


def plot_polarization_ratio(polarization_ratio, plotName, labels,
                            number_of_quantiles):
    """
    Generate a plot to visualize the polarization ratio between A and B
    compartments. It presents how well 2 compartments are seperated.
    """

    for i, r in enumerate(polarization_ratio):
        plt.plot(r, marker="o", label=labels[i])
    plt.axhline(1, c='grey', ls='--', lw=1)
    plt.axvline(number_of_quantiles / 2, c='grey', ls='--', lw=1)
    plt.legend(loc='best')
    plt.xlabel('Quantiles')
    plt.ylabel('signal within comp. / signla between comp.')
    plt.title('compartment polarization ratio')
    plt.savefig(plotName)


def main(args=None):
    """Run the C++ implementation of hicCompartmentalization with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicCompartmentalization')(args)
