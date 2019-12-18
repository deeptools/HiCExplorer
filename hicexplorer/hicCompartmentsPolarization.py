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


def parse_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="""
        Rearrange the average interaction frequencies using the first PC values
        to represent the global compartmentalization signal. To our knowledge
        this has been first introduced and implemented by Wibke Schwarzer et
        al. 2017 (Nature. 2017 Nov 2; 551(7678): 51–56)

        $ hicCompartmentsPolarization --obsexp_matrices obsExpMatrix.h5 \
        --pca pc1.bedgraph -o global_signal.png
        """
                                     )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--obsexp_matrices', '-m',
                                help='HiCExplorer matrices in h5/cool format.',
                                nargs="+",
                                required=True)

    parserRequired.add_argument('--pca',  # TODO bedgraph or bigwig
                                help='a PCA vector as a bedgraph file with '
                                'no header. In case of several matrices with '
                                ' different conditions, ie. control'
                                'treatment, the PCA of control can be used. '
                                'Note that only one PCA can be provided.',
                                required=True)

    parserRequired.add_argument('--outputFileName', '-o',
                                help='Plot to represent the polarization of '
                                'A/B compartments.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--quantile', '-q',
                           help='number of quantiles',
                           default=30,
                           type=int)

    parserOpt.add_argument('--outliers',
                           help='precentage of outlier to remove',
                           default=0,
                           type=float)

    parserOpt.add_argument('--outputMatrix',
                           help='output .npz file includes all the '
                                'generated matrices',
                           default=None)

    parserOpt.add_argument('--offset',
                           help='set nan for the distances mentioned as '
                                'offset from main diagonal, only positive '
                                'values are accepted! Example: if '
                                '--offset 0, then values of main diagonal will'
                                ' set to nan and if --offset 0 1 then on top '
                                'of the main diagonal, +1 and -1 diagonal '
                                'values are also set to nan. ',
                           nargs='+',
                           type=int,
                           default=None)

    parserOpt.add_argument('-h',
                           action='help',
                           help='show the help message and exit.')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def count_interactions(obs_exp, pc1, quantiles_number, offset):
    "Counts the total interaction on obs_exp matrix per quantile and "
    "normalizes it by the number of bins per quantile."
    chromosomes = pc1["chr"].unique()
    sum = np.zeros((quantiles_number, quantiles_number))
    count = np.zeros((quantiles_number, quantiles_number))
    for chrom in chromosomes:
        pc1_chr = pc1.loc[pc1["chr"] == chrom].reset_index(drop=True)
        chr_range = obs_exp.getChrBinRange(chrom)

        chr_submatrix = obs_exp.matrix[chr_range[0]:chr_range[1],
                                       chr_range[0]:chr_range[1]]

        if offset:
            for dist in offset:
                assert(dist >= 0)
                indices = np.arange(0, chr_submatrix.shape[0] - dist)
                chr_submatrix[indices, indices + dist] = np.nan
                chr_submatrix[indices + dist, indices] = np.nan

        for qi in range(0, quantiles_number):
            row_indices = pc1_chr.loc[pc1_chr["quantile"] == qi].index
            if row_indices.empty:
                continue
            for qj in range(0, quantiles_number):
                col_indices = pc1_chr.loc[pc1_chr["quantile"] == qj].index
                if col_indices.empty:
                    continue
                submatrix = chr_submatrix[np.ix_(row_indices, col_indices)]
                submatrix = submatrix.todense()
                submatrix = submatrix[~np.isnan(submatrix)]  # remove nans
                sum[qi, qj] += np.sum(submatrix)
                sum[qj, qi] += np.sum(submatrix)
                count[qi, qj] += submatrix.shape[1]
                count[qj, qi] += submatrix.shape[1]

    return sum / count


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
    """
    Main function to generate the polarization plot.
    """
    args = parse_arguments().parse_args(args)

    pc1_bedgraph = pd.read_table(args.pca, header=None, sep="\t")
    pc1 = pd.DataFrame(pc1_bedgraph.values, columns=["chr", "start", "end",
                                                     "pc1"])
    if args.outliers != 0:
        quantile = [args.outliers / 100, (100 - args.outliers) / 100]
        boundaries = np.nanquantile(pc1['pc1'].values.astype(float), quantile)
        quantiled_bins = np.linspace(boundaries[0], boundaries[1],
                                     args.quantile)
    else:
        quantile = [j / (args.quantile - 1) for j in range(0, args.quantile)]
        quantiled_bins = np.nanquantile(pc1['pc1'].values.astype(float),
                                        quantile)

    pc1["quantile"] = np.searchsorted(quantiled_bins,
                                      pc1['pc1'].values.astype(float),
                                      side="right")  # it does return the bin size instead of -1 for the last bin
    polarization_ratio = []
    output_matrices = []
    labels = []
    for matrix in args.obsexp_matrices:
        obs_exp = hm.hiCMatrix(matrix)
        name = ".".join(matrix.split("/")[-1].split(".")[0:-1])
        labels.append(name)

        normalised_sum_per_quantile = count_interactions(obs_exp, pc1,
                                                         args.quantile,
                                                         args.offset)
        if args.outputMatrix:
            output_matrices.append(normalised_sum_per_quantile)

        polarization_ratio.append(within_vs_between_compartments(
                                  normalised_sum_per_quantile,
                                  args.quantile))
    if args.outputMatrix:
        np.savez(args.outputMatrix, [matrix for matrix in output_matrices])
    plot_polarization_ratio(
        polarization_ratio, args.outputFileName, labels, args.quantile)
