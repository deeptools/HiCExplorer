import numpy as np
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import dia_matrix
import logging
log = logging.getLogger(__name__)

from hicmatrix import HiCMatrix as hm

from hicexplorer._version import __version__

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
        Rearrange the average interaction frequencies using the first PC values
        to represent the global compartmentalization signal. To our knowledge
        this has been first introduced and implemented by Wibke Schwarzer et
        al. 2017 (Nature. 2017 Nov 2; 551(7678): 51–56)

        $ hicCompartmentsPolarization --obsexp_matrices obsExpMatrix.h5 --pca pc1.bedgraph\
        -o global_signal.png
        """
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--obsexp_matrices', '-m',
                                help='HiCExplorer matrices in h5/cool format.',
                                nargs="+",
                                required=True)

    parserRequired.add_argument('--pca',
                                help='a PCA vector as a bedgraph file with '
                                'no header. In case of several matrices with '
                                ' different conditions, ie. control'
                                'treatment, the PCA of control can be '
                                'used. Note that only one PCA can be provided.',
                                required=True)

    parserRequired.add_argument('--outputFileName', '-o',
                                help='Plot to represent the polarization of '
                                'A/B compartments.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--quantile', '-q',
                           help='number of quantiles',
                           default=30)

    parserOpt.add_argument('--outliers',
                           help='precentage of outlier to remove',
                           default=0)

    parserOpt.add_argument('--outputMatrix',
                           help='output .npz file includes all the generated matrices',
                           default=None)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit.')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def count_interactions(obs_exp, pc1, quantiles_number):
    "Counts the total interaction on obs_exp matrix per quantile and "
    "normalizes it by the number of bins per quantile."
    chromosomes = pc1["chr"].unique()
    normalised_sum_per_quantile = np.zeros((quantiles_number, quantiles_number))

    for chrom in chromosomes:
        pc1_chr = pc1.loc[pc1["chr"] == chrom].reset_index(drop=True)
        chr_range = obs_exp.getChrBinRange(chrom)
        chr_submatrix = obs_exp.matrix[chr_range[0]:chr_range[1],
                                       chr_range[0]:chr_range[1]]
        chr_submatrix = chr_submatrix.todense()

        np.fill_diagonal(chr_submatrix, 0)

        for qi in range(0, quantiles_number):
            row_indices = pc1_chr.loc[pc1_chr["quantile"] == qi].index

            for qj in range(0, quantiles_number):
                col_indices = pc1_chr.loc[pc1_chr["quantile"] == qj].index
                data = chr_submatrix[row_indices, :][:, col_indices]

                if data.shape[0] * data.shape[1] != 0:
                    normalised_sum_per_quantile[qi, qj] += (np.sum(data) /
                                                            (data.shape[0] * data.shape[1]))

    return normalised_sum_per_quantile


def within_vs_between_compartments(normalised_sum_per_quantile, quantiles_number):
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
    Generates a plot to visualize the polarization ratio between A and B
    compartments.
    """
    for i, r in enumerate(polarization_ratio):
        plt.plot(r, marker="o", label=labels[i])
    plt.axhline(1, c='grey', ls='--', lw=1)
    plt.axvline(number_of_quantiles / 2, c='grey', ls='--', lw=1)
    plt.legend(loc='upper right')
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
        q0, qn = np.nanquantile(pc1['pc1'].values.astype(float), quantile)
        q_bins = np.linspace(q0, qn, args.quantile)
    else:
        quantile = [j / (args.quantile - 1) for j in range(0, args.quantile)]
        q_bins = np.nanquantile(pc1['pc1'].values.astype(float), quantile)

    pc1["quantile"] = np.searchsorted(q_bins, pc1['pc1'].values.astype(float))
    polarization_ratio = []
    output_matrices = []
    labels = []
    for matrix in args.obsexp_matrices:
        obs_exp = hm.hiCMatrix(matrix)
        name = ".".join(matrix.split("/")[-1].split(".")[0:-1])
        labels.append(name)

        normalised_sum_per_quantile = count_interactions(obs_exp, pc1,
                                                         args.quantile)
        if args.outputMatrix:
            output_matrices.append(normalised_sum_per_quantile)

        polarization_ratio.append(within_vs_between_compartments(normalised_sum_per_quantile,
                                                                 args.quantile))
    if args.outputMatrix:
        np.savez(args.outputMatrix, [matrix for matrix in output_matrices])
    plot_polarization_ratio(polarization_ratio, args.outputFileName, labels, args.quantile)
