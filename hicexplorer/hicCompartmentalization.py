"""hicCompartmentalization, implemented in C++.

main() runs the C++ executable hicCompartmentalization and parse_arguments()
is built from its --help-json (hicexplorer._cpp). get_indices,
count_interactions and within_vs_between_compartments keep the name, the
signature and the behaviour they had in HiCExplorer 3.7.6; their bodies call
hicexplorer._hicx.hicCompartmentalization
(cpp/python/src/api/hicCompartmentalization.cpp). plot_polarization_ratio
only draws the values it is given with matplotlib and stays Python
(cpp/scripts/public_api_allowlist.tsv).
"""
import matplotlib
import matplotlib.pyplot as plt

import logging
matplotlib.use('Agg')
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicCompartmentalization as _impl


def parse_arguments(args=None):
    """The argparse parser of hicCompartmentalization, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicCompartmentalization')


def get_indices(obs_exp, row):
    return _impl.get_indices(obs_exp, row)


def count_interactions(obs_exp, pc1, quantiles_number, offset):
    "Counts the total interaction on obs_exp matrix per quantile and "
    return _impl.count_interactions(obs_exp, pc1, quantiles_number, offset)


def within_vs_between_compartments(normalised_sum_per_quantile,
                                   quantiles_number):
    """
    This function computes the interaction between two compartments and whithin
    each compartment. It adds one quantile at the time to compute these values
    on a larger block of the matrix (normalised_sum_per_quantile)until the
    maximum number of quantiles reached.
    """
    return _impl.within_vs_between_compartments(normalised_sum_per_quantile, quantiles_number)


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
