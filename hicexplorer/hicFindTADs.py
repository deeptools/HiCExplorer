"""hicexplorer.hicFindTADs, implemented in C++.

main() runs the C++ executable hicFindTADs and parse_arguments() is built
from its --help-json (hicexplorer._cpp). Every other function and every
method of HicFindTads keeps its 3.7.6 signature, docstring and behaviour, and
its body calls hicexplorer._hicx.hicFindTADs
(cpp/python/src/api/hicFindTADs.cpp, cpp/core/include/hicx/find_tads.hpp).
The deviations from 3.7.6 are listed in
hicexplorer/test/api/cases_hicFindTADs.py.

The hicFindTADs tool also filters boundaries by CTCF motif orientation
(--ctcfFile, a HiCExplorer 4 extension, with the divergent rule: a reverse
'-' site at or before the boundary, a forward '+' site at or after it). No
name of this module has a CTCF filter, so the API and 3.7.6 agree on it.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import logging
from hicexplorer.utilities import enlarge_bins  # noqa: F401 (the C++ hicFindTADs looks it up in this module)
from hicexplorer.utilities import toString, toBytes  # noqa: F401 (the C++ hicFindTADs looks them up in this module)

log = logging.getLogger(__name__)

from hicexplorer._hicx import hicFindTADs as _impl

# this is a holder vor the
hic_ma = None


def parse_arguments(args=None):
    """The argparse parser of hicFindTADs, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicFindTADs')


def compute_matrix_wrapper(args):
    return _impl.compute_matrix_wrapper(args)


def get_cut_weight_by_bin_id(matrix, cut, depth, return_mean=False):
    """
    like get_cut_weight which is the 'diamond' representing the counts
    from a region -dept to +depth from the given bin position (cut):
    """
    return _impl.get_cut_weight_by_bin_id(matrix, cut, depth, return_mean)


def get_idx_of_bins_at_given_distance(hic_matrix, idx, window_len):
    """
    Returns the right and left indices that are at a given distance in bp
    from `idx`
    Parameters
    ----------
    idx reference index
    window_len distance in bp
    Returns
    -------
    tuple, with left and right bin indices
    """
    return _impl.get_idx_of_bins_at_given_distance(hic_matrix, idx, window_len)


def get_cut_weight(hic_matrix, cut, window_len, return_mean=False):
    """
    'diamond' (when matrix is rotated 45 degrees) representing the counts
    from a region -window_len to +window_len from the given bin position (cut):
    """
    return _impl.get_cut_weight(hic_matrix, cut, window_len, return_mean)


def get_triangle(hic_matrix, cut, window_len, return_mean=False):
    """
    like get_cut_weight which is the 'diamond' representing the counts
    from a region -window_len to +window_len from the given bin position (cut):
    """
    return _impl.get_triangle(hic_matrix, cut, window_len, return_mean)


def get_incremental_step_size(min_win_size, max_win_size, start_step_len):
    """
    generates a list of incremental windows sizes (measured in bins)
    :param min_win_size: starting window size
    :param max_win_size: end window size
    :param start_step_len: start length
    :return: incremental_step list of bin lengths
    """
    return _impl.get_incremental_step_size(min_win_size, max_win_size, start_step_len)


def compute_matrix(bins_list, min_win_size=8, max_win_size=50, step_len=2):
    """
    Receives a number of bins for which the tad-score should be computed.
    Parameters
    ----------
    hic_ma Hi-C matrix object
    bins_list list of bins to process
    min_win_size
    max_win_size
    step_len
    Returns
    -------
    """
    return _impl.compute_matrix(bins_list, min_win_size, max_win_size, step_len)


class HicFindTads(object):
    def __init__(self, matrix, num_processors=1, max_depth=None, min_depth=None, step=None, delta=0.01,
                 min_boundary_distance=None, use_zscore=True, p_correct_for_multiple_testing="fdr", p_threshold_comparisons=0.01,
                 pChromosomes=None):
        """
        Parameters
        ----------
        matrix  Either a filename or a Hi-C Matrix object
        num_processors
        max_depth max window distance to consider (total window length is 2* max depth)
        min_depth min window to consider (total window length is 2* max min)
        step progression step from min_depth to max depth. The value give is for the first step, wich iteratively grows
                exponentially (exponent = 1.5) until it reaches the max_depth value.
        delta
        min_boundary_distance
        use_zscore boolean. By default is true. Set to other option
        pCorrectForMultipleTesting Multiple comparisons method: FDR, Bonferroni or None
        pThresholdComparisons The threshold for the Multiple comparisons. It is used as p-value for Bonferroni or as q-value for FDR.
        pChromosomes The chromomes that should be included for the analysis.
        """
        return _impl.HicFindTads___init__(self, matrix, num_processors, max_depth, min_depth, step, delta, min_boundary_distance, use_zscore, p_correct_for_multiple_testing, p_threshold_comparisons, pChromosomes)

    def set_matrix(self, pMatrix, pChromosomes):
        return _impl.HicFindTads_set_matrix(self, pMatrix, pChromosomes)

    def set_variables(self):
        """
        Checks the value of the max_depth, min_depth and step variables, setting default parameters
        if neccesary.
        Returns
        -------
        None
        """
        return _impl.HicFindTads_set_variables(self)

    @staticmethod
    def peakdetect(y_axis, x_axis=None, lookahead=3, delta=0, chrom=None):
        """
        Based on the MATLAB script at:
        http://billauer.co.il/peakdet.html
        function for detecting local maximum and minimum in a signal.
        Discovers peaks by searching for values which are surrounded by lower
        or larger values for maximum and minimum respectively
        keyword arguments:
        :param: y_axis -- A list containing the signal over which to find peaks
        :param: x_axis -- (optional) A x-axis whose values correspond to the y_axis list
            and is used in the return to specify the position of the peaks. If
            omitted an index of the y_axis is used. (default: None)
        :param: lookahead -- (optional) distance to look ahead from a peak candidate to
            determine if it is the actual peak
        :param: delta -- (optional) this specifies a minimum difference between a peak and
            the following points, before a peak may be considered a peak. Useful
            to hinder the function from picking up false peaks towards to end of
            the signal. To work well delta should be set to delta >= RMSnoise * 5.
        :return: -- two lists [max_peaks, min_peaks] containing the positive and
            negative peaks respectively. Each cell of the lists contains a tuple
            of: (position, peak_value)
            to get the average peak value do: np.mean(max_peaks, 0)[1] on the
            results to unpack one of the lists into x, y coordinates do:
            x, y = zip(*tab)
        """
        return _impl.HicFindTads_peakdetect(y_axis, x_axis, lookahead, delta, chrom)

    @staticmethod
    def delta_wrt_window(min_idx_list, matrix_avg, chrom, window_len=10):
        """
        Computes the local 'delta' for each minima identified. The local
        delta is defined as the difference between the TAD-separation score at the minimum with respect
        to each of the TAD-separation scores of the `window_len` to the left, plus the `window_len` to
        the right of each minimum. The minimum itself is excluded.
        Parameters
        ----------
        matrix TAD score matrix for different window sizes
        chrom list of chrom names for each matrix row
        window_len the length of the window to look for the local TAD-score average
        Returns
        -------
        A dict of each minima delta to the mean of the TAD-separation score surrounding the minima. The keys
        are the min_idx
        """
        return _impl.HicFindTads_delta_wrt_window(min_idx_list, matrix_avg, chrom, window_len)

    @staticmethod
    def find_consensus_minima(tad_score_matrix, lookahead=3, chrom=None):
        """
        Finds the minimum over the average values per column
        :param tad_score_matrix: TAD-separation score matrix
        :return:
        """
        return _impl.HicFindTads_find_consensus_minima(tad_score_matrix, lookahead, chrom)

    def hierarchical_clustering(self, boundary_list, clusters_cutoff=[]):
        """
        :param boundary_list: is a list of tuples each containing
        the location of a boundary. The order should be sorted
        and contain the following values:
            (chrom, start, value)
        :param clusters_cutoff: List of values to separate clusters. The
            clusters found at those value thresholds are returned.
        :return: z_value, clusters
        For z_value, the format used is similar as the scipy.cluster.hierarchy.linkage() function
        which is described as follows:
        A 4 by :math:`(n-1)` matrix ``z_value`` is returned. At the
        :math:`i`-th iteration, clusters with indices ``z_value[i, 0]`` and
        ``z_value[i, 1]`` are combined to form cluster :math:`n + i`. A
        cluster with an index less than :math:`n` corresponds to one of
        the :math:`n` original observations. The distance between
        clusters ``z_value[i, 0]`` and ``z_value[i, 1]`` is given by ``z_value[i, 2]``. The
        fourth value ``z_value[i, 3]`` represents the number of original
        observations in the newly formed cluster.
        The difference is that instead of a 4 times n-1 array, a
        6 times n-1 array is returned. Where positions 4, and 5
        correspond to the genomic coordinates of ``z_value[i, 0]`` and ``z_value[i, 1]``
        """
        return _impl.HicFindTads_hierarchical_clustering(self, boundary_list, clusters_cutoff)

    def save_linkage(Z, file_name):
        """
        :param Z: Z has a format similar to the scipy.cluster.linkage matrix (see function
                    hierarchical_clustering).
        :param file_name: File name to save the results
        :return: None
        """
        return _impl.HicFindTads_save_linkage(Z, file_name)

    def get_domains(boundary_list):
        """
        returns for each boundary a chrom, start,end position
        corresponding to each TAD
        :param boundary_position: list of boundary chromosomal positions
        :return: list of (chrom, start, end, value) tuples.
        """
        return _impl.HicFindTads_get_domains(boundary_list)

    def save_bedgraph_matrix(self, outfile):
        """
        Save matrix as chrom, start, end ,row, values separated by tab
        I call this a bedgraph matrix and the ending is .bm
        Returns
        -------
        None
        """
        return _impl.HicFindTads_save_bedgraph_matrix(self, outfile)

    def save_clusters(clusters, file_prefix):
        """
        :param clusters: is a dictionary whose key is the cut of used to create it.
                         the value is a list of tuples, each representing
                          a genomec interval as ('chr', start, end).
        :param file_prefix: file prefix to save the resulting bed files
        :return: list of file names created
        """
        return _impl.HicFindTads_save_clusters(clusters, file_prefix)

    def save_domains_and_boundaries(self, prefix):
        return _impl.HicFindTads_save_domains_and_boundaries(self, prefix)

    def compute_spectra_matrix(self, perchr=True):
        """
        Uses multiple processors to compute the TAD-score
        Returns
        -------
        bed matrix (chrom start end value1 value2 ... value_n)
        """
        return _impl.HicFindTads_compute_spectra_matrix(self, perchr)

    def load_bedgraph_matrix(self, filename, pChromosomes=None):
        return _impl.HicFindTads_load_bedgraph_matrix(self, filename, pChromosomes)

    def min_pvalue(self, min_idx):
        """
        For each putative local minima, find the -window_len diammond and the +window_len diamond
        and compare with the local minima using wilcoxon rank sum.
        Parameters
        ----------
        min_idx list of local minima
        Returns
        -------
        list of p-values per each local minima
        """
        return _impl.HicFindTads_min_pvalue(self, min_idx)

    def find_boundaries(self):
        return _impl.HicFindTads_find_boundaries(self)


def print_args(args):
    """
    Print to stderr the parameters used
    Parameters
    ----------
    args
    Returns
    -------
    """
    return _impl.print_args(args)


def main(args=None):
    """Run the C++ implementation of hicFindTADs with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicFindTADs')(args)
