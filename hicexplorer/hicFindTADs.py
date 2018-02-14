from __future__ import division
import argparse
import time
import scipy.stats
import scipy.sparse
from hicexplorer import HiCMatrix
import hicexplorer.parserCommon
from hicexplorer._version import __version__
import numpy as np
import logging
log = logging.getLogger(__name__)
from .utilities import convertNansToZeros
from .utilities import remove_outliers


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        conflict_handler='resolve',
        description="""
Uses a measure called TAD-separation score to identify the degree of separation between
the left and right regions at each Hi-C matrix bin. This is done for a
running window of different sizes. Then, TADs are called as those
positions having a local TAD-separation score minimum. The TAD-separation score is
measured using the z-score of the Hi-C matrix and is defined as the mean zscore of all
the matrix contacts between the left and right regions (diamond).

To find the TADs, the program  needs to compute first the
TAD scores at different window sizes. Then, the results of that computation
are used to call the TADs. This is convenient to test different filtering criteria quickly
as the demanding step is the computation of TAD-separation scores.

 A simple example usage is:

$ hicFindTads -m hic_matrix.h5 --outPrefix TADs --correctForMultipleTesting fdr

The bedgraph file produced by this tool can be used to plot the so-called insulation score
along the genome or at specific regions. This score is much more reliable across samples
than the number of TADs or the TADs width that can vary depending on the sequencing depth
because of the lack of information at certain bins, and depending on the parameters used
with this tool.""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Corrected Hi-C matrix to use for the computations',
                                required=True)

    parserRequired.add_argument('--outPrefix',
                                help='File prefix to save the resulting files: 1. <prefix>_tad_separation.bm '
                                'The format of the output file is chrom start end TAD-sep1 TAD-sep2 TAD-sep3 .. etc. '
                                'We call this format a bedgraph matrix and can be plotted using '
                                '`hicPlotTADs`. Each of the TAD-separation scores in the file corresponds to '
                                'a different window length starting from --minDepth to --maxDepth. '
                                '2. <prefix>_zscore_matrix.h5, the zscore matrix used for the computation of '
                                'the TAD-separation score.  3. < prefix > _boundaries.bed, which'
                                'contains the positions of boundaries. The genomic coordinates in this file '
                                'correspond to the resolution used. Thus, for Hi-C bins of '
                                '10.000bp the boundary position is 10.000bp long. For restriction fragment '
                                'matrices the boundary position varies depending on the fragment length '
                                'at the boundary. 4. <prefix>_domains.bed '
                                'contains the TADs positions. This is a non-overlapping set of genomic '
                                'positions. 5. <prefix>_boundaries.gff Similar to the boundaries bed file '
                                'but with extra information (pvalue, delta). 6. <prefix>_score.bedgraph file '
                                'contains the TAD-separation score '
                                'measured at each Hi-C bin coordinate. Is useful to visualize in a genome '
                                'browser. The delta and pvalue settings are saved as part of the name.',
                                required=True)

    parserRequired.add_argument('--correctForMultipleTesting',
                                help='Select the bonferroni or false discovery rate for a multiple comparison. Bonferroni '
                                'controlls the familywise error rate (FWER) and needs a p-value. The false discovery rate '
                                '(FDR) controls the likelyhood of type I errors and needs a q-value. As a third option '
                                'it is possible to not use a multiple comparison method at all.',
                                type=str,
                                default="fdr",
                                choices=['fdr', 'bonferroni', 'None'],
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--minDepth',
                           help='Minimum window length (in bp) to be considered to the left and to the right '
                           'of each Hi-C bin. This number should be at least 3 times '
                           'as large as the bin size of the Hi-C matrix.',
                           metavar='INT bp',
                           type=int)

    parserOpt.add_argument('--maxDepth',
                           help='Maximum window length to be considered to the left and to the right '
                           'of the cut point in bp. This number should around 6-10 times '
                           'as large as the bin size of the Hi-C matrix.',
                           metavar='INT bp',
                           type=int)

    parserOpt.add_argument('--step',
                           help='Step size when moving from --minDepth to --maxDepth. Note, the step size'
                           'grows exponentially as '
                           '`maxDeph + (step * int(x)**1.5) for x in [0, 1, ...]` until  it '
                           'reaches `maxDepth`. For example, selecting  step=10,000, minDepth=20,000 '
                           'and maxDepth=150,000 will compute TAD-scores for window sizes: '
                           '20,000, 30,000, 40,000, 70,000 and 100,000',
                           metavar='INT bp',
                           type=int)

    parserOpt.add_argument('--TAD_sep_score_prefix',
                           help='Sometimes it is useful to change some of the parameters without recomputing the '
                           'z-score matrix and the TAD-separation score. For this case, the prefix containing the '
                           'TAD separation score and the z-score matrix can be given. If this option is given, '
                           'new boundaries will be computed but the values of --minDepth, --maxDepth and --step will '
                           'not be used.',
                           required=False)

    parserOpt.add_argument('--thresholdComparisons',
                           help='P-value threshold for the bonferroni correction / q-value for FDR. '
                           'The probability of a local minima to be a boundary '
                           'is estimated by comparing the distribution (Wilcoxon ranksum) of '
                           'the  zscores between the left and right '
                           'regions (diamond) at the local minimum with the matrix zscores for a '
                           'diamond at --minDepth to the left and a diamond --minDepth to the right. '
                           'If --correctForMultipleTesting is \'None\' the threshold is applied on the '
                           'raw p-values without any multiple testing correction. Set it to \'1\' if no threshold should be used.',
                           type=float,
                           default=0.01)

    parserOpt.add_argument('--delta',
                           help='Minimum threshold of the difference between the TAD-separation score of a '
                           'putative boundary and the mean of the TAD-sep. score of surrounding bins. '
                           'The delta value reduces spurious boundaries that are shallow, which usually '
                           'occur at the center of large TADs when the TAD-sep. score is flat. Higher '
                           'delta threshold values produce more conservative boundary estimations. By '
                           'default a value of 0.01 is used.',
                           type=float,
                           default=0.01)

    parserOpt.add_argument('--minBoundaryDistance',
                           help='Minimum distance between boundaries (in bp). This parameter can be '
                           'used to reduce spurious boundaries caused by noise.',
                           type=int)

    parserOpt.add_argument('--numberOfProcessors', '-p',
                           help='Number of processors to use ',
                           type=int,
                           default=1)

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit.')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser
    

def nbinom_est_dist(size, prob, triu_ma, hic_matrix):
    # compute a mapping from mean to distance
    mean2dist = {'mean': [], 'dist': []}
    for dist in np.sort(list(size)):
        mean = scipy.stats.nbinom.mean(size[dist], prob[dist])
        if not np.isnan(mean):
            mean2dist['mean'].append(mean)
            mean2dist['dist'].append(dist)
    mean2dist['mean'] = np.array(mean2dist['mean'])
    mean2dist['dist'] = np.array(mean2dist['dist'])

    # the values have to be computed for all the
    # matrix excepting inter chromosome
    row, col = np.triu_indices(triu_ma.shape[0])
    dist_list, chrom_list = hic_matrix.getDistList(row, col,
                                                   hic_matrix.cut_intervals)
    triu_ma = triu_ma.tolil()
    transf_ma = np.zeros(len(dist_list))
    for idx, orig_dist in enumerate(dist_list):
        if orig_dist == -1:
            continue
        data = triu_ma[row[idx], col[idx]]
        try:
            size[orig_dist]
        except KeyError:
            continue
        if _nbinomPvalue(data,
                         size[orig_dist],
                         prob[orig_dist]) < 5:
            # pass
            continue

        # get largest closest mean
        if data > mean2dist['mean'].max() or \
                data < mean2dist['mean'].min():
            dist = orig_dist
        else:
            try:
                mean_idx = np.flatnonzero(data < mean2dist['mean'])[-1]
                dist = mean2dist['dist'][mean_idx]
            except IndexError:
                dist = orig_dist

        # min distance should be 1, otherwise the
        # sparse matrix will treat nans as cero distance
        transf_ma[idx] = dist_list[idx] - dist
        transf_ma[idx] = dist

    # set the new values back into the original matrix
    triu_ma = scipy.sparse.coo_matrix((transf_ma, (row, col)), shape=triu_ma.shape)
    # fill the lower triangle
    triu_ma = triu_ma + scipy.sparse.triu(triu_ma, 1).T
    triu_ma = triu_ma.tocsr()
    triu_ma.eliminate_zeros()

    return triu_ma


def transformMatrix(hicma, method, per_chr=False, original_matrix=None, depth_in_bins=None):
    methods_avail = {'residuals': _residuals, 'obs/exp': _obsExp,
                     'z-score': _zscore, 't-score': _tscore,
                     'p-value': _pvalue, 'nbinom-p-value': _nbinomPvalue,
                     'nbinom-expected': _nbinomExpected,
                     'log-norm': _lognormPvalue,
                     'chi-squared': _chi2Pvalue}

    counts_by_distance, cut_intervals = hicma.getCountsByDistance(per_chr=per_chr)
    if method in ['nbinom-p-value', 'nbinom-expected', 'nbinom-est-dist']:
        size, prob = fitNegBinom_Rserve(counts_by_distance, per_chr=per_chr,
                                        plot_distribution=True)
    elif method == 'log-norm':
        mu_, sigma = fitDistribution(counts_by_distance, 'lognorm')
    else:
        if per_chr:
            mu_ = {}
            sigma = {}
            n_value = {}
            for chrom in list(counts_by_distance):
                mu_[chrom] = dict([(x, np.mean(counts_by_distance[chrom][x]))
                                   for x in counts_by_distance[chrom]])
                sigma[chrom] = dict([(x, np.std(counts_by_distance[chrom][x]))
                                     for x in counts_by_distance[chrom]])
                n_value[chrom] = dict([(x, len(counts_by_distance[chrom][x]))
                                       for x in counts_by_distance[chrom]])
        else:
            mu_ = dict([(x, np.mean(counts_by_distance[x]))
                        for x in counts_by_distance])
            sigma = dict([(x, np.std(counts_by_distance[x]))
                          for x in counts_by_distance])
            n_value = dict([(x, len(counts_by_distance[x]))
                            for x in counts_by_distance])

    # use only the upper half of the matrix
    triu_ma = scipy.sparse.triu(hicma.matrix, format='coo')
    if original_matrix:
        orig_ma = original_matrix.matrix
        if per_chr:
            noise_level = {}
            for chrom in list(counts_by_distance):
                chr_range = original_matrix.getChrBinRange(chrom)
                chr_submatrix = orig_ma[chr_range[0]:chr_range[1],
                                        chr_range[0]:chr_range[1]]
                noise_level[chrom] = np.median(chr_submatrix.data)
        else:
            noise_level = np.median(orig_ma.data)

        log.debug('noise error set to {}\n'.format(noise_level))
    else:
        noise_level = None
    log.debug("finish computing fitting parameters")

    ########################
    # after the distributions are fitted
    # now the matrix values are evaluated
    if method == 'nbinom-est-dist':
        triu_ma = nbinom_est_dist(size, prob, triu_ma, hicma)

    else:
        under_noise = 0
        dist_list, chrom_list = hicma.getDistList(triu_ma.row, triu_ma.col, cut_intervals)

        assert len(dist_list) == len(triu_ma.data), "lists not of equal size"
        susprow_list = []
        suspcol_list = []
        transf_ma = np.zeros(len(triu_ma.data))
        start_time = time.time()
        # transform each value  of the data matrix to p-value, obs/exp, correlation etc.
        log.debug("computing transform values\n")
        for idx, data in enumerate(triu_ma.data):
            # skip if original value is less than noise level
            if noise_level:
                if per_chr:
                    if dist_list[idx] == -1:
                        continue
                    elif (orig_ma[triu_ma.row[idx],
                                  triu_ma.col[idx]] <=
                          noise_level[chrom_list[idx]]):
                        under_noise += 1
                        continue
                elif orig_ma[triu_ma.row[idx],
                             triu_ma.col[idx]] <= noise_level:
                    under_noise += 1
                    continue

            if method in ['nbinom-p-value', 'nbinom-expected']:
                if dist_list[idx] == -1:
                    continue
                if per_chr:
                    transf_ma[idx] = methods_avail[method](
                        data,
                        size[chrom_list[idx]][dist_list[idx]],
                        prob[chrom_list[idx]][dist_list[idx]])
                else:
                    transf_ma[idx] = methods_avail[method](data,
                                                           size[dist_list[idx]],
                                                           prob[dist_list[idx]])
                if data > 3 * orig_ma[triu_ma.row[idx], triu_ma.col[idx]]:
                    log.debug("skipping p-value {} for "
                              "value {} at {}, norm-value {}\n".
                              format(transf_ma[idx], chrom_list[idx],
                                     orig_ma[triu_ma.row[idx],
                                             triu_ma.col[idx]], data))
                    continue
                if transf_ma[idx] > 4.5 and \
                        data > 2 * orig_ma[triu_ma.row[idx], triu_ma.col[idx]]:
                    susprow = triu_ma.row[idx]
                    suspcol = triu_ma.col[idx]
                    log.debug("suspicious p-value {} for "
                              "value {} at {}, norm-value {}\n".
                              format(transf_ma[idx], chrom_list[idx],
                                     orig_ma[susprow, suspcol], data))
                    susprow_list.append(susprow)
                    suspcol_list.append(suspcol)

            if method in ['obs/exp', 'residuals']:
                if per_chr:
                    if dist_list[idx] == -1:
                        continue
                    fit_mu = mu_[chrom_list[idx]]
                else:
                    fit_mu = mu_

                transf_ma[idx] = methods_avail[method](data,
                                                       fit_mu[dist_list[idx]])

            else:
                if per_chr:
                    if dist_list[idx] == -1:
                        continue
                    fit_mu = mu_[chrom_list[idx]]
                    fit_sigma = sigma[chrom_list[idx]]
                    fit_n = n_value[chrom_list[idx]]
                else:
                    fit_mu = mu_
                    fit_sigma = sigma
                    fit_n = n_value

                transf_ma[idx] = methods_avail[method](data,
                                                       fit_mu[dist_list[idx]],
                                                       fit_sigma[dist_list[idx]],
                                                       fit_n[dist_list[idx]])

            if idx > 0 and (idx == 10000 or idx % 500000 == 0):
                endtime = time.time()
                estimated = (float(len(transf_ma) - idx) *
                             (endtime - start_time)) / idx
                mmin, sec = divmod(estimated, 60)
                hour, mmin = divmod(mmin, 60)
                log.debug("iteration: {} Estimated remaining time "
                          "{:.0f}:{:.0f}:{:.0f}".format(idx, hour, mmin, sec))

        """
        print "problematic bins:"
        for uniq in np.concatenate([susprow_list, suspcol_list]):
            print hicma.cut_intervals[uniq]
        """
        # set the new values back into the original matrix
        triu_ma.data = transf_ma
        # fill the lower triangle
        triu_ma = triu_ma + scipy.sparse.triu(triu_ma, 1).T
        triu_ma = triu_ma.tocsr()
        triu_ma.eliminate_zeros()

    return triu_ma


def mylog(data):
    return np.log(data + 1)


def getZscores(hicma):
    return transformMatrix(hicma, 'z-score')


def getTscores(hicma):
    return transformMatrix(hicma, 't-score')


def getResiduals(hicma):
    return transformMatrix(hicma, 'residuals')


def getObsExp(hicma):
    return transformMatrix(hicma, 'obs/exp')


def getPearson(matrix):
    matrix = convertNansToZeros(matrix).todense()
    from scipy.stats import pearsonr
    from scipy.sparse import csr_matrix
    numRows, numCols = matrix.shape
    # create matrix to hold computed pval
    pMa = np.zeros(shape=(numCols, numRows))
    pMa[:, :] = 0
    for row in range(numRows):
        if row % 10 == 0:
            log.debug("{} rows processed ({:.2f})\n".format(row, float(row) / numRows))
        for col in range(numCols):
            if not np.isnan(pMa[col, row]):
                pMa[row, col] = pMa[col, row]
                continue
            try:
                # pearsonr returns two values, the first is the
                # correlation, the second is a pvalue.
                pMa[row, col] = pearsonr(np.asarray(matrix[row, :])[0], np.asarray(matrix[:, col].T)[0])[0]
            except Exception:
                continue

    return convertNansToZeros(csr_matrix(pMa)).todense()


def applyFdr(matrix):
    """
    compute false discovery rate
    but only for half the symmetric matrix

    :param: matrix a sparse matrix
    """
    # use only the upper half of the matrix
    mat = scipy.sparse.triu(matrix)
    # exp is required because the matrix contains
    # -log(pvalues)
    mat.data = -np.log(_fdr(np.exp(-mat.data)))
    mat = mat + scipy.sparse.triu(mat, 1).T

    return mat


def _fdr(pvalues):
    """
    Code translated from the R function
    p.adjust(pvalues, method='fdr')

    :param: list of pvalues
    """
    len_pvals = len(pvalues)
    if len_pvals <= 1:
        return pvalues
    i = np.arange(1, len_pvals + 1)[::-1]
    order = np.argsort(pvalues)[::-1]
    rorder = np.argsort(order)

    # clip values bigger than 1 or smaller than 0
    return np.clip(float(len_pvals) / i * pvalues[order], 0, 1)[rorder]


def fitNegBinom_Rserve(countsByDistance, plot_distribution=False,
                       per_chr=False):
    """
    Fits a negative binomial distribution to the counts
    found at each different distance.

    The fitting is attempted first using a python method, and
    if this fails R is used through Rserve.

    For the fitting, the outliers are removed. Outliers are
    defined as those having a z-score higher than 3.4. This
    number as defined after exploring different values
    of z-scores and estimating the best goodness of fit.
    The Hi-C data is expected to contain outliers but they
    are problematic to fit and test the goodness of fit of a
    distribution, that's why there are removed.
    """

    # if the counts are per chromosome,
    # use the function recursively
    if per_chr:
        size = {}
        prob = {}
        for chrom in list(countsByDistance):
            log.info('computing negative binomial for '
                     '{}\n'.format(chrom))
            size[chrom], prob[chrom] = \
                fitNegBinom_Rserve(countsByDistance[chrom],
                                   plot_distribution=plot_distribution)
        return size, prob

    import pyRserve
    import matplotlib.pyplot as plt
    try:
        conn = pyRserve.connect()
        conn.r('library("MASS")')

    except Exception:
        log.exception("Could not connect to Rserve. Check that Rserve is up and running")
        exit(1)
    size = {}
    mu = {}
    prob = {}
    pval = {}
    good = 0
    bad = 0

    for dist in np.sort(list(countsByDistance)):
        if dist == -1:  # skip intra chromosomal counts
            continue
        size[dist] = np.nan
        mu[dist] = np.nan
        prob[dist] = np.nan
        if sum(countsByDistance[dist]) == 0.0:
            log.debug("no counts for bins at distance {}".format(dist))
            continue
        if np.any(np.isnan(countsByDistance[dist])) is True:
            log.error("matrix contains NaN values\n")

        counts = remove_outliers(countsByDistance[dist])
        if len(counts) <= 20:
            continue
        # the values in countsByDistance of a corrected matrix
        # are float values, but integers are needed for
        # the negative binomial fitting in R.
        counts_int = np.round(counts).astype('int')

        # try first using the python fit for the
        # negative binomial
        try:
            size[dist], prob[dist] = fit_nbinom(counts)
        except ValueError:
            # try with R..
            try:
                res = conn.r.fitdistr(counts_int, 'negative binomial')
            except Exception:
                continue
            size[dist] = res[0]['size']
            mu[dist] = res[0]['mu']

            if np.isnan(size[dist]) or np.isnan(mu[dist]):
                log.debug("for dist={}, size={}, mu={}, "
                          "len={}\n".format(dist, size[dist],
                                            mu[dist], len(counts)))
                continue

            # The output from 'fitdistr' are size and mu.
            # but the scipy function that is based on the negative binomial
            # needs size and probability as parameters. However,
            # prob = size / ( size + mu )
            prob[dist] = size[dist] / (size[dist] + mu[dist])

        log.info(".")  # print a . to show progress

        # evaluate fit of the counts distribution with respect to
        # the negative binomial  distribution using the parameters
        # returned by R
        fitted_dist = scipy.stats.nbinom.rvs(size[dist],
                                             prob[dist],
                                             size=len(counts) * 2)
        pval[dist] = scipy.stats.ks_2samp(counts_int,
                                          fitted_dist)[1]

#        pval[dist] = scipy.stats.wilcoxon(counts, fitted_dist)[1]
        if pval[dist] < 0.01:
            bad += 1
            log.debug(
                "\nThe fit p-value {} for {} is too low to consider "
                "the distribution negative binomial".format(
                    pval[dist], dist))
        else:
            good += 1

        if (plot_distribution and
                dist in [50000] + range(0, max(list(countsByDistance)), 1000000)):
            # actual and fitted distributions are plotted
            # next to each other

            diff = counts.max() - counts.min()
            if diff >= 1000:
                nbins = 50
            elif 1000 > diff >= 100:
                nbins = 30
            elif 100 > diff >= 50:
                nbins = diff / 2
            else:
                nbins = (counts.max() - counts.min())
            freq, bins = np.histogram(counts.astype(int), nbins,
                                      normed=True)
            plt.hist(counts, bins, linewidth=0.1, alpha=0.8, normed=True)
            # plt.hist(fitted_dist, bins, histtype='step', linestyle='solid',
            #          linewidth=1.5, color='black', normed=True)
            pdf_fitted = scipy.stats.nbinom.pmf(bins.astype('int'),
                                                size[dist],
                                                prob[dist])
            plt.plot(bins.astype(int), pdf_fitted,
                     label='fitted nbinom gf={:.3f}'.format(pval[dist]))

            fig_name = '/tmp/fitt_{}_{}.png'.format('nbinom', dist)
            plt.title('{} bp; size: {}, prob: {}'.format(dist,
                                                         size[dist],
                                                         prob[dist]))
            plt.ylim(0, np.max(freq) + np.max(freq) * 0.2)
            plt.legend()
            plt.savefig(fig_name, dpi=200)
            plt.close()
            log.debug("check {}".format(fig_name))

    log.debug("good {}, bad {}\n".format(good, bad))

    return size, prob


def fitNegBinom(countsByDistance):
    """
    Replaced by fitNegBinom_Rserve
    Returns a tuple of two dictionaries,
    one contains the size fit and the
    other the probability
    """
    log.debug("fit neg binom")
    size = {}
    mu = {}
    prob = {}
    pval = {}
    for dist in countsByDistance:
        size[dist] = np.nan
        mu[dist] = np.nan
        prob[dist] = np.nan
        if sum(countsByDistance[dist]) == 0.0:
            log.debug("no counts for bins at distance {}".format(dist))
            continue
        if len(countsByDistance[dist]) <= 2:
            continue
        # the values in countsByDistance of a corrected matrix
        # are float values, but integers are needed for
        # the negative binomial.
        counts = np.rint(countsByDistance[dist]).astype('int')
        try:
            size[dist], prob[dist] = fit_nbinom(counts)
        except ValueError as error:
            log.debug("could not compute pval for dist={}. "
                      "Message:\n {}".format(dist, error))

        if np.isnan(size[dist]) or np.isnan(prob[dist]):
            log.debug("for dist={}, size={}, prob={}, len={}".format(dist,
                                                                     size[dist],
                                                                     prob[dist],
                                                                     len(counts)))
        else:
            # evaluate fit of the counts distribution with respect to the
            # negative binomial distribution using the parameters returned by R
            pval[dist] = scipy.stats.wilcoxon(
                counts,
                scipy.stats.nbinom.rvs(size[dist],
                                       prob[dist],
                                       size=len(counts)))[1]
            if pval[dist] < 0.001:
                log.debug("problem with {} when fitting a negative binomial. "
                          "The fit p-value ({}) is too low to consider the "
                          "distribution negative binomial".format(dist, pval[dist]))

    return (size, prob)


def fitDistribution(countsByDistance, distribution, plot_distribution=False):
    """
    Generic method to fit continuous
    distributions to the  Hi-C countsByDistance
    The distribution names are the ones supported
    by scipy.
    """
    mu = {}
    sigma = {}
    pval = {}
    good = 0
    bad = 0
    good_nb = 0
    bad_nb = 0
    import pyRserve
    try:
        conn = pyRserve.connect()
        conn.r('library("MASS")')
    except Exception:
        log.exception("Could not connect to Rserve. Check that Rserve is up and running")
        exit(1)

    for distnc in np.sort(list(countsByDistance)):
        if distnc == -1:  # skip intra chromosomal counts
            continue
        if sum(countsByDistance[distnc]) == 0.0:
            log.debug("no counts for bins at distance {}".format(distnc))
            continue
        if len(countsByDistance[distnc]) <= 2:
            continue
        log.info('.')
        # TEMP code to compare with negative binomial ###

        # the values in countsByDistance of a corrected matrix
        # are float values, but integers are needed for
        # the negative binomial.

        counts_nb = remove_outliers(np.round(countsByDistance[distnc]).astype('int'))

        # try first using the python fit for the
        # negative binomial
        try:
            size, prob = fit_nbinom(remove_outliers(countsByDistance[distnc]))
        except ValueError:
            # try with R..
            res = conn.r.fitdistr(counts_nb, 'negative binomial')
            size = res[0]['size']
            mu_ = res[0]['mu']

            if np.isnan(size) or np.isnan(mu_):
                log.debug("for dist={}, size={}, mu={}, len={}".format(
                    distnc, size, mu_, len(counts_nb)))
                continue

            prob = size / (size + mu_)
        nbin = scipy.stats.nbinom(size, prob)
        #####

        counts = remove_outliers(countsByDistance[distnc])
        counts[counts == 0] = 0.01
        dist = getattr(scipy.stats, distribution)
        param = dist.fit(counts, floc=0)
        if np.any(np.isnan(param)):
            log.debug('\n{} no params computed'.format(distnc))
            import ipdb
            ipdb.set_trace()
        mu[distnc] = param[-1]
        sigma[distnc] = param[0]

        # estimate the goodness of fit pvalue
        fitted_dist = dist.rvs(*param[:-2],
                               loc=param[-2], scale=param[-1],
                               size=len(counts) * 2)
        pval[distnc] = scipy.stats.ks_2samp(counts,
                                            fitted_dist)[1]
        fitted_dist_nb = scipy.stats.nbinom.rvs(size,
                                                prob,
                                                size=len(counts_nb) * 2)

        pval_nb = scipy.stats.ks_2samp(counts_nb,
                                       fitted_dist_nb)[1]
        if pval[distnc] < 0.01:
            bad += 1
        else:
            good += 1

        if pval_nb < 0.01:
            bad_nb += 1
        else:
            good_nb += 1
        if pval[distnc] < 0.01:
            log.warning("\nproblem with {}, p-value for "
                        "{} fit: {} (NB fit: {})".format(distnc,
                                                         distribution,
                                                         pval[distnc],
                                                         pval_nb))

        if (plot_distribution and
                distnc in range(50000, max(countsByDistance.keys()), 500000)):

            import matplotlib.pyplot as plt
            freq, bins = np.histogram(counts, 30,
                                      normed=True)
            plt.close()  # to avoid overlaps
            plt.hist(counts, bins, linewidth=0.1, alpha=0.8, normed=True)
#            plt.hist(fitted_dist, bins, histtype='step', linestyle='solid',
#                      linewidth=1.5, color='black', normed=True)
#            plt.hist(fitted_dist_nb, bins, histtype='step', linestyle='solid',
#                      linewidth=1.5, color='grey', normed=True)
            pdf_fitted = dist.pdf(bins,
                                  *param[:-2],
                                  loc=param[-2],
                                  scale=param[-1])
            ##
            plt.plot(bins.astype(int), nbin.pmf(bins.astype('int')),
                     label='NB {:.2f}'.format(pval_nb))
            ##
            plt.plot(bins, pdf_fitted,
                     label='{} {:.2f}'.format(distribution, pval[distnc]))
            fig_name = '/tmp/fitt_{}_{}.png'.format(distribution, distnc)
            plt.title('{} bp'.format(distnc))
            plt.ylim(0, np.max(freq) + np.max(freq) * 0.2)
            plt.legend()
            plt.savefig(fig_name, dpi=200)
            plt.close()
            log.debug("check {}".format(fig_name))
    log.debug("good {}, bad {}, good_nb {}, bad_nb {}".format(good, bad, good_nb,
                                                              bad_nb))
    return (mu, sigma)


def fitChisquared(countsByDistance):
    shape = {}
    loc = {}
    scale = {}
    pval = {}
    for x in countsByDistance:
        if len(countsByDistance[x]) > 2:
            counts = countsByDistance[x]
            if len(counts) > 10000:
                rand = np.arange(10000)
                np.random.shuffle(rand)
                counts = counts[rand]
            shape[x], loc[x], scale[x] = scipy.stats.chi2.fit(counts, 3, floc=0)
            # estimate the fit pvalue
            pval[x] = scipy.stats.wilcoxon(counts,
                                           scipy.stats.chi2.rvs(shape[x],
                                                                loc=loc[x],
                                                                scale=scale[x],
                                                                size=len(counts)))[1]
            if pval[x] < 0.001:
                log.warning("problem with {}, p-value for log-norm fit: {}".format(x, pval[x]))

    return (shape, loc, scale)


def _residuals(value, mu):
    return value - mu


def _obsExp(value, mu):
    return value / mu


def _zscore(value, mu, sigma, n):
    return (value - mu) / sigma


def _tscore(value, mu, sigma, n):
    return (value - mu) / (sigma / np.sqrt(n))


def _pvalue(value, mu, sigma, n):
    #    tscore = _tscore(value, mu, sigma, n)
    pvalue = -np.log(scipy.stats.t.sf(value, n, loc=mu, scale=sigma))
    return pvalue


def _lognormPvalue(value, mu, sigma, n):
    """
    n is not used, but is generic
    to other calls (see _pvalue).
    """
    return scipy.stats.norm.sf(np.log(value), mu, sigma)


def _chi2Pvalue(value, shape, loc, scale):
    """
    n is not used, but is generic
    to other calls (see _pvalue).
    """
    return -np.log(scipy.stats.chi2.sf(np.log(value), shape, loc=loc, scale=scale))


def _nbinomPvalue(value, size, prob):
    pvalue = -np.log(scipy.stats.nbinom.sf(value, size, prob))
    return pvalue


def _nbinomExpected(value, size, prob):
    mean = scipy.stats.nbinom.mean(size, prob)
    return mean


def fit_nbinom(k):
    """
    from https://gist.github.com/gjx/5987413
    Note: for some cases this function does
    not work. brentq fails with
    error: f(a) and f(b) must have different signs
    Rserve may be more accurate
    """
    from scipy.special import psi
    from scipy.optimize import brentq
    N = len(k)
    n = brentq(lambda r: sum(psi(k + r)) - N * psi(r) +
               N * np.log(r / (r + sum(k / N))),
               np.finfo(np.float128).eps,
               np.max(k))
    p = n / (n + sum(k / N))  # Note: this `p` = 1 - `p` from Wikipedia
    return n, p


def main(args=None):
    log.debug(args)
    args = parse_arguments().parse_args(args)

    hic_ma = HiCMatrix.hiCMatrix(args.matrix)
    try:
        hic_ma.maskBins(hic_ma.nan_bins)
    except AttributeError:
        pass

    if args.skipDiagonal:
        hic_ma.diagflat()

    if args.method == 'obs/exp':
        hic_ma.convert_to_obs_exp_matrix(maxdepth=args.depth, perchr=args.perchr)
    else:
        hic_ma.convert_to_zscore_matrix(maxdepth=args.depth, perchr=args.perchr)

    hic_ma.save(args.outFileName)
