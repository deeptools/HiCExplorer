import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import numpy as np
import time
import logging
log = logging.getLogger(__name__)


def iterativeCorrection(matrix, v=None, M=50, tolerance=1e-5, verbose=False):
    """
    adapted from cytonised version in mirnylab
    original code from: ultracorrectSymmetricWithVector
    https://bitbucket.org/mirnylab/mirnylib/src/924bfdf5ed344df32743f4c03157b0ce49c675e6/mirnylib/numutils_new.pyx?at=default

    Main method for correcting DS and SS read data.
    Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction
    :param matrix: a scipy sparse matrix
    :param tolerance: Tolerance is the maximum allowed relative
                      deviation of the marginals.
    """
    if verbose:
        log.setLevel(logging.INFO)

    total_bias = np.ones(matrix.shape[0], 'float64')

    if np.isnan(matrix.sum()):
        log.warn("[iterative correction] the matrix contains nans, they will be replaced by zeros.")
        matrix.data[np.isnan(matrix.data)] = 0

    matrix = matrix.astype(float)
    W = matrix.tocoo()

    if np.abs(matrix - matrix.T).mean() / (1. * np.abs(matrix.mean())) > 1e-10:
        raise ValueError("Please provide symmetric matrix!")

    start_time = time.time()
    log.info("starting iterative correction")
    for iternum in range(M):
        iternum += 1
        s = np.array(W.sum(axis=1)).flatten()
        mask = (s == 0)
        s = s / np.mean(s[~mask])

        total_bias *= s
        deviation = np.abs(s - 1).max()

        s = 1.0 / s

        # The following code  is an optimization of this
        # for i in range(N):
        #     for j in range(N):
        #         W[i,j] = W[i,j] / (s[i] * s[j])

        W.data *= np.take(s, W.row)
        W.data *= np.take(s, W.col)
        if np.any(W.data > 1e100):
            log.error("*Error* matrix correction is producing extremely large values. "
                      "This is often caused by bins of low counts. Use a more stringent "
                      "filtering of bins.")
            exit(1)
        if verbose:
            if iternum % 5 == 0:
                end_time = time.time()
                estimated = (float(M - iternum) * (end_time - start_time)) / iternum
                m, sec = divmod(estimated, 60)
                h, m = divmod(m, 60)
                log.info("pass {} Estimated time {:.0f}:{:.0f}:{:.0f}".format(iternum, h, m, sec))
                log.info("max delta - 1 = {} ".format(deviation))

        if deviation < tolerance:
            log.info("[iterative correction] {} iterations used\n".format(iternum + 1))
            break

    # scale the total bias such that the sum is 1.0
    corr = total_bias[total_bias != 0].mean()
    total_bias = np.divide(total_bias, corr)
    W.data = W.data * corr * corr
    if np.any(W.data > 1e10):
        log.error("*Error* matrix correction produced extremely large values. "
                  "This is often caused by bins of low counts. Use a more stringent "
                  "filtering of bins.")
        exit(1)

    return W.tocsr(), total_bias
