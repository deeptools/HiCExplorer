import numpy as np
import time
import sys


def iterativeCorrection(matrix, v=None, M=50, diag=-1,
                        tolerance=1e-5, verbose=False):
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

    total_bias = np.ones(matrix.shape[0], 'float64')
    if v is None:
        v = np.zeros(matrix.shape[0], 'float64')  # single-sided reads

    if np.isnan(matrix.sum()):
        sys.stderr.write("[iterative correction] the matrix contains nans, "\
        "they will be replaced by zeros\n")
        matrix.data[np.isnan(matrix.data)] = 0

    matrix = matrix.astype(float)
    W = matrix.tocoo()

    if np.abs(matrix - matrix.T).mean() / (1. * np.abs(matrix.mean())) > 1e-10:
        raise ValueError("Please provide symmetric matrix!")

    v = np.array(v, 'float64')
    start_time = time.time()
    if verbose:
        print "starting iterative correction"
    for iternum in xrange(M):
        iternum += 1
        s = np.array(W.sum(axis=1)).flatten()
        for dd in xrange(diag + 1):  # excluding the diagonal
            if dd == 0:
                s -= W.diagonal()
            else:
                dia = np.diagonal(W.todense(), dd)  # unefficient
                #print dia
                s[dd:] = s[dd:] - dia
                s[:-dd] = s[:-dd] - dia
        mask = (s == 0)
        s = s / np.mean(s[mask == False])
        s[mask] = 1
        total_bias *= s
        deviation = np.abs(s - 1).max()

        s  = 1.0 / s

        # The following code  is an optimization of this
        # for i in range(N):
        #     for j in range(N):
        #         W[i,j] = W[i,j] / (s[i] * s[j])

        W.data *= np.take(s, W.row)
        W.data *= np.take(s, W.col)

        if verbose:
            if iternum % 5 == 0:
                end_time = time.time()
                estimated = \
                    (float(M - iternum ) * (end_time - start_time) ) / iternum
                m, sec = divmod(estimated, 60)
                h, m = divmod(m, 60)
                print "pass {} Estimated time {:.0f}:{:.0f}:{:.0f}".format(
                    iternum, h, m, sec)
                print "max delta - 1 = {} ".format(deviation)

        if deviation < tolerance:
            sys.stderr.write("[iterative correction] {} iterations "
                             "used\n".format(iternum + 1))
            break

    # scale the total bias such that the sum is 1.0
    corr = total_bias[total_bias != 0].mean()
    total_bias /= corr
    W.data = W.data * corr * corr

    return W.tocsr(), total_bias

def iterativeCorrection_dekker(x,v = None,M=50,diag = -1, 
                                    tolerance=1e-5):
    """Main method for correcting DS and SS read data. Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction"""
    if M == None:
        M = 9999
    totalBias = np.ones(len(x),float)    
    if v == None: v = np.zeros(len(x),float)  #single-sided reads    
    x = np.array(x,np.double,order = 'C')
    if not  np.abs(x-x.T).mean() / (1. * np.abs(x.mean())) < 1e-10:
        raise ValueError("Please provide symmetric matrix!")
    _x = x
    v = np.array(v,float,order = "C")        
    N = len(x)       
    for iternum in xrange(M):         
        s0 = np.sum(_x,axis = 1)         
        mask = [s0 == 0]
#        v[mask] = 0   #no SS reads if there are no DS reads here        
#        nv = v / (totalBias * (totalBias[mask==False]).mean())
        s = s0 #+ nv
        for dd in xrange(diag + 1):   #excluding the diagonal 
            if dd == 0:
                s -= np.diagonal(_x)
            else:
                dia = np.diagonal(_x,dd)
                #print dia
                s[dd:] = s[dd:] -  dia
                s[:-dd] = s[:-dd] - dia 
        s = s / np.mean(s[s0!=0])
        s[s0==0] = 1                
        totalBias *= s
          
        for i in range(N):
            for j in range(N):
                _x[i,j] = _x[i,j] / (s[i] * s[j])
        
        if iternum % 5==0 and iternum > 0:
            print "max delta - 1 = {} ".format(np.abs(s-1).max())

        if M == 9999:
            if np.abs(s-1).max() < tolerance:
                print "IC used {} iterations".format(iternum+1)
                break

                         
    corr = totalBias[s0!=0].mean()  #mean correction factor
    x  = x * corr * corr #renormalizing everything
    totalBias /= corr
    return x, v/totalBias, totalBias
