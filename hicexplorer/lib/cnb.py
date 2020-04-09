import numpy as np

import logging
log = logging.getLogger(__name__)
import copy

import sys
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import nbinom
from scipy.special import gammaln
from scipy import special

def pdf(pX, pR, pP):
    """
    PDF for a continuous generalization of NB distribution
    """

    gamma_part = gammaln(pR + pX) - gammaln(pX + 1) - gammaln(pR)
    return np.exp(gamma_part + (pR * np.log(pP)) + special.xlog1py(pX, -pP))

def cdf(pX, pR, pP):
    """
    Cumulative density function of a continuous generalization of NB distribution
    """
    # if pX == 0:
    # return 0
    return special.betainc(pR, pX + 1, pP)