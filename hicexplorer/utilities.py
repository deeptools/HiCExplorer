"""hicexplorer.utilities, implemented in C++.

Every public function keeps the name, the signature and the behaviour it had
in HiCExplorer 3.7.6; its body calls hicexplorer._hicx.utilities
(cpp/python/src/api/utilities.cpp, cpp/core/include/hicx/utilities.hpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from matplotlib import use as mplt_use
mplt_use('Agg')
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import utilities as _impl


def readBed(pBedFile):
    return _impl.readBed(pBedFile)


def readTargetBed(pBedFile):
    return _impl.readTargetBed(pBedFile)


def writableFile(string):
    return _impl.writableFile(string)


def remove_outliers(data, max_deviation=3.5):
    """
    The method is based on the median absolute deviation. See
    Boris Iglewicz and David Hoaglin (1993),
    "Volume 16: How to Detect and Handle Outliers",
    The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    returns the list, without the outliers
    """
    return _impl.remove_outliers(data, max_deviation)


def convertNansToZeros(ma):
    return _impl.convertNansToZeros(ma)


def convertInfsToZeros(ma):
    return _impl.convertInfsToZeros(ma)


def convertInfsToZeros_ArrayFloat(pArray, pToEpsilon=False):
    return _impl.convertInfsToZeros_ArrayFloat(pArray, pToEpsilon)


def convertNansToOnes(pArray):
    return _impl.convertNansToOnes(pArray)


def myAverage(valuesArray, avgType='mean'):
    return _impl.myAverage(valuesArray, avgType)


def enlarge_bins(bin_intervals):
    r"""
    takes a list of consecutive, but not
    directly touching, bin intervals
    and joins them such that the
    end and start of consecutive bins
    is the same.
    >>> bin_intervals = [('chr1', 10, 50, 1), ('chr1', 50, 80, 2),
    ... ('chr2', 10, 60, 3), ('chr2', 70, 90, 4)]
    >>> enlarge_bins(bin_intervals)
    [('chr1', 0, 50, 1), ('chr1', 50, 80, 2), ('chr2', 0, 65, 3), ('chr2', 65, 90, 4)]
    """
    return _impl.enlarge_bins(bin_intervals)


def genomicRegion(string):
    """
    validates and cleans up a string corresponding to a genomic
    region in the form ideally of chromosome:start-end
    but other forms are also possible like start
    and end containing comas.
    This code is intended to be used to validate and
    format a argparse parameter.
    :return: string in the form chrom:start:end
    """
    return _impl.genomicRegion(string)


def getUserRegion(chromSizes, regionString, max_chunk_size=1e6):
    """
    Verifies if a given region argument, given by the user
    is valid. The format of the regionString is chrom:start:end:tileSize
    where start, end and tileSize are optional.
    # this should work in doctest but it does not. So I
    # commented it.
    #>>> data = getUserRegion({'chr2': 1000}, "chr1:10:10")
    #Traceback (most recent call last):
    #    ...
    #NameError: Unknown chromosome: chr1
    #Known chromosomes are: ['chr2']
    >>> getUserRegion({'chr2': 1000}, "chr2:10:1001")
    ([('chr2', 1000)], 10, 1000, 990)

    #Test chunk and regions size reduction to match tile size
    >>> getUserRegion({'chr2': 200000}, "chr2:10:123344:3")
    ([('chr2', 123344)], 9, 123345, 123336)
    """
    return _impl.getUserRegion(chromSizes, regionString, max_chunk_size)


def expected_interactions_in_distance(pLength_chromosome, pChromosome_count, pSubmatrix):
    """
        Computes the function I_chrom(s) for a given chromosome.
    """
    return _impl.expected_interactions_in_distance(pLength_chromosome, pChromosome_count, pSubmatrix)


def expected_interactions_non_zero(pSubmatrix):
    """
        Computes the expected number of interactions per distance
    """
    return _impl.expected_interactions_non_zero(pSubmatrix)


def expected_interactions_thread(pData, pDistances, pMinDistance, pMaxDistance, pSize, pQueue):
    return _impl.expected_interactions_thread(pData, pDistances, pMinDistance, pMaxDistance, pSize, pQueue)


def expected_interactions(pSubmatrix, pThreads=None):
    """
        Computes the expected number of interactions per distance
    """
    return _impl.expected_interactions(pSubmatrix, pThreads)


def compute_zscore(pSubmatrix, pDepth, pThreads):
    return _impl.compute_zscore(pSubmatrix, pDepth, pThreads)


def obs_exp_matrix_lieberman(pSubmatrix, pLength_chromosome, pChromosome_count):
    """
        Creates normalized contact matrix M* by
        dividing each entry by the gnome-wide
        expected contacts for loci at
        that genomic distance. Method: Lieberman-Aiden 2009
    """
    return _impl.obs_exp_matrix_lieberman(pSubmatrix, pLength_chromosome, pChromosome_count)


def obs_exp_matrix_non_zero(pSubmatrix, ligation_factor=False, pInplace=True, pToEpsilon=False, pThreads=None):
    """
        Creates normalized contact matrix M* by
        dividing each entry by the gnome-wide
        expected contacts for loci at
        that genomic distance.
        exp_i,j = sum(interactions at distance abs(i-j)) / number of non-zero
        interactions at abs(i-j). If ligation_factor, then
        exp_i,j = exp_i,j * sum(row(i)) * sum(row(j)) / sum(matrix)
        This factor has been used by Homer software to correct for the effect
        of proximity ligation
    """
    return _impl.obs_exp_matrix_non_zero(pSubmatrix, ligation_factor, pInplace, pToEpsilon, pThreads)


def obs_exp_matrix(pSubmatrix, pInplace=True, pToEpsilon=False, pThreads=None, pDistance=None):
    """
        Creates normalized contact matrix M* by
        dividing each entry by the gnome-wide
        expected contacts for loci at
        that genomic distance.
        exp_i,j = sum(interactions at distance abs(i-j)) / number of non-zero
        interactions at abs(i-j)
    """
    return _impl.obs_exp_matrix(pSubmatrix, pInplace, pToEpsilon, pThreads, pDistance)


def toString(s):
    """
    This takes care of python2/3 differences
    """
    return _impl.toString(s)


def toBytes(s):
    """
    Like toString, but for functions requiring bytes in python3
    """
    return _impl.toBytes(s)


def check_chrom_str_bytes(pIteratableObj, pObj):
    return _impl.check_chrom_str_bytes(pIteratableObj, pObj)


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.
    """
    return _impl.change_chrom_names(chrom)


def opener(filename):
    """
    Determines if a file is compressed or not
    """
    return _impl.opener(filename)


def remove_non_ascii(pText):
    """
    This function converts all non-ascii characters to a most alike representation.
    Code from:
    https://stackoverflow.com/questions/20078816/replace-non-ascii-characters-with-a-single-space/20079244
    """
    return _impl.remove_non_ascii(pText)


def check_cooler(pFileName):
    return _impl.check_cooler(pFileName)


def in_units(pBasePosition):
    return _impl.in_units(pBasePosition)

