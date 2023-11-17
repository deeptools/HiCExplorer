import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import sys
import numpy as np
import argparse
from matplotlib import use as mplt_use
mplt_use('Agg')
from unidecode import unidecode
import cooler
from copy import deepcopy
import time
from multiprocessing import Process, Queue
import traceback
import logging
log = logging.getLogger(__name__)


def readBed(pBedFile):
    viewpoints = []
    with open(pBedFile, 'r') as file:
        for line in file.readlines():
            if line.startswith('#'):
                continue
            _line = line.strip().split('\t')
            if len(line) == 0:
                continue
            try:
                chrom, start, end = _line[:3]
            except ValueError:
                _line = line.strip().split()
                chrom, start, end = _line[:3]
            viewpoints.append((chrom, start, end))

    return viewpoints


def writableFile(string):
    try:
        open(string, 'w').close()
    except IOError:
        msg = "{} file can be opened for writing".format(string)
        log.debug(msg)
        raise argparse.ArgumentTypeError(msg)
    return string


def remove_outliers(data, max_deviation=3.5):
    """
    The method is based on the median absolute deviation. See
    Boris Iglewicz and David Hoaglin (1993),
    "Volume 16: How to Detect and Handle Outliers",
    The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    returns the list, without the outliers
    """
    median = np.median(data)
    b_value = 1.4826  # value for normal distribution
    mad = b_value * np.median(np.abs(data - median))

    if mad > 0:
        deviation = abs(data - median) / mad
        """
        outliers = data[deviation > max_deviation]
        print "outliers removed {}".format(len(outliers))
        print outliers
        """
        data = data[deviation <= max_deviation]
    return data


def convertNansToZeros(ma):
    if len(ma.data) == 0:
        return ma
    nan_elements = np.flatnonzero(np.isnan(ma.data))
    # data_type = type(ma.data[0])
    if len(nan_elements) > 0:
        # if data_type == np.float
        ma.data[nan_elements] = 0.0
    return ma


def convertInfsToZeros(ma):
    if len(ma.data) == 0:
        return ma
    inf_elements = np.flatnonzero(np.isinf(ma.data))
    if len(inf_elements) > 0:
        ma.data[inf_elements] = 0
    return ma


def convertInfsToZeros_ArrayFloat(pArray, pToEpsilon=False):
    if len(pArray) == 0:
        return pArray
    nan_elements = np.flatnonzero(np.isnan(pArray))
    if len(nan_elements) > 0:
        if pToEpsilon:
            pArray[nan_elements] = 0.000001
        else:
            pArray[nan_elements] = 0.0

    inf_elements = np.flatnonzero(np.isinf(pArray))
    if len(inf_elements) > 0:
        if pToEpsilon:
            pArray[inf_elements] = 0.000001
        else:
            pArray[inf_elements] = 0.0
    return pArray


def convertNansToOnes(pArray):
    if len(pArray) == 0:
        return pArray
    nan_elements = np.flatnonzero(np.isnan(pArray))
    if len(nan_elements) > 0:
        pArray[nan_elements] = 1.0
    return pArray


def myAverage(valuesArray, avgType='mean'):

    valuesArray = valuesArray[np.logical_not(np.isnan(valuesArray))]
    if avgType == 'mean':
        mean = np.mean(valuesArray)
    else:
        mean = np.median(valuesArray)

    return mean


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
    # enlarge remaining bins
    chr_start = True
    for idx in range(len(bin_intervals) - 1):
        chrom, start, end, extra = bin_intervals[idx]
        chrom_next, start_next, end_next, extra_next = bin_intervals[idx + 1]

        if chr_start is True:
            start = 0
            chr_start = False
            bin_intervals[idx] = (chrom, start, end, extra)
        if chrom == chrom_next and end != start_next:
            middle = start_next - int((start_next - end) / 2)
            bin_intervals[idx] = (chrom, start, middle, extra)
            bin_intervals[idx + 1] = (chrom, middle, end_next, extra_next)
        if chrom != chrom_next:
            chr_start = True

    chrom, start, end, extra = bin_intervals[-1]
    bin_intervals[-1] = (chrom, start, end, extra)

    return bin_intervals


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
    # remove whitespaces using split,join trick
    region = ''.join(string.split())
    if region == '':
        return None
    # remove undesired characters that may be present and
    # replace - by :
    if sys.version_info[0] == 2:
        region = region.translate(None, ",;|!{}()").replace("-", ":")
    if sys.version_info[0] == 3:
        region = region.translate(str.maketrans('', '', ",;|!{}()")).replace("-", ":")
    if len(region) == 0:
        raise argparse.ArgumentTypeError(
            "{} is not a valid region".format(string))
    return region


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
    region = regionString.split(":")
    chrom = region[0]
    chromSizes = dict(chromSizes)

    try:
        chromSizes[chrom]
    except KeyError:
        raise NameError("Unknown chromosome: %s\nKnown "
                        "chromosomes are: %s " % (chrom, list(chromSizes)))
    try:
        regionStart = int(region[1])
    except IndexError:
        regionStart = 0
    try:
        regionEnd = int(region[2]) if int(region[2]) <= chromSizes[chrom] \
            else chromSizes[chrom]
    except IndexError:
        regionEnd = chromSizes[chrom]
    if regionStart > regionEnd or regionStart < 0:
        raise NameError("%s not valid. The format is chrom:start:end. "
                        "Without comas, dashes or dots. " % (regionString))
    try:
        tilesize = int(region[3])
    except IndexError:
        tilesize = None

    chromSizes = [(chrom, regionEnd)]

    # if tilesize is given, make regionStart and regionEnd
    # multiple of tileSize
    if tilesize:
        regionStart -= regionStart % tilesize
        regionEnd += tilesize - (regionEnd % tilesize)

    chunkSize = int(regionEnd - regionStart)
    if chunkSize > max_chunk_size:
        chunkSize = max_chunk_size
        if tilesize and tilesize < chunkSize:
            chunkSize -= chunkSize % tilesize

    return (chromSizes, regionStart, regionEnd, int(chunkSize))


def expected_interactions_in_distance(pLength_chromosome, pChromosome_count, pSubmatrix):
    """
        Computes the function I_chrom(s) for a given chromosome.
    """
    expected_interactions = np.zeros(pSubmatrix.shape[0])
    row, col = pSubmatrix.nonzero()
    distance = np.absolute(row - col)

    for i, distance_ in enumerate(distance):
        expected_interactions[distance_] += pSubmatrix.data[i]

    count_times_i = np.arange(float(len(expected_interactions)))
    pChromosome_count = np.int32(pChromosome_count)
    pLength_chromosome = np.int32(pLength_chromosome)
    count_times_i *= pChromosome_count
    count_times_i -= pLength_chromosome
    count_times_i *= np.int32(-1)

    expected_interactions = np.divide(expected_interactions, count_times_i)
    # log.debug('exp_obs_matrix_lieberman {}'.format(expected_interactions))

    return expected_interactions


def expected_interactions_non_zero(pSubmatrix):
    """
        Computes the expected number of interactions per distance
    """

    expected_interactions = np.zeros(pSubmatrix.shape[0])
    row, col = pSubmatrix.nonzero()
    distance = np.absolute(row - col)

    occurences = np.zeros(pSubmatrix.shape[0])
    # data_type = type(pSubmatrix.data[0])
    for i, distance_ in enumerate(distance):
        expected_interactions[distance_] += pSubmatrix.data[i]
        occurences[distance_] += 1
    expected_interactions = np.divide(expected_interactions, occurences)

    mask = np.isnan(expected_interactions)
    expected_interactions[mask] = 0
    mask = np.isinf(expected_interactions)
    expected_interactions[mask] = 0

    return expected_interactions


def expected_interactions_thread(pData, pDistances, pMinDistance, pMaxDistance, pSize, pQueue):
    try:
        expected_interactions = np.zeros(pSize)

        for i in range(pMinDistance, pMaxDistance, 1):
            mask = pDistances == i
            expected_interactions[i] = np.sum(pData[mask])

        pQueue.put(expected_interactions)
        return
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return


def expected_interactions(pSubmatrix, pThreads=None):
    """
        Computes the expected number of interactions per distance
    """
    expected_interactions = np.zeros(pSubmatrix.shape[0])

    row, col = pSubmatrix.nonzero()
    distance = np.absolute(row - col)
    occurrences = np.arange(pSubmatrix.shape[0] + 1, 1, -1)
    # occurences = np.zeros(pSubmatrix.shape[0])

    if len(distance) == 0:
        return None
    min_distance = distance.min()
    max_distance = distance.max()
    # time_start = time.time()
    if pThreads is not None and pThreads:
        queue = [None] * pThreads
        process = [None] * pThreads
        distances_per_threads = (max_distance - min_distance) // pThreads
        all_data_collected = False
        thread_done = [False] * pThreads

        for i in range(pThreads):

            if i < pThreads - 1:
                min_distance_thread = min_distance + (i * distances_per_threads)
                max_distance_thread = min_distance + ((i + 1) * distances_per_threads)
            else:
                min_distance_thread = min_distance + (i * distances_per_threads)
                max_distance_thread = max_distance + 1
            queue[i] = Queue()
            process[i] = Process(target=expected_interactions_thread, kwargs=dict(
                pData=pSubmatrix.data,
                pDistances=distance,
                pMinDistance=min_distance_thread,
                pMaxDistance=max_distance_thread,
                pSize=pSubmatrix.shape[0],
                pQueue=queue[i]
            )
            )

            process[i].start()
        fail_flag = False
        fail_message = ''
        while not all_data_collected:
            for i in range(pThreads):
                if queue[i] is not None and not queue[i].empty():
                    expected_interactions_thread_ = queue[i].get()
                    if isinstance(expected_interactions_thread_, str) and 'Fail: ' in expected_interactions_thread_:
                        fail_flag = True
                        fail_message = expected_interactions_thread_
                    else:
                        expected_interactions += expected_interactions_thread_
                    queue[i] = None
                    process[i].join()
                    process[i].terminate()
                    process[i] = None
                    thread_done[i] = True
            all_data_collected = True
            for thread in thread_done:
                if not thread:
                    all_data_collected = False
            time.sleep(1)
        if fail_flag:
            return fail_message
    else:

        for i in range(min_distance, max_distance + 1, 1):
            mask = distance == i
            expected_interactions[i] = np.sum(pSubmatrix.data[mask])
    # log.info('exp inter: {}'.format(time.time() - time_start))
    # for i, distance_ in enumerate(distance):
    #     expected_interactions[distance_] += pSubmatrix.data[i]
        # occurences[distance_] += 1
    expected_interactions = np.divide(expected_interactions, occurrences)

    mask = np.isnan(expected_interactions)
    expected_interactions[mask] = 0
    mask = np.isinf(expected_interactions)
    expected_interactions[mask] = 0

    del row
    del col
    del distance
    del occurrences
    del mask

    return expected_interactions

# def expected_interactions(pSubmatrix):

#     instances, features = pSubmatrix.nonzero()
#     distances = np.absolute(instances - features)
#     sum_per_distance = np.ones(pSubmatrix.shape[0])
#     binary_interactions_per_distance = np.ones(pSubmatrix.shape[0])

#     for i, distance in enumerate(distances):
#         sum_per_distance[distance] += pSubmatrix.data[i]
#         binary_interactions_per_distance[distance] += 1

#     return sum_per_distance / binary_interactions_per_distance


def compute_zscore(pSubmatrix, pDepth, pThreads):
    # depth is to be expected in matrix units
    # indice == distance

    expected_interactions_array = expected_interactions(pSubmatrix, pThreads=pThreads)
    occurrences = np.arange(pSubmatrix.shape[0] + 1, 1, -1)
    row, col = pSubmatrix.nonzero()
    distance = np.absolute(row - col)
    if pDepth is not None:
        max_distance = min(pSubmatrix.shape[0], pDepth)
        mask = distance > pDepth
        pSubmatrix.data[mask] = 0
    else:
        max_distance = pSubmatrix.shape[0]

    x_minus_mu = np.zeros(len(pSubmatrix.data))
    # sum_for_sigma = np.zeros(len(expected_interactions_array))
    for distance_index in range(max_distance):
        mask = distance == distance_index
        x_minus_mu[mask] -= expected_interactions_array[distance_index]
        sum_for_sigma = np.sum(np.square(x_minus_mu[mask]))

        sum_for_sigma += np.square((0 - expected_interactions_array[distance_index]) * (occurrences[distance_index] - np.sum(mask)))
        pSubmatrix.data[mask] = x_minus_mu[mask] / np.sqrt(sum_for_sigma / occurrences[distance_index])

    return pSubmatrix


def obs_exp_matrix_lieberman(pSubmatrix, pLength_chromosome, pChromosome_count):
    """
        Creates normalized contact matrix M* by
        dividing each entry by the gnome-wide
        expected contacts for loci at
        that genomic distance. Method: Lieberman-Aiden 2009
    """

    expected_interactions_in_distance_ = expected_interactions_in_distance(pLength_chromosome, pChromosome_count, pSubmatrix)
    row, col = pSubmatrix.nonzero()
    distance = np.ceil(np.absolute(row - col) / 2).astype(np.int32)

    if len(pSubmatrix.data) > 0:
        data_type = type(pSubmatrix.data[0])

        expected = expected_interactions_in_distance_[distance]
        pSubmatrix.data = pSubmatrix.data.astype(np.float32)
        pSubmatrix.data = np.divide(pSubmatrix.data, expected)
        pSubmatrix.data = convertInfsToZeros_ArrayFloat(pSubmatrix.data).astype(data_type)
    return pSubmatrix


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
    if pInplace:
        submatrix = pSubmatrix
    else:
        submatrix = deepcopy(pSubmatrix)
    expected_interactions_in_distance = expected_interactions_non_zero(submatrix)

    row_sums = np.array(submatrix.sum(axis=1).T).flatten()
    total_interactions = submatrix.sum()

    row, col = submatrix.nonzero()

    submatrix.data = submatrix.data.astype(np.float32)

    for i in range(len(row)):
        expected = expected_interactions_in_distance[np.absolute(row[i] - col[i])]
        if ligation_factor:
            expected *= row_sums[row[i]] * row_sums[col[i]] / total_interactions

        submatrix.data[i] = np.divide(submatrix.data[i], expected)

    if pToEpsilon:
        epsilon = 0.000000001
    else:
        epsilon = 0
    mask = np.isnan(submatrix.data)
    submatrix.data[mask] = epsilon
    mask = np.isinf(submatrix.data)
    submatrix.data[mask] = epsilon
    submatrix.eliminate_zeros()
    return submatrix


def obs_exp_matrix(pSubmatrix, pInplace=True, pToEpsilon=False, pThreads=None, pDistance=None):
    """
        Creates normalized contact matrix M* by
        dividing each entry by the gnome-wide
        expected contacts for loci at
        that genomic distance.
        exp_i,j = sum(interactions at distance abs(i-j)) / number of non-zero
        interactions at abs(i-j)
    """
    # time_start = time.time()
    if pDistance is not None:
        row, col = pSubmatrix.nonzero()
        distance = np.absolute(row - col)
        mask = distance >= pDistance
        pSubmatrix.data[mask] = 0
        pSubmatrix.eliminate_zeros()

    expected_interactions_in_distance_ = expected_interactions(pSubmatrix, pThreads)
    if expected_interactions_in_distance_ is None:
        return None
    # log.info('time exp: {}'.format(time.time() - time_start))
    # time_start = time.time()

    row, col = pSubmatrix.nonzero()
    distance = np.ceil(np.absolute(row - col) / 2).astype(np.int32)
    if not pInplace:
        pSubmatrix_copy = deepcopy(pSubmatrix)
    if len(pSubmatrix.data) > 0:
        data_type = type(pSubmatrix.data[0])

        expected = expected_interactions_in_distance_[distance]
        if pInplace:
            pSubmatrix.data = pSubmatrix.data.astype(np.float32)
            pSubmatrix.data = np.divide(pSubmatrix.data, expected)
            pSubmatrix.data = convertInfsToZeros_ArrayFloat(pSubmatrix.data, pToEpsilon).astype(data_type)
        else:
            pSubmatrix_copy.data = pSubmatrix_copy.data.astype(np.float32)
            pSubmatrix_copy.data = np.divide(pSubmatrix_copy.data, expected)
            pSubmatrix_copy.data = convertInfsToZeros_ArrayFloat(pSubmatrix_copy.data, pToEpsilon).astype(data_type)
        del expected

    del expected_interactions_in_distance_
    del row
    del col
    del distance
    # log.info('test obs/exp: {}'.format(time.time() - time_start))

    if pInplace:
        return pSubmatrix
    else:
        return pSubmatrix_copy


def toString(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):  # or isinstance(s, np.bytes_):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')
    if isinstance(s, list):
        return [toString(x) for x in s]
    if isinstance(s, np.ndarray):
        return s.astype(str)
    return s


def toBytes(s):
    """
    Like toString, but for functions requiring bytes in python3
    """
    if sys.version_info[0] == 2:
        return s
    if isinstance(s, bytes):
        return s
    # if isinstance(s, np.bytes_):
    #     return np.bytes_(s)
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [toBytes(x) for x in s]
    return s


def check_chrom_str_bytes(pIteratableObj, pObj):
    # determine type
    if isinstance(pObj, list) and len(pObj) > 0:
        type_ = type(pObj[0])
    else:
        type_ = type(pObj)
    if not isinstance(type(next(iter(pIteratableObj))), type_):
        if type(next(iter(pIteratableObj))) is str:
            pObj = toString(pObj)
        elif type(next(iter(pIteratableObj))) in [bytes, np.bytes_]:
            pObj = toBytes(pObj)
    return pObj


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.
    """
    # TODO: mapping from chromosome names like mithocondria is missing

    # python 2 / 3 issue with string, bytes and np.bytes_
    # if chrom.startswith('chr'):

    chrom = toString(chrom)
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom


def opener(filename):
    """
    Determines if a file is compressed or not
    """
    import gzip
    f = open(filename, 'rb')
    # print("gzip or not?", f.read(2))

    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f


# def check_chrom_str_bytes(pChrom, pInstanceToCompare):
#     """
#     Checks and changes pChroms to str or bytes depending on datatype
#     of pInstanceToCompare
#     """
#     if type(next(iter(pInstanceToCompare))) is str:
#         pChrom = toString(pChrom)
#     elif type(next(iter(pInstanceToCompare))) in [bytes, np.bytes_]:
#         pChrom = toBytes(pChrom)
#     return pChrom
def remove_non_ascii(pText):
    """
    This function converts all non-ascii characters to a most alike representation.
    Code from:
    https://stackoverflow.com/questions/20078816/replace-non-ascii-characters-with-a-single-space/20079244
    """
    return unidecode(pText)


def check_cooler(pFileName):
    if pFileName.endswith('.cool') or cooler.fileops.is_cooler(pFileName) or '.mcool' in pFileName:
        return True
    return False


def in_units(pBasePosition):
    pBasePosition = float(pBasePosition)
    # log.debug("pBasePosition {}".format(pBasePosition))
    if pBasePosition > 1.5e6:
        labels = "{:.2f} ".format(pBasePosition / 1e6)
        labels += " Mbp"
    elif pBasePosition > 1500:
        labels = "{:.0f}".format(pBasePosition / 1e3)
        labels += " Kbp"
    else:
        labels = "{:.2f} ".format((pBasePosition))
        labels += " bp"
    return labels
