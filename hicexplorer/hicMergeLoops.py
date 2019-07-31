import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import pandas as pd
from pybedtools import BedTool
import logging
log = logging.getLogger(__name__)

from intervaltree import Interval, IntervalTree
import hicmatrix.HiCMatrix as hm

from hicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
        Merges the loop locations which were detected on different resolutions.
        Loops need to have format as follows:
        chr start end chr start end

        A merge happens if x and y position of a loop overlap in both cases, all loops are considered with bin sizes of the lowest resolution.
        This means for a loop with coordinates x and y, the overlap to all other loops is search for (x - lowest resolution) and (y + lowest resolution).
        If two or more locations should be merged, the one with the lowest resolution is taken as the merged loop.
       
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--inputFiles', '-i',
                                help='The loop files from hicDetectLoops. To use files from other sources, '
                                'please follow \'chr start end chr start end\' format and remove any header.',
                                required=True,
                                nargs='+')
    parserRequired.add_argument('--outFileName', '-o',
                                help='The name of the merged loop file.',
                                required=True)
    parserRequired.add_argument('--lowestResolution', '-r',  # loop or domain
                                help='The lowest resolution of all loop files, i.e. 5kb, 10kb and 25kb, please use 25000.',
                                required=True,
                                type=int)

    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def mergeLoops(pDataFrame, pLowestResolution, pTupleX, pTupleY):
    hic = hm.hiCMatrix()
    target_regions_intervaltree_x = hic.intervalListToIntervalTree(pTupleX)[0]
    target_regions_intervaltree_y = hic.intervalListToIntervalTree(pTupleY)[0]
    new_loop_list = []

    for i, loop in enumerate(pDataFrame.values):
        # neighborhood factor to extent the search range. This allows to consider the smaller bin sizes
        # like they would be bins of the lowest resolution
        neighborhood_factor_x = int(
            pLowestResolution) - abs(int(loop[2]) - int(loop[1]))
        neighborhood_factor_y = int(
            pLowestResolution) - abs(int(loop[5]) - int(loop[4]))

        if loop[0] in target_regions_intervaltree_x:
            x_interval = target_regions_intervaltree_x[loop[0]].overlap(
                loop[1]-neighborhood_factor_x-1, loop[2]+neighborhood_factor_x+1)
        if loop[3] in target_regions_intervaltree_y:
            y_interval = target_regions_intervaltree_y[loop[0]].overlap(
                loop[4]-neighborhood_factor_y-1, loop[5]+neighborhood_factor_y+1)

        if len(x_interval) <= 1 or len(y_interval) <= 1:
            continue

        dict_of_interest_x = {}
        list_of_interest = []
        for data in x_interval:
            dict_of_interest_x[data[2]] = [data[0], data[1]]
        for data in y_interval:
            if data[2] in dict_of_interest_x:
                list_of_interest.append(data)

        max_index = 0
        max_distance = 0
        all_id_list = []
        for data in list_of_interest:
            if abs(data[0] - data[1]) > max_distance:
                max_distance = abs(data[0] - data[1])
                max_index = data[2]
            all_id_list.append(data[2])
        for data in x_interval:
            if data[2] == max_index:
                continue
            if data[2] not in all_id_list:
                continue
            target_regions_intervaltree_x[loop[0]].remove(data)

        for data in y_interval:
            if data[2] == max_index:
                continue
            if data[2] not in all_id_list:
                continue
            target_regions_intervaltree_y[loop[0]].remove(data)

    result_list = []
    result_list_index = []
    dict_x = {}
    dict_y = {}
    for chromosome_x, chromosome_y in zip(target_regions_intervaltree_x, target_regions_intervaltree_y):
        target_regions_intervaltree_x[chromosome_x] = sorted(
            target_regions_intervaltree_x[chromosome_x])
        target_regions_intervaltree_y[chromosome_y] = sorted(
            target_regions_intervaltree_y[chromosome_y])

        for x in target_regions_intervaltree_x[chromosome_x]:
            dict_x[x[2]] = (x[0], x[1])
        for y in target_regions_intervaltree_y[chromosome_y]:
            dict_y[y[2]] = (y[0], y[1])
        for x in dict_x:
            if x in dict_y:
                result_list_index.append(x)

        dict_x = None
        dict_x = {}
        dict_y = None
        dict_y = {}
    return result_list_index


def readFile(pFile):
    return pd.read_csv(pFile, sep='\t', header=(-1))


def main(args=None):

    lowest_resolution = args.lowestResolution

    files = args.inputFiles
    outfile_name = args.outFileName

    dataframe = None

    for file in files:
        if dataframe is None:
            dataframe = readFile(file)
        else:
            dataframe = dataframe.append(readFile(file), ignore_index=True)

    dataframe_bedtool = BedTool.from_dataframe(dataframe)
    dataframe = dataframe_bedtool.sort().to_dataframe(
        disable_auto_names=True, header=None)

    dataframe.drop_duplicates(keep=False, inplace=True)

#     dataframe.sort_values([0,1,2,3,4,5], ascending=[True, True, True, True, True, True])
    tuples_x = [tuple(x) for x in dataframe[[0, 1, 2]].values]
    tuples_y = [tuple(x) for x in dataframe[[3, 4, 5]].values]

    result_list_index = mergeLoops(
        dataframe, lowest_resolution, tuples_x, tuples_y)
    result_dataframe = dataframe.iloc[sorted(result_list_index), :]
    result_dataframe.to_csv(outfile_name, sep='\t', header=False, index=False)
