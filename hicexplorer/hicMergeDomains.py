import argparse
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib import pyplot as plt
from graphviz import Digraph
import os
import logging
log = logging.getLogger(__name__)
from hicexplorer._version import __version__

def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""
hicMergeDomains takes as input multiple TAD domain files from hicFindTads. It merges TADs from different resolutions to one TAD domains file,
considers protein peaks from known TAD binding sites and computes a dependency graph of the TADs.

Two TADs are considered as one if they don't overlap at x bins given by `--value`; TAD borders need to match the protein peaks given by `--proteinFile`;
a relation between two TADs is given by their overlap of area in percent, parameter `--percent`. The protein peaks are only considered if in one bin at least `--minPeak`.

An example usage is:

`$ hicMergeDomains --domainFiles 10kbtad_domains.bed 50kbtad_domains.bed --proteinFile ctcf_sorted.bed --outputMergedList two_files_ctcf --outputRelationList two_files_relation_ctcf --outputTreePlotPrefix two_files_plot_ctcf --outputTreePlotFormat pdf`

""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--domainFiles', '-d',
                                help='The domain files of the different resolutions is required',
                                required=True,
                                nargs='+')

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--proteinFile', '-p',
                           help='In order to be able to better assess the relationship between TADs, the associated '
                                'protein file (e.g. CTCF for mammals) can be included. The protein file is required in broadpeak format')

    parserOpt.add_argument('--minimumNumberOfPeaks', '-m',
                           help='Optional parameter to adjust the number of protein peaks when adapting the resolution to '
                                'the domain files. At least minimumNumberOfPeaks of unique peaks must be in a bin to considered. Otherwise the bin is treated like it has no peaks.',
                           type=int, default=1)

    parserOpt.add_argument('--value', '-v',
                           help='Determine a value by how much the boundaries of two TADs must at least differ to consider '
                                'them as two separate TADs.',
                           type=int, default=5000)

    parserOpt.add_argument('--percent', '-pe',
                           help='For the relationship determination, a percentage is required from which area coverage '
                                'the TADs are related to each other.'
                           'For example, a relationship should be entered from 5 percent area coverage -p 0.05',
                           type=float, default=0.5)
    parserOpt.add_argument('--outputMergedList', '-om',
                           help='File name for the merged domains list',
                           default='mergedDomains.bed',
                           required=False)
    parserOpt.add_argument('--outputRelationList', '-or',
                           help='File name for the relationship list of the TADs.',
                           default='relationList.txt',
                           required=False)
    parserOpt.add_argument('--outputTreePlotPrefix', '-ot',
                           help='File name prefix for the relationship tree of the TADs',
                           default='relationship_tree_',
                           required=False)
    parserOpt.add_argument('--outputTreePlotFormat', '-of',
                           help='File format of the relationship tree. Supported formats are listed on: https://www.graphviz.org/doc/info/output.html',
                           default='pdf',
                           required=False)
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser


def create_list_of_file(file, binSizeList):
    """Creates a list from a file and recognizes its BinSize at the same time"""
    with open(file) as f:
        newList = [line.rstrip() for line in f]
    splittedList = []
    binSize = 10000000
    """ It is not necessary to iterate over the entire length of the list
        The first 20 positions are sufficient to determine the BinSize """
    pos = 0
    for line in newList:
        x = line.split("\t")
        splittedList.append(x)
        if pos < 20:
            pos += 1
            zeros, num = 0, x[1]
            while num[len(num) - 1:] == '0':
                zeros += 1
                num = num[0:(len(num) - 1)]
            num = num[len(num) - 1:]
            while zeros != 0:
                num += '0'
                zeros -= 1
            if int(num) < binSize:
                binSize = int(num)
    log.debug('binSize {}'.format(binSize))
    return splittedList, binSize


def merge_list(d1, d2, pValue=5000):
    """Combine two lists into one, taking into account the value,
       with which one can decide whether two very similar TADs should
       be seen as one or as two """
    pos1, pos2 = 0, 0
    merged_list = []
    while True:
        if pos1 == len(d1):
            break
        while d1[pos1][0] == d2[pos2][0]:
            """Checks whether the left boundary of TAD1 is less than or equal to the left boundary of TAD2"""
            if int(d1[pos1][1]) <= int(d2[pos2][1]):
                """ Checks whether either the left or right border of the two TADs exceed the minimum distance value"""
                if (abs(int(d1[pos1][1]) - int(d2[pos2][1])) > pValue) or (
                        abs(int(d1[pos1][2]) - int(d2[pos2][2])) > pValue):
                    merged_list.append(d1[pos1])
                if pos1 + 1 != len(d1):
                    pos1 += 1
                else:
                    """As long as List2 has not reached its end and still share the same chromosome, the remaining
                    TADs are added """
                    while (pos2 < len(d2)) and (d1[pos1][0] == d2[pos2][0]):
                        merged_list.append(d2[pos2])
                        pos2 += 1
                    break
                """Checks whether the left boundary of TAD2 is less than or equal to the left boundary of TAD1"""
            elif (int(d1[pos1][1]) > int(d2[pos2][1])):
                """ Checks whether either the left or right border of the two TADs exceed the minimum distance value"""
                if (abs(int(d1[pos1][1]) - int(d2[pos2][1])) > pValue) or (
                        abs(int(d1[pos1][2]) - int(d2[pos2][2])) > pValue):
                    merged_list.append(d2[pos2])
                if pos2 + 1 != len(d2):
                    pos2 += 1
                else:
                    """As long as List1 has not reached its end and still share the same chromosome, the remaining
                    TADs are added"""
                    while (pos1 < len(d1)) and (d1[pos1][0] == d2[pos2][0]):
                        merged_list.append(d1[pos1])
                        pos1 += 1
                    break
        if pos1 == len(d1):
            break
        old_pos2 = pos2
        """ check if the current chromosome at pos1 is in d2 """
        while pos2 < len(d2):
            if not (d1[pos1][0] == d2[pos2][0]):
                pos2 += 1
            else:
                break
        """ if the chromosome is not found in d2, add all TADs with the current chromosome from d1 to the list """
        if pos2 == len(d2):
            pos2 = old_pos2
            chrom = d1[pos1][0]
            while d1[pos1][0] == chrom:
                merged_list.append(d1[pos1])
                pos1 += 1
                if pos1 == len(d1):
                    break
    """ add all TADs from d2 to the list that have not yet been added """
    pos2 = 0
    while pos2 < len(d2):
        if d2[pos2] not in merged_list:
            merged_list.append(d2[pos2])
        pos2 += 1
    return merged_list


def add_id(domainList):
    """Give each element of a list its individual number """
    id_number = 1
    for tad in domainList:
        tad[3] = "ID_" + str(id_number)
        id_number += 1
    return domainList


def create_relationsship_list(domainList, pPercent=0.5):
    """Give each element of a list its individual number """
    relationList = []
    tad1, tad2 = 0, 1
    while tad1 < len(domainList):
        while tad2 < len(domainList):
            """ Since the list is sorted, the second loop is terminated if either tad1 is to the
                right of tad2 or if they are no longer on the same chromosome """
            if (int(domainList[tad1][2]) < int(domainList[tad2][1])) or (domainList[tad1][0] != domainList[tad2][0]):
                break
            """ The minimum area of TAD1 is calculated using the "percent" parameter in order
            to detect overlapping TADs """
            minArea = (float(domainList[tad1][2]) - float(domainList[tad1][1])) * pPercent
            """ Checks whether the two TADs are overlapping TADs """
            if ((float(domainList[tad1][2]) - minArea) > float(domainList[tad2][1])) and \
                    ((float(domainList[tad1][2]) + minArea) <= float(domainList[tad2][2])):
                if (float(domainList[tad1][2]) - float(domainList[tad1][1])) > \
                        (float(domainList[tad2][2]) - int(domainList[tad2][1])):
                    add_relation_to_list(relationList, domainList[tad1][0], domainList[tad1][3], domainList[tad2][3])
                else:
                    add_relation_to_list(relationList, domainList[tad2][0], domainList[tad2][3], domainList[tad1][3])
                """ Checks whether TAD1 is completely in TAD2 """
            elif ((not float(domainList[tad2][1]) < (float(domainList[tad1][1]))) and
                  (not int(domainList[tad2][2]) > float(domainList[tad1][2]))):
                add_relation_to_list(relationList, domainList[tad1][0], domainList[tad1][3], domainList[tad2][3])
            tad2 += 1
        tad2 = tad1 + 2
        tad1 += 1
    return relationList


def add_relation_to_list(rList, chromosom, parent, child):
    """ Add a discovered relationship to a list """
    if len(rList) == 0:
        rList.append([chromosom, parent, [child]])
        return rList
    pos = len(rList) - 1
    while True:
        if rList[pos][1] == parent:
            rList[pos][2].append(child)
            return rList
        elif int(rList[pos][1][3:]) < int(parent[3:]):
            rList.append([chromosom, parent, [child]])
            return rList
        else:
            pos -= 1


def write_in_file(pList, name, relation=False):
    """ Create and save various lists in their individual files """
    filename = name
    with open(filename, 'w') as myfile:
        i = 0
        if relation:
            while i < len(pList):
                element = 0
                while element < len(pList[i][2]):
                    string = pList[i][0] + '\t' + pList[i][1] + '\t' + pList[i][2][element] + '\n'
                    myfile.write(string)
                    element += 1
                i += 1
        else:
            while i < len(pList):
                element = 0
                string = ""
                while element < len(pList[i]) - 1:
                    string += pList[i][element] + '\t'
                    element += 1
                string += pList[i][element] + '\n'
                myfile.write(string)
                i += 1


def create_tree(rList, dList, pFileNameSuffix, pFormat):
    """ Create a tree diagram from a list and
       thus represent the relationships of individual TADs per chromosome """
    name = pFileNameSuffix + '_' + rList[0][0]
    g = Digraph(filename=name, strict=True, format=pFormat)
    sList = create_small_list(dList)
    chrom = rList[0][0]
    posRList = 0
    posSList = 0
    while sList[posSList][0] != chrom:
        posSList += 1
    while posRList < len(rList):
        if chrom == rList[posRList][0]:
            while int(sList[posSList][1][0][3:]) < int(rList[posRList][1][3:]):
                g.node(sList[posSList][1][0])
                sList[posSList][1].remove(sList[posSList][1][0])
            if rList[posRList][1] in sList[posSList][1]:
                sList[posSList][1].remove(rList[posRList][1])
            for child in rList[posRList][2]:
                g.edge(rList[posRList][1], child)
                if child in sList[posSList][1]:
                    sList[posSList][1].remove(child)
        else:
            while len(sList[posSList][1]) != 0:
                g.node(sList[posSList][1][0])
                sList[posSList][1].remove(sList[posSList][1][0])
            print("Saved relation tree of " + chrom)
            g.render(name, cleanup=True)
            chrom = rList[posRList][0]
            name = pFileNameSuffix + '_' + chrom
            g = Digraph(filename=name, format=pFormat)
            posSList = 0
            while sList[posSList][0] != chrom:
                posSList += 1
        posRList += 1
    print("Saved relation tree of " + chrom)
    g.render(cleanup=True)


def create_small_list(dList):
    """ Creates a nested list from a domain file, in which the TADs
    are sorted according to their chromosomes """
    sList = [[dList[0][0], []]]
    chrom = dList[0][0]
    pos = 0
    while pos < len(dList):
        if dList[pos][0] == chrom:
            sList[len(sList) - 1][1].append(dList[pos][3])
        else:
            chrom = dList[pos][0]
            sList.append([dList[pos][0], [dList[pos][3]]])
        pos += 1
    return sList


def read_protein(file):
    """ Reads a protein file and saves it in a list """
    with open(file) as f:
        newList = [line.rstrip() for line in f]
    splittedList = [[]]
    actualChr = newList[0].split("\t")[0]
    for line in newList:
        x = line.split("\t")
        if x[0] == actualChr:
            splittedList[len(splittedList) - 1].append(x[0:3])
        else:
            actualChr = x[0]
            splittedList.append([x[0:3]])
    return splittedList


def merge_protein(pProteinList, binSize, minPeak):
    """ Adapts the protein list to the corresponding BinSize, taking into
       account the minimum number of peaks per section """
    newProteinList = []
    for chromosome in pProteinList:
        newProteinList.append([])
        currentBoundaryLeft = 0
        currentBoundaryRight = binSize
        count = 0
        for peak in chromosome:
            if int(peak[1]) <= currentBoundaryRight:
                count += 1
            else:
                if count >= minPeak:
                    newProteinList[len(newProteinList) - 1].append([peak[0], currentBoundaryLeft, currentBoundaryRight, count])
                currentBoundaryLeft = currentBoundaryRight
                currentBoundaryRight = currentBoundaryLeft + binSize
                count = 0
                if int(peak[1]) < currentBoundaryRight:
                    count += 1
                else:
                    while int(peak[1]) > currentBoundaryLeft:
                        currentBoundaryLeft += binSize
                    currentBoundaryRight = currentBoundaryLeft + binSize
                    count = 1
    return newProteinList


def compare_boundaries_protein(bList, cList, paraScore=0.2):
    """ check the boundaries for a protein-peak """
    posTad = 0
    posPeak = 0
    chromPosition = 0
    removedTads = []

    while posTad < len(bList):
        if posTad < len(bList) and chromPosition < len(cList) and bList[posTad][0] != cList[chromPosition][0][0]:
            chromPosition = 0
            while posTad < len(bList) and chromPosition < len(cList) and bList[posTad][0] != cList[chromPosition][0][0]:
                chromPosition += 1
            posPeak = 0
        elif posTad < len(bList) and chromPosition < len(cList):
            while ((posPeak + 1) < len(cList[chromPosition])) and \
                    (int(bList[posTad][1]) > int(cList[chromPosition][posPeak][2])):
                posPeak += 1
            if int(bList[posTad][1]) < int(cList[chromPosition][posPeak][1]):
                if (posTad != 0) and (abs(float(bList[posTad - 1][4]) - float(bList[posTad][4])) < paraScore):
                    if (posTad != (len(bList) - 1)) and \
                            (abs(float(bList[posTad + 1][4]) - float(bList[posTad][4])) < paraScore):
                        removedTads.append(bList[posTad])
                        bList.remove(bList[posTad])

        posTad += 1
    return bList


def create_list_with_protein(bList, minPeak, cList=None):
    binSize = 0
    bList, binSize = create_list_of_file(bList, binSize)
    if cList is not None:
        cList = merge_protein(cList, binSize, minPeak)
        bList = compare_boundaries_protein(bList, cList)
    return bList


def main(args=None):
    args = parse_arguments().parse_args(args)
    pValue = args.value
    listOfDomains = []
    proteinList = None
    if len(args.domainFiles) == 1 and args.proteinFile is None:
        log.error('Please use multiple or domain files or at least one domain file and one protein file.')
        # raise Exception(')
        exit(1)
    if args.proteinFile is not None:
        proteinList = read_protein(args.proteinFile)

    mergedList = create_list_with_protein(args.domainFiles[0], args.minimumNumberOfPeaks, proteinList)

    if len(args.domainFiles) > 1:
        for domain in args.domainFiles[1:]:
            listOfDomains.append(create_list_with_protein(domain, args.minimumNumberOfPeaks, proteinList))
        for domain in listOfDomains:
            mergedList = merge_list(mergedList, domain, pValue)

    mergedListWithId = add_id(mergedList)
    write_in_file(mergedListWithId, args.outputMergedList)
    if len(args.domainFiles) > 1:
        relationList = create_relationsship_list(mergedListWithId, args.percent)
        write_in_file(relationList, args.outputRelationList, True)
        create_tree(relationList, mergedListWithId, args.outputTreePlotPrefix, args.outputTreePlotFormat)
