#!/usr/bin/env python
import argparse
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from matplotlib import pyplot as plt
from graphviz import Digraph
import os


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--domain1', '-d',
                                help='The domains.bed file of the first matrix is required',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--domainList', '-l',
                                help='The second and all other Domain.bed files are required',
                           nargs='+')

    parserOpt.add_argument('--ctcfFile', '-c',
                           help='In order to be able to better assess the relationship between TADs, the associated '
                                'CTCF file can be included. The CTCF-file is required in broadpeak format')
    

    parserOpt.add_argument('--minPeak', '-m',
                           help='Optional parameter to adjust the number of CTCF peaks when adapting the resolution to '
                                'the domain files',
                           type=int, default=1)


    parserOpt.add_argument('--value', '-v',
                           help='Determine a value by how much the boundaries of two TADs must at least differ to view '
                                'them as two separate TADs.',
                           type=int, default=5000)

    parserOpt.add_argument('--percent', '-p',
                           help='For the relationship determination, a percentage is required from which area coverage '
                                'the TADs are related to each other.'
			   'For example, a relationship should be entered from 5 percent area coverage -p 0.05',
                           type=float, default=0.5)

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
                num = num[0:(len(num)-1)]
            num = num[len(num)-1:]
            while zeros != 0:
                num += '0'
                zeros -= 1
            if int(num) < binSize:
                binSize = int(num)
    return splittedList, binSize


def merge_list(d1, d2, pValue= 5000):
    """Combine two lists into one, taking into account the value,
       with which one can decide whether two very similar TADs should
       be seen as one or as two """
    pos1, pos2 = 0, 0
    merged_list = []
    visited_chrom = []
    while True:
        if pos1 == len(d1):
            break
        while d1[pos1][0] == d2[pos2][0]:
            """ Checks whether the left boundary of TAD1 is less than or equal to the left boundary of TAD2 """
            if int(d1[pos1][1]) <= int(d2[pos2][1]):
                """ Checks whether either the left or right border of the two TADs exceed the minimum distance value """
                if (abs(int(d1[pos1][1]) - int(d2[pos2][1])) > pValue) or (
                        abs(int(d1[pos1][2]) - int(d2[pos2][2])) > pValue):
                    merged_list.append(d1[pos1])
                if pos1+1 != len(d1):
                    pos1 += 1
                else:
                    """ As long as List2 has not reached its end and still share the same chromosome, the remaining 
                    TADs are added """
                    while (pos2 < len(d2)) and (d1[pos1][0] == d2[pos2][0]):
                        merged_list.append(d2[pos2])
                        pos2 += 1
                    break
                """ Checks whether the left boundary of TAD2 is less than or equal to the left boundary of TAD1 """
            elif (int(d1[pos1][1]) > int(d2[pos2][1])):
                """ Checks whether either the left or right border of the two TADs exceed the minimum distance value"""
                if (abs(int(d1[pos1][1]) - int(d2[pos2][1])) > pValue) or (
                        abs(int(d1[pos1][2]) - int(d2[pos2][2])) > pValue):
                    merged_list.append(d2[pos2])
                if pos2+1 != len(d2):
                    pos2 += 1
                else:
                    """ As long as List1 has not reached its end and still share the same chromosome, the remaining 
                    TADs are added """
                    while (pos1 < len(d1))and (d1[pos1][0] == d2[pos2][0]):
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


def create_relationsship_list(domainList, pProzent = 0.5):
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
            minArea = (float(domainList[tad1][2])-float(domainList[tad1][1]))*pProzent
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
    pos = len(rList)-1
    while True:
        if rList[pos][1] == parent:
            rList[pos][2].append(child)
            return rList
        elif int(rList[pos][1][3:]) < int(parent[3:]):
            rList.append([chromosom, parent, [child]])
            return rList
        else:
            pos -= 1
        

def write_in_file(l, name, relation = False):
    """ Create and save various lists in their individual files """
    filename = name
    myfile = open(filename, 'w')
    i = 0
    if relation:
        while i < len(l):
            element = 0
            while element < len(l[i][2]):
                string = l[i][0] + '\t' + l[i][1] + '\t' + l[i][2][element] + '\n'
                myfile.write(string)
                element += 1
            i += 1
        myfile.close()
    else:
        while i < len(l):
            element = 0
            string = ""
            while element < len(l[i])-1:
                string += l[i][element] + '\t'
                element += 1
            string += l[i][element] + '\n'
            myfile.write(string)
            i += 1
        myfile.close()


def create_tree(rList, dList):
    """ Create a tree diagram from a list and
       thus represent the relationships of individual TADs per chromosome """
    name = "./trees/" + rList[0][0] + '_relations'
    g = Digraph(filename=name, strict=True)
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
            while len(sList[posSList][1])!= 0:
                g.node(sList[posSList][1][0])
                sList[posSList][1].remove(sList[posSList][1][0])
            print("Saved relation tree of " + chrom)
            if not os.path.exists("trees"):
                os.makedirs("trees")
            g.render(name)
            chrom = rList[posRList][0]
            name = "./trees/" + chrom + '_relations'
            g = Digraph(filename=name)
            posSList = 0
            while sList[posSList][0] != chrom:
                posSList += 1
        posRList += 1
    print("Saved relation tree of " + chrom)
    g.render()


def create_small_list(dList):
    """ Creates a nested list from a domain file, in which the TADs
    are sorted according to their chromosomes """
    sList = [[dList[0][0], []]]
    chrom = dList[0][0]
    pos = 0
    while pos < len(dList):
        if dList[pos][0] == chrom:
            sList[len(sList)-1][1].append(dList[pos][3])
        else:
            chrom = dList[pos][0]
            sList.append([dList[pos][0], [dList[pos][3]]])
        pos += 1
    return sList


def read_ctcf(file):
    """ Reads a CTCF file and saves it in a list """
    with open(file) as f:
        newList = [line.rstrip() for line in f]
    splittedList = [[]]
    actualChr = newList[0].split("\t")[0]
    for line in newList:
        x = line.split("\t")
        if x[0] == actualChr:
            splittedList[len(splittedList)-1].append(x[0:3])
        else:
            actualChr = x[0]
            splittedList.append([x[0:3]])
    return splittedList


def merge_ctcf(ctcfList, binSize, minPeak):
    """ Adapts the CTCF list to the corresponding BinSize, taking into
       account the minimum number of peaks per section """
    newCtcfList = []
    for chromosome in ctcfList:
        newCtcfList.append([])
        deletedCtcf = []
        currentBoundaryLeft = 0
        currentBoundaryRight = binSize
        count = 0
        for peak in chromosome:
            if int(peak[1]) <= currentBoundaryRight:
                count += 1
            else:
                if count >= minPeak:
                    newCtcfList[len(newCtcfList)-1].append([peak[0], currentBoundaryLeft, currentBoundaryRight, count])
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
    return newCtcfList


def compare_boundaries_ctcf(bList, cList, paraScore = 0.2):
    """ check the boundaries for a ctcf-peak """
    posTad = 0
    posPeak = 0
    chromPosition = 0
    removedTads = []
    while posTad < len(bList):
        if bList[posTad][0] != cList[chromPosition][0][0]:
            chromPosition = 0
            while bList[posTad][0] != cList[chromPosition][0][0]:
                chromPosition += 1
            posPeak = 0
        else:
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


def create_list_with_ctcf(bList, minPeak, cList = None):
    binSize = 0
    bList, binSize = create_list_of_file(bList, binSize)
    if cList is not None:
        cList = merge_ctcf(cList, binSize, minPeak)
        bList = compare_boundaries_ctcf(bList, cList)
    return bList


def main(args=None):
    args = parse_arguments().parse_args(args)
    pValue = args.value
    listOfDomains = []
    if args.ctcfFile is not None:
        ctcfList = read_ctcf(args.ctcfFile)
        mergedList= create_list_with_ctcf(args.domain1, args.minPeak, ctcfList)
        for domain in args.domainList:
            listOfDomains.append(create_list_with_ctcf(domain, args.minPeak, ctcfList))
    else:
        mergedList= create_list_with_ctcf(args.domain1, args.minPeak)
        for domain in args.domainList:
            listOfDomains.append(create_list_with_ctcf(domain, args.minPeak))
    for domain in listOfDomains:
        mergedList = merge_list(mergedList, domain, pValue)
    mergedListWithId = add_id(mergedList)
    write_in_file(mergedListWithId, "mergedDomains.bed")
    relationList = create_relationsship_list(mergedListWithId, args.percent)
    write_in_file(relationList, "relationList.bed", True)
    create_tree(relationList, mergedListWithId)



