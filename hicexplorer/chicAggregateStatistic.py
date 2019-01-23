import argparse
import sys
import numpy as np
import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities

from hicexplorer._version import __version__
from .lib import Viewpoint
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os

import math
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Aggregates the statistics of interaction files and prepares them for chicDifferentialTest')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for aggregation of the statistics.',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--targetFile', '-tf',
                                help='path to the target files which contains the target regions to prepare data for differential analysis.',
                                required=True,
                                nargs='+')
   
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileNameSuffix', '-o',
                            help='File name suffix to save the result.',
                            required=False,
                            default='_aggregate_target.bed')
    
    
    parserOpt.add_argument("--mergeBins", "-mb", action='store_true', help="Merge neighboring significant interactions to one. The value is averaged.")


    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def filter_scores(pScoresDictionary, pTargetRegions):

    accepted_scores = {}
    for target in pTargetRegions:
        start = int(target[1])
        end = int(target[2])

        for key in pScoresDictionary:
            if int(pScoresDictionary[key][5]) >= start and int(pScoresDictionary[key][6]) <= end:
                accepted_scores[key] = pScoresDictionary[key]
                break
    return accepted_scores


def merge_neighbors(pScoresDictionary, pMergeThreshold = 1000):

    key_list = list(pScoresDictionary.keys())

    # [[start, ..., end]]
    neighborhoods = []
    neighborhoods.append([key_list[0], key_list[0]])
    scores = [pScoresDictionary[key_list[0]]]
    
    for key in key_list[1:]:
        
        if np.absolute(key - neighborhoods[-1][0]) <= pMergeThreshold or np.absolute(key - neighborhoods[-1][1]) <= pMergeThreshold:
            neighborhoods[-1][-1] = key
            scores[-1] += pScoresDictionary[key]
        else:
            neighborhoods.append([key, key])
            scores.append(pScoresDictionary[key])

    for i in range(len(neighborhoods)):
        scores[i] /= len(neighborhoods[i])

    return neighborhoods, scores

def write (pOutFileName, pNeighborhoods, pScores, pInteractionLines, pThreshold):

    with open(pOutFileName, 'w') as file:
        file.write('# Significant regions with rbz-score higher as ' + str(pThreshold) + '\n')
        file.write('#ChrViewpoint\tStart\tEnd\tGene\tChrInteraction\tStart\tEnd\tRel Inter viewpoint\trbz-score viewpoint\tRaw viewpoint\tRel Inter target\trbz-score target\tRaw target')
        file.write('\n')
        for i in range(len(pNeighborhoods)):
            start = pNeighborhoods[i][0]
            end = pNeighborhoods[i][1]
            pInteractionLines[start]
            pInteractionLines[end][-3:]
            new_end = pInteractionLines[end][6]

            pInteractionLines[0][-3:]

            new_line = '\t'.join(pInteractionLines[start][:6])
            new_line += '\t' + new_end
            new_line += '\t' + '\t'.join(format(float(x), "10.5f") for x in pInteractionLines[0][8:])

            new_line += '\t' + format(pScores[i][0], '10.5f') + '\t' + format(pScores[i][1], '10.5f') + '\t' + format(pScores[i][2], '10.5f')

            new_line += '\n'
            file.write(new_line)


def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    background_data = None
    relative_interaction = False
    rbz_score = True
    # read all interaction files.
    for interactionFile, targetFile in zip(args.interactionFile, args.targetFile):
        header, interaction_data, interaction_file_data = viewpointObj.readInteractionFileForAggregateStatistics(interactionFile)
        log.debug('header {}'.format(header))
        # log.debug('interaction_data {}'.format(interaction_data))
        log.debug('interaction_file_data {}'.format(list(interaction_file_data.items())[1990:2010]))

        target_regions = utilities.readBed(targetFile)
        log.debug('target_regions {}'.format(target_regions))
        accepted_scores = filter_scores(interaction_file_data, target_regions)

        if len(accepted_scores) == 0:
            log.error('No target regions found')
            sys.exit(0)
        log.debug('accepted_scores {}'.format(accepted_scores))
        merged_neighborhood = merge_neighbors(accepted_scores)
        outFileName = interactionFile.split('.')[0] + '_' + args.outFileNameSuffix
        write(outFileName, merged_neighborhood[0], merged_neighborhood[1], interaction_file_data, args.acceptThreshold)


    