import numpy as np

import logging
log = logging.getLogger(__name__)
import copy

import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import nbinom
from scipy.special import gammaln
from scipy import special
import h5py

from hicexplorer.lib import cnb


class Viewpoint():

    def __init__(self, pHiCMatrix=None):
        self.hicMatrix = pHiCMatrix

    def readReferencePointFile(self, pBedFile, pGene=True):
        '''
        This function reads a text file which contains reference points. Reference points are

        Reference points need to be tab seperated
        and contain the chromosome, the reference start point and/or the reference end point, and the gene:

        chr start gene

        or

        chr start end gene

        Per line one reference point.

        Returns a list of tuples with (chr, start, end), if only a start index is given, end = start and a list of genes.

        '''
        viewpoints = []
        gene_list = []
        with open(pBedFile, 'r') as file:
            for line in file.readlines():
                _line = line.strip().split('\t')
                if len(line) == 0 or len(_line) == 0 or line == '':
                    continue
                if len(_line) == 3:
                    chrom, start, end = _line[0], _line[1], _line[1]
                    if pGene:
                        gene_list.append(_line[2])
                elif len(_line) > 3:
                    chrom, start, end = _line[:3]
                    if pGene:
                        gene_list.append(_line[3])
                else:
                    continue

                viewpoints.append((chrom, start, end))
        if pGene:
            return viewpoints, gene_list
        else:
            return viewpoints

    # def readInteractionFile(self, pBedFile):
    #     '''
    #     Reads an interaction file produced by chicViewpoint. Contains header information, these lines
    #     start with '#'.
    #     Interactions files contain:
    #     Chromosome Viewpoint, Start, End, Gene, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
    #     Relative number of interactions, p-score based on relative interactions.

    #     This function returns:
    #     - header as  a string
    #     - interaction data in relation to relative position as a dict e.g. {-1000:0.1, -1500:0.2}
    #     - p-score in relation to relative position as a dict (same format as interaction data)
    #     - interaction_file_data: the raw line in relation to the relative position. Needed for additional output file.
    #     '''
    #     # use header info to store reference point, and based matrix
    #     interaction_data = {}
    #     p_score = {}
    #     interaction_file_data = {}
    #     genomic_coordinates = {}
    #     with open(pBedFile) as fh:
    #         fh.readline()
    #         header = fh.readline()
    #         for line in fh.readlines():
    #             # Addition header information for end users
    #             if line.strip().startswith('#'):
    #                 continue

    #             _line = line.strip().split('\t')
    #             # relative postion and relative interactions
    #             interaction_data[int(_line[-5])] = float(_line[-4])
    #             p_score[int(_line[-5])] = float(_line[-3])
    #             interaction_file_data[int(_line[-5])] = _line
    #             genomic_coordinates[int(_line[-5])] = [_line[0], _line[1], _line[2]]

    #     return header, interaction_data, p_score, interaction_file_data, genomic_coordinates

    def readInteractionFile(self, pFilePath, pInternalIdentifierTriplet):
        '''
        Reads an interaction file produced by chicViewpoint. Contains header information, these lines
        start with '#'.
        Interactions files contain:
        Chromosome Viewpoint, Start, End, Gene, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-value based on negative binomial distribution for the relative position.

        This function returns:
        - header as  a string
        - interaction data in relation to relative position as a dict together with p-value, raw and x-fold: e.g. {-1000:[0.1, 0.01, 2.3, 5]}
        - interaction_file_data: the raw line in relation to the relative position. Needed for additional output file.
        '''

        # interaction_data[int(_line[-5])] = np.array([float(_line[-4]), float(_line[-3]), float(_line[-1]), float(_line[-2])])
        # interaction_file_data[int(_line[-5])] = _line

        interactionFileHDF5Object = h5py.File(pFilePath, 'r')
        arrays_to_retrieve = ['relative_position_list', 'interaction_data_list', 'pvalue', 'raw', 'xfold']
        internal_path = '/'.join(pInternalIdentifierTriplet)
        data = []

        if internal_path not in interactionFileHDF5Object:
            log.debug('internal path not found: {}'.format(internal_path))
            return {}, {}, []
        if list(interactionFileHDF5Object[internal_path].keys()) == 0:
            log.debug('interactionFileHDF5Object[internal_path] {}'.format(interactionFileHDF5Object[internal_path]))
            return {}, {}, []

        for array_name in arrays_to_retrieve:
            try:
                if internal_path + '/' + array_name in interactionFileHDF5Object:
                    data.append(np.array(interactionFileHDF5Object[internal_path + '/' + array_name][:]))
                # log.debug('datatype of {} {} '.format(array_name, data[-1].dtype))
                # data.append(interactionFileHDF5Object.get( internal_path + '/' + array_name).value)
            except Exception as exp:
                log.debug(internal_path + '/' + array_name)
                log.debug(str(exp))

        for array_name in ['start_list', 'end_list']:
            try:
                if internal_path + '/' + array_name in interactionFileHDF5Object:
                    data.append(np.array(interactionFileHDF5Object[internal_path + '/' + array_name][:]))
                # data.append(interactionFileHDF5Object.get( internal_path + '/' + array_name).value)
            except Exception as exp:
                log.debug(internal_path + '/' + array_name)
                log.debug(str(exp))

        if internal_path + '/' + 'chromosome' in interactionFileHDF5Object:
            chromosome = interactionFileHDF5Object.get(internal_path + '/' + 'chromosome')[()].decode("utf-8")
        else:
            log.debug('internal path not found: {}'.format(internal_path + '/' + 'chromosome'))

        if internal_path + '/' + 'gene' in interactionFileHDF5Object:
            gene = interactionFileHDF5Object.get(internal_path + '/' + 'gene')[()].decode("utf-8")
        else:
            log.debug('internal path not found: {}'.format(internal_path + '/' + 'gene'))

        if internal_path + '/' + 'sum_of_interactions' in interactionFileHDF5Object:
            sum_of_interactions = interactionFileHDF5Object.get(internal_path + '/' + 'sum_of_interactions')[()]
        else:
            log.debug('internal path not found: {}'.format(internal_path + '/' + 'sum_of_interactions'))

        reference_point_start = None
        reference_point_end = None
        if 'reference_point_start' in interactionFileHDF5Object[internal_path]:
            reference_point_start = interactionFileHDF5Object.get(internal_path + '/' + 'reference_point_start')[()]
        if 'reference_point_end' in interactionFileHDF5Object[internal_path]:
            reference_point_end = interactionFileHDF5Object.get(internal_path + '/' + 'reference_point_end')[()]

        interactionFileHDF5Object.close()
        # log.debug('chromosomes: {}'.format(chromosome))

        # log.debug(data)
        interaction_data = {}
        interaction_file_data = {}

        try:
            for i in range(len(data[0])):
                interaction_data[data[0][i]] = list([float(data[1][i]), float(data[2][i]), float(data[3][i]), float(data[4][i])])
                interaction_file_data[data[0][i]] = list([str(chromosome), int(float(data[5][i])), int(float(data[6][i])), str(gene), float(sum_of_interactions), int(float(data[0][i])), float(data[1][i]), float(data[2][i]), float(data[4][i]), float(data[3][i])])
        except Exception:
            return {}, {}, []
        # log.debug('153')
        # log.debug('interaction_data {}'.format(interaction_data[0]))
        # log.debug('interaction_file_data {}'.format(interaction_file_data[0]))

        return interaction_data, interaction_file_data, [reference_point_start, reference_point_end]

        # -5: relative position relative_position_list
        # -4: relative interaction: interaction_data_list
        # -3: p-value: pvalue
        # -1: raw: raw
        # -2: x-fold: xfold

        # # use header info to store reference point, and based matrix
        # interaction_data = {}
        # interaction_file_data = {}
        # with open(pBedFile) as fh:
        #     fh.readline()
        #     header = fh.readline()
        #     fh.readline()

        #     for line in fh.readlines():
        #         # Addition header information for end users
        #         if line.strip().startswith('#'):
        #             continue
        #         if not line:
        #             continue
        #         _line = line.strip().split('\t')
        #         # relative postion and relative interactions

        #         interaction_data[int(
        #             _line[-5])] = np.array([float(_line[-4]), float(_line[-3]), float(_line[-1]), float(_line[-2])])
        #         interaction_file_data[int(_line[-5])] = _line
        # return header, interaction_data, interaction_file_data

    def interactionDataForPlot(self, pFilePath, pInternalIdentifierTriplet):

        # header, interaction_data, p_value_data, #_interaction_file_data_raw#, genomic_coordinates = self.readInteractionFile(pInteractionFile)
        # matrix_name, viewpoint, upstream_range, downstream_range, gene, _ = header.split('\t')
        interaction_data, interaction_file_data, reference_point = self.readInteractionFile(pFilePath, pInternalIdentifierTriplet)

        #             _line = line.strip().split('\t')
        #             # relative postion and relative interactions
        #             interaction_data[int(_line[-5])] = float(_line[-4])
        #             p_score[int(_line[-5])] = float(_line[-3])
        #             interaction_file_data[int(_line[-5])] = _line
        #             genomic_coordinates[int(_line[-5])] = [_line[0], _line[1], _line[2]]

        # arrays_to_retrieve = ['relative_position_list', 'interaction_data_list', 'pvalue', 'raw', 'xfold']
        interaction_data = {}
        p_score = {}
        # interaction_file_data = {}
        genomic_coordinates = {}

        for key, value in interaction_file_data.items():
            interaction_data[key] = value[6]
            p_score[key] = value[7]
            genomic_coordinates[key] = value[:3]

        return interaction_data, p_score, genomic_coordinates, reference_point

    def readBackgroundDataFile(self, pBedFile, pRange, pFixateRange, pMean=False):
        '''
        Reads a background data file, containing per line a tab delimited content:
        Relative position to viewpoint, relative interaction count to total number of interactions of all viewpoints over all samples, SEM value of this data point.

        '''
        distance = {}

        with open(pBedFile) as fh:
            _ = fh.readline()
            for line in fh.readlines():
                _line = line.split('\t')
                if not pMean:
                    distance[int(_line[0])] = [float(_line[1]), float(_line[2]), float(_line[3])]
                else:
                    distance[int(_line[0])] = [float(_line[-1])]

        max_key = max(distance)
        min_key = min(distance)

        # check if max is really pFixateRange or it needs to be changed
        if max_key > pFixateRange:
            max_key = pFixateRange

        if min_key < -pFixateRange:
            min_key = -pFixateRange

        keys = list(distance.keys())
        inc = np.absolute(np.absolute(keys[0]) - np.absolute(keys[1]))
        if max_key < pRange[1]:
            i = max_key
            while i < pRange[1]:
                i += inc
                distance[i] = distance[max_key]

        if min_key > -pRange[0]:
            i = min_key
            while i > -pRange[0]:
                i -= inc
                distance[i] = distance[min_key]
        return distance

    def writeInteractionFile(self, pBedFile, pData, pHeader, pPValueData, pXfold, pDecimalPlaces=12):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-values based on negative binomial distribution per relative distance, raw interaction data
        '''

        with open((pBedFile + '.txt').strip(), 'w') as fh:
            fh.write('{}\n'.format(pHeader))
            for j, interaction in enumerate(pData):
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.{decimal_places}f}\t{:.{decimal_places}f}\t{:.{decimal_places}f}\t{:.{decimal_places}f}\n".
                         format(interaction[0], interaction[1],
                                interaction[2], interaction[3], interaction[4], interaction[5], interaction[6], pPValueData[j], pXfold[j], interaction[7], decimal_places=pDecimalPlaces))
        return

    def writeInteractionFileHDF5(self, pInteractionFileGroupH5Object, pFileName, pData, pReferencePoint):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-values based on negative binomial distribution per relative distance, raw interaction data
        '''

        # chrom, start_list, end_list, gene, sum_of_interactions, relative_position_list, interaction_data_list, pPValueList, xFoldList, raw_data_list
        # log.debug(pFileName)

        groupObject, file_name = self.createUniqueHDFGroup(pInteractionFileGroupH5Object, pFileName)
        # if counter != 0:
        #     pass
        # log.debug('Gene name {} occurred {} times! Stored as {}_{}'.format(pFileName, counter, pFileName, counter))
        # groupObject.create_dataset("header", data=pHeader)
        # groupObject.create_dataset("chromosome", data=pData[0].decode("utf-8"))
        groupObject["chromosome"] = str(pData[0])

        groupObject.create_dataset("start_list", data=pData[1], compression="gzip", compression_opts=9)
        groupObject.create_dataset("end_list", data=pData[2], compression="gzip", compression_opts=9)
        groupObject.create_dataset("gene", data=pData[3])
        groupObject.create_dataset("sum_of_interactions", data=pData[4])
        groupObject.create_dataset("relative_position_list", data=pData[5], compression="gzip", compression_opts=9)
        groupObject.create_dataset("interaction_data_list", data=pData[6], compression="gzip", compression_opts=9)
        groupObject.create_dataset("pvalue", data=pData[7], compression="gzip", compression_opts=9)
        groupObject.create_dataset("xfold", data=pData[8], compression="gzip", compression_opts=9)
        groupObject.create_dataset("raw", data=pData[9], compression="gzip", compression_opts=9)

        # groupObject.create_dataset("start_list", data=pData[1])
        # groupObject.create_dataset("end_list", data=pData[2])
        # groupObject.create_dataset("gene", data=pData[3])
        # groupObject.create_dataset("sum_of_interactions", data=pData[4])
        # groupObject.create_dataset("relative_position_list", data=pData[5])
        # groupObject.create_dataset("interaction_data_list", data=pData[6])
        # groupObject.create_dataset("pvalue", data=pData[7])
        # groupObject.create_dataset("xfold", data=pData[8])
        # groupObject.create_dataset("raw", data=pData[9])

        groupObject.create_dataset("reference_point_start", data=int(pReferencePoint[0]))
        groupObject.create_dataset("reference_point_end", data=int(pReferencePoint[1]))

        # with open((pBedFile + '.txt').strip(), 'w') as fh:
        #     fh.write('{}\n'.format(pHeader))
        #     for j, interaction in enumerate(pData):
        #         fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.{decimal_places}f}\t{:.{decimal_places}f}\t{:.{decimal_places}f}\t{:.{decimal_places}f}\n".
        #                  format(interaction[0], interaction[1],
        #                         interaction[2], interaction[3], interaction[4], interaction[5], interaction[6], pPValueData[j], pXfold[j], interaction[7], decimal_places=pDecimalPlaces))
        return file_name

    def createUniqueHDFGroup(self, pGroupObject, pAdditionalGroupName):

        success = False
        counter = 0
        # log.debug('createUniqueHDFGroup')
        additional_group_name_raw = pAdditionalGroupName
        while not success:
            try:
                if counter != 0:
                    pAdditionalGroupName = additional_group_name_raw + '_' + str(counter)
                    file_name = pAdditionalGroupName
                else:
                    file_name = pAdditionalGroupName
                groupObject = pGroupObject.create_group(file_name)
                success = True
            except ValueError:
                counter += 1
                log.debug('counter +1: {}'.format(counter))
        # log.debug('createUniqueHDFGroup EXIT')

        return groupObject, file_name

    def computeViewpoint(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        This function computes a viewpoint for a given sample and a given pReferencePoint within the
        range of pRegion_start and pRegion_end.

        All interactions with the reference point of one relative distance to it are summed up,
        if the reference point is larger than one bin of the Hi-C matrix, it is considered as one bin and the values are summed together.
        '''
        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(
            pReferencePoint)

        view_point_range = self.getViewpointRangeAsMatrixIndices(
            pChromViewpoint, pRegion_start, pRegion_end)
        view_point_range = list(view_point_range)
        view_point_range[1] += 1

        elements_of_viewpoint = (view_point_range[1] - view_point_range[0])
        data_list = np.zeros(elements_of_viewpoint)
        _view_point_start = view_point_start

        while _view_point_start <= view_point_end:
            chrom, start, end, _ = self.hicMatrix.getBinPos(_view_point_start)
            for j, idx in zip(range(elements_of_viewpoint), range(view_point_range[0], view_point_range[1], 1)):
                data_list[j] += self.hicMatrix.matrix[_view_point_start, idx]

            _view_point_start += 1

        elements_of_viewpoint = elements_of_viewpoint - \
            (view_point_end - view_point_start)
        data_list_new = np.zeros(elements_of_viewpoint)

        index_before_viewpoint = view_point_start - view_point_range[0]

        # elements before the viewpoint
        data_list_new[0:index_before_viewpoint] = data_list[0:index_before_viewpoint]

        # summation because the viewpoint can not be only one bin but can contain multiple
        data_list_new[index_before_viewpoint] = np.sum(
            data_list[index_before_viewpoint: index_before_viewpoint + view_point_end - view_point_start + 1])

        # elements after the viewpoint
        data_list_new[index_before_viewpoint + 1:] = data_list[index_before_viewpoint +
                                                               view_point_end - view_point_start + 1:]
        return data_list_new, index_before_viewpoint

    def createInteractionFileData(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene, pSumOfInteractions):
        '''
        Creates out of internal information a list of tuples which can be written to an interaction file.
        Tuple contains:
        Chromosome viewpoint, start, end, chromosome interaction, start, end, relative_position, interaction data
        '''
        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(
            pReferencePoint)
        view_point_range = self.getViewpointRangeAsMatrixIndices(
            pChromViewpoint, pRegion_start, pRegion_end)
        view_point_range = list(view_point_range)
        view_point_range[1] += 1
        # elements_of_viewpoint = view_point_range[1] - view_point_range[0] - (
        #     view_point_end - view_point_start) + 1

        interactions_list = []
        chrom, start, _, _ = self.hicMatrix.getBinPos(view_point_start)
        _, _, end, _ = self.hicMatrix.getBinPos(view_point_end)
        interaction_positions = list(
            range(view_point_range[0], view_point_start, 1))
        interaction_positions.extend([view_point_start])
        interaction_positions.extend(
            list(range(view_point_end + 1, view_point_range[1], 1)))
        relative_position = -1
        for j, idx in zip(range(len(pInteractionData)), interaction_positions):
            try:

                chrom_second, start_second, end_second, _ = self.hicMatrix.getBinPos(
                    idx)
                if relative_position < 0:
                    relative_position = int(start_second) - int(start)
                else:
                    relative_position = int(end_second) - int(end)

                interactions_list.append((chrom_second, start_second, end_second, pGene, str(
                    pSumOfInteractions), relative_position, float(pInteractionData[j]), float(pInteractionDataRaw[j])))
            except Exception:
                log.error('Failed to get bin position of index {}'.format(idx))
                exit(1)
        return interactions_list

    def createInteractionFileDataHDF5(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene, pSumOfInteractions, pPValueList, xFoldList):
        '''
        Creates out of internal information a list of tuples which can be written to an interaction file.
        Tuple contains:
        Chromosome viewpoint, start, end, chromosome interaction, start, end, relative_position, interaction data
        '''
        chrom = ''
        start_list = []
        end_list = []
        gene = pGene
        # sum_of_interactions_list = pSumOfInteractions
        relative_position_list = []
        interaction_data_list = []
        raw_data_list = []

        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(
            pReferencePoint)
        view_point_range = self.getViewpointRangeAsMatrixIndices(
            pChromViewpoint, pRegion_start, pRegion_end)
        view_point_range = list(view_point_range)
        view_point_range[1] += 1

        chrom, start, _, _ = self.hicMatrix.getBinPos(view_point_start)
        _, _, end, _ = self.hicMatrix.getBinPos(view_point_end)
        interaction_positions = list(
            range(view_point_range[0], view_point_start, 1))
        interaction_positions.extend([view_point_start])
        interaction_positions.extend(
            list(range(view_point_end + 1, view_point_range[1], 1)))
        relative_position = -1
        for j, idx in zip(range(len(pInteractionData)), interaction_positions):
            try:

                chrom_second, start_second, end_second, _ = self.hicMatrix.getBinPos(
                    idx)
                if relative_position < 0:
                    relative_position = int(start_second) - int(start)
                else:
                    relative_position = int(end_second) - int(end)

                chrom = chrom_second
                start_list.append(start_second)
                end_list.append(end_second)
                relative_position_list.append(relative_position)
                interaction_data_list.append(float(pInteractionData[j]))
                raw_data_list.append(float(pInteractionDataRaw[j]))

            except Exception:
                log.error('Failed to get bin position of index {}'.format(idx))
                exit(1)
        return [chrom, start_list, end_list, gene, pSumOfInteractions, relative_position_list, interaction_data_list, pPValueList, xFoldList, raw_data_list]

    def getViewpointRangeAsMatrixIndices(self, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        Returns the matrix indices of a chromosome and a specific position.
        '''

        _range = self.hicMatrix.getRegionBinRange(
            pChromViewpoint, pRegion_start, pRegion_end)

        return _range

    def getReferencePointAsMatrixIndices(self, pReferencePoint):
        '''
        Returns the correct matrix indices of a given reference point.
        '''
        if len(pReferencePoint) == 2:
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(
                pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[1]))
        elif len(pReferencePoint) == 3:
            # log.debug('pReferencePoint: {}'.format(pReferencePoint))
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(
                pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
        else:
            log.error("No valid reference point given. {}".format(
                pReferencePoint))
            exit(1)
        # log.debug('view_point_start: {} view_point_end {}'.format(view_point_start, view_point_end))

        return view_point_start, view_point_end

    def smoothInteractionValues(self, pData, pWindowSize):
        '''
        Adds -pWindowsSize/2 and +pWindowsSize/2 around pData[i] and averages pData[i] by pWindowSize to
        smooth the interaction values.
        '''
        window_size = np.int(np.floor(pWindowSize / 2))
        window_size_upstream = window_size
        if pWindowSize % 2 == 0:
            window_size_upstream -= 1

        average_contacts = np.zeros(len(pData))

        # add upstream and downstream, handle regular case
        for i in range(window_size_upstream, len(pData) - window_size):
            start = i - window_size_upstream
            end = i + window_size + 1
            average_contacts[i] = np.mean(pData[start:end])

        # handle border conditions
        for i in range(window_size):
            start = i - window_size_upstream
            if start < 0:
                start = 0
            end = i + window_size + 1

            average_contacts[i] = np.mean(pData[start:end])
            average_contacts[-(i + 1)] = np.mean(pData[-end:])

        # average_contacts = average_contacts
        return average_contacts

    def computeRelativeValues(self, pData, pDenominator=None):
        '''
        Computes the relative values of pData by adding all data points together and divding all data points by the result.
        pData[i] = pData[i] / sum(pData)
        '''
        pOutput = np.array(pData, dtype=float)
        if pDenominator:
            pOutput /= pDenominator
        else:
            sumValue = np.sum(pOutput)
            pOutput /= sumValue
        return pOutput

    def calculateViewpointRange(self, pViewpoint, pRange):
        '''
        This function computes the correct start and end position of a viewpoint given the viewpoint and the range.
        '''

        chr_bin_range = self.hicMatrix.getChrBinRange(pViewpoint[0])
        if chr_bin_range is None:
            return None, None, None
        max_length = self.hicMatrix.getBinPos(
            chr_bin_range[1] - 1)[2]
        # log.debug('max_length {}'.format(max_length))
        bin_size = self.hicMatrix.getBinSize()
        _range = [pRange[0], pRange[1]]
        region_start = int(pViewpoint[1]) - pRange[0]
        if region_start < 0:
            region_start = 0
            _range[0] = int(pViewpoint[1])

        region_end = int(pViewpoint[2]) + pRange[1]
        if region_end > max_length:
            # -1 is important, otherwise self.hicMatrix.getRegionBinRange will crash
            region_end = max_length - 1
            _range[1] = (max_length - int(pViewpoint[2])) + bin_size
        return region_start, region_end, _range

    def getDataForPlotting(self, pFilePath, pInteractionFile, pRange, pBackgroundModel, pResolution):
        # header, interaction_data, p_value_data, _interaction_file_data_raw, genomic_coordinates = self.readInteractionFile(pInteractionFile)
        # matrix_name, viewpoint, upstream_range, downstream_range, gene, _ = header.split('\t')

        # data = pViewpointObj.readInteractionFile(pFilePath, sample)
        interaction_data, p_value_data, genomic_coordinates, viewpoint = self.interactionDataForPlot(pFilePath, pInteractionFile)

        data = []
        p_value = []
        data_background = None
        viewpoint_index_end = None
        if pRange:

            interaction_data_keys = copy.deepcopy(
                list(interaction_data.keys()))
            for key in interaction_data_keys:
                if key >= -pRange[0] and key <= pRange[1]:
                    continue
                interaction_data.pop(key, None)
            if pBackgroundModel:
                background_data_keys_sorted = sorted(pBackgroundModel)

                for key in background_data_keys_sorted:
                    if key >= -pRange[0] and key <= pRange[1]:
                        continue
                    pBackgroundModel.pop(key, None)
                background_data_keys_sorted = sorted(pBackgroundModel)

        if pBackgroundModel:
            # viewpoint_index = background_data_keys_sorted.index(0)
            viewpoint_index_start = background_data_keys_sorted.index(0)

            data_background = []

            for key in background_data_keys_sorted:
                # log.debug('key {}'.format(key))
                if key in interaction_data:
                    if key == 0:
                        chromosome, start, end = genomic_coordinates[key]
                        # log.debug('peak width {}'.format(np.abs(int(start) - int(end))))
                        # log.debug('pResolution {}'.format(pResolution))

                        if np.abs(int(start) - int(end)) > pResolution:

                            peak_width = np.abs(int(start) - int(end)) // pResolution
                            viewpoint_index_end = peak_width
                            i = 0
                            while i < peak_width:
                                data.append(interaction_data[key])
                                i += 1
                    else:
                        data.append(interaction_data[key])

                    if key in p_value_data:
                        if key == 0:
                            chromosome, start, end = genomic_coordinates[key]
                            if np.abs(int(start) - int(end)) > pResolution:
                                peak_width = np.abs(int(start) - int(end)) // pResolution
                                i = 0
                                while i < peak_width:
                                    p_value.append(p_value_data[key])
                                    i += 1
                                # log.debug('peak width: {}'.format(peak_width))
                        else:
                            p_value.append(p_value_data[key])
                    if key == 0:
                        chromosome, start, end = genomic_coordinates[key]
                        if np.abs(int(start) - int(end)) > pResolution:
                            peak_width = np.abs(int(start) - int(end)) // pResolution
                            i = 0
                            while i < peak_width:
                                data_background.append(pBackgroundModel[key][0])
                                log.debug('pBackgroundModel[key][0] {}'.format(pBackgroundModel[key][0]))
                                log.debug('peak_width {}'.format(peak_width))
                                log.debug('data_background[-1] {}'.format(data_background[-1]))

                                i += 1
                            # log.debug('peak width: {}'.format(peak_width))
                    else:
                        data_background.append(pBackgroundModel[key][0])

            if viewpoint_index_end is None:
                viewpoint_index_end = viewpoint_index_start
            else:
                viewpoint_index_end += viewpoint_index_start

        else:
            data = []
            interaction_key = sorted(interaction_data)
            for key in interaction_key:
                # log.debug('key {}'.format(key))
                if key == 0:
                    chromosome, start, end = genomic_coordinates[key]
                    # log.debug('peak width {}'.format(np.abs(int(start) - int(end))))
                    # log.debug('pResolution {}'.format(pResolution))

                    if np.abs(int(start) - int(end)) > pResolution:

                        peak_width = np.abs(int(start) - int(end)) // pResolution
                        viewpoint_index_end = peak_width
                        i = 0
                        while i < peak_width:
                            data.append(interaction_data[key])
                            i += 1
                else:
                    data.append(interaction_data[key])
                if key in p_value_data:
                    if key == 0:
                        chromosome, start, end = genomic_coordinates[key]
                        if np.abs(int(start) - int(end)) > pResolution:
                            peak_width = np.abs(int(start) - int(end)) // pResolution
                            i = 0
                            while i < peak_width:
                                p_value.append(p_value_data[key])
                                i += 1
                            # log.debug('peak width: {}'.format(peak_width))
                    else:
                        p_value.append(p_value_data[key])

            viewpoint_index_start = interaction_key.index(0)
            if viewpoint_index_end is None:
                viewpoint_index_end = viewpoint_index_start
            else:
                viewpoint_index_end += viewpoint_index_start
        return data, data_background, p_value, viewpoint_index_start, viewpoint_index_end, viewpoint

    def plotViewpoint(self, pAxis, pData, pColor, pLabelName, pHighlightRegion=None, pHighlightSignificantRegion=None):
        data_plot_label = pAxis.plot(
            range(len(pData)), pData, '-' + pColor, alpha=0.9, label=pLabelName, linewidth=1)
        if pHighlightRegion:
            for region in pHighlightRegion:
                pAxis.axvspan(region[0], region[1], color='red', alpha=0.3)
        if pHighlightSignificantRegion:
            for region in pHighlightSignificantRegion:
                pAxis.axvspan(region[0], region[1], color=pColor, alpha=0.3)
        return data_plot_label

    def plotBackgroundModel(self, pAxis, pBackgroundData, pXFold=None):
        pBackgroundData = np.array(pBackgroundData)
        data_plot_label = pAxis.plot(range(len(pBackgroundData)), pBackgroundData, '-r', alpha=0.5, label='background model', linewidth=1)
        if pXFold:
            upper_values = pBackgroundData * pXFold
            lower_values = pBackgroundData
            pAxis.fill_between(range(len(pBackgroundData)), upper_values, lower_values, facecolor='r', alpha=0.5)
        return data_plot_label

    def plotPValue(self, pAxis, pAxisLabel, pPValueData, pLabelText, pCmap, pFigure, pValueSignificanceLevels):

        _z_score = np.empty([2, len(pPValueData)])
        _z_score[:, :] = pPValueData
        pAxis.xaxis.set_visible(False)
        pAxis.yaxis.set_visible(False)
        divider = make_axes_locatable(pAxisLabel)
        cax = divider.append_axes("left", size="20%", pad=0.09)

        if pPValueData is not None:
            # log.debug('pValueSignificanceLevels {}'.format(pValueSignificanceLevels))
            img = pAxis.contourf(_z_score, cmap=pCmap)
            colorbar = pFigure.colorbar(
                img, cax=cax, ticks=[min(pPValueData), max(pPValueData)])
            colorbar.ax.set_ylabel('p-value', size=6)

        elif pValueSignificanceLevels:
            pValueSignificanceLevels.insert(0, -1)
            pValueSignificanceLevels.append(1)

            img = pAxis.contourf(_z_score, levels=pValueSignificanceLevels, colors=['#CC0000', '#FFD43B', '#306998', '#FFFFFF'])
            colorbar = pFigure.colorbar(img, cax=cax, ticks=[pValueSignificanceLevels[1], pValueSignificanceLevels[2], pValueSignificanceLevels[3]])
            colorbar.ax.tick_params(labelsize=6)
            colorbar.ax.set_ylabel('p-value', size=6)

        pAxisLabel.text(0.45, 0, pLabelText, size=7)
        pAxisLabel.xaxis.set_visible(False)
        pAxisLabel.yaxis.set_visible(False)
        pAxisLabel.set_frame_on(False)

    def writePlotData(self, pInteractionFileDataRaw, pFileName, pBackgroundModel):
        interaction_file_data_raw_sorted = sorted(pInteractionFileDataRaw)
        with open(pFileName + '.txt', 'w') as output_file:
            output_file.write(
                '#ChrViewpoint\tStart\tEnd\tGene\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\trbz-score\tRaw')

            if pBackgroundModel:
                output_file.write('\tbackground model\tbackground model SEM\n')
            else:
                output_file.write('\n')
            for key in interaction_file_data_raw_sorted:

                _array = pInteractionFileDataRaw[key]

                output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(_array[0], _array[1], _array[2],
                                                                                      _array[3], _array[4], _array[5],
                                                                                      _array[6], _array[7], _array[8], _array[9], _array[10]))
                if pBackgroundModel:
                    output_file.write(
                        '\t{}\t{}\n'.format(_array[11], _array[12]))
                else:
                    output_file.write('\n')

    def interactionBackgroundData(self, pBackground, pRange):

        background_model = []
        background_data_keys_sorted = sorted(pBackground)
        for key in background_data_keys_sorted:
            if key >= -pRange[0] and key <= pRange[1]:
                background_model.append(pBackground[key])
        return np.array(background_model)

    def rbz_score(self, pRelativeInteractions, pBackgroundModel, pBackgroundModelSEM):
        _rbz_score = np.empty(len(pRelativeInteractions))
        if len(pRelativeInteractions) != len(pBackgroundModel) or \
                len(pRelativeInteractions) != len(pBackgroundModelSEM):
            log.info('Computing of rbz-score failed, data is having different size. ' +
                     '\nrelative interactions {} background model {} background model SEM {}'.format(len(pRelativeInteractions),
                                                                                                     len(pBackgroundModel), len(pBackgroundModelSEM)))
            return None
        _rbz_score = pRelativeInteractions - pBackgroundModel
        _rbz_score /= pBackgroundModelSEM

        mask = np.isnan(_rbz_score)
        _rbz_score[mask] = -1
        mask = np.isinf(_rbz_score)
        _rbz_score[mask] = -1
        return _rbz_score

    def readRejectedFile(self, pFilePath, pDifferentialHighlightTriplet, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint):
        # list of start and end point of regions to highlight

        if pDifferentialHighlightTriplet is None or len(pDifferentialHighlightTriplet) == 0:
            return None

        differentialFileHDF5Object = h5py.File(pFilePath, 'r')

        # arrays_to_retrieve = ['relative_position_list', 'interaction_data_list', 'pvalue', 'raw', 'xfold']
        internal_path = '/'.join(pDifferentialHighlightTriplet)

        if internal_path not in differentialFileHDF5Object:
            log.debug('internal path not found: {}'.format(internal_path))
            return None
        if list(differentialFileHDF5Object[internal_path].keys()) == 0:
            log.debug('differentialFileHDF5Object[internal_path] {}'.format(differentialFileHDF5Object[internal_path]))
            return None

        try:
            start_list = np.array(differentialFileHDF5Object[internal_path + '/' + 'start_list'][:])
            end_list = np.array(differentialFileHDF5Object[internal_path + '/' + 'end_list'][:])
            relative_distance_list = np.array(differentialFileHDF5Object[internal_path + '/' + 'relative_distance_list'][:])
        except Exception:
            return None
        if len(start_list) == 0:
            return None

        highlight_areas_list = []
        # for bed_file in pDifferentialHighlightFiles:
        reference_point_start, reference_point_end = pViewpoint

        for start, end, relative_distance in zip(start_list, end_list, relative_distance_list):

            if int(relative_distance) >= -pRange[0] and int(relative_distance) <= pRange[1]:

                width = (end - start) / pResolution
                if int(relative_distance) < 0:
                    relative_position_genomic_coordinates = start - int(reference_point_start)
                    viewpointIndex = pViewpointIndexStart
                else:
                    relative_position_genomic_coordinates = start - int(reference_point_end)
                    viewpointIndex = pViewpointIndexEnd

                # relative_position_genomic_coordinates  = reference_point_position
                relative_position = viewpointIndex + \
                    (relative_position_genomic_coordinates / pResolution)
                highlight_areas_list.append(
                    [relative_position, relative_position + width])
        if len(highlight_areas_list) == 0:
            return None
        log.debug('highlight_areas_list {}'.format(highlight_areas_list))
        return highlight_areas_list

    def pvalues(self, pBackgroundModel, pDataList, pIndexReferencePoint):
        # cdf
        p_value_list = []

        for i, data_element in enumerate(pDataList):
            relative_distance = i - pIndexReferencePoint
            # check if relative distance is available.
            # if not, use either min or max key for the distribution
            if relative_distance not in pBackgroundModel:
                if relative_distance < 0:
                    relative_distance = min(pBackgroundModel.keys())
                else:
                    relative_distance = max(pBackgroundModel.keys())
            if data_element == 0.0:
                p_value_list.append(0.0)
            else:
                p_value_list.append(cnb.cdf(data_element, pBackgroundModel[relative_distance][0], pBackgroundModel[relative_distance][1]))

        # for reason I do not understand atm the values needs to be inverted again, it seems it is not enough to do this in try/catch region
        p_value_list = np.array(p_value_list, dtype=np.float64)

        p_value_list = 1 - p_value_list

        # remove possible occuring nan with a p-value of 1
        mask = np.isnan(p_value_list)
        mask_inf = np.isinf(p_value_list)
        p_value_list = np.array(p_value_list)
        mask = np.logical_or(mask, mask_inf)
        p_value_list[mask] = 1.0
        return p_value_list

    def merge_neighbors(self, pScoresDictionary, pMergeThreshold=1000):
        if pScoresDictionary is None or len(pScoresDictionary) == 0:
            log.debug('dict is None or empty')
            return None
        key_list = list(pScoresDictionary.keys())

        merge_ids = []
        non_merge = []
        for i, (key_pre, key_suc) in enumerate(zip(key_list[:-1], key_list[1:])):

            if np.absolute(int(pScoresDictionary[key_pre][5]) - int(pScoresDictionary[key_suc][5])) <= pMergeThreshold:
                if len(merge_ids) > 0 and merge_ids[-1][-1] == key_pre:
                    merge_ids[-1].append(key_suc)
                else:
                    merge_ids.append([key_pre, key_suc])
            else:
                if i == len(key_list) - 1:
                    non_merge.append(key_suc)
                if merge_ids is not None and len(merge_ids) > 0 and merge_ids[-1][-1] != key_pre:
                    non_merge.append(key_pre)
                elif merge_ids is not None and len(merge_ids) == 0:
                    non_merge.append(key_pre)

        scores_dict = {}
        merged_lines_dict = {}
        for element in merge_ids:
            lines = []
            lines.append(pScoresDictionary[element[0]])
            index_maximum_element = 0
            base_element = pScoresDictionary[element[0]]
            values = np.array(list(map(float, base_element[-4:])))
            max_value = float(base_element[-1])
            for i, key in enumerate(element[1:]):
                base_element[-6] = pScoresDictionary[key][-6]
                values += np.array(list(map(float, pScoresDictionary[key][-4:])))
                lines.append(pScoresDictionary[key])
                if max_value < float(pScoresDictionary[key][-1]):
                    max_value = float(pScoresDictionary[key][-1])
                    index_maximum_element = i + 1
            base_element = pScoresDictionary[element[index_maximum_element]]
            base_element[-4] = values[0]
            base_element[-3] = values[1]
            base_element[-2] = values[2]
            base_element[-1] = values[3]

            base_element[2] = pScoresDictionary[element[-1]][2]
            base_element[1] = pScoresDictionary[element[0]][1]

            scores_dict[element[index_maximum_element]] = base_element
            merged_lines_dict[element[index_maximum_element]] = lines

        for key in non_merge:
            scores_dict[key] = pScoresDictionary[key]
            merged_lines_dict[key] = [pScoresDictionary[key]]

        return scores_dict, merged_lines_dict

    def readSignificantRegionsFile(self, pFilePath, pInternalIdentifierTriplet, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint):
        # list of start and end point of regions to highlight

        _, interaction_file_data, _ = self.readInteractionFile(pFilePath, pInternalIdentifierTriplet)
        highlight_areas_list = []
        p_values = []
        if len(pViewpoint) == 2:
            reference_point_start, reference_point_end = pViewpoint
        else:
            log.debug('viewpoint_split {}, file: {} {}'.format(pViewpoint, pFilePath, pInternalIdentifierTriplet))
            return None, None

        for _line in interaction_file_data.values():
            start = int(_line[1])
            end = int(_line[2])

            if int(_line[5]) >= -pRange[0] and int(_line[5]) <= pRange[1]:

                width = (end - start) / pResolution
                if int(_line[5]) < 0:
                    relative_position_genomic_coordinates = start - int(reference_point_start)
                    viewpointIndex = pViewpointIndexStart
                else:
                    relative_position_genomic_coordinates = start - int(reference_point_end)
                    viewpointIndex = pViewpointIndexEnd

                relative_position = viewpointIndex + \
                    (relative_position_genomic_coordinates / pResolution)
                highlight_areas_list.append(
                    [relative_position, relative_position + width])
                p_values.append([int(relative_position), int(relative_position + width), float(_line[-3])])

        if len(highlight_areas_list) == 0:
            return None, None

        return highlight_areas_list, p_values

    def readTargetHDFFile(self, pFile):
        present_genes = {}
        targetDict = {}
        targetFileHDF5Object = h5py.File(pFile, 'r')
        keys_targetFile = list(targetFileHDF5Object.keys())
        log.debug('keys_interactionFile {}'.format(keys_targetFile))
        for outer_matrix in keys_targetFile:
            if outer_matrix not in present_genes:
                present_genes[outer_matrix] = {}
            inner_matrix_object = targetFileHDF5Object[outer_matrix]
            keys_inner_matrices = list(inner_matrix_object.keys())
            for inner_matrix in keys_inner_matrices:
                if inner_matrix not in present_genes[outer_matrix]:
                    present_genes[outer_matrix][inner_matrix] = []
                inner_object = inner_matrix_object[inner_matrix]
                gene_object = inner_object['genes']
                keys_genes = list(gene_object.keys())
                for gene in keys_genes:
                    targetDict[gene] = [outer_matrix, inner_matrix, 'genes', gene]
                    present_genes[outer_matrix][inner_matrix].append(gene)
        targetFileHDF5Object.close()
        return targetDict, present_genes

    def readAggregatedFileHDF(self, pAggregatedFileName, pInternalPath):
        aggregatedFileHDF5Object = h5py.File(pAggregatedFileName, 'r')

        internal_path = '/'.join(pInternalPath)
        chromosome = aggregatedFileHDF5Object.get(internal_path + '/' + 'chromosome')[()].decode("utf-8")

        gene_name = aggregatedFileHDF5Object.get(internal_path + '/' + 'gene_name')[()].decode("utf-8")
        reference_point_start = aggregatedFileHDF5Object.get(internal_path + '/' + 'reference_point_start')[()]
        reference_point_end = aggregatedFileHDF5Object.get(internal_path + '/' + 'reference_point_end')[()]

        # log.debug('reference_point_start {}'.format(reference_point_start))
        # log.debug('reference_point_end {}'.format(reference_point_end))


        start_list = np.array(aggregatedFileHDF5Object[internal_path + '/' + 'start_list'][:])
        end_list = np.array(aggregatedFileHDF5Object[internal_path + '/' + 'end_list'][:])
        relative_distance_list = np.array(aggregatedFileHDF5Object[internal_path + '/' + 'relative_distance_list'][:])
        raw_target_list = np.array(aggregatedFileHDF5Object[internal_path + '/' + 'raw_target_list'][:])
        sum_of_interactions = float(aggregatedFileHDF5Object.get(internal_path + '/' + 'sum_of_interactions')[()])
        aggregatedFileHDF5Object.close()
        line_content = []
        data = []
        for i in range(len(start_list)):
            line_content.append([chromosome, start_list[i], end_list[i], gene_name, sum_of_interactions, relative_distance_list[i], raw_target_list[i], reference_point_start, reference_point_end])
            data.append([sum_of_interactions, raw_target_list[i]])

        return line_content, data

    def readDifferentialFile(self, pFilePath, pQuadruple):

        differentialFileHDF5Object = h5py.File(pFilePath, 'r')

        internal_path = '/'.join(pQuadruple)
        output_list = []
        for item in ['accepted', 'all', 'rejected']:

            try:
                item_object = differentialFileHDF5Object[internal_path + '/' + item]
                chromosome = item_object.get('chromosome')[()].decode("utf-8")
                reference_point_start = int(item_object.get('reference_point_start')[()])
                reference_point_end = int(item_object.get('reference_point_end')[()])
                
                # log.debug('reference_point_start {}'.format(int(reference_point_start)))
                # log.debug('reference_point_end {}'.format(int(reference_point_end)))

                start_list = np.array(item_object['start_list'][:])
                end_list = np.array(item_object['end_list'][:])
                relative_distance_list = np.array(item_object['relative_distance_list'][:])
                gene = item_object.get('gene')[()].decode("utf-8")
                pvalue_list = np.array(item_object['pvalue_list'][:])

                raw_target_list_1 = np.array(item_object['raw_target_list_1'][:])
                raw_target_list_2 = np.array(item_object['raw_target_list_2'][:])
                sum_of_interactions_1 = float(item_object.get('sum_of_interactions_1')[()])
                sum_of_interactions_2 = float(item_object.get('sum_of_interactions_2')[()])

                chromosome = [chromosome] * len(start_list)
                gene = [gene] * len(start_list)
                reference_point_start_list = [reference_point_start]  * len(start_list)
                reference_point_end_list = [reference_point_end]  * len(start_list)
                sum_of_interactions_1 = [sum_of_interactions_1] * len(start_list)
                sum_of_interactions_2 = [sum_of_interactions_2] * len(start_list)
                output_list.append(zip(chromosome, start_list, end_list, gene, relative_distance_list, sum_of_interactions_1, raw_target_list_1, sum_of_interactions_2, raw_target_list_2, pvalue_list, reference_point_start_list, reference_point_end_list))
                # log.debug('new item {}'.format(list(output_list[-1])))
            except Exception as exp:
                # log.debug('readDifferentialFile exception {}'.format(str(exp)))
                output_list.append([])

        return output_list
