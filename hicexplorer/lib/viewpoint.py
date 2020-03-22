import numpy as np

import logging
log = logging.getLogger(__name__)
import copy

import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import nbinom
from scipy.special import gammaln
from scipy import special


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

    def readInteractionFile(self, pBedFile):
        '''
        Reads an interaction file produced by chicViewpoint. Contains header information, these lines
        start with '#'.
        Interactions files contain:
        Chromosome Viewpoint, Start, End, Gene, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-score based on relative interactions.

        This function returns:
        - header as  a string
        - interaction data in relation to relative position as a dict e.g. {-1000:0.1, -1500:0.2}
        - p-score in relation to relative position as a dict (same format as interaction data)
        - interaction_file_data: the raw line in relation to the relative position. Needed for additional output file.
        '''
        # use header info to store reference point, and based matrix
        interaction_data = {}
        p_score = {}
        interaction_file_data = {}
        genomic_coordinates = {}
        with open(pBedFile) as fh:
            fh.readline()
            header = fh.readline()
            for line in fh.readlines():
                # Addition header information for end users
                if line.strip().startswith('#'):
                    continue

                _line = line.strip().split('\t')
                # relative postion and relative interactions
                interaction_data[int(_line[-5])] = float(_line[-4])
                p_score[int(_line[-5])] = float(_line[-3])
                interaction_file_data[int(_line[-5])] = _line
                genomic_coordinates[int(_line[-5])] = [_line[0], _line[1], _line[2]]

        return header, interaction_data, p_score, interaction_file_data, genomic_coordinates

    def readInteractionFileForAggregateStatistics(self, pBedFile):
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
        # use header info to store reference point, and based matrix
        interaction_data = {}
        interaction_file_data = {}
        with open(pBedFile) as fh:
            fh.readline()
            header = fh.readline()
            fh.readline()

            for line in fh.readlines():
                # Addition header information for end users
                if line.strip().startswith('#'):
                    continue
                if not line:
                    continue
                _line = line.strip().split('\t')
                # relative postion and relative interactions
                interaction_data[int(
                    _line[-5])] = np.array([float(_line[-4]), float(_line[-3]), float(_line[-1]), float(_line[-2])])
                interaction_file_data[int(_line[-5])] = _line
        return header, interaction_data, interaction_file_data

    def readBackgroundDataFile(self, pBedFile, pRange, pMean=False):
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

    def writeInteractionFile(self, pBedFile, pData, pHeader, pPValueData, pXfold):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-values based on negative binomial distribution per relative distance, raw interaction data
        '''

        with open((pBedFile + '.txt').strip(), 'w') as fh:
            fh.write('{}\n'.format(pHeader))
            for j, interaction in enumerate(pData):
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\t{:.12f}\t{:.12f}\t{:.12f}\n".
                         format(interaction[0], interaction[1],
                                interaction[2], interaction[3], interaction[4], interaction[5], interaction[6], pPValueData[j], pXfold[j], interaction[7]))
        return

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
            log.debug('pReferencePoint: {}'.format(pReferencePoint))
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(
                pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
        else:
            log.error("No valid reference point given. {}".format(
                pReferencePoint))
            exit(1)
        log.debug('view_point_start: {} view_point_end {}'.format(view_point_start, view_point_end))

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

    def getDataForPlotting(self, pInteractionFile, pRange, pBackgroundModel, pResolution):
        header, interaction_data, p_value_data, _interaction_file_data_raw, genomic_coordinates = self.readInteractionFile(pInteractionFile)
        matrix_name, viewpoint, upstream_range, downstream_range, gene, _ = header.split('\t')

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
        return header, data, data_background, p_value, viewpoint_index_start, viewpoint_index_end

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

    def readRejectedFile(self, pDifferentialHighlightFiles, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint):
        # list of start and end point of regions to highlight
        # [[start, end], [start, end]]
        highlight_areas_list = []
        # for bed_file in pDifferentialHighlightFiles:
        _, reference_point_start, reference_point_end = pViewpoint.split('_')

        with open(pDifferentialHighlightFiles) as fh:
            # skip header
            for line in fh.readlines():
                if line.startswith('#'):
                    continue
                _line = line.split('\t')

                start = int(_line[1])
                end = int(_line[2])

                if int(_line[4]) >= -pRange[0] and int(_line[4]) <= pRange[1]:

                    width = (end - start) / pResolution
                    if int(_line[4]) < 0:
                        # reference_point_position = reference_point_start
                        relative_position_genomic_coordinates = start - int(reference_point_start)
                        viewpointIndex = pViewpointIndexStart
                    else:
                        relative_position_genomic_coordinates = start - int(reference_point_end)
                        viewpointIndex = pViewpointIndexEnd

                    log.debug('_line[4] {}'.format(_line[4]))
                    log.debug('relative_position_genomic_coordinates {}'.format(relative_position_genomic_coordinates))
                    log.debug('start {} end {}'.format(start, end))
                    log.debug('reference_point_start {} reference_point_end {}'.format(reference_point_start, reference_point_end))

                    # start

                    # relative_position_genomic_coordinates  = reference_point_position
                    relative_position = viewpointIndex + \
                        (relative_position_genomic_coordinates / pResolution)
                    highlight_areas_list.append(
                        [relative_position, relative_position + width])
        if len(highlight_areas_list) == 0:
            return None
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
                p_value_list.append(self.cdf(data_element, pBackgroundModel[relative_distance][0], pBackgroundModel[relative_distance][1]))

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

    # def computeSumOfDensities(self, pBackgroundModel, pArgs, pXfoldMaxValue=None):
    #     background_nbinom = {}
    #     background_sum_of_densities_dict = {}
    #     max_value = 0

    #     fixateRange = int(pArgs.fixateRange)
    #     for distance in pBackgroundModel:
    #         max_value_distance = int(pBackgroundModel[distance][2])
    #         if max_value < int(pBackgroundModel[distance][2]):
    #             max_value = int(pBackgroundModel[distance][2])

    #         if pXfoldMaxValue is not None:
    #             if max_value_distance == 0:
    #                 max_value_distance = 1
    #             if pXfoldMaxValue == 0:
    #                 pXfoldMaxValue = 1
    #             max_value_distance *= pXfoldMaxValue

    #         if -int(pArgs.fixateRange) < distance and int(pArgs.fixateRange) > distance:
    #             # background_nbinom[distance] = nbinom(pBackgroundModel[distance][0], pBackgroundModel[distance][1])
    #             background_nbinom[distance] = (pBackgroundModel[distance][0], pBackgroundModel[distance][1])

    #             sum_of_densities = np.zeros(max_value_distance)
    #             for j in range(max_value_distance):
    #                 if j >= 1:
    #                     sum_of_densities[j] += sum_of_densities[j - 1]
    #                 # sum_of_densities[j] += background_nbinom[distance].pmf(j)
    #                 sum_of_densities[j] += pmf(j, background_nbinom[distance][0], background_nbinom[distance][1])

    #             background_sum_of_densities_dict[distance] = sum_of_densities

    #     # background_nbinom[fixateRange] = nbinom(pBackgroundModel[fixateRange][0], pBackgroundModel[fixateRange][1])
    #     background_nbinom[fixateRange] = (pBackgroundModel[fixateRange][0], pBackgroundModel[fixateRange][1])

    #     sum_of_densities = np.zeros(max_value)
    #     for j in range(max_value):
    #         if j >= 1:
    #             sum_of_densities[j] += sum_of_densities[j - 1]
    #         # sum_of_densities[j] += background_nbinom[fixateRange].pmf(j)
    #         sum_of_densities[j] += pmf(j, background_nbinom[fixateRange][0], background_nbinom[fixateRange][1])

    #     background_sum_of_densities_dict[fixateRange] = sum_of_densities
    #     # background_nbinom[-fixateRange] = nbinom(pBackgroundModel[-fixateRange][0], pBackgroundModel[-fixateRange][1])
    #     background_nbinom[-fixateRange] = (pBackgroundModel[-fixateRange][0], pBackgroundModel[-fixateRange][1])

    #     sum_of_densities = np.zeros(max_value)
    #     for j in range(max_value):
    #         if j >= 1:
    #             sum_of_densities[j] += sum_of_densities[j - 1]
    #         # sum_of_densities[j] += background_nbinom[-fixateRange].pmf(j)
    #         sum_of_densities[j] += pmf(j, background_nbinom[-fixateRange][0], background_nbinom[-fixateRange][1])

    #     background_sum_of_densities_dict[-fixateRange] = sum_of_densities

    #     min_key = min(background_sum_of_densities_dict)
    #     max_key = max(background_sum_of_densities_dict)

    #     for key in pBackgroundModel.keys():
    #         if key in background_sum_of_densities_dict:
    #             continue
    #         if key < min_key:
    #             background_sum_of_densities_dict[key] = background_sum_of_densities_dict[min_key]
    #         elif key > max_key:
    #             background_sum_of_densities_dict[key] = background_sum_of_densities_dict[max_key]

    #     return background_sum_of_densities_dict

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

    def readSignificantRegionsFile(self, pSignificantFile, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint):
        # list of start and end point of regions to highlight
        # [[start, end], [start, end]]
        highlight_areas_list = []
        p_values = []
        viewpoint_split = pViewpoint.split('_')
        if len(viewpoint_split) == 3:
            _, reference_point_start, reference_point_end = viewpoint_split
        else:
            log.debug('viewpoint_split {}, file: {}'.format(viewpoint_split, pSignificantFile))
            return None, None
        with open(pSignificantFile) as fh:
            # skip header
            for line in fh.readlines():
                if line.startswith('#'):
                    continue
                _line = line.split('\t')

                start = int(_line[1])
                end = int(_line[2])

                if int(_line[5]) >= -pRange[0] and int(_line[5]) <= pRange[1]:

                    width = (end - start) / pResolution
                    if int(_line[5]) < 0:
                        # reference_point_position = reference_point_start
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

    def pdf(self, pX, pR, pP):
        """
        PDF for a continuous generalization of NB distribution
        """

        gamma_part = gammaln(pR + pX) - gammaln(pX + 1) - gammaln(pR)
        return np.exp(gamma_part + (pR * log(pP)) + special.xlog1py(pX, -pP))

    def cdf(self, pX, pR, pP):
        """
        Cumulative density function of a continuous generalization of NB distribution
        """
        # if pX == 0:
        # return 0
        return special.betainc(pR, pX + 1, pP)
