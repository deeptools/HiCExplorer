import numpy as np

import logging
log = logging.getLogger(__name__)
import copy


class Viewpoint():

    def __init__(self, pHiCMatrix=None):
        self.hicMatrix = pHiCMatrix

    def readReferencePointFile(self, pBedFile):
        '''
        This function reads a text file which contains reference points. Reference points are 

        Reference points need to be tab seperated
        and contain the chromosome, the reference start point and/or the reference end point:

        chr start

        or 

        chr start end

        Per line one reference point.

        Returns a list of tuples with (chr, start, end), if only a start index is given, end = start. 
        '''
        viewpoints = []
        gene_list = []
        with open(pBedFile, 'r') as file:
            for line in file.readlines():
                _line = line.strip().split('\t')
                if len(line) == 0:
                    continue
                if len(_line) == 2:
                    chrom, start, end = _line[0], _line[1], _line[1]
                else:
                    chrom, start, end = _line[:3]
                viewpoints.append((chrom, start, end))
                gene_list.append(_line[3])
        return viewpoints, gene_list

    def readInteractionFile(self, pBedFile):
        '''
        Reads an interaction file produced by chicViewpoint. Contains header information, these lines
        start with '#'. 
        Interactions files contain:
        Chromosome Viewpoint, Start, End, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, z-score based on relative interactions.

        This function returns:
        - header as  a string
        - interaction data in relation to relative position as a dict e.g. {-1000:0.1, -1500:0.2}  
        - z-score in relation to relative position as a dict (same format as interaction data)
        - interaction_file_data: the raw line in relation to the relative position. Needed for additional output file.
        '''
        # use header info to store reference point, and based matrix
        interaction_data = {}
        z_score = {}
        interaction_file_data = {}
        with open(pBedFile) as fh:
            header = fh.readline()
            for line in fh.readlines():
                # Addition header information for end users
                if line.strip().startswith('#'):
                    continue

                _line = line.strip().split('\t')
                # relative postion and relative interactions
                interaction_data[int(_line[-4])] = float(_line[-3])
                z_score[int(_line[-4])] = float(_line[-2])
                interaction_file_data[int(_line[-4])] = _line
        return header, interaction_data, z_score, interaction_file_data

    def readBackgroundDataFile(self, pBedFile):
        '''
        Reads a background data file, containing per line a tab delimited content:
        Relative position to viewpoint, relative interaction count to total number of interactions of all viewpoints over all samples, SEM value of this data point.
        '''
        distance = {}
        with open(pBedFile) as fh:
            for line in fh.readlines():
                _line = line.split('\t')
                distance[int(_line[0])] = [float(_line[1]), float(_line[2])]

        return distance

    def writeInteractionFile(self, pBedFile, pData, pHeader, pZscoreData):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Viewpoint, Start, End, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, z-score based on relative interactions, raw interaction data
        '''
        with open(pBedFile + '.bed', 'w') as fh:
            fh.write('#{}\n'.format(pHeader))
            for j, interaction in enumerate(pData):
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\t{:.12f}\t{:.12f}\n".
                         format(interaction[0], interaction[1], interaction[2],
                                interaction[3], interaction[4], interaction[5],
                                interaction[6], interaction[7], pZscoreData[j], interaction[8]))
        return

    def computeViewpoint(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        This function computes a viewpoint for a given sample and a given pReferencePoint within the  
        range of pRegion_start and pRegion_end. 

        All interactions with the reference point of one relative distance to it are summed up,
        if the reference point is larger than one bin of the Hi-C matrix, it is considered as one bin and the values are summed together.
        '''
        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(pReferencePoint)

        view_point_range = self.getViewpointRangeAsMatrixIndices(pChromViewpoint, pRegion_start, pRegion_end)

        elements_of_viewpoint = (view_point_range[1] - view_point_range[0])

        data_list = np.zeros(elements_of_viewpoint)
        _view_point_start = view_point_start

        while _view_point_start <= view_point_end:
            chrom, start, end, _ = self.hicMatrix.getBinPos(_view_point_start)
            for j, idx in zip(range(elements_of_viewpoint), range(view_point_range[0], view_point_range[1], 1)):
                data_list[j] += self.hicMatrix.matrix[_view_point_start, idx]

            _view_point_start += 1

        elements_of_viewpoint = elements_of_viewpoint - (view_point_end - view_point_start)
        data_list_new = np.zeros(elements_of_viewpoint)

        index_before_viewpoint = view_point_start - view_point_range[0]

        # elements before the viewpoint
        data_list_new[0:index_before_viewpoint] = data_list[0:index_before_viewpoint]

        # summation because the viewpoint can not be only one bin but can contain multiple
        data_list_new[index_before_viewpoint] = np.sum(data_list[index_before_viewpoint: index_before_viewpoint + view_point_end - view_point_start + 1])

        # elements after the viewpoint
        data_list_new[index_before_viewpoint + 1:] = data_list[index_before_viewpoint + view_point_end - view_point_start + 1:]
        return data_list_new

    def createInteractionFileData(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw):
        '''
        Creates out of internal information a list of tuples which can be written to an interaction file.
        Tuple contains:
        Chromosome viewpoint, start, end, chromosome interaction, start, end, relative_position, interaction data
        '''
        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(pReferencePoint)
        view_point_range = self.getViewpointRangeAsMatrixIndices(pChromViewpoint, pRegion_start, pRegion_end)
        elements_of_viewpoint = view_point_range[1] - view_point_range[0] - (view_point_end - view_point_start) + 1

        interactions_list = []
        chrom, start, _, _ = self.hicMatrix.getBinPos(view_point_start)
        _, _, end, _ = self.hicMatrix.getBinPos(view_point_end)
        interaction_positions = list(range(view_point_range[0], view_point_start, 1))
        interaction_positions.extend([view_point_start])
        interaction_positions.extend(list(range(view_point_end + 1, view_point_range[1], 1)))
        relative_position = -1
        for j, idx in zip(range(elements_of_viewpoint), interaction_positions):

            chrom_second, start_second, end_second, _ = self.hicMatrix.getBinPos(idx)
            if relative_position < 0:
                relative_position = int(start_second) - int(start)
            else:
                relative_position = int(end_second) - int(end)

            interactions_list.append((chrom, start, end, chrom_second, start_second, end_second, relative_position, float(pInteractionData[j]), float(pInteractionDataRaw[j])))

        return interactions_list

    def getViewpointRangeAsMatrixIndices(self, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        Returns the matrix indices of a chromosome and a specific position.
        '''
        return self.hicMatrix.getRegionBinRange(pChromViewpoint, pRegion_start, pRegion_end)

    def getReferencePointAsMatrixIndices(self, pReferencePoint):
        '''
        Returns the correct matrix indices of a given reference point.
        '''
        if len(pReferencePoint) == 2:
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[1]))
        elif len(pReferencePoint) == 3:
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
        else:
            log.error("No valid reference point given. {}".format(pReferencePoint))
            exit(1)
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
        return average_contacts

    def computeRelativeValues(self, pData):
        '''
        Computes the relative values of pData by adding all data points together and divding all data points by the result.
        pData[i] = pData[i] / sum(pData)
        '''
        sumValue = np.sum(pData)
        pData /= sumValue
        return pData

    def calculateViewpointRange(self, pViewpoint, pRange):
        '''
        This function computes the correct start and end position of a viewpoint given the viewpoint and the range.
        '''
        max_length = self.hicMatrix.getBinPos(self.hicMatrix.getChrBinRange(pViewpoint[0])[1] - 1)[2]

        region_start = int(pViewpoint[1]) - pRange[0]
        if region_start < 0:
            region_start = 0

        region_end = int(pViewpoint[2]) + pRange[1]
        if region_end > max_length:
            region_end = max_length - 1

        return region_start, region_end

    def getDataForPlotting(self, pInteractionFile, pRange, pBackgroundModel):
        header, interaction_data, z_score_data, _interaction_file_data_raw = self.readInteractionFile(pInteractionFile)
        matrix_name, viewpoint, upstream_range, downstream_range, gene = header.split('\t')

        data = []
        z_score = []
        interaction_file_data_raw = {}
        if pRange:
            background_data_keys_sorted = sorted(pBackgroundModel)

            interaction_data_keys = copy.deepcopy(list(interaction_data.keys()))
            for key in interaction_data_keys:
                if key >= -pRange[0] and key <= pRange[1]:
                    continue
                interaction_data.pop(key, None)
            if pBackgroundModel:
                for key in background_data_keys_sorted:
                    if key >= -pRange[0] and key <= pRange[1]:
                        continue
                    pBackgroundModel.pop(key, None)
                background_data_keys_sorted = sorted(pBackgroundModel)

        if pBackgroundModel:
            viewpoint_index = background_data_keys_sorted.index(0)

            data_background = []
            data_background_mean = []

            for key in background_data_keys_sorted:
                if key in interaction_data:
                    data.append(interaction_data[key])

                    if key in z_score_data:
                        z_score.append(z_score_data[key])

                    data_background.append(pBackgroundModel[key][0])
                    data_background_mean.append(pBackgroundModel[key][1])

                    if key in _interaction_file_data_raw:
                        line_data_raw = _interaction_file_data_raw[key]
                        line_data_raw.append(pBackgroundModel[key][0])
                        line_data_raw.append(pBackgroundModel[key][1])

                        interaction_file_data_raw[key] = line_data_raw

        else:
            data = []
            interaction_key = sorted(interaction_data)
            for key in interaction_key:
                data.append(interaction_data[key])
            viewpoint_index = interaction_key.index(0)

        return header, data, data_background, data_background_mean, z_score, interaction_file_data_raw, viewpoint_index

    def plotViewpoint(self, pAxis, pData, pColor, pLabelName):
        data_plot_label = pAxis.plot(range(len(pData)), pData, '--' + pColor, alpha=0.9, label=pLabelName)

        return data_plot_label

    def plotBackgroundModel(self, pAxis, pBackgroundData, pBackgroundDataMean):
        pBackgroundData = np.array(pBackgroundData)
        pBackgroundDataMean = np.array(pBackgroundDataMean)
        data_plot_label = pAxis.plot(range(len(pBackgroundData)), pBackgroundData, '--r', alpha=0.5, label='background model')
        pAxis.fill_between(range(len(pBackgroundData)), pBackgroundData + pBackgroundDataMean, pBackgroundData - pBackgroundDataMean, facecolor='red', alpha=0.3)
        return data_plot_label

    def plotZscore(self, pAxis, pAxisLabel, pZscoreData, pLabelText, pCmap):

        _z_score = np.empty([2, len(pZscoreData)])
        _z_score[:, :] = pZscoreData
        pAxis.xaxis.set_visible(False)
        pAxis.yaxis.set_visible(False)
        pAxis.contourf(_z_score, cmap=pCmap)
        pAxisLabel.text(0, 0, pLabelText)
        pAxisLabel.xaxis.set_visible(False)
        pAxisLabel.yaxis.set_visible(False)
        pAxisLabel.set_frame_on(False)

    def writePlotData(self, pInteractionFileDataRaw, pFileName, pBackgroundModel):
        interaction_file_data_raw_sorted = sorted(pInteractionFileDataRaw)
        with open(pFileName + '.bed', 'w') as output_file:
            output_file.write('#ChrViewpoint\tStart\tEnd\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\tZ-score\tRaw')

            if pBackgroundModel:
                output_file.write('\tbackground model\tbackground model SEM\n')
            else:
                output_file.write('\n')
            for key in interaction_file_data_raw_sorted:

                _array = pInteractionFileDataRaw[key]

                output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(_array[0], _array[1], _array[2],
                                                                                  _array[3], _array[4], _array[5],
                                                                                  _array[6], _array[7], _array[8], _array[9]))
                if pBackgroundModel:
                    output_file.write('\t{}\t{}\n'.format(_array[10], _array[11]))
                else:
                    output_file.write('\n')
