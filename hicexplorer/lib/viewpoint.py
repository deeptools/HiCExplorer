"""hicexplorer.lib.viewpoint, implemented in C++.

Viewpoint keeps the attribute (hicMatrix), the methods and the behaviour it
had in HiCExplorer 3.7.6; every method body calls hicexplorer._hicx.lib_viewpoint
with self (cpp/python/src/api/lib_viewpoint.cpp, cpp/core/include/hicx/chic_viewpoint.hpp
and chic_hdf5.hpp). plotViewpoint, plotBackgroundModel and plotPValue draw
with matplotlib; the data they draw is computed by the extension.
"""
import logging
log = logging.getLogger(__name__)

from mpl_toolkits.axes_grid1 import make_axes_locatable

from hicexplorer._hicx import lib_viewpoint as _impl


class Viewpoint():
    def __init__(self, pHiCMatrix=None):
        _impl.Viewpoint___init__(self, pHiCMatrix)

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
        return _impl.Viewpoint_readReferencePointFile(self, pBedFile, pGene)

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
        return _impl.Viewpoint_readInteractionFile(self, pFilePath, pInternalIdentifierTriplet)

    def interactionDataForPlot(self, pFilePath, pInternalIdentifierTriplet):
        return _impl.Viewpoint_interactionDataForPlot(self, pFilePath, pInternalIdentifierTriplet)

    def readBackgroundDataFile(self, pBedFile, pRange, pFixateRange, pMean=False):
        '''
        Reads a background data file, containing per line a tab delimited content:
        Relative position to viewpoint, relative interaction count to total number of interactions of all viewpoints over all samples, SEM value of this data point.

        '''
        return _impl.Viewpoint_readBackgroundDataFile(self, pBedFile, pRange, pFixateRange, pMean)

    def writeInteractionFile(self, pBedFile, pData, pHeader, pPValueData, pXfold, pDecimalPlaces=12):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-values based on negative binomial distribution per relative distance, raw interaction data
        '''
        return _impl.Viewpoint_writeInteractionFile(self, pBedFile, pData, pHeader, pPValueData, pXfold, pDecimalPlaces)

    def writeInteractionFileHDF5(self, pInteractionFileGroupH5Object, pFileName, pData, pReferencePoint):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, p-values based on negative binomial distribution per relative distance, raw interaction data
        '''
        return _impl.Viewpoint_writeInteractionFileHDF5(self, pInteractionFileGroupH5Object, pFileName, pData, pReferencePoint)

    def createUniqueHDFGroup(self, pGroupObject, pAdditionalGroupName):
        return _impl.Viewpoint_createUniqueHDFGroup(self, pGroupObject, pAdditionalGroupName)

    def computeViewpoint(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        This function computes a viewpoint for a given sample and a given pReferencePoint within the
        range of pRegion_start and pRegion_end.

        All interactions with the reference point of one relative distance to it are summed up,
        if the reference point is larger than one bin of the Hi-C matrix, it is considered as one bin and the values are summed together.
        '''
        return _impl.Viewpoint_computeViewpoint(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end)

    def createInteractionFileData(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene, pSumOfInteractions):
        '''
        Creates out of internal information a list of tuples which can be written to an interaction file.
        Tuple contains:
        Chromosome viewpoint, start, end, chromosome interaction, start, end, relative_position, interaction data
        '''
        return _impl.Viewpoint_createInteractionFileData(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene, pSumOfInteractions)

    def createInteractionFileDataHDF5(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene, pSumOfInteractions, pPValueList, xFoldList):
        '''
        Creates out of internal information a list of tuples which can be written to an interaction file.
        Tuple contains:
        Chromosome viewpoint, start, end, chromosome interaction, start, end, relative_position, interaction data
        '''
        return _impl.Viewpoint_createInteractionFileDataHDF5(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene, pSumOfInteractions, pPValueList, xFoldList)

    def getViewpointRangeAsMatrixIndices(self, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        Returns the matrix indices of a chromosome and a specific position.
        '''
        return _impl.Viewpoint_getViewpointRangeAsMatrixIndices(self, pChromViewpoint, pRegion_start, pRegion_end)

    def getReferencePointAsMatrixIndices(self, pReferencePoint):
        '''
        Returns the correct matrix indices of a given reference point.
        '''
        return _impl.Viewpoint_getReferencePointAsMatrixIndices(self, pReferencePoint)

    def smoothInteractionValues(self, pData, pWindowSize):
        '''
        Adds -pWindowsSize/2 and +pWindowsSize/2 around pData[i] and averages pData[i] by pWindowSize to
        smooth the interaction values.
        '''
        return _impl.Viewpoint_smoothInteractionValues(self, pData, pWindowSize)

    def computeRelativeValues(self, pData, pDenominator=None):
        '''
        Computes the relative values of pData by adding all data points together and divding all data points by the result.
        pData[i] = pData[i] / sum(pData)
        '''
        return _impl.Viewpoint_computeRelativeValues(self, pData, pDenominator)

    def calculateViewpointRange(self, pViewpoint, pRange):
        '''
        This function computes the correct start and end position of a viewpoint given the viewpoint and the range.
        '''
        return _impl.Viewpoint_calculateViewpointRange(self, pViewpoint, pRange)

    def getDataForPlotting(self, pFilePath, pInteractionFile, pRange, pBackgroundModel, pResolution):
        return _impl.Viewpoint_getDataForPlotting(self, pFilePath, pInteractionFile, pRange, pBackgroundModel, pResolution)

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
        pBackgroundData, upper_values = _impl.Viewpoint_backgroundModelBand(pBackgroundData, pXFold)
        data_plot_label = pAxis.plot(range(len(pBackgroundData)), pBackgroundData, '-r', alpha=0.5, label='background model', linewidth=1)
        if upper_values is not None:
            pAxis.fill_between(range(len(pBackgroundData)), upper_values, pBackgroundData, facecolor='r', alpha=0.5)
        return data_plot_label

    def plotPValue(self, pAxis, pAxisLabel, pPValueData, pLabelText, pCmap, pFigure, pValueSignificanceLevels):
        _z_score = _impl.Viewpoint_pValueRows(pPValueData)
        pAxis.xaxis.set_visible(False)
        pAxis.yaxis.set_visible(False)
        divider = make_axes_locatable(pAxisLabel)
        cax = divider.append_axes("left", size="20%", pad=0.09)

        if pPValueData is not None:
            img = pAxis.contourf(_z_score, cmap=pCmap)
            colorbar = pFigure.colorbar(
                img, cax=cax, ticks=_impl.Viewpoint_pValueTicks(pPValueData))
            colorbar.ax.set_ylabel('p-value', size=6)

        elif pValueSignificanceLevels:
            ticks = _impl.Viewpoint_pValueLevels(pValueSignificanceLevels)
            img = pAxis.contourf(_z_score, levels=pValueSignificanceLevels, colors=['#CC0000', '#FFD43B', '#306998', '#FFFFFF'])
            colorbar = pFigure.colorbar(img, cax=cax, ticks=ticks)
            colorbar.ax.tick_params(labelsize=6)
            colorbar.ax.set_ylabel('p-value', size=6)

        pAxisLabel.text(0.45, 0, pLabelText, size=7)
        pAxisLabel.xaxis.set_visible(False)
        pAxisLabel.yaxis.set_visible(False)
        pAxisLabel.set_frame_on(False)

    def writePlotData(self, pInteractionFileDataRaw, pFileName, pBackgroundModel):
        return _impl.Viewpoint_writePlotData(self, pInteractionFileDataRaw, pFileName, pBackgroundModel)

    def interactionBackgroundData(self, pBackground, pRange):
        return _impl.Viewpoint_interactionBackgroundData(self, pBackground, pRange)

    def rbz_score(self, pRelativeInteractions, pBackgroundModel, pBackgroundModelSEM):
        return _impl.Viewpoint_rbz_score(self, pRelativeInteractions, pBackgroundModel, pBackgroundModelSEM)

    def readRejectedFile(self, pFilePath, pDifferentialHighlightTriplet, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint):
        return _impl.Viewpoint_readRejectedFile(self, pFilePath, pDifferentialHighlightTriplet, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint)

    def pvalues(self, pBackgroundModel, pDataList, pIndexReferencePoint):
        return _impl.Viewpoint_pvalues(self, pBackgroundModel, pDataList, pIndexReferencePoint)

    def merge_neighbors(self, pScoresDictionary, pMergeThreshold=1000):
        return _impl.Viewpoint_merge_neighbors(self, pScoresDictionary, pMergeThreshold)

    def readSignificantRegionsFile(self, pFilePath, pInternalIdentifierTriplet, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint):
        return _impl.Viewpoint_readSignificantRegionsFile(self, pFilePath, pInternalIdentifierTriplet, pViewpointIndexStart, pViewpointIndexEnd, pResolution, pRange, pViewpoint)

    def readTargetHDFFile(self, pFile):
        return _impl.Viewpoint_readTargetHDFFile(self, pFile)

    def readAggregatedFileHDF(self, pAggregatedFileName, pInternalPath):
        return _impl.Viewpoint_readAggregatedFileHDF(self, pAggregatedFileName, pInternalPath)

    def readDifferentialFile(self, pFilePath, pQuadruple):
        return _impl.Viewpoint_readDifferentialFile(self, pFilePath, pQuadruple)
