import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.textpath
import matplotlib.colors
import matplotlib.gridspec
import matplotlib.cm
import mpl_toolkits.axisartist as axisartist
from matplotlib.patches import Rectangle

import scipy.sparse
from collections import OrderedDict

from bx.intervals.intersection import IntervalTree, Interval

import hicexplorer.HiCMatrix as HiCMatrix
import hicexplorer.utilities


DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 3  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
DEFAULT_WIDTH_RATIOS = (0.95, 0.05)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.12, 'top': 0.9}

class MultiDict(OrderedDict):
    """
    Class to allow identically named
    sections in configuration file
    by appending the section number as
    for example:
    1. section name
    """
    _unique = 0

    def __setitem__(self, key, val):
        if isinstance(val, OrderedDict):
            self._unique += 1
            key = "{}. [{}]".format(str(self._unique), key)
        OrderedDict.__setitem__(self, key, val)

class PlotTracks(object):

    def __init__(self, tracks_file, fig_width=DEFAULT_FIGURE_WIDTH,
                 fig_height=None, fontsize=None, dpi=None):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.vlines_intval_tree = None
        start = self.print_elapsed(None)
        self.parse_tracks(tracks_file)
        if fontsize:
            fontsize = fontsize
        else:
            fontsize = float(fig_width) * 0.4

        font = {'size': fontsize}
        matplotlib.rc('font', **font)
        #import ipdb;ipdb.set_trace()
        # initialize each track
        self.track_obj_list = []
        for idx, properties in enumerate(self.track_list):
            if 'spacer' in properties:
                continue
            if 'x-axis' in properties:
                self.track_obj_list.append(PlotXAxis(properties))
                continue
            if properties['file_type'] == 'bedgraph':
                self.track_obj_list.append(PlotBedGraph(properties))

            elif properties['file_type'] == 'bigwig':
                self.track_obj_list.append(PlotBigWig(properties))

            elif properties['file_type'] == 'bedgraph_matrix':
                self.track_obj_list.append(PlotBedGraphMatrix(properties))

            elif properties['file_type'] == 'hic_matrix':
                self.track_obj_list.append(PlotHiCMatrix(properties))

            elif properties['file_type'] == 'bed':
                self.track_obj_list.append(PlotBed(properties))

        print "time initializing tracks"
        start = self.print_elapsed(start)

    def get_tracks_height(self, start, end):
        # prepare layout based on the tracks given.
        # The main purpose of the following loop is
        # to get the height of each of the tracks
        track_height = []
        for track_dict in self.track_list:
            height = DEFAULT_TRACK_HEIGHT

            if 'width' in track_dict:
                height = track_dict['width']

            # compute the height of a Hi-C track
            # based on the depth such that the
            # resulting plot appears proportional
            #
            #      /|\
            #     / | \
            #    /  |d \   d is the depth that we want to be proportional
            #   /   |   \  when plotted in the figure
            # ------------------
            #   region len
            #
            # d (in cm) =  depth (in bp) * width (in cm) / region len (in bp)

            elif 'depth' in track_dict and track_dict['file_type'] == 'hic_matrix':
                # to compute the actual width of the figure the margins and the region
                # set for the legends have to be considered
                # DEFAULT_MARGINS[1] - DEFAULT_MARGINS[0] is the proportion of plotting area

                hic_width = self.fig_width * (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) * DEFAULT_WIDTH_RATIOS[0]
                scale_factor = 0.6  # the scale factor is to obtain a pleasing result.
                height = scale_factor * track_dict['depth'] * hic_width / (end - start)

            track_height.append(height)

        return track_height

    def plot(self, file_name, chrom, start, end, title=None):
        track_height = self.get_tracks_height(start, end)
        if self.fig_height:
            fig_height = self.fig_height
        else:
            fig_height = sum(track_height)

        sys.stderr.write("Figure size in cm is {} x {}. Dpi is set to {}\n".format(self.fig_width,
                                                                                   fig_height, self.dpi))
        fig = plt.figure(figsize=self.cm2inch(self.fig_width, fig_height))
        fig.suptitle(title)

        grids = matplotlib.gridspec.GridSpec(len(track_height), 2,
                                             height_ratios=track_height,
                                             width_ratios=[1, 0.05])
        axis_list = []
        for idx, track in enumerate(self.track_obj_list):
            axis = axisartist.Subplot(fig, grids[idx, 0])
            fig.add_subplot(axis)
            axis.axis[:].set_visible(False)
            # to make the background transparent
            axis.patch.set_visible(False)
            label_axis = plt.subplot(grids[idx, 1])
            label_axis.set_axis_off()
            track.plot(axis, label_axis, chrom, start, end)
            axis_list.append(axis)

        if self.vlines_intval_tree:
            self.plot_vlines(axis_list, chrom, start, end)

        fig.subplots_adjust(wspace=0, hspace=0.1,
                            left=DEFAULT_MARGINS['left'],
                            right=DEFAULT_MARGINS['right'],
                            bottom=DEFAULT_MARGINS['bottom'],
                            top=DEFAULT_MARGINS['top'])

        print "time before saving"
        start = self.print_elapsed(start)
        print "saving {}".format(file_name)

        fig.savefig(file_name, dpi=self.dpi)
        print "time saving "
        start = self.print_elapsed(start)

    def plot_vlines(self, axis_list, chrom_region, start_region, end_region):
        """
        Plots dotted lines from the top of the first plot to the bottom
        of the last plot at the specified positions.

        :param axis_list: list of plotted axis
        :return: None
        """
        vlines_list = []
        for region in self.vlines_intval_tree[chrom_region].find(start_region - 10000,
                                                                 end_region + 10000):
            vlines_list.append(region.start)

        from matplotlib.patches import ConnectionPatch
        a_ymax = axis_list[0].get_ylim()[1]
        b_ymin = axis_list[-1].get_ylim()[0]

        for start_pos in vlines_list:
            con = ConnectionPatch(xyA=(start_pos, a_ymax),
                                  xyB=(start_pos, b_ymin),
                                  coordsA="data", coordsB="data",
                                  axesA=axis_list[0],
                                  axesB=axis_list[-1],
                                  arrowstyle="-",
                                  linestyle='dashed',
                                  linewidth=0.5,
                                  zorder=100)
            axis_list[0].add_artist(con)

    def parse_tracks(self, tracks_file):
        """
        Parses a configuration file

        :param tracks_file: file path containing the track configuration
        :return: array of dictionaries and vlines_file. One dictionary per track
        """
        from ConfigParser import SafeConfigParser
        from ast import literal_eval
        parser = SafeConfigParser(None, MultiDict)
        parser.readfp(open(tracks_file, 'r'))

        track_list = []
        vlines_file = None
        for section_name in parser.sections():
            track_options = dict({"section_name": section_name})
            if section_name.endswith('[spacer]'):
                track_options['spacer'] = True
            elif section_name.endswith('[x-axis]'):
                track_options['x-axis'] = True
            for name, value in parser.items(section_name):
                if name in ['max_value', 'min_value', 'depth', 'width'] and value != 'auto':
                    track_options[name] = literal_eval(value)
                else:
                    track_options[name] = value

            if 'type' in track_options and track_options['type'] == 'vlines':
                vlines_file = track_options['file']
            else:
                track_list.append(track_options)

        for track_dict in track_list:
            warn = None
            if 'file' in track_dict and track_dict['file'] != '':
                self.check_file_exists(track_dict)
                if 'file_type' not in track_dict:
                    track_dict['file_type'] = self.guess_filetype(track_dict)

                #  set some default values
                if 'title' not in track_dict:
                    warn = "\ntitle not set for 'section {}'\n".format(track_dict['section_name'])
                    track_dict['title'] = ''
                if warn:
                    sys.stderr.write(warn)

        self.track_list = track_list
        if vlines_file:
            self.vlines_intval_tree = file_to_intervaltree(vlines_file)


    def check_file_exists(self, track_dict):
        """
        Checks if a file or list of files exists
        :param track_dict: dictionary of track values. Should contain
                            a 'file' key containing the path of the file
                            or files to be checked separated by space
                            For example: file1 file2 file3
        :return: None
        """
        file_names = [x for x in track_dict['file'].split(" ") if x != '']
        for file_name in file_names:
            try:
                open(file_name, 'r').close()
            except IOError:
                sys.stderr.write("\n*ERROR*\nFile in section [{}] "
                                 "not found:\n{}\n\n".format(track_dict['section_name'],
                                                             file_name))
                exit(1)

    @staticmethod
    def guess_filetype(track_dict):
        """

        :param track_dict: dictionary of track values with the 'file' key
                    containing a string path of the file or files. Only the ending
                     of the last file is used in case when there are more files
        :return: string file type detected
        """
        file = track_dict['file'].strip()
        file_type = None

        if file.endswith(".bed"):
            file_type = 'bed'
        elif file.endswith(".npz"):
            file_type = 'hic_matrix'
        elif file.endswith(".bw"):
            file_type = 'bigwig'
        elif file.endswith(".bg"):
            file_type = 'bedgraph'
        elif file.endswith(".bm"):
            file_type = 'bedgraph_matrix'
        else:
            exit("Section [{}]: can not identify file type. Please specify "
                 "the file_type for {}".format(track_dict['section_name'],
                                               file))
        return file_type

    @staticmethod
    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i/inch for i in tupl[0])
        else:
            return tuple(i/inch for i in tupl)


    @staticmethod
    def print_elapsed(start):
        import time
        if start:
            print time.time() - start
        return time.time()


def file_to_intervaltree(file_name):
    """
    converts a BED like file into a bx python interval tree
    :param file_name: string file name
    :return: interval tree dictionary. They key is the chromosome/contig name and the
    value is an IntervalTree. Each of the intervals have as 'value' the fields[3:] if any.
    """
    # iterate over a BED like file
    # saving the data into an interval tree
    # for quick retrieval
    file_h = open(file_name, 'r')
    line_number = 0
    valid_intervals = 0
    prev_chrom = None
    prev_start = -1
    prev_line = None
    interval_tree = {}
    for line in file_h.readlines():
        line_number += 1
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = None
        start = None
        end = None
        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(line_number, detail)
            exit(msg)

        try:
            start = int(start)
        except ValueError as detail:
            msg = "Error reading line: {}. The start field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            exit(msg)

        try:
            end = int(end)
        except ValueError as detail:
            msg = "Error reading line: {}. The end field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            exit(msg)

        if prev_chrom == chrom:
            assert prev_start <= start, \
                "Bed file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

        if chrom not in interval_tree:
            interval_tree[chrom] = IntervalTree()

        value = None
        if fields > 3:
            value = fields[3:]

        assert end > start, \
                "Start position larger or equal than end for line\n{} ".format(line)

        interval_tree[chrom].insert_interval(Interval(start, end, value=value))
        valid_intervals += 1

    if valid_intervals == 0:
        sys.stderr.write("No valid intervals were found in file {}".format(file_name))

    file_h.close()
    return interval_tree


class TrackPlot(object):

    def __init__(self, properties_dict):
        self.properties = properties_dict

class PlotBedGraph(TrackPlot):

    def __init__(self, properties_dict):
        #super(self.__class__, self).__init__(*args, **kwargs)
        self.properties = properties_dict
        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BEDGRAPH_COLOR
        self.interval_tree = file_to_intervaltree(self.properties['file'])

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        score_list = []
        pos_list = []

        for region in self.interval_tree[chrom_region].find(start_region - 10000,
                                                           end_region + 10000):
            score_list.append(float(region.value[0]))
            pos_list.append(region.start + (region.end - region.start)/2)

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BEDGRAPH_COLOR

        if 'extra' in self.properties and \
                self.properties['extra'][0] == '4C':
            # draw a vertical line for each fragment region center
            self.ax.fill_between(pos_list, score_list,
                            facecolor=self.properties['color'],
                            edgecolor='none')
            self.ax.vlines(pos_list, [0], score_list, color='olive', linewidth=0.5)
            self.ax.plot(pos_list, score_list, '-', color='slateblue', linewidth=0.7)
        else:
            try:
                self.ax.fill_between(pos_list, score_list, facecolor=self.properties['color'])
            except ValueError:
                exit("Invalid color {} for {}".format(self.properties['color'], self.properties['file']))

        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.ax.set_xlim(start_region, end_region)

        ymin, ymax = self.ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if float(ymax) % 1 == 0:
            ymax_print = int(ymax)
        else:
            ymax_print = "{:.1f}".format(ymax)
        self.ax.set_ylim(ymin, ymax)
        ydelta = ymax - ymin
        small_x = 0.01 * (end_region - start_region)

        if 'show data range' in self.properties and \
                self.properties['show data range'] == 'no':
            pass
        else:
            # by default show the data range
            self.ax.text(start_region-small_x, ymax - ydelta * 0.2,
                    "[{}-{}]".format(ymin, ymax_print),
                    horizontalalignment='left', size='small',
                    verticalalignment='bottom')

        self.label_ax.text(0.15, 0, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='bottom', transform=self.label_ax.transAxes)

class PlotBedGraphMatrix(PlotBedGraph):

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """

        start_pos = []
        matrix_rows = []
        for region in self.interval_tree[chrom_region].find(start_region - 10000,
                                                           end_region + 10000):
            start_pos.append(region.start)
            values = map(float, region.values)
            matrix_rows.append(values)

        matrix = np.vstack(matrix_rows).T
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            matrix = np.flipud(matrix)

        vmin = None
        vmax = None
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            vmax = self.properties['max_value']

        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            vmin = self.properties['min_value']

        if 'type' in self.properties and self.properties['type'] == 'lines':
            for row in matrix:
                self.ax.plot(start_pos, row)
            self.ax.plot(start_pos, matrix.mean(axis=0), "--")
        else:
            x, y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
            img = self.ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax, shading='gouraud')
            img.set_rasterized(True)
        self.ax.set_xlim(start_region, end_region)
        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)

        self.label_ax.text(0.15, 0, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='bottom', transform=self.label_ax.transAxes)


class PlotBigWig(TrackPlot):

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        from bx.bbi.bigwig_file import BigWigFile
        self.bw = BigWigFile(open(self.properties['file'], 'r'))
        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BIGWIG_COLOR


    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        formated_region = "{}:{}-{}".format(chrom_region, start_region, end_region)
        print "plotting {}".format(self.properties['file'])
        # compute the score in bins of 10000 SLOW
    #    score = np.array([self.bw.query(region[0], x, x+10000,1)[0]['mean']
    #                      for x in range(region[1], region[2], 10000)])

        num_bins = 700
        if 'number of bins' in self.properties:
            try:
                num_bins = int(self.properties['number of bins'])
            except TypeError:
                exit("'number of bins' value: {} for bigwig file {} "
                     "is not valid.".format(self.properties['number of bins'],
                                            self.properties['file']))

        if end_region - start_region < 2e6:
            scores = self.bw.get_as_array(chrom_region, start_region, end_region)
            if scores is None:
                # usually bw returns None when the chromosome
                # is not found. So, we try either
                # removing or appending 'chr'
                if chrom_region.startswith('chr'):
                    scores = self.bw.get_as_array(chrom_region[3:] + chrom_region, start_region, end_region)
                else:
                    scores = self.bw.get_as_array('chr' + chrom_region, start_region, end_region)

                if scores is None:
                    exit("Can not read region {} from bigwig file:\n\n"
                         "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                         "and that the region is valid".format(formated_region, self.properties['file']))

            if 'nans to zeros' in self.properties and self.properties['nans to zeros'] is True:
                scores[np.isnan(scores)] = 0

            scores = np.ma.masked_invalid(scores)

            lins = np.linspace(0, len(scores), num_bins).astype(int)
            scores_per_bin = [np.mean(scores[lins[x]:lins[x+1]]) for x in range(len(lins)-1)]
            _x = lins + start_region
            x_values = [float(_x[x] + _x[x+1])/2 for x in range(len(lins)-1)]
            self.ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                                 color=self.properties['color'],
                                 facecolor=self.properties['color'])

        else:
            # this method produces shifted regions. It is not clear to me why this happens.
            # Thus I only activate the faster but shifted method for large regions
            # when the previous method would be to slow
            scores = self.bw.query(chrom_region, start_region, end_region, num_bins)
            if scores is None:

                if chrom_region.startswith('chr'):
                    scores = self.bw.query(chrom_region[3:] + chrom_region, start_region, end_region)
                else:
                    scores = self.bw.query('chr' + chrom_region, start_region, end_region)

                if scores is None:
                    exit("Can not read region {} from bigwig file:\n\n"
                         "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                         "and that the region is valid".format(formated_region, self.properties['file']))

            scores = [x['mean'] for x in scores]
            x_values = np.linspace(start_region, end_region, num_bins)
            self.ax.fill_between(x_values, scores, linewidth=0.1,
                                 color=self.properties['color'],
                                 facecolor=self.properties['color'], zorder=1)

        self.ax.set_xlim(start_region, end_region)
        ymin, ymax = self.ax.get_ylim()
        if 'max_value' in self.properties and ['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.ax.set_ylim(ymax, ymin)
        else:
            self.ax.set_ylim(ymin, ymax)

    #    self.ax.set_yticks([ymax])
        ydelta = ymax - ymin

    #    self.ax.set_yticklabels(["{}-{}".format(int(ymin), int(ymax))], size='large')
        # set min max
        if float(ymax) % 1 == 0:
            ymax_print = int(ymax)
        else:
            ymax_print = "{:.1f}".format(ymax)
        small_x = 0.01 * (end_region - start_region)
        if 'show data range' in self.properties and self.properties['show data range'] == 'no':
            pass
        else:
            # by default show the data range
            self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                         "[{}-{}]".format(int(ymin), ymax_print),
                         horizontalalignment='left', size='small',
                         verticalalignment='bottom')

        """
        self.ax.text(region_end, ymax - ydelta * 0.2, self.properties['title'],
                horizontalalignment='right', size='large',
                verticalalignment='bottom')

        """
        self.label_ax.text(0.15, 0, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='bottom')

        return self.ax

class PlotHiCMatrix(TrackPlot):

    def __init__(self, properties_dict):
        # to avoid the color bar to span all the
        # width of the axis I pass two axes
        # to plot_matrix
        self.properties = properties_dict

        self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'])

        if 'show_masked_bins' in self.properties and self.properties['show_masked_bins'] == 'yes':
            pass
        else:
            self.hic_ma.maskBins(self.hic_ma.nan_bins)

        new_intervals = hicexplorer.utilities.enlarge_bins(self.hic_ma.cut_intervals)
        self.hic_ma.interval_trees, self.hic_ma.chrBinBoundaries = \
            self.hic_ma.intervalListToIntervalTree(new_intervals)

        self.hic_ma.cut_intervals = new_intervals

        # select only the upper triangle of the matrix
        self.hic_ma.matrix = scipy.sparse.triu(self.hic_ma.matrix, k=0, format='csr')

        # fill the main diagonal, otherwise it looks
        # not so good. The main diagonal is filled
        # with an array containing the max value found
        # in the matrix
        if sum(self.hic_ma.matrix.diagonal()) == 0:
            sys.stderr.write("Filling main diagonal because is empty and " \
                "otherwise it looks bad...\n")
            max_value = self.hic_ma.matrix.data.max()
            main_diagonal = scipy.sparse.dia_matrix(([max_value]*self.hic_ma.matrix.shape[0], [0]),
                                                    shape=self.hic_ma.matrix.shape)
            self.hic_ma.matrix = self.hic_ma.matrix + main_diagonal

        self.plot_inverted = False
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.plot_inverted = True

        self.norm = None

        if 'colormap' not in self.properties:
            self.properties['colormap'] = DEFAULT_MATRIX_COLORMAP

        self.cmap = matplotlib.cm.get_cmap(self.properties['colormap'])
        self.cmap.set_bad('white')

        self.cmap.set_bad('black')

        if 'boundaries_file' in self.properties:
            self.boundaries_obj = PlotBoundaries({'file':self.properties['boundaries_file']})

    def plot(self, ax, label_ax, chrom, region_start, region_end):
        import copy
        self.cbar_ax = copy.copy(label_ax)
        self.label_ax = label_ax
        self.label_ax.set_axis_off()
        self.ax = ax

        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges
    
        # get bin id of start and end of region in given chromosome
        chr_start_id, chr_end_id = self.hic_ma.getChrBinRange(chrom)
        chr_start = self.hic_ma.cut_intervals[chr_start_id][1]
        chr_end = self.hic_ma.cut_intervals[chr_end_id-1][1]
        start_bp = max(chr_start, region_start - self.properties['depth'])
        end_bp = min(chr_end, region_end + self.properties['depth'])
    
        idx, start_pos = zip(*[(idx, x[1]) for idx, x in
                               enumerate(self.hic_ma.cut_intervals)
                               if x[0] == chrom and x[1] >= start_bp and x[2] <= end_bp])
    
        idx = idx[0:-1]
        # select only relevant matrix part
        matrix = self.hic_ma.matrix[idx, :][:, idx]
        matrix = np.asarray(matrix.todense().astype(float))

        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                matrix += 1
                self.norm = matplotlib.colors.LogNorm()

            elif self.properties['transform'] == '-log':
                mask = matrix == 0
                matrix[mask] = matrix[mask == False].min()
                matrix = -1 * np.log(matrix)

        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            vmax = self.properties['max_value']
    
        else:
            # try to use a 'aesthetically pleasant' max value
            vmax = np.percentile(matrix.diagonal(1), 80)
    
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            vmin = self.properties['min_value']
        else:
            bin_size = self.hic_ma.getBinSize()
            depth_bins = int(self.properties['depth'] / bin_size)
            vmin = np.median(matrix.diagonal(depth_bins))
    
        sys.stderr.write("setting min, max values for track {} to: {}, {}\n".format(self.properties['section_name'],
                                                                                    vmin, vmax))
        img = self.pcolormesh_45deg(matrix, start_pos, vmax=vmax, vmin=vmin)
        img.set_rasterized(True)
        if self.plot_inverted:
            self.ax.set_ylim(self.properties['depth'], 0)
        else:
            self.ax.set_ylim(0, self.properties['depth'])
    
        # ##plot boundaries
        # if a boundaries file is given, plot the
        # tad boundaries as line delineating the TAD triangles
        if 'boundaries_file' in self.properties:
            self.boundaries_obj.plot(ax, label_ax, chrom, region_start, region_end)

        self.ax.set_xlim(region_start, region_end)
        if 'x labels' in self.properties and self.properties['x labels'] != 'no':
            ticks = self.ax.get_xticks()
            labels = ["{:.2f}".format((x / 1e6))
                      for x in ticks]
            labels[-1] += "Mbp"
            self.ax.get_xaxis().set_tick_params(
                which='both',
                bottom='on',
                top='off',
                direction='out')
    
            self.ax.set_xticklabels(labels)
        else:
            self.ax.get_xaxis().set_tick_params(
                which='both',
                bottom='off',
                top='off',
                direction='out')
            self.ax.axes.get_xaxis().set_visible(False)

        self.ax.set_frame_on(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.cbar_ax.patch.set_alpha(0.0)
        cobar = plt.colorbar(img, ax=self.cbar_ax, fraction=0.95)
        cobar.solids.set_edgecolor("face")
        self.label_ax.text(0.3, 0.0, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='bottom', transform=self.label_ax.transAxes)

    def pcolormesh_45deg(self, matrix_c, start_pos_vector, vmin=None,
                         vmax=None):
        """
        Turns the matrix 45 degrees and adjusts the
        bins to match the actual start end positions.
        """
        import itertools
        # code for rotating the image 45 degrees
        n = matrix_c.shape[0]
        # create rotation/scaling matrix
        t = np.array([[1, 0.5], [-1, 0.5]])
        # create coordinate matrix and transform it
        matrix_a = np.dot(np.array([(i[1], i[0])
                                    for i in itertools.product(start_pos_vector[::-1],
                                                               start_pos_vector)]), t)
        # this is to convert the indices into bp ranges
        x = matrix_a[:, 1].reshape(n+1, n+1)
        y = matrix_a[:, 0].reshape(n+1, n+1)
        # plot
        im = self.ax.pcolormesh(x, y, np.flipud(matrix_c),
                                vmin=vmin, vmax=vmax, cmap=self.cmap, norm=self.norm)
        return im

class PlotXAxis(TrackPlot):

    def plot(self, ax, label_axis, chrom_region, region_start, region_end):
        ax.set_xlim(region_start, region_end)
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f} kb".format((x / 1e3))
                      for x in ticks]

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f} kb".format((x / 1e3))
                      for x in ticks]
        else:
            labels = ["{:,.1f} Mbp".format((x / 1e6))
                      for x in ticks]
            # labels[-1] += "Mbp"

        ax.axis["x"] = ax.new_floating_axis(0, 0.5)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        if 'fontsize' in self.properties:
            ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"].set_axis_direction("top")

class PlotBoundaries(PlotBedGraph):
    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        intervals = []
        for region in self.interval_tree[chrom_region].find(start_region - 10000,
                                                            end_region + 10000):
            intervals.append((region.start, region.end))

        start, end = zip(*intervals)
        start_list = np.array(start)
        end = np.array(end)
        # check if intervals are consecutive or 1bp positions demarcating the boundaries
        if np.any(end - start_list > 1):
            # intervals are consecutive, but only the boundaries are need.
            start_list = end

        prev_start = -1
        x = []
        y = []
        for start in start_list:
            if prev_start is None:
                # draw only half a triangle
                length = start - prev_start
                x1 = prev_start
                y1 = length
                x2 = start
                y2 = 0
                x.extend([x1, x2])
                y.extend([y1, y2])
            else:
                x1 = prev_start
                x2 = x1 + (start - prev_start) / 2
                x3 = start
                y1 = 0
                y2 = (start - prev_start)
                y3 = 0
                x.extend([x1, x2, x3])
                y.extend([y1, y2, y3])

            prev_start = start
        ax.plot(x, y,  color='black')

class PlotBed(TrackPlot):

    def __init__(self, *args, **kwarg):
        super(PlotBed, self).__init__(*args, **kwarg)
        import readBed
        bed_file_h = readBed.ReadBed(open(self.properties['file'], 'r'))
        self.bed_type = bed_file_h.file_type
        valid_intervals = 0
        self.interval_tree = {}
        for bed in bed_file_h:
            if bed.chromosome not in self.interval_tree:
                self.interval_tree[bed.chromosome] = IntervalTree()

            self.interval_tree[bed.chromosome].insert_interval(Interval(bed.start, bed.end, value=bed))
            valid_intervals += 1

        if valid_intervals == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(self.properties['file_name']))

        from matplotlib import font_manager
        if 'fontsize' in self.properties:
            self.fp = font_manager.FontProperties(size=self.properties['fontsize'])
        else:
            self.fp = font_manager.FontProperties()

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BED_COLOR
        if 'border_color' not in self.properties:
            self.properties['border_color'] = 'black'
        if 'labels' not in self.properties:
            self.properties['labels'] = 'on'
        if 'style' not in self.properties:
            self.properties['style'] = 'flybase'

        # to improve the visualization of the genes
        # it is good to have an estimation of the label
        # length. In the following code I try to get
        # the length of a 'w'.
        if self.properties['labels'] == 'on':
            text_path = matplotlib.textpath.TextPath((0, 0), 'w', prop=self.fp)
            self.len_w = text_path.get_extents().width * 1000
            #self.len_w = text_path.get_extents().width * 300
        else:
            self.len_w = 1

        self.colormap = None
        # check if the color given is a color map
        color_options = [m for m in matplotlib.cm.datad]
        if self.properties['color'] in color_options:
            norm = matplotlib.colors.Normalize(vmin=self.properties['min_value'],
                                               vmax=self.properties['max_value'])
            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def get_y_pos(self, bed):
        """
        The y_pos is set such that regions to be plotted do not overlap. To override this
        the properties['collapsed'] needs to be set.
        :return: int y position
        """

        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if 'display' in self.properties and self.properties['display'] == 'interlaced':
            ypos = 100 if self.counter % 2 == 0 else 1

        elif 'display' in self.properties and  self.properties['display'] == 'collapsed':
            ypos = 0

        else:
            # 1. check for the number of other intervals that overlap
            #    with the given interval
            #
            #  1=========       4=========
            #       2=========
            #         3============
            #
            # for 1, min_free_row is 0
            # for 2, min_free_row is 1
            # for 3, min_free_row is 2
            # for 4, min_free_row is 0

            # check for overlapping features
            match = self.region_intervals.find(bed.start, bed.end)
            if len(match) == 0:
                min_free_row = 0
            else:
                rows_used = np.zeros(self.max_num_row + 2)
                for x in match:
                    # in the Interval, the 'value' field stores the row
                    # in which the bed region was printed
                    rows_used[x.value] = 1
                min_free_row = min(np.flatnonzero(rows_used == 0))

            # check if the label may be larger than the interval, if this is the case
            # set the interval to match the expected label length
            print min_free_row, len(bed.name) * self.len_w, bed.end - bed.start
            if self.properties['labels'] == 'on' and \
                        bed.end - bed.start < len(bed.name) * self.len_w:
                print "long label"
                self.region_intervals.add_interval(Interval(bed.start,
                                                            bed.start + (len(bed.name) * self.len_w),
                                                            min_free_row))
            else:
                self.region_intervals.add_interval(Interval(bed.start,
                                                            bed.end + 2 * self.small_relative,
                                                            min_free_row))
            if min_free_row > self.max_num_row:
                self.max_num_row = min_free_row

            if self.properties['labels'] == 'off':
                scale_factor = 230
            else:
                #scale_factor = self.len_w * 2
                scale_factor = 330

            ypos = min_free_row * scale_factor

        print ypos
        return ypos

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        from matplotlib.patches import Rectangle
        self.counter = 0
        self.small_relative = 0.005 * (end_region-start_region)
        self.max_num_row = 1
        self.region_intervals = IntervalTree()

        ax.set_frame_on(False)
        for region in self.interval_tree[chrom_region].find(start_region, end_region):
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED10
            bed with rbg at end 
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44            

            BED6
            bed with rbg at end
            chr2L   0       70000   ID_5    0.26864549832   .
            """
            self.counter += 1
            bed = region.value
            if self.colormap:
                # translate value field (in the example above is 0 or 0.2686...) into a color
                rgb = self.colormap.to_rgba(bed.score)
                edgecolor = self.colormap.to_rgba(bed.score)
            else:
                rgb = self.properties['color']
                edgecolor = self.properties['border_color']

            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type == 'bed10' and len(bed.rgb) == 3:
                try:
                    rgb = [float(x)/255 for x in rgb]
                    edgecolor = self.properties['color']
                except IndexError:
                    pass

            ypos = self.get_y_pos(bed)
            if self.bed_type == 'bed12':
                if self.properties['style'] == 'flybase':
                    self.draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor)
                else:
                    self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor)


            else:
                self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if 'labels' in self.properties and self.properties['labels'] == 'off':
                pass
            else:
                ax.text(bed.start, ypos + 125, bed.name, horizontalalignment='left',
                        verticalalignment='top', fontproperties=self.fp)

        if self.counter == 0:
            sys.stderr.write("*Warning* No intervals were found for file {} \n"
                             "in section '{}' for the interval plotted ({}:{}-{}).\n"
                             "".format(self.properties['file'],
                                       self.properties['section_name'],
                                       chrom_region, start_region, end_region))

        ax.set_ylim((self.max_num_row + 1) * 330, -25)

        if 'display' in self.properties:
            if self.properties['display'] == 'domain':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

        ax.set_xlim(start_region, end_region)

        label_ax.text(0.15, 1.0, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='top', transform=label_ax.transAxes)

    def draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws an interval with direction (if given)
        """

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end-bed.start, 100, edgecolor=edgecolor,
                                   facecolor=rgb, linewidth=0.5))
        else:
            vertices = self._draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=0.5))

    def draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene using different styles
        """
        from matplotlib.patches import Polygon
        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos+50, ypos+50], 'black', linewidth=0.5, zorder=-1)

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x0 < bed.thick_start < x1:
                positions.append((x0, bed.thick_start, 'UTR'))
                positions.append((bed.thick_start, x1, 'coding'))

            elif x0 < bed.thick_end < x1:
                positions.append((x0, bed.thick_end, 'coding'))
                positions.append((bed.thick_end, x1, 'UTR'))

            else:
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    type = 'UTR'
                else:
                    type = 'coding'

                positions.append((x0, x1, type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = 'grey'
        else:
            _rgb = rgb

        vertices = self._draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=0.5))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = 'grey'
            else:
                _rgb = rgb
            vertices = [(start_pos, ypos), (start_pos, ypos + 100),
                        (end_pos, ypos + 100), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=0.5))

#            ax.add_patch(Rectangle((start_pos, ypos), end_pos-start_pos, 100, edgecolor=edgecolor,
#                                   facecolor=_rgb, linewidth=0.5))

    def _draw_arrow(self, ax, start, end, strand, ypos):
        """
        Draws a filled arrow
        :param ax:
        :param start:
        :param end:
        :param strand:
        :param ypos:
        :param rgb:
        :return: None
        """
        from matplotlib.patches import Polygon
        if strand == '+':
            x0 = start
            x1 = end #- self.small_relative
            y0 = ypos
            y1 = ypos + 100
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            -----------------\
            ---------------- /

            """

            vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + self.small_relative, y0 + 50), (x1, y0)]

        else:
            x0 = start #+ self.small_relative
            x1 = end
            y0 = ypos
            y1 = ypos + 100
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            /-----------------
            \-----------------
            """
            vertices = [(x0, y0), (x0 - self.small_relative, y0 + 50), (x0, y1), (x1, y1), (x1, y0)]

        return vertices


    def draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor):
            """
            draws a gene like in flybase gbrowse.
            """
            from matplotlib.patches import Polygon

            if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
                self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
                return

            # draw 'backbone', a line from the start until the end of the gene
            ax.plot([bed.start, bed.end], [ypos+50, ypos+50], 'black', linewidth=0.5, zorder=-1)

            for idx in range(0, bed.block_count):
                x0 = bed.start + bed.block_starts[idx]
                x1 = x0 + bed.block_sizes[idx]
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    y0 = ypos + 25
                    y1 = ypos + 75
                else:
                    y0 = ypos
                    y1 = ypos + 100

                if x0 < bed.thick_start < x1:
                    vertices = ([(x0, ypos+25), (x0, ypos+75), (bed.thick_start, ypos+75), (bed.thick_start, ypos+100),
                                 (bed.thick_start, ypos+100), (x1, ypos+100), (x1, ypos),
                                 (bed.thick_start, ypos), (bed.thick_start, ypos+25)])

                elif x0 < bed.thick_end < x1:
                    vertices = ([(x0, ypos), (x0, ypos+100), (bed.thick_end, ypos+100), (bed.thick_end, ypos+75),
                                 (x1, ypos+75), (x1, ypos+25), (bed.thick_end, ypos+25), (bed.thick_end, ypos)])
                else:
                    vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

                ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                     linewidth=0.1,
                                     edgecolor='none',
                                     facecolor=rgb))

                if idx < bed.block_count - 1:
                    # plot small arrows using the character '<' or '>' over the back bone
                    intron_length = bed.block_starts[idx+1] - (bed.block_starts[idx] + bed.block_sizes[idx])
                    marker = 5 if bed.strand == '+' else 4
                    if intron_length > 3 * self.small_relative:
                        pos = np.arange(x1 + 1 * self.small_relative,
                                        x1 + intron_length + self.small_relative, int(2 * self.small_relative))
                        ax.plot(pos, np.zeros(len(pos)) + ypos + 50, '.', marker=marker,
                                fillstyle='none', color='blue', markersize=3)

                    elif intron_length > self.small_relative:
                        intron_center = x1 + int(intron_length)/2
                        ax.plot([intron_center], [ypos+50], '.', marker=5,
                                fillstyle='none', color='blue', markersize=3)

