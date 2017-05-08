# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.textpath
import matplotlib.colors
import matplotlib.gridspec
import matplotlib.cm
import mpl_toolkits.axisartist as axisartist
from matplotlib.patches import Rectangle
import textwrap
import os.path


import scipy.sparse
from collections import OrderedDict

from intervaltree import IntervalTree, Interval

import hicexplorer.HiCMatrix as HiCMatrix
import hicexplorer.utilities

reload(sys)
sys.setdefaultencoding('utf-8')

DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 3  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
DEFAULT_WIDTH_RATIOS = (0.93, 0.07)
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
                 fig_height=None, fontsize=None, dpi=None, track_label_width=None):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.vlines_intval_tree = None
        self.vlines_properties = None
        start = self.print_elapsed(None)
        self.parse_tracks(tracks_file)
        if fontsize:
            fontsize = fontsize
        else:
            fontsize = float(fig_width) * 0.3
        if track_label_width is None:
            self.width_ratios = DEFAULT_WIDTH_RATIOS
        else:
            self.width_ratios = (1 - track_label_width, track_label_width)

        font = {'size': fontsize}
        matplotlib.rc('font', **font)
        # initialize each track
        self.track_obj_list = []
        for idx, properties in enumerate(self.track_list):
            if 'spacer' in properties:
                self.track_obj_list.append(PlotSpacer(properties))
                continue
            elif 'x-axis' in properties:
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

            elif properties['file_type'] == 'links':
                self.track_obj_list.append(PlotArcs(properties))

            elif properties['file_type'] == 'boundaries':
                self.track_obj_list.append(PlotBoundaries(properties))

            if 'title' in properties:
                # adjust titles that are too long
                # if the track label space is small
                if track_label_width < 0.1:
                    properties['title'] = textwrap.fill(properties['title'].encode("UTF-8"), 12)
                else:
                    properties['title'] = textwrap.fill(properties['title'].encode("UTF-8"), 50)

        print "time initializing track(s):"
        self.print_elapsed(start)

    def get_tracks_height(self, start, end):
        # prepare layout based on the tracks given.
        # The main purpose of the following loop is
        # to get the height of each of the tracks
        # because for the hi-C the height is variable with respect
        # to the range being plotted, the function is called
        # when each plot is going to be printed.
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

                hic_width = \
                    self.fig_width * (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) * self.width_ratios[0]
                scale_factor = 0.6  # the scale factor is to obtain a 'pleasing' result.
                depth = min(track_dict['depth'], (end - start))

                height = scale_factor * depth * hic_width / (end - start)
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
        if title:
            fig.suptitle(title)

        grids = matplotlib.gridspec.GridSpec(len(track_height), 2,
                                             height_ratios=track_height,
                                             width_ratios=self.width_ratios)
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

        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=DEFAULT_MARGINS['left'],
                            right=DEFAULT_MARGINS['right'],
                            bottom=DEFAULT_MARGINS['bottom'],
                            top=DEFAULT_MARGINS['top'])

        fig.savefig(file_name, dpi=self.dpi, transparent=False)
        return fig.get_size_inches()

    def plot_vlines(self, axis_list, chrom_region, start_region, end_region):
        """
        Plots dotted lines from the top of the first plot to the bottom
        of the last plot at the specified positions.

        :param axis_list: list of plotted axis
        :param chrom_region chromosome name
        :param start_region start position
        :param end_region end position

        :return: None
        """
        vlines_list = []
        if 'line width' in self.vlines_properties:
            line_width = self.vlines_properties['line width']
        else:
            line_width = 0.5

        if chrom_region not in self.vlines_intval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)
        for region in sorted(self.vlines_intval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            vlines_list.append(region.begin)

        for idx, track in enumerate(self.track_obj_list):
            ymin, ymax = axis_list[idx].get_ylim()

            axis_list[idx].vlines(vlines_list, ymin, ymax, linestyle='dashed', zorder=10,
                                  linewidth=line_width,
                                  color=(0, 0, 0, 0.7), alpha=0.5)

        return

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

        tracks_file_path = os.path.dirname(tracks_file)

        track_list = []
        for section_name in parser.sections():
            track_options = dict({"section_name": section_name})
            if section_name.endswith('[spacer]'):
                track_options['spacer'] = True
            elif section_name.endswith('[x-axis]'):
                track_options['x-axis'] = True
            for name, value in parser.items(section_name):
                if name in ['max_value', 'min_value', 'depth', 'width', 'line width', 'fontsize', 'scale factor', 'number of bins'] and value != 'auto':
                    track_options[name] = literal_eval(value)
                else:
                    track_options[name] = value

            if 'type' in track_options and track_options['type'] == 'vlines':
                self.vlines_properties = self.check_file_exists(track_options, tracks_file_path)
            elif 'skip' in track_options and track_options['skip'] != 'no':
                pass
            else:
                track_list.append(track_options)

        updated_track_list = []
        for track_dict in track_list:
            warn = None
            if 'file' in track_dict and track_dict['file'] != '':
                track_dict = self.check_file_exists(track_dict, tracks_file_path)
                if 'file_type' not in track_dict:
                    track_dict['file_type'] = self.guess_filetype(track_dict)

                #  set some default values
                if 'title' not in track_dict:
                    warn = "\ntitle not set for 'section {}'\n".format(track_dict['section_name'])
                    track_dict['title'] = ''
                if warn:
                    sys.stderr.write(warn)
            updated_track_list.append(track_dict)
        self.track_list = updated_track_list
        if self.vlines_properties:
            self.vlines_intval_tree, __, __ = file_to_intervaltree(self.vlines_properties['file'])

    @staticmethod
    def check_file_exists(track_dict, tracks_path):
        """
        Checks if a file or list of files exists. If the file does not exists
        tries to check if the file may be relative to the track_file path, in such case
        the path is updated.
        :param track_dict: dictionary of track values. Should contain
                            a 'file' key containing the path of the file
                            or files to be checked separated by space
                            For example: file1 file2 file3
        :param tracks_path: path of the tracks file
        :return: None
        """
        for file_type in ['boundaries_file', 'file']:
            if file_type in track_dict:
                file_names = [x for x in track_dict[file_type].split(" ") if x != '']
                full_path_file_names = []
                for file_name in file_names:
                    try:
                        open(file_name, 'r').close()
                        full_path_file_names.append(file_name)
                    except IOError:
                        try:
                            # try to find the file in the same path as the
                            # the path of the
                            name_with_tracks_path = tracks_path + "/" + file_name
                            open(name_with_tracks_path, 'r').close()
                            full_path_file_names.append(name_with_tracks_path)
                        except IOError:
                            sys.stderr.write("\n*ERROR*\nFile in section [{}] "
                                             "not found:\n{}\n\n".format(track_dict['section_name'],
                                                                         file_name))
                            sys.exit(1)

                track_dict[file_type] = " ".join(full_path_file_names)
        return track_dict

    @staticmethod
    def guess_filetype(track_dict):
        """

        :param track_dict: dictionary of track values with the 'file' key
                    containing a string path of the file or files. Only the ending
                     of the last file is used in case when there are more files
        :return: string file type detected
        """
        file = track_dict['file'].strip()
        if file.endswith(".bed") or file.endswith(".bed.gz"):
            file_type = 'bed'
        elif file.endswith(".npz") or file.endswith(".h5"):
            file_type = 'hic_matrix'
        elif file.endswith(".bw"):
            file_type = 'bigwig'
        elif file.endswith(".bg") or file.endswith(".bg.gz"):
            file_type = 'bedgraph'
        elif file.endswith(".bm") or file.endswith(".bm.gz"):
            file_type = 'bedgraph_matrix'
        else:
            sys.exit("Section {}: can not identify file type. Please specify "
                     "the file_type for {}".format(track_dict['section_name'], file))
        return file_type

    @staticmethod
    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    @staticmethod
    def print_elapsed(start):
        import time
        if start:
            print time.time() - start
        return time.time()


def opener(filename):
    """
    Determines if a file is compressed or not
    """
    import gzip
    f = open(filename, 'rb')
    if f.read(2) == '\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f


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
    file_h = opener(file_name)
    line_number = 0
    valid_intervals = 0
    prev_chrom = None
    prev_start = -1
    prev_line = None
    interval_tree = {}
    min_value = np.inf
    max_value = -np.inf

    for line in file_h.readlines():
        line_number += 1
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        try:
            start = int(start)
        except ValueError as detail:
            msg = "Error reading line: {}. The start field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        try:
            end = int(end)
        except ValueError as detail:
            msg = "Error reading line: {}. The end field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        if prev_chrom == chrom:
            assert prev_start <= start, \
                "Bed file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

        if chrom not in interval_tree:
            interval_tree[chrom] = IntervalTree()

        value = None

        if fields > 3:
            value = fields[3:]
            try:
                line_min = min(map(float, value))
                if line_min < min_value:
                    min_value = line_min

                line_max = max(map(float, value))
                if line_max > max_value:
                    max_value = line_max
            except ValueError:
                pass

        assert end > start, "Start position larger or equal than end for line\n{} ".format(line)

        interval_tree[chrom].add(Interval(start, end, value))
        valid_intervals += 1

    if valid_intervals == 0:
        sys.stderr.write("No valid intervals were found in file {}".format(file_name))
    file_h.close()

    return interval_tree, min_value, max_value


class TrackPlot(object):
    """
    The TrackPlot object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends TrackPlot
    should be created.

    It is expected that all TrackPlot objects have a plot method.

    """

    def __init__(self, properties_dict):
        self.properties = properties_dict


class PlotSpacer(TrackPlot):

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        ax.set_xlim(start_region, end_region)
        pass


class PlotBedGraph(TrackPlot):

    def __init__(self, properties_dict):
        # super(self.__class__, self).__init__(*args, **kwargs)
        self.properties = properties_dict
        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BEDGRAPH_COLOR
        self.interval_tree, ymin, ymax = file_to_intervaltree(self.properties['file'])

        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = ymax

        if 'min_value' not in self.properties or self.properties['min_value'] == 'auto':
            self.properties['min_value'] = ymin

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        score_list = []
        pos_list = []

        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)

        for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            score_list.append(float(region.data[0]))
            pos_list.append(region.begin + (region.end - region.begin) / 2)

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BEDGRAPH_COLOR

        if 'extra' in self.properties and self.properties['extra'][0] == '4C':
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
                sys.stderr.write("Invalid color {} for {}. "
                                 "Using gray instead.".format(self.properties['color'], self.properties['file']))
                self.ax.fill_between(pos_list, score_list, facecolor='gray')

        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.ax.set_xlim(start_region, end_region)

        ymax = self.properties['max_value']
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
            self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                         "[{}-{}]".format(ymin, ymax_print),
                         horizontalalignment='left', size='small',
                         verticalalignment='bottom')

        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)


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
        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)
        for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            start_pos.append(region.begin)
            values = map(float, region.data)
            matrix_rows.append(values)

        matrix = np.vstack(matrix_rows).T
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            matrix = np.flipud(matrix)

        if 'type' in self.properties and self.properties['type'] == 'lines':
            for row in matrix:
                self.ax.plot(start_pos, row, color='grey', linewidth=0.5)

            self.ax.plot(start_pos, matrix.mean(axis=0), linestyle="--", marker="|")
            ymax = self.properties['max_value']
            ymin = self.properties['min_value']
            self.ax.set_ylim(ymin, ymax)

            if 'show data range' in self.properties and self.properties['show data range'] == 'no':
                pass
            else:
                if float(ymax) % 1 == 0:
                    ymax_print = int(ymax)
                else:
                    ymax_print = "{:.1f}".format(ymax)

                if float(ymin) % 1 == 0:
                    ymin_print = int(ymin)
                else:
                    ymin_print = "{:.1f}".format(ymin)

                ydelta = ymax - ymin
                small_x = 0.01 * (end_region - start_region)
                # by default show the data range
                self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                             "[{}-{}]".format(ymin_print, ymax_print),
                             horizontalalignment='left',
                             verticalalignment='bottom')
            if 'plot horizontal lines' in self.properties and self.properties['horizontal lines']:
                # plot horizontal lines to compare values
                self.ax.hlines(np.arange(0, 1.1, 0.1), start_region, end_region, linestyle="--",
                               zorder=0, color='grey')

        else:
            x, y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
            shading = 'gouraud'
            vmax = self.properties['max_value']
            vmin = self.properties['min_value']

            img = self.ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax, shading=shading)
            img.set_rasterized(True)

        self.ax.set_xlim(start_region, end_region)
        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.properties['title'].replace('\\\\', '\\')
        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)


class PlotBigWig(TrackPlot):

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])
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
                num_bins = 700
                sys.stderr.write("'number of bins' value: {} for bigwig file {} "
                                 "is not valid. Using default value (700)".format(self.properties['number of bins'],
                                                                                  self.properties['file']))

        if chrom_region not in self.bw.chroms().keys():
            chrom_region = change_chrom_names(chrom_region)

        if chrom_region not in self.bw.chroms().keys():
            sys.stderr.write("Can not read region {} from bigwig file:\n\n"
                             "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                             "and that the region is valid".format(formated_region, self.properties['file']))

        if end_region - start_region < 2e6:
            scores = self.bw.values(chrom_region, start_region, end_region)

            if 'nans to zeros' in self.properties and self.properties['nans to zeros'] is True:
                scores[np.isnan(scores)] = 0

            scores = np.ma.masked_invalid(scores)

            lins = np.linspace(0, len(scores), num_bins).astype(int)
            scores_per_bin = [np.mean(scores[lins[x]:lins[x + 1]]) for x in range(len(lins) - 1)]
            _x = lins + start_region
            x_values = [float(_x[x] + _x[x + 1]) / 2 for x in range(len(lins) - 1)]
            if 'type' in self.properties and self.properties['type'] == 'line':
                self.ax.plot(x_values, scores_per_bin, linewidth=0.1, color=self.properties['color'])

            else:
                self.ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                                     color=self.properties['color'],
                                     facecolor=self.properties['color'])

        else:
            # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
            # of the memory. This only occurs when calling trackPlot from different
            # processors. Reloading the file solves the problem.
            num_tries = 0
            while num_tries < 5:
                num_tries += 1
                try:
                    scores = np.array(self.bw.stats(chrom_region, start_region, end_region, nBins=num_bins)).astype(float)
                except Exception as e:
                    import pyBigWig
                    self.bw = pyBigWig.open(self.properties['file'])

                    print "error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".format(e, num_tries)
                    pass
                else:
                    if num_tries > 1:
                        print "After {} the scores could be computed".format(num_tries)
                    break

            x_values = np.linspace(start_region, end_region, num_bins)
            self.ax.fill_between(x_values, scores, linewidth=0.1,
                                 color=self.properties['color'],
                                 facecolor=self.properties['color'], zorder=1)
        self.ax.set_xlim(start_region, end_region)
        ymin, ymax = self.ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
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
                         horizontalalignment='left',
                         verticalalignment='top')

        """
        self.ax.text(region_end, ymax - ydelta * 0.2, self.properties['title'],
                horizontalalignment='right', size='large',
                verticalalignment='bottom')

        """
        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center')

        return self.ax


class PlotHiCMatrix(TrackPlot):

    def __init__(self, properties_dict):
        # to avoid the color bar to span all the
        # width of the axis I pass two axes
        # to plot_matrix
        self.properties = properties_dict

        self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'])
        if len(self.hic_ma.matrix.data) == 0:
            exit("Matrix {} is empty".format(self.properties['file']))
        if 'show_masked_bins' in self.properties and self.properties['show_masked_bins'] == 'yes':
            pass
        else:
            self.hic_ma.maskBins(self.hic_ma.nan_bins)

        # check that the matrix can be log transformed
        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                if self.hic_ma.matrix.data.min() + 1 < 0:
                    exit("\n*ERROR*\nMatrix contains negative values.\n"
                         "log1p transformation can not be applied to \n"
                         "values in matrix: {}".format(self.properties['file']))

            elif self.properties['transform'] == '-log':
                if self.hic_ma.matrix.data.min() < 0:
                    exit("\n*ERROR*\nMatrix contains negative values.\n"
                         "log(-1 * <values>) transformation can not be applied to \n"
                         "values in matrix: {}".format(self.properties['file']))

            elif self.properties['transform'] == 'log':
                if self.hic_ma.matrix.data.min() < 0:
                    exit("\n*ERROR*\nMatrix contains negative values.\n"
                         "log transformation can not be applied to \n"
                         "values in matrix: {}".format(self.properties['file']))

        new_intervals = hicexplorer.utilities.enlarge_bins(self.hic_ma.cut_intervals)
        self.hic_ma.interval_trees, self.hic_ma.chrBinBoundaries = \
            self.hic_ma.intervalListToIntervalTree(new_intervals)

        self.hic_ma.cut_intervals = new_intervals
        binsize = self.hic_ma.getBinSize()
        max_depth_in_bins = int(self.properties['depth'] / binsize)

        # work only with the lower matrix
        # and remove all pixels that are beyond
        # 2 * max_depth_in_bis which are not required
        # (this is done by subtracting a second sparse matrix
        # that contains only the lower matrix that wants to be removed.
        limit = 2 * max_depth_in_bins
        self.hic_ma.matrix = scipy.sparse.triu(self.hic_ma.matrix, k=0, format='csr') - \
            scipy.sparse.triu(self.hic_ma.matrix, k=limit, format='csr')
        self.hic_ma.matrix.eliminate_zeros()

        # fill the main diagonal, otherwise it looks
        # not so good. The main diagonal is filled
        # with an array containing the max value found
        # in the matrix
        if sum(self.hic_ma.matrix.diagonal()) == 0:
            sys.stderr.write("Filling main diagonal with max value "
                             "because it empty and looks bad...\n")
            max_value = self.hic_ma.matrix.data.max()
            main_diagonal = scipy.sparse.dia_matrix(([max_value] * self.hic_ma.matrix.shape[0], [0]),
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
            self.boundaries_obj = PlotBoundaries({'file': self.properties['boundaries_file']})

    def plot(self, ax, label_ax, chrom, region_start, region_end):
        import copy
        self.cbar_ax = copy.copy(label_ax)
        self.label_ax = label_ax
        self.label_ax.set_axis_off()
        self.ax = ax

        chrom_sizes = self.hic_ma.get_chromosome_sizes()
        if chrom not in chrom_sizes.keys():
            chrom = change_chrom_names(chrom)
        if region_end > chrom_sizes[chrom]:
            sys.stderr.write("*Error*\nThe region to plot extends beyond the chromosome size. Please check.\n")
            sys.stderr.write("{} size: {}. Region to plot {}-{}\n".format(chrom, chrom_sizes[chrom],
                                                                          region_start, region_end))

        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges

        # get bin id of start and end of region in given chromosome
        chr_start_id, chr_end_id = self.hic_ma.getChrBinRange(chrom)
        chr_start = self.hic_ma.cut_intervals[chr_start_id][1]
        chr_end = self.hic_ma.cut_intervals[chr_end_id - 1][1]
        start_bp = max(chr_start, region_start - self.properties['depth'])
        end_bp = min(chr_end, region_end + self.properties['depth'])

        idx, start_pos = zip(*[(idx, x[1]) for idx, x in
                               enumerate(self.hic_ma.cut_intervals)
                               if x[0] == chrom and x[1] >= start_bp and x[2] <= end_bp])

        idx = idx[0:-1]
        # select only relevant matrix part
        matrix = self.hic_ma.matrix[idx, :][:, idx]
        # limit the 'depth' based on the length of the region being viewed

        region_len = region_end - region_start
        depth = min(self.properties['depth'], int(region_len * 1.25))
        depth_in_bins = int(1.5 * region_len / self.hic_ma.getBinSize())

        if depth < self.properties['depth']:
            # remove from matrix all data points that are not visible.
            matrix = matrix - scipy.sparse.triu(matrix, k=depth_in_bins, format='csr')
        matrix = np.asarray(matrix.todense().astype(float))
        if 'scale factor' in self.properties:
            matrix = matrix * self.properties['scale factor']

        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                matrix += 1
                self.norm = matplotlib.colors.LogNorm()

            elif self.properties['transform'] == '-log':
                mask = matrix == 0
                matrix[mask] = matrix[mask is False].min()
                matrix = -1 * np.log(matrix)

            elif self.properties['transform'] == 'log':
                mask = matrix == 0
                matrix[mask] = matrix[mask is False].min()
                matrix = np.log(matrix)

        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            vmax = self.properties['max_value']

        else:
            # try to use a 'aesthetically pleasant' max value
            vmax = np.percentile(matrix.diagonal(1), 80)

        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            vmin = self.properties['min_value']
        else:
            if depth_in_bins > matrix.shape[0]:
                depth_in_bins = matrix.shape[0] - 5
            vmin = np.median(matrix.diagonal(int(region_len / self.hic_ma.getBinSize())))

        sys.stderr.write("setting min, max values for track {} to: {}, {}\n".format(self.properties['section_name'],
                                                                                    vmin, vmax))
        img = self.pcolormesh_45deg(matrix, start_pos, vmax=vmax, vmin=vmin)
        img.set_rasterized(True)
        if self.plot_inverted:
            self.ax.set_ylim(depth, 0)
        else:
            self.ax.set_ylim(0, depth)

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
        try:
            if 'transform' in self.properties and \
                    self.properties['transform'] in ['log', 'log1p']:
                # get a useful log scale
                # that looks like [1, 2, 5, 10, 20, 50, 100, ... etc]
                """
                The following code is problematic with some versions of matplotlib.
                Should be uncommented once the problem is clarified
                from matplotlib.ticker import LogFormatter
                formatter = LogFormatter(10, labelOnlyBase=False)
                aa = np.array([1, 2, 5])
                tick_values = np.concatenate([aa * 10**x for x in range(10)])
                cobar = plt.colorbar(img, ticks=tick_values, format=formatter, ax=self.cbar_ax, fraction=0.95)
                """
                aa = np.array([0, 1, 2, 3, 4, 5])
                tick_values = set(np.concatenate([aa * 10**x for x in range(10)]))
                cobar = plt.colorbar(img, ticks=list(tick_values), ax=self.cbar_ax, fraction=0.95)
            else:
                cobar = plt.colorbar(img, ax=self.cbar_ax, fraction=0.95)
            cobar.solids.set_edgecolor("face")
            cobar.ax.set_ylabel(self.properties['title'])

            # adjust the labels of the colorbar
            labels = cobar.ax.get_yticklabels()
            ticks = cobar.ax.get_yticks()
            if ticks[0] == 0:
                # if the label is at the start of the colobar
                # move it above avoid being cut or overlapping with other track
                labels[0].set_verticalalignment('bottom')
            if ticks[-1] == 1:
                # if the label is at the end of the colobar
                # move it a bit inside to avoid overlapping
                # with other labels
                labels[-1].set_verticalalignment('top')
            cobar.ax.set_yticklabels(labels)

        except ValueError:
            pass
        # self.label_ax.text(0.3, 0.0, self.properties['title'], rotation=90,
        #                   horizontalalignment='left', size='large',
        #                   verticalalignment='bottom', transform=self.label_ax.transAxes)

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
        x = matrix_a[:, 1].reshape(n + 1, n + 1)
        y = matrix_a[:, 0].reshape(n + 1, n + 1)
        # plot
        im = self.ax.pcolormesh(x, y, np.flipud(matrix_c),
                                vmin=vmin, vmax=vmax, cmap=self.cmap, norm=self.norm)
        return im


class PlotXAxis(TrackPlot):

    def __init__(self, *args, **kwargs):
        super(PlotXAxis, self).__init__(*args, **kwargs)
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 15

    def plot(self, ax, label_axis, chrom_region, region_start, region_end):
        ax.set_xlim(region_start, region_end)
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"

        print ticks
        print labels
        ax.axis["x"] = ax.new_floating_axis(0, 0.5)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"].set_axis_direction("top")


class PlotBoundaries(TrackPlot):

    def __init__(self, *args, **kwargs):
        super(PlotBoundaries, self).__init__(*args, **kwargs)

        line_number = 0
        interval_tree = {}
        intervals = []
        prev_chrom = None
        valid_intervals = 0

        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                    continue
                try:
                    chrom, start, end = line.strip().split('\t')[0:3]
                except Exception as detail:
                    msg = 'Could not read line\n{}\n. {}'.format(line, detail)
                    sys.exit(msg)

                try:
                    start = int(start)
                    end = int(end)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                assert start <= end, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)

                if prev_chrom and chrom != prev_chrom:
                    start_array, end_array = zip(*intervals)
                    start_array = np.array(start_array)
                    end_array = np.array(end_array)
                    # check if intervals are consecutive or 1bp positions demarcating the boundaries
                    if np.any(end_array - start_array == 1):
                        # The file contains only boundaries at 1bp position.
                        end_array = start_array[1:]
                        start_array = start_array[:-1]
                    interval_tree[prev_chrom] = IntervalTree()
                    for idx in range(len(start_array)):
                        interval_tree[prev_chrom].add(Interval(start_array[idx], end_array[idx]))
                        valid_intervals += 1
                    intervals = []

                intervals.append((start, end))

                # each interval spans from the smallest start to the largest end
                prev_chrom = chrom

        start, end = zip(*intervals)
        start = np.array(start)
        end = np.array(end)
        # check if intervals are consecutive or 1bp positions demarcating the boundaries
        if np.any(end - start == 1):
            # The file contains only boundaries at 1bp position.
            end = start[1:]
            start = start[:-1]
        interval_tree[chrom] = IntervalTree()
        for idx in range(len(start)):
            interval_tree[chrom].add(Interval(start[idx], end[idx]))
            valid_intervals += 1

        if valid_intervals == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        x = []
        y = []
        if chrom_region not in self.interval_tree:
            orig = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            print 'changing {} to {}'.format(orig, chrom_region)
        for region in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            """
                  /\
                 /  \
                /    \
            _____________________
               x1 x2 x3
            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)
            x.extend([x1, x2, x3])
            y.extend([y1, y2, y1])

        ax.plot(x, y, color='black')
        ax.set_xlim(start_region, end_region)


class PlotBed(TrackPlot):

    def __init__(self, *args, **kwarg):
        super(PlotBed, self).__init__(*args, **kwarg)

        from matplotlib import font_manager
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 12
        else:
            self.properties['fontsize'] = float(self.properties['fontsize'])

        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BED_COLOR
        if 'border_color' not in self.properties:
            self.properties['border_color'] = 'black'
        if 'labels' not in self.properties:
            self.properties['labels'] = 'on'
        if 'style' not in self.properties:
            self.properties['style'] = 'flybase'
        if 'display' not in self.properties:
            self.properties['display'] = 'stacked'
        if 'interval height' not in self.properties:
            self.properties['interval_height'] = 100

        # to set the distance between rows
        self.row_scale = self.properties['interval_height'] * 2.3

        self.colormap = None
        # check if the color given is a color map
        color_options = [m for m in matplotlib.cm.datad]
        if self.properties['color'] in color_options and 'min_value' in self.properties and 'max_value' in self.properties:
            norm = matplotlib.colors.Normalize(vmin=self.properties['min_value'],
                                               vmax=self.properties['max_value'])
            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

        # self.process_bed()

    def process_bed(self, fig_width, region_start, region_end):

        # to improve the visualization of the genes
        # it is good to have an estimation of the label
        # length. In the following code I try to get
        # the length of a 'W'.
        if self.properties['labels'] == 'on':
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties['fontsize'] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
            print "len of w set to: {} bp".format(self.len_w)
        else:
            self.len_w = 1

        import readBed

        bed_file_h = readBed.ReadBed(opener(self.properties['file']))
        self.bed_type = bed_file_h.file_type
        valid_intervals = 0
        self.max_num_row = {}
        self.interval_tree = {}
        prev = None
        # check for the number of other intervals that overlap
        #    with the given interval
        #            1         2
        #  012345678901234567890123456
        #  1=========       4=========
        #       2=========
        #         3============
        #
        # for 1 row_last_position = [9]
        # for 2 row_last_position = [9, 14]
        # for 3 row_last_position = [9, 14, 19]
        # for 4 row_last_position = [26, 14, 19]

        row_last_position = []  # each entry in this list contains the end position
        # of genomic interval. The list index is the row
        # in which the genomic interval was plotted.
        # Any new genomic interval that wants to be plotted,
        # knows the row to use by finding the list index that
        # is larger than its start
        for bed in bed_file_h:
            if prev is not None:
                if prev.chromosome == bed.chromosome:
                    if bed.start < prev.start:
                        import ipdb
                        ipdb.set_trace()
                    assert bed.start >= prev.start, "*Error* Bed file not sorted. Please sort using sort -k1,1 -k2,2n"
                if prev.chromosome != bed.chromosome:
                    # init var
                    row_last_position = []
                    self.max_num_row[bed.chromosome] = 0
            else:
                self.max_num_row[bed.chromosome] = 0

            # check for overlapping genes including
            # label size (if plotted)

            if self.properties['labels'] == 'on':
                bed_extended_end = int(bed.end + (len(bed.name) * self.len_w))
            else:
                bed_extended_end = (bed.end + 2 * self.small_relative)

            # get smallest free row
            if len(row_last_position) == 0:
                free_row = 0
                row_last_position.append(bed_extended_end)
            else:
                # get list of rows that are less than bed.start, then take the min
                idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                if len(idx_list):
                    free_row = min(idx_list)
                    row_last_position[free_row] = bed_extended_end
                else:
                    free_row = len(row_last_position)
                    row_last_position.append(bed_extended_end)

            if bed.chromosome not in self.interval_tree:
                self.interval_tree[bed.chromosome] = IntervalTree()

            self.interval_tree[bed.chromosome].add(Interval(bed.start, bed.end, (bed, free_row)))
            valid_intervals += 1
            prev = bed
            if free_row > self.max_num_row[bed.chromosome]:
                self.max_num_row[bed.chromosome] = free_row

        print self.max_num_row

        if valid_intervals == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(self.properties['file_name']))

    def get_y_pos(self, bed, free_row):
        """
        The y_pos is set such that regions to be plotted do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.

        The algorithm uses a interval tree (self.region_interval) to check the overlaps
        and a sort of coverage vector 'rows used' to identify the row in which to plot
        :return: int y position
        """

        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if self.properties['display'] == 'interlaced':
            ypos = self.properties['interval_height'] if self.counter % 2 == 0 else 1

        elif self.properties['display'] == 'collapsed':
            ypos = 0

        else:
            ypos = free_row * self.row_scale
        return ypos

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.counter = 0
        self.small_relative = 0.004 * (end_region - start_region)
        self.process_bed(ax.get_figure().get_figwidth(), start_region, end_region)

        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)

        genes_overlap = sorted(self.interval_tree[chrom_region][start_region:end_region])

        # turn labels off when too many intervals are visible.
        if self.properties['labels'] != 'off' and len(genes_overlap) > 60:
            self.properties['labels'] = 'off'

        max_num_row_local = 1
        max_ypos = 0
        for region in genes_overlap:
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED9
            bed with rbg at end
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

            BED6
            bed with rbg at end
            chr2L   0       70000   ID_5    0.26864549832   .
            """
            self.counter += 1
            # an namedTuple object is stored in the region.value together with the free value for the gene
            bed, free_row = region.data
            if self.colormap:
                # translate value field (in the example above is 0 or 0.2686...) into a color
                rgb = self.colormap.to_rgba(bed.score)
                edgecolor = self.colormap.to_rgba(bed.score)
            else:
                rgb = self.properties['color']
                edgecolor = self.properties['border_color']

            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                    if 'border_color' in self.properties:
                        edgecolor = self.properties['border_color']
                    else:
                        edgecolor = self.properties['color']
                except IndexError:
                    pass

            ypos = self.get_y_pos(bed, free_row)

            # do not plot if the maximum interval rows to plot is reached
            if 'gene rows' in self.properties and free_row >= int(self.properties['gene rows']):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if self.bed_type == 'bed12':
                if self.properties['style'] == 'flybase':
                    self.draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor)
                else:
                    self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor)
            else:
                self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if self.properties['labels'] == 'off':
                pass
            elif bed.start > start_region and bed.end < end_region:
                ax.text(bed.end + self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                        bed.name, horizontalalignment='left',
                        verticalalignment='center', fontproperties=self.fp)

        if self.counter == 0:
            sys.stderr.write("*Warning* No intervals were found for file {} \n"
                             "in section '{}' for the interval plotted ({}:{}-{}).\n"
                             "".format(self.properties['file'],
                                       self.properties['section_name'],
                                       chrom_region, start_region, end_region))

        ymax = 0

        if 'global max row' in self.properties and self.properties['global max row'] == 'yes':
            ymin = self.max_num_row[chrom_region] * self.row_scale

        elif 'gene rows' in self.properties:
            ymin = int(self.properties['gene rows']) * self.row_scale
        else:
            ymin = max_ypos + self.properties['interval_height']

        print("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        """
        ax.text(start_region, (self.max_num_row[bed.chromosome] + 1) * 330 * 0.5,
                     "{}".format(self.properties['title']),
                     horizontalalignment='left', size='small',
                     verticalalignment='bottom')
        """

        if 'display' in self.properties:
            if self.properties['display'] == 'domain':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

        ax.set_xlim(start_region, end_region)

        label_ax.text(0.15, 1, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='top', transform=label_ax.transAxes)

    def draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws an interval with direction (if given)
        """
        from matplotlib.patches import Polygon

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, self.properties['interval_height'],
                                   edgecolor=edgecolor, facecolor=rgb, linewidth=0.5))
        else:
            vertices = self._draw_arrow(ax, bed.start, bed.end, bed.strand, ypos)
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
        half_height = float(self.properties['interval_height']) / 2
        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

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
            vertices = [(start_pos, ypos), (start_pos, ypos + self.properties['interval_height']),
                        (end_pos, ypos + self.properties['interval_height']), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=0.5))

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
        half_height = float(self.properties['interval_height']) / 2
        if strand == '+':
            x0 = start
            x1 = end  # - self.small_relative
            y0 = ypos
            y1 = ypos + self.properties['interval_height']
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            -----------------\
            ---------------- /

            """

            vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + self.small_relative, y0 + half_height), (x1, y0)]

        else:
            x0 = start  # + self.small_relative
            x1 = end
            y0 = ypos
            y1 = ypos + self.properties['interval_height']
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            /-----------------
            \-----------------
            """
            vertices = [(x0, y0), (x0 - self.small_relative, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

        return vertices

    def draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene like in flybase gbrowse.
        """
        from matplotlib.patches import Polygon

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
        quarter_height = float(self.properties['interval_height']) / 4
        three_quarter_height = quarter_height * 3

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end:
                y0 = ypos + quarter_height
                y1 = ypos + three_quarter_height
            else:
                y0 = ypos
                y1 = ypos + self.properties['interval_height']

            if x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + quarter_height), (x0, ypos + three_quarter_height),
                             (bed.thick_start, ypos + three_quarter_height),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (x1, ypos + self.properties['interval_height']), (x1, ypos),
                             (bed.thick_start, ypos), (bed.thick_start, ypos + quarter_height)])

            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + three_quarter_height),
                             (x1, ypos + three_quarter_height),
                             (x1, ypos + quarter_height),
                             (bed.thick_end, ypos + quarter_height),
                             (bed.thick_end, ypos)])
            else:
                vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 linewidth=0.1,
                                 edgecolor='none',
                                 facecolor=rgb))

            if idx < bed.block_count - 1:
                # plot small arrows using the character '<' or '>' over the back bone
                intron_length = bed.block_starts[idx + 1] - (bed.block_starts[idx] + bed.block_sizes[idx])
                marker = 5 if bed.strand == '+' else 4
                if intron_length > 3 * self.small_relative:
                    pos = np.arange(x1 + 1 * self.small_relative,
                                    x1 + intron_length + self.small_relative, int(2 * self.small_relative))
                    ax.plot(pos, np.zeros(len(pos)) + ypos + half_height, '.', marker=marker,
                            fillstyle='none', color='blue', markersize=3)

                elif intron_length > self.small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    ax.plot([intron_center], [ypos + half_height], '.', marker=5,
                            fillstyle='none', color='blue', markersize=3)


class PlotArcs(TrackPlot):

    def __init__(self, *args, **kwarg):
        import collections
        super(PlotArcs, self).__init__(*args, **kwarg)
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is an score.

        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        Arc = collections.namedtuple("Arc", ['chrom1', 'start1', 'end1', 'start2', 'end2', 'score'])
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2, score = line.strip().split('\t')
                except Exception as detail:
                    msg = 'File not valid. The format is chrom1 start1, end1, ' \
                          'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line)
                    sys.exit(msg)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)
                try:
                    score = float(score)
                except ValueError as detail:
                    msg = "Error reading line: {}. The score is not valid {}. " \
                          "\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                if chrom1 != chrom2:
                    sys.stderr.write("Only links in same chromosome are considere. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                value = Arc._make([chrom1, start1, end1, start2, end2, score])

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, value))
                valid_intervals += 1

        if valid_intervals == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        from matplotlib.patches import Arc
        max_diameter = 0
        count = 0

        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)
        arcs_in_region = sorted(self.interval_tree[chrom_region][region_start:region_end])

        for idx, interval in enumerate(arcs_in_region):
            # skip arcs whose start and end are outside the ploted region
            if interval.start < region_start and interval.end > region_end:
                continue

            diameter = interval.end - interval.start
            center = interval.start + float(diameter) / 2
            if diameter > max_diameter:
                max_diameter = diameter
            count += 1
            ax.plot([center], [diameter])
            if 'line width' in self.properties:
                line_width = float(self.properties['line width'])
            else:
                line_width = 0.5 * np.sqrt(interval.value.score)
            ax.add_patch(Arc((center, 0), diameter,
                             diameter * 2, 0, 0, 180, color=self.properties['color'], lw=line_width))

        print "{} arcs plotted".format(count)
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(max_diameter, -1)
        else:
            ax.set_ylim(-1, max_diameter)

        ax.set_xlim(region_start, region_end)
        label_ax.text(0.3, 0.0, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='bottom', transform=label_ax.transAxes)


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.
    """
    # TODO: mapping from chromosome names like mithocondria is missing
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom
