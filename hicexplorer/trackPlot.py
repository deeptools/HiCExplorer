from bx.intervals.intersection import IntervalTree, Interval

DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 3  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
DEFAULT_WIDTH_RATIOS = (0.95, 0.05)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.12, 'top': 0.9}


class IntervalFile(object):
    def __init__(self, file_name, type=None):
        # iterate over the file contents and save them for later usage
        file_h = open(file_name, 'r')
        intervals = []
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

            interval_tree[chrom].insert_interval(Interval(start, end, value=value))
            valid_intervals += 1
        if valid_intervals == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(file_name))


        self.interval_tree = interval_tree


class TrackPlot(object):

    def __init__(self, file_handler, ax, label_ax, properties_dict, region):
        self.file_handler = file_handler
        self.ax = ax
        self.label_ax = label_ax
        self.properties = properties_dict
        chrom_region, start_region, end_region = region
        self.chrom_region = chrom_region
        self.start_region = start_region
        self.end_region = end_region


class PlotBedGraph(TrackPlot):

    def plot(self):
    #def plot_bedgraph(ax, label_ax, self.properties, region):
        score_list = []
        pos_list = []

        for region in self.file_handler[self.chrom_region].find(self.start_region - 10000,
                                                                self.end_region + 10000):
            score_list.append(float(region.value[0]))
            pos_list.append(start + (end - start)/2)

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
                self.ax.fill_between(pos_list, score_list,
                                facecolor=self.properties['color'])
            except ValueError:
                exit("Invalid color {} for {}".format(self.properties['color'], self.properties['file']))

        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.ax.set_xlim(region[1], region[2])

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


class PlotBigWig(TrackPlot):

    def plot(self):
        bw = self.file_handler
        # compute the score in bins of 10000 SLOW
    #    score = np.array([bw.query(region[0], x, x+10000,1)[0]['mean']
    #                      for x in range(region[1], region[2], 10000)])

        num_bins = 700
        if 'number of bins' in self.properties:
            try:
                num_bins = int(self.properties['number of bins'])
            except TypeError:
                exit("'number of bins' value: {} for bigwig file {} "
                     "is not valid.".format(self.properties['number of bins'],
                                            self.properties['file']))

        scores = bw.get_as_array(self.chrom_region, self.start_region, self.end_region)
        if scores is None:
            # usually bw returns None when the chromosome
            # is not found. So, we try either
            # removing or appending 'chr'
            if self.chrom_region.startswith('chr'):
                scores = bw.get_as_array(self.chrom_region[3:] + self.chrom_region, self.start_region, self.end_region)
            else:
                scores = bw.get_as_array('chr' + self.chrom_region, self.start_region, self.end_region)

            if scores is None:
                exit("Can not read region {}:{}-{} from bigwig file:\n\n"
                     "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                     "and that the region is valid".format(self.chrom_region, self.start_region, self.end_region,
                                                           self.properties['file']))

        if 'nans to zeros' in self.properties and self.properties['nans to zeros'] is True:
            scores[np.isnan(scores)] = 0

        scores = np.ma.masked_invalid(scores)

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BIGWIG_COLOR

        if self.end_region - self.start_region < 2e6:
            if scores is None:
                self.chrom_region = self.chrom_region.replace('chr', '')
                scores = np.ma.masked_invalid(bw.get_as_array(self.chrom_region, self.start_region, self.end_region))
            if scores is None:
                sys.stderr.write("could not find values for region {}\n".format(region))

            else:
                lins = np.linspace(0, len(scores), num_bins).astype(int)
                scores_per_bin = [np.mean(scores[lins[x]:lins[x+1]]) for x in range(len(lins)-1)]
                _x = lins + self.start_region
                x_values = [float(_x[x] + _x[x+1])/2 for x in range(len(lins)-1)]
                ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                                color=self.properties['color'],
                                facecolor=self.properties['color'])

        else:
            # this method produces shifted regions. It is not clear to me why this happens.
            # Thus I only activate the faster but shifted method for large regions
            # when the previous method would be to slow
            score = bw.query(self.chrom_region, self.start_region, self.end_region, num_bins)
            if score is None:
                self.chrom_region = self.chrom_region.replace('chr', '')
                score = bw.query(self.chrom_region, self.start_region, self.end_region, num_bins)
            if score is None:
                sys.stderr.write("could not find values for region {}\n".format(region))
            else:
                score = [x['mean'] for x in score]
                x_values = np.linspace(self.start_region, self.end_region, num_bins)
                ax.fill_between(x_values, score, linewidth=0.1,
                                color=self.properties['color'],
                                facecolor=self.properties['color'], zorder=1)

        ax.set_xlim(region[1], region[2])
        ymin, ymax = ax.get_ylim()
        if 'max_value' in self.properties and ['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':

            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)

    #    ax.set_yticks([ymax])
        ydelta = ymax - ymin

    #    ax.set_yticklabels(["{}-{}".format(int(ymin), int(ymax))], size='large')
        # set min max
        if float(ymax) % 1 == 0:
            ymax_print = int(ymax)
        else:
            ymax_print = "{:.1f}".format(ymax)
        small_x = 0.01 * (self.end_region - self.start_region)
        if 'show data range' in self.properties and self.properties['show data range'] == 'no':
            pass
        else:
            # by default show the data range
            ax.text(self.start_region - small_x, ymax - ydelta * 0.2,
                    "[{}-{}]".format(int(ymin), ymax_print),
                    horizontalalignment='left', size='small',
                    verticalalignment='bottom')

        """
        ax.text(region_end, ymax - ydelta * 0.2, self.properties['title'],
                horizontalalignment='right', size='large',
                verticalalignment='bottom')

        """
        label_ax.text(0.15, 0, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='bottom')





file = "/data/manke/group/ramirez/HiC-ThomasLing/data/bedfiles/dm3.genes.chrX.bed"
ii = IntervalFile(file)
ii.interval_tree['chrX'].fetch(10000, 2000)