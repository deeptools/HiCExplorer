
"""

[hic]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/hic_data/s2/hic-norm/RF/merge.npz
title =
colormap = RdYlBu_r
depth = 100000
#min_value =2.8
#max_value = 3.0
transform = log1p
boundaries_file = conductance_vs_hic/boundaries_all.bed
x labels = yes
type = arcplot
type = interaction
#optional in case it can not be guessed by the file ending
file_type = hic_matrix
# show masked bins plots as white lines
# those bins that were not used during the correction
# the default is to extend neighboring bins to
# obtain an aesthetically pleasant output
show_masked_bins = yes

[x-axis]
#optional
fontsize=20
# default is bottom meaning below the axis line
where=top

[spacer]

[bigwig]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/external/Graveley_mRNA-seq/GSM390060_Kc167-4_spa.bw
title = Kc RNA-seq (Cherbas et al.)
color = black
min_value = 0
#max_value = auto
width = 1.5
number of bins = 500
nans to zeros = True
# optional. If type is not specified, then a fill curve is plotted.
type = line
#optional in case it can not be guessed by the file ending
file_type = bigwig

[simple bed]
file = file.bed
title = peaks
color = red
# optional boder color. Set to none for no border color
border_color = black
width = 0.5
# optional. If not given is guessed from the file ending
file_type = bed

[bed genes]
# example of a genes track
# has the same options as a simple
# bed but if the type=genes is given
# the the file is interpreted as gene
# file. If the bed file contains the exon
# structure then this is plotted. Otherwise
# a region **with direction** is plotted.
file = genes.bed
title = genes
color = darkblue
#if color is a valid colormap, then the score is mapped
# to the colormap
#color = RdBlGn
width = 5
# to turn off/on printing of labels
labels = off
# optional. If not given is guessed from the file ending
file_type = bed
# optional: font size can be given if default are not good
fontsize = 10
# optional
# display = collapsed
# display = interleaved
#optional, default is black
#border_color = black
# style to plot the genes when they have exon information
#style = UCSC
#style = flybase
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size
# genes in all images
#gene rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global max row = yes


[chrom states]
# this is a case of a bed file that is plotted 'collapsed'
# useful to plot chromatin states if the bed file contains
# the color to plot
file = chromatinStates.bed
title = chromatin states
# color is replaced by the color in the bed file
# in this case
color = black
# optional boder color. Set to none for no border color
border_color = black
# default behaviour when plotting intervals from a
# bed file is to 'expand' them such that they
# do not overlap. The display = collapsed
# directive overlaps the intervals.
display = collapsed
width = 0.3

[bedgraph]
file = file.bg
title = bedgraph track
color = green
width = 0.2
# optional, otherwise guseed from file ending
file_type = bedgraph


[bedgraph matrix]
# a bedgraph matrix file is like a bedgraph, except that per bin there
# are more than one value separated by tab: E.g.
# chrX	18279	40131	0.399113	0.364118	0.320857	0.274307
# chrX	40132	54262	0.479340	0.425471	0.366541	0.324736
file = spectra_conductance.bm
title = conductance spectra
width = 1.5
orientation = inverted
min_value = 0.10
max_value = 0.70
# if type is set as lines, then the TAD score lines are drawn instead
# of the matrix
# set to lines if a heatmap representing the matrix
# is not wanted
type = lines
file_type = bedgraph_matrix
plot horizontal lines=False

[dendrogram]
# a dendrogram is saved in a bed file that contains
# further fields defining the linkage.
file =  linkage.bed
title = dendrogram
width = 2
orientation = inverted
# required on this case. Otherwise
# the bed file will be plotted as regions
file_type = dendrogram
# y positions to draw horizontal lines
hlines = 0.2 0.3

[vlines]
# vertical dotted lines from the top to the bottom of the figure
# can be drawn. For this a bed file is required
# but only the first two columns (chromosome name and start
# are used.
# vlines can also be given at the command line as a list
# of genomic positions. However, sometimes to give a file
# is more convenient in case many lines want to be plotted.
file = regions.bed
type = vlines


"""

from __future__ import division
import sys
import argparse
import matplotlib
matplotlib.use('Agg')

import hicexplorer.trackPlot
from hicexplorer._version import __version__

DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 3  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
# DEFAULT_WIDTH_RATIOS = (0.95, 0.05)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.12, 'top': 0.9}


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description='Plots the diagonal,  and some values close to '
        'the diagonal of a  HiC matrix. The diagonal of the matrix is '
        'plotted horizontally for a region. I will not draw the diagonal '
        'for the whole chromosome',
        usage="%(prog)s --tracks tracks.ini --region chr1:1000000-4000000 -o image.png")

    parser.add_argument('--tracks',
                        help='File containing the instructions to plot the tracks ',
                        type=argparse.FileType('r'),
                        required=True,
                        )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--region',
                       help='Region to plot, the format is chr:start-end')

    group.add_argument('--BED',
                       help='Instead of a region, a file containing the regions to plot, in BED format, '
                       'can be given. If this is the case, multiple files will be created using a prefix '
                       'the value of --outFileName',
                       type=argparse.FileType('r')
                       )

    parser.add_argument('--width',
                        help='figure width in centimeters',
                        type=float,
                        default=DEFAULT_FIGURE_WIDTH)

    parser.add_argument('--height',
                        help='Figure height in centimeters. If not given, the figure height is computed '
                             'based on the widths of the tracks and on the depth '
                             'of the Hi-C tracks, if any. Setting this '
                             'this parameter can cause the Hi-C tracks to look collapsed or expanded.',
                        type=float)

    parser.add_argument('--title', '-t',
                        help='Plot title',
                        required=False)

    parser.add_argument('--scoreName', '-s',
                        help='Score name',
                        required=False)

    parser.add_argument('--outFileName', '-out',
                        help='File name to save the image, file prefix in case multiple images '
                             'are stored',
                        required=True)

    parser.add_argument('--zMax',
                        help='zMax',
                        type=float,
                        default=None)

    parser.add_argument('--vlines',
                        help='Genomic cooordindates separated by space. E.g. '
                        '--vlines 150000 3000000 124838433 ',
                        type=int,
                        nargs='+'
                        )

    parser.add_argument('--fontSize',
                        help='Font size for the labels of the plot',
                        type=float,
                        )

    parser.add_argument('--dpi',
                        help='Resolution for the image in case the'
                             'ouput is a raster graphics image (e.g png, jpg)',
                        type=int,
                        default=72
                        )

    parser.add_argument('--trackLabelFraction',
                        help='By default the space dedicated to the track labels is 0.05 of the'
                             'plot width. This fraction can be changed with this parameter if needed.',
                        default=0.05,
                        type=float)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def get_region(region_string):
    """
    splits a region string into
    a chrom, start_region, end_region tuple
    The region_string format is chr:start-end
    """
    if region_string:
        region_string = region_string.translate(
            None, ",.;|!{}()").replace("-", ":")
        region = region_string.split(":")
        chrom = region[0]
        try:
            region_start = int(region[1])
        except IndexError:
            region_start = 0
        try:
            region_end = int(region[2])
        except IndexError:
            region_end = 1e15  # a huge number

        if region_end <= region_start:
            exit("Please check that the region end is larger than the region start.\n"
                 "Values given:\nstart: {}\nend: {}\n".format(region_start, region_end))

        return chrom, region_start, region_end


def main(args=None):

    args = parse_arguments().parse_args(args)
    trp = hicexplorer.trackPlot.PlotTracks(args.tracks.name, args.width, fig_height=args.height,
                                           fontsize=args.fontSize, dpi=args.dpi,
                                           track_label_width=args.trackLabelFraction)

    if args.BED:
        count = 0
        for line in args.BED.readlines():
            count += 1
            try:
                chrom, start, end = line.strip().split('\t')[0:3]
            except ValueError:
                continue
            try:
                start, end = map(int, [start, end])
            except ValueError as detail:
                sys.stderr.write("Invalid value found at line\t{}\t. {}\n".format(line, detail))
            file_name = "{}_{}:{}-{}".format(args.outFileName, chrom, start, end)
            if end - start < 200000:
                start -= 100000
                end += 100000
            sys.stderr.write("saving {}'\n".format(file_name))
            trp.plot(file_name, chrom, start, end, title=args.title)
    else:
        region = get_region(args.region)
        trp.plot(args.outFileName, *region, title=args.title)
