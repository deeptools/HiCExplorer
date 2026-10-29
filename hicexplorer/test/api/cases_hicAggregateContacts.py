"""Equivalence cases for hicexplorer.hicAggregateContacts (apieq.py).

The drivers rebuild what hicAggregateContacts.main hands to the functions:
the matrix (Li_et_al_2015 as h5 and cool, small_test_matrix_50kb_res.h5 for
several chromosomes) loaded with hicmatrix and prepared as main prepares it
(maskBins, NaN to 0, enlarge_bins), agg_info with its chromosome
coordinates, and the BED files of test_data/hicAggregateContacts. The steps a
case does not test are taken from the module under test in each
interpreter, the step it tests is the target; the drivers return agg_info,
the clustering and the files and figure contents the plot functions write
and draw (plt.savefig is intercepted, so no image is written).

Deviations from 3.7.6 (cpp/python/src/api/hicAggregateContacts.cpp):
method='spectral' raises NotImplementedError (main never reaches it);
largeRegionsOperation other than first, last and center raises TypeError;
the matrix must be csr, compute_avg and the clustering take float64
submatrices of one shape, BED lines must be str; a window past the matrix
raises ValueError.
"""
import argparse
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicAggregateContacts.'
LI_H5 = os.path.join(TEST_DATA, 'Li_et_al_2015.h5')
LI_COOL = os.path.join(TEST_DATA, 'Li_et_al_2015.cool')
SMALL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
BED = os.path.join(TEST_DATA, 'hicAggregateContacts', 'test_regions.bed')
BED_STRAND = os.path.join(TEST_DATA, 'hicAggregateContacts', 'test_regions_strand.bed')
BED_LARGE = os.path.join(TEST_DATA, 'hicAggregateContacts', 'test_regions_large.bed')
ROW1 = os.path.join(TEST_DATA, 'hicAggregateContacts', 'bed1_row-wise.bed')
ROW2 = os.path.join(TEST_DATA, 'hicAggregateContacts', 'bed2_row-wise.bed')
ROW1_STRAND = os.path.join(TEST_DATA, 'hicAggregateContacts', 'bed1_row-wise_strand.bed')
ROW2_STRAND = os.path.join(TEST_DATA, 'hicAggregateContacts', 'bed2_row-wise_strand.bed')


def module():
    import importlib
    return importlib.import_module('hicexplorer.hicAggregateContacts')


def prepared(path):
    """The matrix as main prepares it."""
    from hicmatrix import HiCMatrix as hm
    import importlib
    utilities = importlib.import_module('hicexplorer.utilities')
    ma = hm.hiCMatrix(path)
    ma.maskBins(ma.nan_bins)
    ma.matrix.data[np.isnan(ma.matrix.data)] = 0
    ma.maskBins(ma.nan_bins)
    ma.setCutIntervals(utilities.enlarge_bins(ma.cut_intervals))
    return ma


def new_agg_info(ma):
    chrom_coord = dict()
    chrom_list = ma.getChrNames()
    for chrom in chrom_list:
        first, last = ma.getChrBinRange(chrom)
        chrom_coord[chrom] = (ma.getBinPos(first)[1], ma.getBinPos(last - 1)[2])
    return {"chrom_coord": chrom_coord, "seen": [], "agg_matrix": {c: {} for c in chrom_list},
            "agg_total": {c: {} for c in chrom_list}, "agg_diagonals": {c: {} for c in chrom_list},
            "agg_contact_position": {c: {} for c in chrom_list}, "agg_center_values": {c: {} for c in chrom_list},
            "counter": 0, "used_counter": 0, "empty_mat": 0}


def plain(value):
    """`value` with every numpy.matrix replaced by ('numpy.matrix', its
    ndarray): apieq.canon reduces arrays through tolist(), which nests a
    matrix's rows."""
    if isinstance(value, np.matrix):
        return ('numpy.matrix', np.asarray(value))
    if isinstance(value, dict):
        return type(value)((k, plain(v)) for k, v in value.items())
    if isinstance(value, list):
        return [plain(v) for v in value]
    if isinstance(value, tuple):
        return tuple(plain(v) for v in value)
    return value


def guarded(function):
    try:
        return plain(function())
    except SystemExit as error:
        return ('SystemExit', error.code)


def read_bed(target, path, chrom_list, strand=None):
    with open(path) as fh:
        return guarded(lambda: target(fh, chrom_list, strand))


def aggregate(target, matrix, bed, bed2=None, M_half=5, operation='first', range='50000:900000',
              transform='none', mode='intra-chr', strand=None):
    ma = prepared(matrix)
    agg_info = new_agg_info(ma)
    chrom_list = ma.getChrNames()
    with open(bed) as fh:
        bed1 = module().read_bed_per_chrom(fh, chrom_list, strand)
    if bed2 is None:
        intervals2 = bed1
    else:
        with open(bed2) as fh:
            intervals2 = module().read_bed_per_chrom(fh, chrom_list, strand)
    result = guarded(lambda: target(bed1, intervals2, agg_info, ma, M_half, operation, range=range,
                                    transform=transform, mode=mode, pConsiderStrandDirection=strand))
    return result, plain(agg_info), ma.bin_size


def aggregate_rows(target, matrix, bed1, bed2, M_half=5, operation='first', range='50000:900000',
                   transform='none', mode='intra-chr', strand=None):
    ma = prepared(matrix)
    agg_info = new_agg_info(ma)
    with open(bed1) as f1, open(bed2) as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
    result = guarded(lambda: target(lines1, lines2, agg_info, ma, ma.getChrNames(), M_half, operation,
                                    range=range, transform=transform, mode=mode, perChr=False,
                                    pConsiderStrandDirection=strand))
    return result, plain(agg_info)


def contacts(target, matrix, intervals, M_half=3, mode='all', operation='center', range='50000:900000',
             transform=None):
    """count_contacts on a sequence of intervals sharing agg_info and seen_chrs."""
    ma = prepared(matrix)
    agg_info = new_agg_info(ma)
    seen_chrs = []
    returned = []
    for interval in intervals:
        returned.append(target(interval, ma, M_half, mode, agg_info, operation, seen_chrs, range, transform))
    return returned, plain(agg_info), seen_chrs


def aggregated(matrix=LI_H5, bed=BED, M_half=5, transform='none', mode='intra-chr', range='50000:900000'):
    """agg_info after aggregate_contacts, from the module under test."""
    ma = prepared(matrix)
    agg_info = new_agg_info(ma)
    with open(bed) as fh:
        intervals = module().read_bed_per_chrom(fh, ma.getChrNames())
    module().aggregate_contacts(intervals, intervals, agg_info, ma, M_half, 'first', range=range,
                                transform=transform, mode=mode)
    return agg_info


def genome_info(agg_info):
    info = {"coords": [], "centers": [], "submatrices": [], "clustered_dict": [], "diagonal": []}
    for chrom1 in agg_info["agg_matrix"]:
        for chrom2, matrices in agg_info["agg_matrix"][chrom1].items():
            info["submatrices"] += matrices
            for idx in range(len(matrices)):
                s1, e1, s2, e2 = agg_info["agg_contact_position"][chrom1][chrom2][idx]
                info["coords"].append((chrom1, s1, e1, chrom2, s2, e2))
                info["centers"].append(agg_info["agg_center_values"][chrom1][chrom2][idx])
                info["diagonal"].append(agg_info["agg_diagonals"][chrom1][chrom2][idx])
    return info


def clusters(target, k, method, how='full', max_deviation=2, keep_outlier=False, transform='none', M_half=5,
             matrix=LI_H5):
    info = genome_info(aggregated(matrix=matrix, transform=transform, M_half=M_half))
    try:
        result = target(info, k, method, how, max_deviation, keep_outlier)
    except SystemExit as error:
        return ('SystemExit', error.code)
    identity = {key: result[key] is info[key] for key in ('coords', 'centers', 'submatrices', 'diagonal')}
    return plain(result), identity


def outliers_of_matrix(target, seed, max_deviation):
    return target(np.asmatrix(random_stack((5, 3), seed)), max_deviation)


def average_of_matrices(target, operation):
    return target([np.asmatrix(random_stack((3, 3), s)) for s in range(5)], operation)


def clusters_of(target, info, *args):
    return guarded(lambda: target(info, *args))


def matrices(target, k, method, how='full', perChr=False, max_deviation=2, keep_outlier=False, matrix=LI_H5,
             mode='intra-chr', bed=BED, transform='none'):
    agg_info = aggregated(matrix=matrix, mode=mode, bed=bed, transform=transform)
    return guarded(lambda: target(agg_info, k, method, how, perChr, max_deviation, keep_outlier))


def random_stack(shape, seed=0, zeros=False):
    values = np.random.RandomState(seed).standard_normal(shape)
    if zeros:
        values[...] = 0
        values[0] = 1
    return values


class Handle(object):
    def __init__(self, name):
        self.name = name


def figure_content(fig):
    out = []
    for ax in fig.axes:
        images = [np.asarray(im.get_array()) for im in ax.get_images()]
        lines = [(line.get_xdata().tolist(), np.asarray(line.get_ydata()).tolist()) for line in ax.get_lines()]
        out.append((ax.get_title(), ax.get_xlabel(), ax.get_ylabel(), images, lines))
    return out


def tie_canonical_lines(text):
    """The lines of a contact pairs file with the order of tied lines removed.

    The lines are written in np.argsort order of their last column, the value
    at the centre of the submatrix. np.argsort is not stable, and numpy 1.26
    orders equal keys differently with and without its AVX-512 argsort, so the
    order of lines with an equal last column depends on the CPU numpy runs on.
    Each run of consecutive lines with the same last column is sorted; the
    order of the runs, the number of lines and the content of every line are
    kept.
    """
    lines = text.split('\n')
    trailing = lines and lines[-1] == ''
    if trailing:
        lines = lines[:-1]
    out = []
    run = []
    for line in lines:
        if run and line.split('\t')[-1] != run[0].split('\t')[-1]:
            out.extend(sorted(run))
            run = []
        run.append(line)
    out.extend(sorted(run))
    return '\n'.join(out) + ('\n' if trailing else '')


def tie_canonical_rows(image, center):
    """The rows of a heatmap with the order of rows tied in `center` removed.

    The rows of the diagnostic heatmap are sorted by their value at the centre
    column (np.argsort, see tie_canonical_lines). Each run of consecutive rows
    with the same centre value is sorted lexicographically.
    """
    image = np.asarray(image)
    if image.ndim != 2 or image.shape[0] == 0:
        return image
    rows = image.tolist()
    out = []
    run = []
    for row in rows:
        if run and row[center] != run[0][center]:
            out.extend(sorted(run))
            run = []
        run.append(row)
    out.extend(sorted(run))
    return np.asarray(out, dtype=image.dtype)


def plotting(target, directory, k=1, method='no_clust', perChr=False, prefix=True, pairs=True, plot_type='2d',
             operation='mean', diagnostic=False, matrix=LI_H5, mode='intra-chr', how='full', transform='none',
             ties=False):
    import matplotlib.pyplot as plt
    agg_info = aggregated(matrix=matrix, mode=mode, transform=transform)
    clustered = module().cluster_matrices(agg_info, k=k, method=method, how=how, perChr=perChr)
    args = argparse.Namespace(vMin=None, vMax=None, colorMap='RdYlBu_r', operationType=operation,
                              plotType=plot_type, outFilePrefixMatrix=os.path.join(directory, 'mat') if prefix else None,
                              outFileContactPairs=os.path.join(directory, 'pairs') if pairs else None,
                              disable_bbox_tight=False, outFileName=Handle(os.path.join(directory, 'agg.png')),
                              dpi=50, diagnosticHeatmapFile=Handle(os.path.join(directory, 'heat.png')))
    captured = []
    original = plt.savefig
    plt.savefig = lambda *a, **kw: captured.append(figure_content(plt.gcf()))
    try:
        if diagnostic:
            result = guarded(lambda: target(clustered, 5, args))
        else:
            result = guarded(lambda: target(clustered, k, 5, args))
    finally:
        plt.savefig = original
        plt.close('all')
    files = {}
    for name in sorted(os.listdir(directory)):
        with open(os.path.join(directory, name)) as handle:
            files[name] = handle.read()
    if ties:
        # ties=True: the order of tied lines and rows is not compared, see
        # tie_canonical_lines. Only the contact pairs files and the diagnostic
        # heatmap images (centre column 5, the M_half the driver passes) change.
        for name in files:
            if name.startswith('pairs'):
                files[name] = tie_canonical_lines(files[name])
        if diagnostic:
            captured = [[(title, xlabel, ylabel, [tie_canonical_rows(image, 5) if np.ndim(image) == 2 and np.shape(image)[1] == 11 else image
                                                  for image in images], lines)
                         for (title, xlabel, ylabel, images, lines) in figure] for figure in captured]
    return result, captured, files


BED_SMALL = b"# header\nX\t100000\t100001\tp1\t1\t+\nchrZ\t1\t2\tp2\t1\t-\nX\t500000\t500001\tp3\t1\t-\n" \
            b"2L\t1000\t2000\tp4\t1\t+\n"


TIES_NOTE = ('the lines of the contact pairs files and the rows of the diagnostic heatmaps are in np.argsort order of the '
             'centre value, and numpy 1.26 orders tied keys differently with and without its AVX-512 argsort; the order of '
             'tied lines and rows is not compared (tie_canonical_lines), every other order and every value is')


HEAT_NOTE = TIES_NOTE + '; the summary lines hold column means summed over the rows in that order, so they agree to a relative 1e-12'


CASES = [
    # read_bed_per_chrom
    Case('read_bed/li', M + 'read_bed_per_chrom', driver='read_bed', args=(BED, ['X'])),
    Case('read_bed/dm3', M + 'read_bed_per_chrom', driver='read_bed',
         args=(BED, ['chr2L', 'chr3R', 'chrX', 'chr2R', 'chr3L'])),
    Case('read_bed/strand', M + 'read_bed_per_chrom', driver='read_bed', args=(BED_STRAND, ['X'], True)),
    Case('read_bed/strand_missing', M + 'read_bed_per_chrom', driver='read_bed', args=(BED, ['X'], True)),
    Case('read_bed/small', M + 'read_bed_per_chrom', driver='read_bed', args=(File('a.bed'), ['X', 'chr2L'], True),
         files={'a.bed': BED_SMALL}),
    Case('read_bed/blank_line', M + 'read_bed_per_chrom', driver='read_bed', args=(File('a.bed'), ['X']),
         files={'a.bed': b"X\t1\t2\n\nX\t3\t4\n"}),
    Case('read_bed/one_field', M + 'read_bed_per_chrom', driver='read_bed', args=(File('a.bed'), ['X']),
         files={'a.bed': b"X\n"}),
    Case('read_bed/not_a_number', M + 'read_bed_per_chrom', driver='read_bed', args=(File('a.bed'), ['X']),
         files={'a.bed': b"X\t1.5\t2\n"}),

    # aggregate_contacts
    Case('aggregate/li_intra', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED)),
    Case('aggregate/li_cool_intra', M + 'aggregate_contacts', driver='aggregate', args=(LI_COOL, BED)),
    Case('aggregate/li_total_counts', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED),
         kwargs={'transform': 'total-counts', 'operation': 'center'}),
    Case('aggregate/li_all', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED),
         kwargs={'mode': 'all', 'M_half': 2}),
    Case('aggregate/li_strand', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED_STRAND),
         kwargs={'strand': True, 'M_half': 3}),
    Case('aggregate/li_large_first', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED_LARGE),
         kwargs={'operation': 'first'}),
    Case('aggregate/li_large_last', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED_LARGE),
         kwargs={'operation': 'last'}),
    Case('aggregate/li_large_center', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED_LARGE),
         kwargs={'operation': 'center', 'range': '10000:2000000'}),
    Case('aggregate/li_inter_one_chromosome', M + 'aggregate_contacts', driver='aggregate', args=(LI_H5, BED),
         kwargs={'mode': 'inter-chr'}),
    Case('aggregate/small_inter', M + 'aggregate_contacts', driver='aggregate', args=(SMALL, BED),
         kwargs={'mode': 'inter-chr', 'M_half': 2}),
    Case('aggregate/small_all', M + 'aggregate_contacts', driver='aggregate', args=(SMALL, BED),
         kwargs={'mode': 'all', 'M_half': 1, 'transform': 'total-counts'}),
    Case('aggregate/small_intra', M + 'aggregate_contacts', driver='aggregate', args=(SMALL, BED),
         kwargs={'mode': 'intra-chr', 'M_half': 2, 'range': '100000:5000000'}),
    Case('aggregate/small_two_beds', M + 'aggregate_contacts', driver='aggregate', args=(SMALL, BED, BED_LARGE),
         kwargs={'mode': 'all', 'M_half': 1}),

    # aggregate_contacts_per_row
    Case('rows/li_intra', M + 'aggregate_contacts_per_row', driver='aggregate_rows', args=(LI_H5, ROW1, ROW2)),
    Case('rows/li_all', M + 'aggregate_contacts_per_row', driver='aggregate_rows', args=(LI_H5, ROW1, ROW2),
         kwargs={'mode': 'all'}),
    Case('rows/li_strand', M + 'aggregate_contacts_per_row', driver='aggregate_rows',
         args=(LI_H5, ROW1_STRAND, ROW2_STRAND), kwargs={'strand': True, 'mode': 'all'}),
    Case('rows/small_inter', M + 'aggregate_contacts_per_row', driver='aggregate_rows', args=(SMALL, ROW1, ROW2),
         kwargs={'mode': 'inter-chr', 'M_half': 1}),
    Case('rows/small_intra', M + 'aggregate_contacts_per_row', driver='aggregate_rows', args=(SMALL, ROW1, ROW2),
         kwargs={'mode': 'intra-chr', 'M_half': 1, 'range': '0:100000000'}),

    # count_contacts
    Case('count/sequence', M + 'count_contacts', driver='contacts', args=(LI_H5, [
        [('X', 100000, 100001), ('X', 400000, 400001)],
        [('X', 400000, 400001), ('X', 100000, 100001)],
        [('X', 100000, 100001), ('X', 100000, 100500)],
        [('X', 100000, 100001), ('X', 10 ** 9, 10 ** 9 + 1)],
        [('chrZ', 1, 2), ('X', 1, 2)],
        [('X', 50, 51), ('X', 900000, 900001)],
        [('X', '300000', '301000'), ('X', '700000', '701000')],
        [('X', 1000000, 1000001, '-'), ('X', 1200000, 1200001, '-')],
        [('X', 1300000, 1300001, '+'), ('X', 1600000, 1600001, '-')],
        [('X', 1700000, 1700001, '-'), ('X', 1900000, 1900001, '+')],
        [('X', 2100000, 2100001, '.'), ('X', 2400000, 2400001, '.')],
    ])),
    Case('count/intra_range', M + 'count_contacts', driver='contacts', args=(LI_H5, [
        [('X', 100000, 100001), ('X', 400000, 400001)],
        [('X', 100000, 100001), ('X', 120000, 120001)],
        [('X', 100000, 100001), ('X', 4000000, 4000001)],
    ]), kwargs={'mode': 'intra-chr', 'range': '50000:900000', 'transform': 'total-counts'}),
    Case('count/small_two_chromosomes', M + 'count_contacts', driver='contacts', args=(SMALL, [
        [('chrX', 1000000, 1000001), ('chr2L', 1000000, 1000001)],
        [('chr2L', 2000000, 2000001), ('chrX', 3000000, 3000001)],
        [('chrX', 1000000, 1000001), ('chr2L', 1000000, 1000001)],
    ]), kwargs={'M_half': 1, 'operation': 'last'}),
    Case('count/three_and_four_fields', M + 'count_contacts', driver='contacts',
         args=(LI_H5, [[('X', 1, 2), ('X', 3, 4, '+')]])),
    Case('count/bad_range', M + 'count_contacts', driver='contacts',
         args=(LI_H5, [[('X', 100000, 100001), ('X', 400000, 400001)]]),
         kwargs={'mode': 'intra-chr', 'range': '50000'}),

    # get_outlier_indices
    Case('outliers/random', M + 'get_outlier_indices', args=lambda: (random_stack((20, 9), 1), 2)),
    Case('outliers/default_deviation', M + 'get_outlier_indices', args=lambda: (random_stack((20, 9), 2) * 100,)),
    Case('outliers/with_outlier', M + 'get_outlier_indices',
         args=lambda: (np.vstack([random_stack((10, 4), 3), np.full((1, 4), 50.0)]), 5)),
    Case('outliers/mostly_zero', M + 'get_outlier_indices', args=lambda: (random_stack((6, 5), 4, zeros=True), 2)),
    Case('outliers/nan', M + 'get_outlier_indices',
         args=lambda: (np.where(random_stack((6, 5), 5) > 1, np.nan, 1.0), 2)),
    Case('outliers/list', M + 'get_outlier_indices', args=([[1.0, 2.0], [3.0, 400.0], [2.0, 2.0]], 3)),
    Case('outliers/one_dimensional', M + 'get_outlier_indices', args=lambda: (random_stack((7,), 6), 2)),
    Case('outliers/matrix', M + 'get_outlier_indices', driver='outliers_of_matrix', args=(7, 1)),

    # compute_clusters
    Case('clusters/kmeans_full', M + 'compute_clusters', driver='clusters', args=(2, 'kmeans')),
    Case('clusters/kmeans_three', M + 'compute_clusters', driver='clusters', args=(3, 'kmeans', 'full', 2, True)),
    Case('clusters/kmeans_diagonal', M + 'compute_clusters', driver='clusters', args=(3, 'kmeans', 'diagonal')),
    Case('clusters/kmeans_center', M + 'compute_clusters', driver='clusters', args=(2, 'kmeans', 'center')),
    Case('clusters/hierarchical', M + 'compute_clusters', driver='clusters', args=(3, 'hierarchical')),
    Case('clusters/hierarchical_diagonal', M + 'compute_clusters', driver='clusters',
         args=(4, 'hierarchical', 'diagonal', 2, False, 'total-counts')),
    Case('clusters/outliers_removed', M + 'compute_clusters', driver='clusters',
         args=(2, 'kmeans', 'full', 0.5, False, 'total-counts')),
    Case('clusters/no_clust', M + 'compute_clusters', driver='clusters', args=(1, 'no_clust')),
    Case('clusters/no_clust_two', M + 'compute_clusters', driver='clusters', args=(2, 'no_clust')),
    Case('clusters/center_too_small', M + 'compute_clusters', driver='clusters', args=(2, 'kmeans', 'center'),
         kwargs={'M_half': 0}),
    Case('clusters/empty', M + 'compute_clusters', driver='clusters_of',
         args=({'coords': [], 'centers': [], 'submatrices': [], 'clustered_dict': [], 'diagonal': []}, 2)),
    Case('clusters/too_few_samples', M + 'compute_clusters', driver='clusters_of',
         args=(lambda: ({'coords': [1], 'centers': [2.0], 'submatrices': [np.asmatrix(np.eye(3))],
                         'clustered_dict': [], 'diagonal': [np.asmatrix(np.ones(3))]}, 2)), compare_args=False,
         note='the argument holds numpy.matrix values, which apieq cannot canonicalise; the result is compared'),

    # cluster_matrices
    Case('cluster_matrices/genome_kmeans', M + 'cluster_matrices', driver='matrices', args=(2, 'kmeans')),
    Case('cluster_matrices/genome_no_clust', M + 'cluster_matrices', driver='matrices', args=(1, 'no_clust')),
    Case('cluster_matrices/per_chr', M + 'cluster_matrices', driver='matrices', args=(2, 'hierarchical', 'full', True)),
    Case('cluster_matrices/per_chr_small', M + 'cluster_matrices', driver='matrices',
         args=(3, 'kmeans', 'diagonal', True, 2, False, SMALL, 'intra-chr')),
    Case('cluster_matrices/genome_small_all', M + 'cluster_matrices', driver='matrices',
         args=(2, 'kmeans', 'center', False, 2, False, SMALL, 'all')),

    # compute_avg
    Case('avg/mean', M + 'compute_avg', args=lambda: (random_stack((7, 5, 5), 8), 'mean')),
    Case('avg/median_odd', M + 'compute_avg', args=lambda: (random_stack((7, 5, 5), 9), 'median')),
    Case('avg/median_even', M + 'compute_avg', args=lambda: (random_stack((8, 3, 3), 10), 'median')),
    Case('avg/sum', M + 'compute_avg', args=lambda: (random_stack((4, 3, 3), 11), 'sum')),
    Case('avg/other_is_sum', M + 'compute_avg', args=lambda: (random_stack((4, 3, 3), 12), 'max')),
    Case('avg/median_zero', M + 'compute_avg', args=lambda: (np.array([np.zeros((3, 3)), np.zeros((3, 3)),
                                                                        np.ones((3, 3))]), 'median')),
    Case('avg/all_zero', M + 'compute_avg', args=lambda: (np.zeros((3, 2, 2)), 'median')),
    Case('avg/nan', M + 'compute_avg', args=lambda: (np.array([[[np.nan, 1.0]], [[2.0, 3.0]]]), 'median')),
    Case('avg/list_of_matrices', M + 'compute_avg', driver='average_of_matrices', args=('mean',)),
    Case('avg/empty', M + 'compute_avg', args=([], 'median')),
    Case('avg/two_dimensional', M + 'compute_avg', args=lambda: (random_stack((5, 4), 13), 'median')),

    # plot_aggregated_contacts and plot_diagnostic_heatmaps (drawing on C++ data)
    Case('plot/one_cluster', M + 'plot_aggregated_contacts', driver='plotting', args=(File(),),
         kwargs={'ties': True}, note=TIES_NOTE),
    Case('plot/kmeans', M + 'plot_aggregated_contacts', driver='plotting', args=(File(), 2, 'kmeans'),
         kwargs={'ties': True}, note=TIES_NOTE),
    Case('plot/median_per_chr', M + 'plot_aggregated_contacts', driver='plotting',
         args=(File(), 1, 'no_clust', True), kwargs={'operation': 'median', 'ties': True}, note=TIES_NOTE),
    Case('plot/sum_3d', M + 'plot_aggregated_contacts', driver='plotting', args=(File(), 2, 'hierarchical'),
         kwargs={'operation': 'sum', 'plot_type': '3d', 'pairs': False}),
    Case('plot/small_all', M + 'plot_aggregated_contacts', driver='plotting', args=(File(), 2, 'kmeans'),
         kwargs={'matrix': SMALL, 'mode': 'all', 'transform': 'total-counts', 'ties': True}, note=TIES_NOTE),
    Case('heatmap/one_cluster', M + 'plot_diagnostic_heatmaps', driver='plotting', args=(File(),),
         kwargs={'diagnostic': True, 'prefix': False, 'pairs': False, 'ties': True}, rtol=1e-12, note=HEAT_NOTE),
    Case('heatmap/kmeans', M + 'plot_diagnostic_heatmaps', driver='plotting', args=(File(), 3, 'kmeans'),
         kwargs={'diagnostic': True, 'prefix': False, 'pairs': False, 'ties': True}, rtol=1e-12, note=HEAT_NOTE),
    Case('heatmap/per_chr', M + 'plot_diagnostic_heatmaps', driver='plotting', args=(File(), 2, 'kmeans', True),
         kwargs={'diagnostic': True, 'prefix': False, 'pairs': False, 'matrix': SMALL, 'ties': True},
         rtol=1e-12, note=HEAT_NOTE),
]
