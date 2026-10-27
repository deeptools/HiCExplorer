"""Equivalence cases for hicexplorer.lib.viewpoint.Viewpoint (apieq.py).

The inputs are the repository's cHi-C test data: the matrices
FL-E13-5_chr1.cool and MB-E10-5_chr1.cool (read with hicmatrix in each
interpreter), the reference point files, the background model, and the HDF5
files the cHi-C tools wrote (interaction, significant, target, aggregate and
differential files). Constructed inputs cover what the real data does not
reach: missing datasets, scalar datasets, other dtypes, empty files.

Each case calls the method on a Viewpoint built inside the driver (the
hiCMatrix object itself is not comparable between interpreters), and returns
what the method returned together with the arguments it changes. Files the
methods write are compared by content (text) or by an h5py dump of every
attribute and dataset with dtype, shape, compression and values.

Behaviour that is not reproduced (cpp/python/src/api/lib_viewpoint.cpp):
none of the cases below differs. createUniqueHDFGroup with an existing name
and readBackgroundDataFile with an increment of 0 loop forever in 3.7.6 and
in the C++, so they are not run.
"""
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
CHIC = os.path.join(TEST_DATA, 'cHi-C')
V = 'hicexplorer.lib.viewpoint.Viewpoint.'

FL = os.path.join(CHIC, 'FL-E13-5_chr1.cool')
MB = os.path.join(CHIC, 'MB-E10-5_chr1.cool')
INTERACTIONS = os.path.join(CHIC, 'chicViewpoint', 'two_matrices.hdf5')
SIGNIFICANT = os.path.join(CHIC, 'chicSignificantInteractions', 'significantInteractions_dual.hdf5')
TARGET_DUAL = os.path.join(CHIC, 'chicSignificantInteractions', 'targetFile_dual.hdf5')
TARGET_SINGLE = os.path.join(CHIC, 'chicSignificantInteractions', 'targetFile_single.hdf5')
AGGREGATE = os.path.join(CHIC, 'chicAggregateStatistic', 'aggregate.hdf5')
DIFFERENTIAL = os.path.join(CHIC, 'chicDifferentialTest', 'differential.hdf5')
BACKGROUND = os.path.join(CHIC, 'background.txt')

# a HiCExplorer h5 matrix: numpy.str_, numpy.int64 and numpy.float64 bins and
# int64 counts
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')

EYA1 = ('chr1', '14300280', '14300280')
SOX17 = ('chr1', '4487435', '4487435')
WIDE = ('chr1', '24202042', '24205042')
NEAR_END = ('chr1', '197145200', '197145200')
NEAR_START = ('chr1', '150500', '150500')


# ---------------------------------------------------------------------------
# drivers

_MATRICES = {}


def _matrix(path):
    """hiCMatrix(path), loaded once per interpreter: loading a 1 kb chr1 cool
    builds 197,000 intervals and takes seconds, and no Viewpoint method changes
    the matrix (getBinSize caches the same bin size in both runs)."""
    import hicmatrix.HiCMatrix as hm
    if path not in _MATRICES:
        _MATRICES[path] = hm.hiCMatrix(path)
    return _MATRICES[path]


def _viewpoint(matrix=None):
    from hicexplorer.lib.viewpoint import Viewpoint
    return Viewpoint(_matrix(matrix) if matrix is not None else None)


def exits(function):
    """exit(1) raises SystemExit, which apieq does not catch: return it."""
    def wrapper(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except SystemExit as error:
            return 'SystemExit', error.code
    wrapper.__name__ = function.__name__
    return wrapper


@exits
def on_matrix(target, matrix, *args, **kwargs):
    """target(Viewpoint(hiCMatrix(matrix)), *args)"""
    return target(_viewpoint(matrix), *args, **kwargs)


@exits
def no_matrix(target, *args, **kwargs):
    """target(Viewpoint(), *args), for the methods that do not use the matrix"""
    return target(_viewpoint(), *args, **kwargs)


def with_state(target, matrix, *args):
    """The result and the arguments, for methods that change them."""
    result = target(_viewpoint(matrix), *args)
    return result, list(args)


def construct(target, value):
    from hicexplorer.lib.viewpoint import Viewpoint
    viewpoint = Viewpoint(value)
    return type(viewpoint).__name__, viewpoint.hicMatrix, target is Viewpoint.__init__


def construct_default(target):
    from hicexplorer.lib.viewpoint import Viewpoint
    return Viewpoint().hicMatrix


def write_text(target, path, *args):
    """A method writing path + '.txt': its return value and the file."""
    target(_viewpoint(), path, *args)
    with open(path + '.txt', 'rb') as handle:
        return handle.read()


def write_plot_data(target, data, path, background):
    target(_viewpoint(), data, path, background)
    with open(path + '.txt', 'rb') as handle:
        return handle.read()


def dump_h5(path):
    """Every group, attribute and dataset of an HDF5 file with h5py."""
    import h5py
    out = []

    def attributes(obj):
        return sorted((k, repr(v), type(v).__name__) for k, v in obj.attrs.items())

    with h5py.File(path, 'r') as f:
        out.append(('/', attributes(f)))

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                value = obj[()]
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                out.append((name, 'dataset', obj.dtype.str, obj.shape, obj.compression, obj.compression_opts,
                            obj.chunks, value, attributes(obj)))
            else:
                out.append((name, 'group', attributes(obj)))
        f.visititems(visit)
    return out


def interaction_hdf5(target, path, name, data, reference_point):
    import h5py
    with h5py.File(path, 'w') as f:
        group = f.create_group('matrix/chr1')
        result = target(_viewpoint(), group, name, data, reference_point)
    return result, dump_h5(path)


def unique_group(target, path, names):
    import h5py
    results = []
    with h5py.File(path, 'w') as f:
        for name in names:
            group, file_name = target(_viewpoint(), f, name)
            results.append((group.name, file_name))
    return results, dump_h5(path)


def make_h5(path, layout):
    """An HDF5 file with the datasets of `layout`: {path: value}."""
    import h5py
    with h5py.File(path, 'w') as f:
        for name, value in layout.items():
            if value is None:
                f.create_group(name)
            else:
                f[name] = value


def on_constructed(target, path, layout, *args):
    make_h5(path, layout)
    return target(_viewpoint(), path, *args)


def differential(target, path, quadruple):
    result = target(_viewpoint(), path, quadruple)
    return [list(item) for item in result]


def differential_constructed(target, path, layout, quadruple):
    make_h5(path, layout)
    return differential(target, path, quadruple)


def computed_viewpoint(target, matrix, reference_point, gene, background=None):
    """createInteractionFileData(HDF5) on the viewpoint computeViewpoint gives."""
    viewpoint = _viewpoint(matrix)
    region_start, region_end, _ = viewpoint.calculateViewpointRange(reference_point, (20000, 20000))
    data, _ = viewpoint.computeViewpoint(reference_point, reference_point[0], region_start, region_end)
    relative = viewpoint.computeRelativeValues(data)
    if target.__name__ == 'createInteractionFileData':
        return target(viewpoint, reference_point, reference_point[0], region_start, region_end, relative, data,
                      gene, 123.5)
    pvalues = list(np.linspace(0.0, 1.0, len(data)))
    return target(viewpoint, reference_point, reference_point[0], region_start, region_end, relative, data,
                  gene, 123.5, pvalues, list(relative * 2))


def p_values_on_matrix(target, matrix, reference_point, dtype, smooth):
    viewpoint = _viewpoint(matrix)
    background = viewpoint.readBackgroundDataFile(BACKGROUND, (200000, 200000), 500000)
    region_start, region_end, _ = viewpoint.calculateViewpointRange(reference_point, (200000, 200000))
    data, index = viewpoint.computeViewpoint(reference_point, reference_point[0], region_start, region_end)
    if smooth:
        data = viewpoint.smoothInteractionValues(data, 5)
    return target(viewpoint, background, np.asarray(data, dtype=dtype), index)


def background_band(target, data, xfold):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axis = plt.subplots()
    lines = target(_viewpoint(), axis, data, xfold)
    out = [line.get_xydata().tolist() for line in lines]
    out.append([collection.get_paths()[0].vertices.tolist() for collection in axis.collections])
    plt.close(figure)
    return out


def viewpoint_lines(target, data, color, label, highlight, significant):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axis = plt.subplots()
    lines = target(_viewpoint(), axis, data, color, label, highlight, significant)
    out = [line.get_xydata().tolist() for line in lines], [line.get_label() for line in lines], len(axis.patches)
    plt.close(figure)
    return out


def p_value_plot(target, data, levels):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, (axis, label) = plt.subplots(1, 2)
    target(_viewpoint(), axis, label, data, 'sample', 'RdYlBu', figure, levels)
    colorbar_axes = [a for a in figure.axes if a not in (axis, label)]
    ticks = [list(a.get_yticks()) for a in colorbar_axes]
    collections = axis.collections
    arrays = [np.asarray(c.get_array()).tolist() for c in collections if c.get_array() is not None]
    plt.close(figure)
    return ticks, levels, arrays, [t.get_text() for t in label.texts]


def data_for_plotting(target, path, triplet, prange, background_range, resolution, mean=True):
    viewpoint = _viewpoint()
    background = None
    if background_range is not None:
        background = viewpoint.readBackgroundDataFile(BACKGROUND, background_range, 500000, pMean=mean)
    result = target(viewpoint, path, triplet, prange, background, resolution)
    return result, background


def merge(target, path, triplet, xfold, threshold):
    viewpoint = _viewpoint()
    data = viewpoint.readInteractionFile(path, triplet)
    accepted = {key: line for key, line in data[1].items() if data[0][key][-1] >= xfold}
    result = target(viewpoint, accepted, threshold)
    return result, accepted


def merge_dict(target, scores, *args):
    result = target(_viewpoint(), scores, *args)
    return result, scores


def background_data(target, prange, mean):
    viewpoint = _viewpoint()
    background = viewpoint.readBackgroundDataFile(BACKGROUND, (200000, 200000), 500000, pMean=mean)
    return target(viewpoint, background, prange)


# ---------------------------------------------------------------------------
# inputs

REFERENCE_POINTS = (b"chr1\t100\t200\tgeneA\n"
                    b"\n"
                    b"chr1\t300\tthree\n"
                    b"chr2\t1\n"
                    b"  chr3\t5\t6\tgeneB\textra\t \n"
                    b"chrX\t7\t8\t\n")

BACKGROUND_SMALL = (b"Relative position\tsize nbinom\tprob nbinom\tmax value\tmean value\n"
                    b"-2000\t1.5\t0.5\t3.0\t0.1\n"
                    b"-1000\t2.5\t0.25\t4.0\t0.2\n"
                    b"0\t3.5\t0.125\t5.0\t0.3\n"
                    b"1000\t4.5\t0.0625\t6.0\t0.4\n"
                    b"2000\t5.5\t0.03125\t7.0\t0.5\n")

INTERACTION_LINES = [
    ('chr1', 14290000, 14291000, 'Eya1', '123.5', -10000, 0.125, 3.0),
    ('chr1', 14300000, 14301000, 'Eya1', '123.5', 0, 0.5, 12.0),
    ('chr1', np.int64(14310000), np.int64(14311000), 'Eya1', 123.5, 10000, np.float64(0.25), np.float32(6.0)),
]

PLOT_DATA = {
    1000: ['chr1', 1, 2, 'g', 'chr1', 3, 4, 1000, 0.5, -1.0, 7.0, 0.25, 0.125],
    -1000: ['chr1', 1, 2, 'g', 'chr1', 0, 1, -1000, 0.25, 1.5, 3.0, np.float64(0.5), 0.0625],
}


def scores_dictionary():
    """merge_neighbors input with runs, gaps and an unmerged tail."""
    rows = {}
    positions = [-5000, -4000, -3000, 0, 5000, 6000, 20000, 40000]
    for i, position in enumerate(positions):
        rows[np.int64(position)] = ['chr1', 1000 * i, 1000 * i + 1000, 'g', 10.0 + i, position,
                                    0.1 * i, 0.01 * i, 1.5 + i, float((7 * i) % 5)]
    return (rows,)


def interaction_layout(**overrides):
    base = 'm/chr1/g'
    layout = {
        base + '/relative_position_list': np.array([-1000, 0, 1000], dtype=np.int64),
        base + '/interaction_data_list': np.array([0.1, 0.5, 0.25]),
        base + '/pvalue': np.array([0.5, 0.01, 0.2]),
        base + '/raw': np.array([3.0, 12.0, 6.0]),
        base + '/xfold': np.array([1.0, 2.5, 1.5]),
        base + '/start_list': np.array([9000, 10000, 11000], dtype=np.int64),
        base + '/end_list': np.array([10000, 11000, 12000], dtype=np.int64),
        base + '/chromosome': 'chr1',
        base + '/gene': 'g',
        base + '/sum_of_interactions': 42.0,
        base + '/reference_point_start': np.int64(10500),
        base + '/reference_point_end': np.int64(10500),
    }
    for key, value in overrides.items():
        path = base + '/' + key
        if value is None:
            layout.pop(path, None)
        else:
            layout[path] = value
    return layout


DIFFERENTIAL_EYA1 = ['FL-E13-5_chr1', 'MB-E10-5_chr1', 'chr1', 'Eya1']


CASES = [
    # __init__ ----------------------------------------------------------------
    Case('__init__/value', V + '__init__', args=(5,), driver='construct'),
    Case('__init__/default', V + '__init__', driver='construct_default'),

    # readReferencePointFile ---------------------------------------------------
    Case('readReferencePointFile/real', V + 'readReferencePointFile', driver='no_matrix',
         args=(os.path.join(CHIC, 'referencePoints.bed'),)),
    Case('readReferencePointFile/edge_cases', V + 'readReferencePointFile', driver='no_matrix',
         args=(os.path.join(CHIC, 'referencePoints_edge_cases.bed'),)),
    Case('readReferencePointFile/no_gene', V + 'readReferencePointFile', driver='no_matrix',
         args=(os.path.join(CHIC, 'referencePoints_chicViewpoint.bed'), False)),
    Case('readReferencePointFile/constructed', V + 'readReferencePointFile', driver='no_matrix',
         args=(File('r.bed'),), files={'r.bed': REFERENCE_POINTS}),
    Case('readReferencePointFile/constructed_no_gene', V + 'readReferencePointFile', driver='no_matrix',
         args=(File('r.bed'), 0), files={'r.bed': REFERENCE_POINTS}),
    Case('readReferencePointFile/crlf', V + 'readReferencePointFile', driver='no_matrix',
         args=(File('r.bed'),), files={'r.bed': b"chr1\t1\t2\tg\r\nchr1\t3\tx\r"}),
    Case('readReferencePointFile/empty', V + 'readReferencePointFile', driver='no_matrix',
         args=(File('r.bed'),), files={'r.bed': b""}),
    Case('readReferencePointFile/missing', V + 'readReferencePointFile', driver='no_matrix',
         args=(File('missing.bed'),)),

    # readInteractionFile ------------------------------------------------------
    Case('readInteractionFile/real', V + 'readInteractionFile', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Eya1'])),
    Case('readInteractionFile/real_tuple', V + 'readInteractionFile', driver='no_matrix',
         args=(INTERACTIONS, ('MB-E10-5_chr1', 'chr1', 'Tfap2d'))),
    Case('readInteractionFile/significant', V + 'readInteractionFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'chr1', 'Sox17'])),
    Case('readInteractionFile/through_genes', V + 'readInteractionFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'genes', 'Sox17'])),
    Case('readInteractionFile/missing_path', V + 'readInteractionFile', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Nope'])),
    Case('readInteractionFile/target_file', V + 'readInteractionFile', driver='no_matrix',
         args=(TARGET_DUAL, ['FL-E13-5_chr1', 'MB-E10-5_chr1', 'chr1', 'Eya1'])),
    Case('readInteractionFile/dataset_path', V + 'readInteractionFile', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Eya1', 'pvalue'])),
    Case('readInteractionFile/missing_file', V + 'readInteractionFile', driver='no_matrix',
         args=(File('missing.hdf5'), ['a'])),
    Case('readInteractionFile/not_hdf5', V + 'readInteractionFile', driver='no_matrix',
         args=(File('x.hdf5'), ['a']), files={'x.hdf5': b'not an hdf5 file'}),
    Case('readInteractionFile/constructed', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/float32_int32', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(
             relative_position_list=np.array([-1000, 0, 1000], dtype=np.int32),
             raw=np.array([3.0, 12.0, 6.0], dtype=np.float32)), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/repeated_position', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(
             relative_position_list=np.array([0, -1000, 0], dtype=np.int64)), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/no_chromosome', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(chromosome=None), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/no_gene_no_sum', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(gene=None, sum_of_interactions=None), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/no_raw', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(raw=None), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/short_xfold', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(xfold=np.array([1.0, 2.0])), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/scalar_pvalue', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(pvalue=0.5), ['m', 'chr1', 'g'])),
    Case('readInteractionFile/nan_start', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(start_list=np.array([9000.0, np.nan, 11000.0])),
               ['m', 'chr1', 'g'])),
    Case('readInteractionFile/empty_arrays', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(relative_position_list=np.array([], dtype=np.int64)),
               ['m', 'chr1', 'g'])),
    Case('readInteractionFile/no_reference_point', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), interaction_layout(reference_point_start=None, reference_point_end=None),
               ['m', 'chr1', 'g'])),
    Case('readInteractionFile/empty_group', V + 'readInteractionFile', driver='on_constructed',
         args=(File('i.hdf5'), {'m/chr1/g': None}, ['m', 'chr1', 'g'])),
    Case('readInteractionFile/not_str', V + 'readInteractionFile', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 1])),

    # interactionDataForPlot ---------------------------------------------------
    Case('interactionDataForPlot/real', V + 'interactionDataForPlot', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Sox17'])),
    Case('interactionDataForPlot/missing', V + 'interactionDataForPlot', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chrX', 'Sox17'])),

    # readBackgroundDataFile ---------------------------------------------------
    Case('readBackgroundDataFile/real', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(BACKGROUND, (200000, 200000), 500000)),
    Case('readBackgroundDataFile/real_mean', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(BACKGROUND, [200000, 200000], 500000, True)),
    Case('readBackgroundDataFile/extended', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(BACKGROUND, (650000, 700000), 500000)),
    Case('readBackgroundDataFile/fixate_inside', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(BACKGROUND, (300000, 300000), 250000, True),
         note='the extension overwrites the stored positions beyond fixateRange'),
    Case('readBackgroundDataFile/fixate_not_a_key', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(BACKGROUND, (300000, 300000), 250500)),
    Case('readBackgroundDataFile/small', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('b.txt'), (5000, 3500), 10000), files={'b.txt': BACKGROUND_SMALL}),
    Case('readBackgroundDataFile/small_mean', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('b.txt'), (5000, 3500), 10000, 1), files={'b.txt': BACKGROUND_SMALL}),
    Case('readBackgroundDataFile/header_only', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('b.txt'), (5000, 5000), 10000), files={'b.txt': BACKGROUND_SMALL.split(b'\n')[0] + b'\n'}),
    Case('readBackgroundDataFile/one_line', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('b.txt'), (5000, 5000), 10000), files={'b.txt': b'h\n0\t1\t2\t3\t4\n'}),
    Case('readBackgroundDataFile/short_line', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('b.txt'), (5000, 5000), 10000), files={'b.txt': b'h\n0\t1\n'}),
    Case('readBackgroundDataFile/bad_number', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('b.txt'), (5000, 5000), 10000), files={'b.txt': b'h\nzero\t1\t2\t3\n'}),
    Case('readBackgroundDataFile/missing', V + 'readBackgroundDataFile', driver='no_matrix',
         args=(File('missing.txt'), (5000, 5000), 10000)),

    # writeInteractionFile -----------------------------------------------------
    Case('writeInteractionFile/default', V + 'writeInteractionFile', driver='write_text',
         args=(File('out'), INTERACTION_LINES, '# header\tline', [0.5, 0.25, np.float32(0.125)], [1.0, 2.0, 3.0])),
    Case('writeInteractionFile/decimals', V + 'writeInteractionFile', driver='write_text',
         args=(File('out'), INTERACTION_LINES, 'h', [0.5, 0.25, 0.125], [1.0, 2.0, 3.0], 3)),
    Case('writeInteractionFile/short_pvalues', V + 'writeInteractionFile', driver='write_text',
         args=(File('out'), INTERACTION_LINES, 'h', [0.5], [1.0, 2.0, 3.0])),

    # writeInteractionFileHDF5 / createUniqueHDFGroup --------------------------
    Case('writeInteractionFileHDF5/lists', V + 'writeInteractionFileHDF5', driver='interaction_hdf5',
         args=(File('w.hdf5'), 'Eya1', ['chr1', [1, 2], [2, 3], 'Eya1', 12.5, [-1000, 0], [0.5, 0.25],
                                        [0.1, 0.2], [1.5, 2.5], [3.0, 4.0]], ['14300280', 14300280])),
    Case('writeInteractionFileHDF5/arrays', V + 'writeInteractionFileHDF5', driver='interaction_hdf5',
         args=lambda: (File('w.hdf5'), 'g', [np.str_('chr2'), np.arange(5), np.arange(5) + 1, 'g', np.float64(3.5),
                                            np.arange(-2, 3), np.linspace(0, 1, 5), np.zeros(5),
                                            np.ones(5, dtype=np.float32), np.arange(5.0)], (7, 8))),
    Case('writeInteractionFileHDF5/empty_lists', V + 'writeInteractionFileHDF5', driver='interaction_hdf5',
         args=(File('w.hdf5'), 'g', ['chr1', [], [], 'g', 0.0, [], [], [], [], []], [1, 1])),
    Case('createUniqueHDFGroup/new', V + 'createUniqueHDFGroup', driver='unique_group',
         args=(File('u.hdf5'), ['a', 'b/c', 'b/d'])),

    # computeViewpoint ---------------------------------------------------------
    Case('computeViewpoint/eya1', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, EYA1, 'chr1', 14100280, 14500280)),
    Case('computeViewpoint/wide', V + 'computeViewpoint', driver='on_matrix',
         args=(MB, WIDE, 'chr1', 24102042, 24305042)),
    Case('computeViewpoint/two_fields', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, ('chr1', 4487435), 'chr1', 4400000, 4600000)),
    Case('computeViewpoint/near_start', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, NEAR_START, 'chr1', 0, 250500)),
    Case('computeViewpoint/reference_outside_region', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, EYA1, 'chr1', 14310000, 14500280)),
    Case('computeViewpoint/reversed_reference', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, ('chr1', '14305000', '14300000'), 'chr1', 14200000, 14400000)),
    Case('computeViewpoint/position_outside', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, ('chr1', '300000000', '300000000'), 'chr1', 14200000, 14400000)),
    Case('computeViewpoint/unknown_chromosome', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, ('chrZ', '100', '100'), 'chrZ', 0, 1000)),
    Case('computeViewpoint/bad_reference', V + 'computeViewpoint', driver='on_matrix',
         args=(FL, ('chr1',), 'chr1', 0, 1000)),

    Case('computeViewpoint/h5_matrix', V + 'computeViewpoint', driver='on_matrix',
         args=(H5, ('chr2L', '1000000', '1100000'), 'chr2L', 500000, 1600000)),

    # createInteractionFileData(HDF5) ------------------------------------------
    Case('createInteractionFileData/eya1', V + 'createInteractionFileData', driver='computed_viewpoint',
         args=(FL, EYA1, 'Eya1')),
    Case('createInteractionFileData/wide', V + 'createInteractionFileData', driver='computed_viewpoint',
         args=(MB, WIDE, 'Wide')),
    Case('createInteractionFileDataHDF5/eya1', V + 'createInteractionFileDataHDF5', driver='computed_viewpoint',
         args=(FL, EYA1, 'Eya1')),
    Case('createInteractionFileDataHDF5/wide', V + 'createInteractionFileDataHDF5', driver='computed_viewpoint',
         args=(MB, WIDE, 'Wide')),
    Case('createInteractionFileDataHDF5/near_start', V + 'createInteractionFileDataHDF5',
         driver='computed_viewpoint', args=(FL, NEAR_START, 'NearStart')),
    Case('createInteractionFileDataHDF5/short_data', V + 'createInteractionFileDataHDF5', driver='on_matrix',
         args=(FL, EYA1, 'chr1', 14290000, 14310000, [0.5, 0.25], [2, 1], 'g', 7, [0.1], [1.0])),
    Case('createInteractionFileDataHDF5/h5_matrix', V + 'createInteractionFileDataHDF5', driver='on_matrix',
         args=(H5, ('chr2L', '1000000', '1100000'), 'chr2L', 900000, 1200000, [0.5, 0.25, 0.1, 0.1, 0.2],
               [2, 1, 1, 1, 3], 'g', 7, [0.1], [1.0])),
    Case('createInteractionFileData/h5_matrix', V + 'createInteractionFileData', driver='on_matrix',
         args=(H5, ('chr2L', '1000000', '1100000'), 'chr2L', 900000, 1200000, [0.5, 0.25, 0.1, 0.1, 0.2],
               [2, 1, 1, 1, 3], 'g', np.float64(7))),
    Case('createInteractionFileData/past_the_end', V + 'createInteractionFileData', driver='on_matrix',
         args=(FL, NEAR_END, 'chr1', 197100000, 197195000, list(range(200)), list(range(200)), 'g', 7)),

    # getViewpointRangeAsMatrixIndices / getReferencePointAsMatrixIndices ------
    Case('getViewpointRangeAsMatrixIndices/real', V + 'getViewpointRangeAsMatrixIndices', driver='on_matrix',
         args=(FL, 'chr1', 14100280, 14500280)),
    Case('getViewpointRangeAsMatrixIndices/outside', V + 'getViewpointRangeAsMatrixIndices', driver='on_matrix',
         args=(FL, 'chr1', 14100280, 400000000)),
    Case('getReferencePointAsMatrixIndices/three', V + 'getReferencePointAsMatrixIndices', driver='on_matrix',
         args=(FL, WIDE)),
    Case('getReferencePointAsMatrixIndices/two', V + 'getReferencePointAsMatrixIndices', driver='on_matrix',
         args=(FL, ['chr1', 14300280])),
    Case('getReferencePointAsMatrixIndices/four', V + 'getReferencePointAsMatrixIndices', driver='on_matrix',
         args=(FL, ['chr1', 1, 2, 3])),
    Case('getReferencePointAsMatrixIndices/not_a_number', V + 'getReferencePointAsMatrixIndices',
         driver='on_matrix', args=(FL, ['chr1', 'x', '2'])),
    Case('getReferencePointAsMatrixIndices/outside', V + 'getReferencePointAsMatrixIndices', driver='on_matrix',
         args=(FL, ['chr1', '400000000'])),

    # smoothInteractionValues --------------------------------------------------
    Case('smoothInteractionValues/odd', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.random.default_rng(3).random(41), 5)),
    Case('smoothInteractionValues/even', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.random.default_rng(4).random(40) * 100, 4)),
    Case('smoothInteractionValues/float32', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.random.default_rng(5).random(33).astype(np.float32), 5)),
    Case('smoothInteractionValues/list', V + 'smoothInteractionValues', driver='no_matrix',
         args=([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0], 3)),
    Case('smoothInteractionValues/ints', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.arange(12), 2)),
    Case('smoothInteractionValues/window_one', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.arange(6.0), 1)),
    Case('smoothInteractionValues/window_float', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.arange(9.0) ** 2, 5.0)),
    Case('smoothInteractionValues/longer_than_data', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.arange(3.0), 9)),
    Case('smoothInteractionValues/empty', V + 'smoothInteractionValues', driver='no_matrix',
         args=lambda: (np.array([]), 5)),

    # computeRelativeValues ----------------------------------------------------
    Case('computeRelativeValues/sum', V + 'computeRelativeValues', driver='no_matrix',
         args=lambda: (np.random.default_rng(6).random(50),)),
    Case('computeRelativeValues/denominator', V + 'computeRelativeValues', driver='no_matrix',
         args=lambda: (np.arange(10), 7.5)),
    Case('computeRelativeValues/zero_denominator', V + 'computeRelativeValues', driver='no_matrix',
         args=([1, 2, 3], 0)),
    Case('computeRelativeValues/nan_denominator', V + 'computeRelativeValues', driver='no_matrix',
         args=([1.0, 2.0], float('nan'))),
    Case('computeRelativeValues/float32_in_place', V + 'computeRelativeValues', driver='no_matrix',
         args=lambda: (np.array([1.0, 3.0], dtype=np.float32), np.float32(2.0))),
    Case('computeRelativeValues/zeros', V + 'computeRelativeValues', driver='no_matrix', args=([0.0, 0.0],)),

    # calculateViewpointRange --------------------------------------------------
    Case('calculateViewpointRange/eya1', V + 'calculateViewpointRange', driver='on_matrix',
         args=(FL, EYA1, (200000, 200000))),
    Case('calculateViewpointRange/near_start', V + 'calculateViewpointRange', driver='on_matrix',
         args=(FL, NEAR_START, [500000, 500000])),
    Case('calculateViewpointRange/near_end', V + 'calculateViewpointRange', driver='on_matrix',
         args=(FL, NEAR_END, (500000, 500000))),
    Case('calculateViewpointRange/h5_near_end', V + 'calculateViewpointRange', driver='on_matrix',
         args=(H5, ('chr2L', '23000000', '23000000'), (500000, 500000))),
    Case('calculateViewpointRange/h5_near_start', V + 'calculateViewpointRange', driver='on_matrix',
         args=(H5, ('chr2L', '20000', '30000'), (500000, 500000))),
    Case('calculateViewpointRange/unknown_chromosome', V + 'calculateViewpointRange', driver='on_matrix',
         args=(FL, ('chr7', '1', '1'), (5, 5))),

    # getDataForPlotting -------------------------------------------------------
    Case('getDataForPlotting/no_background', V + 'getDataForPlotting', driver='data_for_plotting',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Eya1'], [200000, 200000], None, 1000)),
    Case('getDataForPlotting/background', V + 'getDataForPlotting', driver='data_for_plotting',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Eya1'], [100000, 150000], (200000, 200000), 1000)),
    Case('getDataForPlotting/narrow_resolution', V + 'getDataForPlotting', driver='data_for_plotting',
         args=(INTERACTIONS, ['MB-E10-5_chr1', 'chr1', 'Sox17'], [5000, 5000], (200000, 200000), 250)),
    Case('getDataForPlotting/narrow_resolution_no_background', V + 'getDataForPlotting',
         driver='data_for_plotting',
         args=(INTERACTIONS, ['MB-E10-5_chr1', 'chr1', 'Sox17'], [5000, 3000], None, 300)),
    Case('getDataForPlotting/no_range', V + 'getDataForPlotting', driver='data_for_plotting',
         args=(INTERACTIONS, ['MB-E10-5_chr1', 'chr1', 'Sox17'], None, None, 1000)),
    Case('getDataForPlotting/no_range_background', V + 'getDataForPlotting', driver='data_for_plotting',
         args=(INTERACTIONS, ['MB-E10-5_chr1', 'chr1', 'Sox17'], None, (200000, 200000), 1000)),
    Case('getDataForPlotting/missing', V + 'getDataForPlotting', driver='data_for_plotting',
         args=(INTERACTIONS, ['MB-E10-5_chr1', 'chr1', 'Nope'], [5000, 5000], None, 1000)),

    # plotViewpoint / plotBackgroundModel / plotPValue (drawing) ---------------
    Case('plotViewpoint/regions', V + 'plotViewpoint', driver='viewpoint_lines',
         args=([0.1, 0.4, 0.2, 0.3], 'g', 'sample', [[0, 1]], [[2.0, 3.5], [1, 2]])),
    Case('plotViewpoint/plain', V + 'plotViewpoint', driver='viewpoint_lines',
         args=(np.arange(5.0), 'b', 'x', None, None)),
    Case('plotBackgroundModel/xfold', V + 'plotBackgroundModel', driver='background_band',
         args=([0.1, 0.4, 0.2, 0.3], 1.5)),
    Case('plotBackgroundModel/float32', V + 'plotBackgroundModel', driver='background_band',
         args=lambda: (np.array([0.1, 0.4, 0.2], dtype=np.float32), 2.5)),
    Case('plotBackgroundModel/ints', V + 'plotBackgroundModel', driver='background_band',
         args=([1, 2, 3], 2)),
    Case('plotBackgroundModel/no_xfold', V + 'plotBackgroundModel', driver='background_band',
         args=([0.1, 0.4, 0.2], None)),
    Case('plotPValue/values', V + 'plotPValue', driver='p_value_plot',
         args=lambda: (np.array([0.5, 0.01, 0.2, 1.0], dtype=np.float32), None)),
    Case('plotPValue/levels', V + 'plotPValue', driver='p_value_plot',
         args=lambda: (None, [0.01, 0.05, 0.1])),

    # writePlotData ------------------------------------------------------------
    Case('writePlotData/background', V + 'writePlotData', driver='write_plot_data',
         args=lambda: (dict(PLOT_DATA), File('plot'), {'x': 1})),
    Case('writePlotData/no_background', V + 'writePlotData', driver='write_plot_data',
         args=lambda: (dict(PLOT_DATA), File('plot'), None)),

    # interactionBackgroundData ------------------------------------------------
    Case('interactionBackgroundData/mean', V + 'interactionBackgroundData', driver='background_data',
         args=((20000, 10000), True)),
    Case('interactionBackgroundData/nbinom', V + 'interactionBackgroundData', driver='background_data',
         args=([5000, 5000], False)),

    # rbz_score ----------------------------------------------------------------
    Case('rbz_score/float64', V + 'rbz_score', driver='no_matrix',
         args=lambda: (np.array([0.5, 0.2, 0.0, 1.0]), np.array([0.25, 0.2, 0.0, 0.5]), np.array([0.1, 0.0, 0.0, 2.0]))),
    Case('rbz_score/float32', V + 'rbz_score', driver='no_matrix',
         args=lambda: (np.array([0.5, 0.3], dtype=np.float32), np.array([0.1, 0.2], dtype=np.float32),
                       np.array([0.3, 0.7], dtype=np.float32))),
    Case('rbz_score/float32_sem64', V + 'rbz_score', driver='no_matrix',
         args=lambda: (np.array([0.5, 0.3], dtype=np.float32), np.array([0.1, 0.2], dtype=np.float32),
                       np.array([0.3, 0.7]))),
    Case('rbz_score/array_and_list', V + 'rbz_score', driver='no_matrix',
         args=lambda: (np.array([0.5, 0.3]), [0.1, 0.2], [0.3, 0.0])),
    Case('rbz_score/lists', V + 'rbz_score', driver='no_matrix', args=([0.5], [0.1], [0.3])),
    Case('rbz_score/ints', V + 'rbz_score', driver='no_matrix',
         args=lambda: (np.array([5, 3]), np.array([1, 2]), np.array([3, 7]))),
    Case('rbz_score/length_mismatch', V + 'rbz_score', driver='no_matrix',
         args=lambda: (np.array([0.5, 0.3]), np.array([0.1]), np.array([0.3, 0.7]))),

    # readRejectedFile ---------------------------------------------------------
    Case('readRejectedFile/eya1', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, DIFFERENTIAL_EYA1 + ['rejected'], 200, 200, 1000, [200000, 200000],
               [np.int64(14300280), np.int64(14300280)])),
    Case('readRejectedFile/narrow_range', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, DIFFERENTIAL_EYA1 + ['all'], 100, 101, 1000, [100000, 30000],
               ['14300280', '14301280'])),
    Case('readRejectedFile/numpy_index', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, DIFFERENTIAL_EYA1 + ['accepted'], np.int64(200), np.int64(201), 500,
               [200000, 200000], [14300280, 14300280])),
    Case('readRejectedFile/out_of_range', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, DIFFERENTIAL_EYA1 + ['rejected'], 200, 200, 1000, [10, 10], [1, 1])),
    Case('readRejectedFile/none', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, None, 200, 200, 1000, [1, 1], [1, 1])),
    Case('readRejectedFile/empty', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, [], 200, 200, 1000, [1, 1], [1, 1])),
    Case('readRejectedFile/missing_path', V + 'readRejectedFile', driver='no_matrix',
         args=(DIFFERENTIAL, ['x', 'y'], 200, 200, 1000, [1, 1], [1, 1])),
    Case('readRejectedFile/not_differential', V + 'readRejectedFile', driver='no_matrix',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Eya1'], 200, 200, 1000, [1, 1], [1, 1])),

    # pvalues ------------------------------------------------------------------
    Case('pvalues/eya1', V + 'pvalues', driver='p_values_on_matrix', args=(FL, EYA1, np.float64, False)),
    Case('pvalues/eya1_float32', V + 'pvalues', driver='p_values_on_matrix', args=(FL, EYA1, np.float32, False)),
    Case('pvalues/sox17_smoothed', V + 'pvalues', driver='p_values_on_matrix', args=(MB, SOX17, np.float64, True)),
    Case('pvalues/wide', V + 'pvalues', driver='p_values_on_matrix', args=(MB, WIDE, np.float64, True)),
    Case('pvalues/list', V + 'pvalues', driver='no_matrix',
         args=({-1: [2.5, 0.25], 0: [3.5, 0.5], 1: [1.5, 0.75]}, [0.0, 1.0, 2, 7.5, 0.25], 2)),
    Case('pvalues/nan_parameters', V + 'pvalues', driver='no_matrix',
         args=({0: [float('nan'), 0.5], 3: [-1.0, 0.5]}, [1.0, 2.0, 3.0, 4.0], 0)),
    Case('pvalues/empty_model', V + 'pvalues', driver='no_matrix', args=({}, [1.0], 0)),

    # merge_neighbors ----------------------------------------------------------
    Case('merge_neighbors/real', V + 'merge_neighbors', driver='merge',
         args=(INTERACTIONS, ['FL-E13-5_chr1', 'chr1', 'Eya1'], 1.5, 1000)),
    Case('merge_neighbors/real_wide_threshold', V + 'merge_neighbors', driver='merge',
         args=(INTERACTIONS, ['MB-E10-5_chr1', 'chr1', 'Tfap2d'], 1.2, 5000)),
    Case('merge_neighbors/constructed', V + 'merge_neighbors', driver='merge_dict', args=scores_dictionary),
    Case('merge_neighbors/constructed_threshold', V + 'merge_neighbors', driver='merge_dict',
         args=lambda: scores_dictionary() + (1500.5,)),
    Case('merge_neighbors/single', V + 'merge_neighbors', driver='merge_dict',
         args=lambda: ({0: ['c', 1, 2, 'g', 1.0, 0, 0.5, 0.1, 2.0, 3.0]},)),
    Case('merge_neighbors/empty', V + 'merge_neighbors', driver='merge_dict', args=lambda: ({},)),
    Case('merge_neighbors/none', V + 'merge_neighbors', driver='merge_dict', args=(None,)),

    # readSignificantRegionsFile -----------------------------------------------
    Case('readSignificantRegionsFile/eya1', V + 'readSignificantRegionsFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'chr1', 'Eya1'], 200, 200, 1000, [200000, 200000],
               [np.int64(14300280), np.int64(14300280)])),
    Case('readSignificantRegionsFile/numpy_index', V + 'readSignificantRegionsFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'chr1', 'Tfap2d'], np.int64(200), np.int64(201), 1000,
               [200000, 200000], [19093103, 19093103])),
    Case('readSignificantRegionsFile/narrow', V + 'readSignificantRegionsFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'chr1', 'Sox17'], 50, 50, 1000, [10000, 10000], [4487435, 4487435])),
    Case('readSignificantRegionsFile/bad_viewpoint', V + 'readSignificantRegionsFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'chr1', 'Sox17'], 50, 50, 1000, [10000, 10000], [None, None, None])),
    Case('readSignificantRegionsFile/missing', V + 'readSignificantRegionsFile', driver='no_matrix',
         args=(SIGNIFICANT, ['FL-E13-5_chr1', 'chr1', 'Nope'], 50, 50, 1000, [10000, 10000], [1, 1])),

    # readTargetHDFFile --------------------------------------------------------
    Case('readTargetHDFFile/dual', V + 'readTargetHDFFile', driver='no_matrix', args=(TARGET_DUAL,)),
    Case('readTargetHDFFile/single', V + 'readTargetHDFFile', driver='no_matrix', args=(TARGET_SINGLE,)),
    Case('readTargetHDFFile/not_a_target_file', V + 'readTargetHDFFile', driver='no_matrix', args=(INTERACTIONS,)),
    Case('readTargetHDFFile/missing', V + 'readTargetHDFFile', driver='no_matrix', args=(File('missing.hdf5'),)),

    # readAggregatedFileHDF ----------------------------------------------------
    Case('readAggregatedFileHDF/eya1', V + 'readAggregatedFileHDF', driver='no_matrix',
         args=(AGGREGATE, ['FL-E13-5_chr1_MB-E10-5_chr1', 'FL-E13-5_chr1', 'chr1', 'Eya1'])),
    Case('readAggregatedFileHDF/through_genes', V + 'readAggregatedFileHDF', driver='no_matrix',
         args=(AGGREGATE, ['FL-E13-5_chr1_MB-E10-5_chr1', 'MB-E10-5_chr1', 'genes', 'Tfap2d'])),
    Case('readAggregatedFileHDF/missing', V + 'readAggregatedFileHDF', driver='no_matrix',
         args=(AGGREGATE, ['FL-E13-5_chr1_MB-E10-5_chr1', 'FL-E13-5_chr1', 'chr1', 'Nope'])),
    Case('readAggregatedFileHDF/missing_list', V + 'readAggregatedFileHDF', driver='on_constructed',
         args=(File('a.hdf5'), {'p/chromosome': 'chr1', 'p/gene_name': 'g', 'p/start_list': [1, 2]}, ['p'])),

    # readDifferentialFile -----------------------------------------------------
    Case('readDifferentialFile/eya1', V + 'readDifferentialFile', driver='differential',
         args=(DIFFERENTIAL, DIFFERENTIAL_EYA1)),
    Case('readDifferentialFile/tfap2d', V + 'readDifferentialFile', driver='differential',
         args=(DIFFERENTIAL, ('FL-E13-5_chr1', 'MB-E10-5_chr1', 'chr1', 'Tfap2d'))),
    Case('readDifferentialFile/missing', V + 'readDifferentialFile', driver='differential',
         args=(DIFFERENTIAL, ['FL-E13-5_chr1', 'x', 'chr1', 'Eya1'])),
    Case('readDifferentialFile/partial', V + 'readDifferentialFile', driver='differential_constructed',
         args=(File('d.hdf5'), {'q/accepted/chromosome': 'chr1', 'q/all': None,
                                'q/rejected/chromosome': 'chr1', 'q/rejected/start_list': [1]}, ['q'])),
]
