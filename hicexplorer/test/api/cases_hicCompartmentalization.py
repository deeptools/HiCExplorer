"""Equivalence cases for hicexplorer.hicCompartmentalization (apieq.py).

The drivers rebuild what hicCompartmentalization.main hands to the
functions: the obs/exp matrix hicPCA/obsexp_norm.h5 as a hiCMatrix, the
first principal component hicCompartmentalization/pca1.bedgraph read by
pandas with main's dtypes, its quantiles (np.nanquantile, np.searchsorted)
and the bin ids of every row (getRegionBinRange, so that count_interactions
is tested on its own). count_interactions also returns the obs/exp matrix
after the call, which --offset changes in place.

Deviations from 3.7.6 (cpp/python/src/api/hicCompartmentalization.cpp):
count_interactions reads float64 and integer matrices only (float32 raises
TypeError naming HiCExplorer 4), within_vs_between_compartments takes a two
dimensional float64 array.
"""
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicCompartmentalization.'
OBS_EXP = os.path.join(TEST_DATA, 'hicPCA', 'obsexp_norm.h5')
PCA = os.path.join(TEST_DATA, 'hicCompartmentalization', 'pca1.bedgraph')


def load(path=OBS_EXP):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(path)


def read_pca(path=PCA, quantiles=30):
    import pandas as pd
    pc1 = pd.read_table(path, header=None, sep="\t", dtype={0: "object", 1: "Int64", 2: "Int64", 3: "float32"})
    pc1 = pc1.rename(columns={0: "chr", 1: "start", 2: "end", 3: "pc1"})
    quantile = [j / (quantiles - 1) for j in range(0, quantiles)]
    quantiled_bins = np.nanquantile(pc1['pc1'].values.astype(float), quantile)
    pc1["quantile"] = np.searchsorted(quantiled_bins, pc1['pc1'].values.astype(float), side="right")
    return pc1


def with_bin_ids(hic, pc1):
    def ids(row):
        a = hic.getRegionBinRange(row['chr'], row['start'], row['end'] - 1)
        return np.arange(a[0], a[1] + 1)
    if len(pc1) == 0:
        pc1["bin_id"] = []
    else:
        pc1["bin_id"] = pc1.apply(ids, axis=1)
    return pc1


def matrix_state(hic):
    m = hic.matrix
    return (m.format, m.shape, m.dtype.str, m.data, m.indices, m.indptr)


# ---------------------------------------------------------------------------
# drivers

def indices_of_rows(target, rows):
    """get_indices on the given rows of pca1.bedgraph as pandas rows."""
    hic = load()
    pc1 = read_pca()
    out = []
    for k in rows:
        try:
            out.append(target(hic, pc1.iloc[k]))
        except Exception as error:  # compared as part of the result
            out.append((type(error).__name__, str(error)))
    return out


def indices_of_dict(target, row):
    return target(load(), row)


def count(target, quantiles=30, offset=None, rows=None, matrix=OBS_EXP, dtype=None, pca_quantiles=None):
    hic = load(matrix)
    if dtype is not None:
        hic.matrix.data = hic.matrix.data.astype(dtype)
    pc1 = read_pca(quantiles=pca_quantiles or (quantiles if quantiles > 1 else 30))
    if rows is not None:
        pc1 = pc1.iloc[rows].reset_index(drop=True)
    if matrix == OBS_EXP:
        pc1 = with_bin_ids(hic, pc1)
    else:
        pc1["bin_id"] = [np.arange(k, k + 2) for k in range(len(pc1))]
    result = target(hic, pc1, quantiles, offset)
    return result, matrix_state(hic)


def count_without_columns(target, columns):
    hic = load()
    pc1 = with_bin_ids(hic, read_pca())
    return target(hic, pc1[columns], 30, None)


def ratios_of_real(target, quantiles=30):
    """within_vs_between_compartments on main's normalised sums."""
    import importlib
    module = importlib.import_module('hicexplorer.hicCompartmentalization')
    hic = load()
    pc1 = with_bin_ids(hic, read_pca(quantiles=quantiles))
    normalised = np.nan_to_num(module.count_interactions(hic, pc1, quantiles, None))
    return target(normalised, quantiles), normalised


def random_values(shape, seed=0):
    return np.random.RandomState(seed).rand(*shape)


CASES = [
    # get_indices
    Case('get_indices/rows', M + 'get_indices', driver='indices_of_rows', args=([0, 1, 2, 100, 448],)),
    Case('get_indices/dict', M + 'get_indices', driver='indices_of_dict',
         args=({'chr': 'chrX', 'start': 100000, 'end': 160000},)),
    Case('get_indices/one_bin', M + 'get_indices', driver='indices_of_dict',
         args=({'chr': 'chrX', 'start': 100, 'end': 200},)),
    Case('get_indices/float_end', M + 'get_indices', driver='indices_of_dict',
         args=({'chr': 'chrX', 'start': 5000, 'end': 20000.5},)),
    Case('get_indices/outside', M + 'get_indices', driver='indices_of_dict',
         args=({'chr': 'chrX', 'start': 100, 'end': 10 ** 10},)),
    Case('get_indices/unknown_chromosome', M + 'get_indices', driver='indices_of_dict',
         args=({'chr': 'chr2', 'start': 100, 'end': 200},)),
    Case('get_indices/missing_key', M + 'get_indices', driver='indices_of_dict',
         args=({'chr': 'chrX', 'start': 100},)),

    # count_interactions
    Case('count/main', M + 'count_interactions', driver='count'),
    Case('count/ten_quantiles', M + 'count_interactions', driver='count', args=(10,)),
    Case('count/offset_main_diagonal', M + 'count_interactions', driver='count', args=(30, [0])),
    Case('count/offsets', M + 'count_interactions', driver='count', args=(12, [0, 1, 3])),
    Case('count/offset_beyond', M + 'count_interactions', driver='count', args=(5, [10 ** 6])),
    Case('count/negative_offset', M + 'count_interactions', driver='count', args=(5, [1, -1])),
    Case('count/some_rows', M + 'count_interactions', driver='count', args=(30, None, list(range(40, 90)))),
    Case('count/no_rows', M + 'count_interactions', driver='count', args=(30, None, [])),
    Case('count/one_quantile', M + 'count_interactions', driver='count', args=(1, None, None, OBS_EXP, None, 4)),
    Case('count/zero_quantiles', M + 'count_interactions', driver='count', args=(0,)),
    Case('count/negative_quantiles', M + 'count_interactions', driver='count', args=(-2,)),
    Case('count/float_quantiles', M + 'count_interactions', driver='count', args=(4.0, None, None, OBS_EXP, None, 4)),
    Case('count/more_quantiles_than_pca', M + 'count_interactions', driver='count',
         args=(40, None, None, OBS_EXP, None, 30)),
    Case('count/int_matrix', M + 'count_interactions', driver='count',
         args=(6, None, list(range(30)), os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5'))),
    Case('count/int_matrix_offset', M + 'count_interactions', driver='count',
         args=(6, [0], list(range(30)), os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5'))),
    Case('count/no_quantile_column', M + 'count_interactions', driver='count_without_columns',
         args=(['chr', 'start', 'end', 'bin_id'],)),
    Case('count/no_chr_column', M + 'count_interactions', driver='count_without_columns',
         args=(['start', 'end', 'quantile', 'bin_id'],)),

    # within_vs_between_compartments
    Case('ratios/main', M + 'within_vs_between_compartments', driver='ratios_of_real'),
    Case('ratios/eight', M + 'within_vs_between_compartments', driver='ratios_of_real', args=(8,)),
    Case('ratios/random', M + 'within_vs_between_compartments', args=lambda: (random_values((30, 30)), 30)),
    Case('ratios/random_large', M + 'within_vs_between_compartments', args=lambda: (random_values((150, 150), 3), 150)),
    Case('ratios/fewer_quantiles', M + 'within_vs_between_compartments', args=lambda: (random_values((10, 10)), 6)),
    Case('ratios/more_quantiles', M + 'within_vs_between_compartments', args=lambda: (random_values((5, 5)), 8)),
    Case('ratios/rectangular', M + 'within_vs_between_compartments', args=lambda: (random_values((4, 7)), 6)),
    Case('ratios/zeros', M + 'within_vs_between_compartments', args=lambda: (np.zeros((4, 4)), 4)),
    Case('ratios/one', M + 'within_vs_between_compartments', args=lambda: (np.ones((3, 3)), 1)),
    Case('ratios/one_dimensional', M + 'within_vs_between_compartments', args=lambda: (np.ones(3), 3)),
]
