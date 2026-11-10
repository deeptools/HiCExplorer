"""The hiCMatrix conversion helpers of hicexplorer._hicx (support/convert.hpp)
on real matrices: reading a hiCMatrix into hicx::MatrixData and writing it
back gives the same matrix, bins and attributes, types included."""
import os
import types

import numpy as np
import pytest

hicmatrix = pytest.importorskip('hicmatrix.HiCMatrix')
from hicexplorer import _hicx  # noqa: E402

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test_data')


@pytest.mark.parametrize('name', ['small_test_matrix_50kb_res.h5', 'small_test_matrix.cool', 'Li_et_al_2015.h5'])
def test_matrix_data_round_trip(name):
    hic = hicmatrix.hiCMatrix(os.path.join(TEST_DATA, name))
    target = types.SimpleNamespace(cut_intervals=hic.cut_intervals)
    _hicx._roundtrip_matrix_data(hic, target)
    assert target.matrix.format == 'csr'
    assert target.matrix.dtype == hic.matrix.dtype
    assert target.matrix.shape == hic.matrix.shape
    assert np.array_equal(target.matrix.indptr, hic.matrix.indptr)
    assert np.array_equal(target.matrix.indices, hic.matrix.indices)
    assert np.array_equal(target.matrix.data, hic.matrix.data)
    assert target.cut_intervals == hic.cut_intervals
    assert [tuple(type(v) for v in t) for t in target.cut_intervals] == \
        [tuple(type(v) for v in t) for t in hic.cut_intervals]
    assert np.array_equal(target.nan_bins, np.asarray(hic.nan_bins))
    for attribute in ('correction_factors', 'distance_counts'):
        original = getattr(hic, attribute)
        if original is None:
            assert getattr(target, attribute) is None
        else:
            assert np.array_equal(getattr(target, attribute), np.asarray(original, dtype=float))
