"""pytest over the cool comparator (comparators/cool.py).

    python -m pytest cpp/scripts/test_cool_comparator.py

Writes two small cool-shaped HDF5 files that differ in exactly one place and
checks the verdict per class: float attributes and float bin columns follow
the declared class (bit identical at E0, E1, E2; relative tolerance at ED, E3,
E4), while integer and string content stays exact at every class.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from comparators import cool  # noqa: E402

SUM = 123456.78901234567
WEIGHT = np.array([0.5, 0.25, np.nan, 1.0 / 3.0])


def _write(path, sum_value=SUM, weight=WEIGHT, nbins_attr=4, counts=None):
    counts = np.array([1.5, 2.0, 3.25]) if counts is None else counts
    with h5py.File(path, "w") as handle:
        handle.attrs["sum"] = sum_value
        handle.attrs["nbins"] = nbins_attr
        handle.attrs["format"] = "HDF5::Cooler"
        handle.create_dataset("chroms/name", data=np.array([b"chr1"]))
        handle.create_dataset("chroms/length", data=np.array([400], dtype=np.int32))
        handle.create_dataset("bins/chrom", data=np.zeros(4, dtype=np.int32))
        handle.create_dataset("bins/start", data=np.arange(0, 400, 100, dtype=np.int32))
        handle.create_dataset("bins/end", data=np.arange(100, 500, 100, dtype=np.int32))
        handle.create_dataset("bins/weight", data=np.asarray(weight, dtype=float))
        handle.create_dataset("pixels/bin1_id", data=np.array([0, 0, 1], dtype=np.int64))
        handle.create_dataset("pixels/bin2_id", data=np.array([0, 1, 1], dtype=np.int64))
        handle.create_dataset("pixels/count", data=counts)
        handle.create_dataset("indexes/bin1_offset", data=np.array([0, 2, 3, 3, 3]))
        handle.create_dataset("indexes/chrom_offset", data=np.array([0, 4]))
    return str(path)


def _pair(tmp_path, **changed):
    return _write(tmp_path / "py.cool"), _write(tmp_path / "cpp.cool", **changed)


LAST_DIGIT_SUM = np.nextafter(SUM, np.inf)


@pytest.mark.parametrize("cls", ["E0", "E1", "E2"])
def test_a_last_digit_sum_difference_fails_the_bit_exact_classes(tmp_path, cls):
    result = cool.compare(*_pair(tmp_path, sum_value=LAST_DIGIT_SUM), cls)
    assert not result.passed
    assert any("attribute sum" in diff for diff in result.diffs)


@pytest.mark.parametrize("cls", ["ED", "E3", "E4"])
def test_a_last_digit_sum_difference_passes_the_tolerant_classes(tmp_path, cls):
    assert cool.compare(*_pair(tmp_path, sum_value=LAST_DIGIT_SUM), cls).passed


def test_a_sum_beyond_three_significant_digits_fails_ed(tmp_path):
    result = cool.compare(*_pair(tmp_path, sum_value=SUM * 1.002), "ED")
    assert not result.passed
    assert any("attribute sum" in diff for diff in result.diffs)


def test_a_sum_beyond_e3_fails_e3_but_passes_ed(tmp_path):
    left, right = _pair(tmp_path, sum_value=SUM * (1 + 1e-9))
    assert not cool.compare(left, right, "E3").passed
    assert cool.compare(left, right, "ED").passed


LAST_DIGIT_WEIGHT = WEIGHT.copy()
LAST_DIGIT_WEIGHT[3] = np.nextafter(WEIGHT[3], 1.0)


@pytest.mark.parametrize("cls", ["E0", "E1", "E2"])
def test_a_last_digit_weight_difference_fails_the_bit_exact_classes(tmp_path, cls):
    result = cool.compare(*_pair(tmp_path, weight=LAST_DIGIT_WEIGHT), cls)
    assert not result.passed
    assert any("bins/weight" in diff for diff in result.diffs)


@pytest.mark.parametrize("cls", ["ED", "E3", "E4"])
def test_a_last_digit_weight_difference_passes_the_tolerant_classes(tmp_path, cls):
    assert cool.compare(*_pair(tmp_path, weight=LAST_DIGIT_WEIGHT), cls).passed


def test_a_weight_turned_nan_fails_ed(tmp_path):
    weight = WEIGHT.copy()
    weight[0] = np.nan
    result = cool.compare(*_pair(tmp_path, weight=weight), "ED")
    assert not result.passed


def test_a_weight_beyond_three_significant_digits_fails_ed(tmp_path):
    weight = WEIGHT.copy()
    weight[1] *= 1.01
    assert not cool.compare(*_pair(tmp_path, weight=weight), "ED").passed


@pytest.mark.parametrize("cls", ["E0", "ED"])
def test_an_integer_attribute_stays_exact_at_every_class(tmp_path, cls):
    result = cool.compare(*_pair(tmp_path, nbins_attr=5), cls)
    assert not result.passed
    assert any("attribute nbins" in diff for diff in result.diffs)


@pytest.mark.parametrize("cls", ["E0", "E1", "E2", "ED"])
def test_identical_files_pass_with_nan_weights(tmp_path, cls):
    assert cool.compare(*_pair(tmp_path), cls).passed
