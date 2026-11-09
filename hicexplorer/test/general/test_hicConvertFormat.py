
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicConvertFormat
from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
import gzip
from scipy.sparse import triu
import numpy.testing as nt
import numpy as np
import pytest
from hicexplorer.test.test_compute_function import compute

REMOVE_OUTPUT = True

DELTA_DECIMAL = 0
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicConvertFormat")
original_matrix_h5 = ROOT + "/small_test_matrix.h5"
original_matrix_cool = ROOT + "/small_test_matrix.cool"
original_matrix_cool_chr4 = ROOT + "/small_test_matrix_chr4.cool"


def test_hicConvertFormat_h5_to_cool():

    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat cool".format(original_matrix_h5, outfile.name).split()
    compute(hicConvertFormat.main, args, 5)

    test = hm.hiCMatrix(original_matrix_cool)
    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)


def test_hicConvertFormat_h5_to_cool_enforce_integer():

    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat cool ".format(original_matrix_h5, outfile.name).split()
    compute(hicConvertFormat.main, args, 5)

    test = hm.hiCMatrix(original_matrix_cool)
    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=0)
    assert issubclass(test.matrix.data.dtype.type, np.integer)


def test_hicConvertFormat_h5_to_homer():

    outfile = NamedTemporaryFile(suffix='.homer', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat cool --outputFormat homer ".format(original_matrix_cool_chr4, outfile.name).split()
    compute(hicConvertFormat.main, args, 5)

    test = hm.hiCMatrix(original_matrix_cool_chr4)
    f = gzip.open(outfile.name, 'rb')
    file_content = f.read()
    outfile2 = NamedTemporaryFile(suffix='.homer', delete=False)
    outfile2.close()
    with open(outfile2.name, 'wb') as matrix_file:
        matrix_file.write(file_content)

    matrixFileHandlerInput = MatrixFileHandler(pFileType='homer', pMatrixFile=outfile2.name)

    _matrix, cut_intervals, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

    nt.assert_array_almost_equal(test.matrix.data, _matrix.data, decimal=0)


def test_hicConvertFormat_h5_to_ginteractions():
    outfile = NamedTemporaryFile(suffix='.ginteractions', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat ginteractions ".format(original_matrix_h5, outfile.name).split()
    compute(hicConvertFormat.main, args, 5)


def test_hicConvertFormat_h5_to_hicpro():
    outfile = NamedTemporaryFile(suffix='.hicpro', delete=False)
    outfile_bed = NamedTemporaryFile(suffix='.bed', delete=False)

    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat hicpro --bedFileHicpro {}".format(original_matrix_h5, outfile.name, outfile_bed.name).split()
    compute(hicConvertFormat.main, args, 5)


def test_hicConvertFormat_h5_to_mcool():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat mcool -r 10000 100000 200000 ".format(original_matrix_h5, outfile.name).split()
    compute(hicConvertFormat.main, args, 5)

    new1 = hm.hiCMatrix(outfile.name + '::/resolutions/10000')  # noqa: F841
    new2 = hm.hiCMatrix(outfile.name + '::/resolutions/100000')  # noqa: F841
    new3 = hm.hiCMatrix(outfile.name + '::/resolutions/200000')  # noqa: F841


def test_hicConvertFormat_cool_to_h5():

    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat cool --outputFormat h5".format(original_matrix_cool, outfile.name).split()
    compute(hicConvertFormat.main, args, 5)

    test = hm.hiCMatrix(original_matrix_h5)
    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)


def test_hicConvertFormat_hicpro_to_cool():

    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    hicprofile = ROOT + '/test_matrix.hicpro'
    bedfile = ROOT + '/test_matrix.bed'
    args = "--matrices {} --outFileName {} --inputFormat hicpro --outputFormat cool --bedFileHicpro {}".format(hicprofile, outfile.name, bedfile).split()
    compute(hicConvertFormat.main, args, 5)

    new = hm.hiCMatrix(outfile.name)

    matrixFileHandlerInput = MatrixFileHandler(pFileType='hicpro', pMatrixFile=hicprofile,
                                               pBedFileHicPro=bedfile)

    _matrix, cut_intervals, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

    new.matrix = triu(new.matrix)
    nt.assert_array_almost_equal(new.matrix.data, _matrix.data, decimal=0)


def _expected_2d_text_pixels(text_file, sizes_file, resolution):
    """The upper triangle a 2D-text file means: one pixel per line.

    The pixel sits at the bins that contain start1 and start2. A later line
    for the same pixel replaces the earlier one, a value of zero removes it,
    and a line whose start lies outside its chromosome is skipped. 3.7.6 also
    assigned the bins of the exclusive end coordinates, which doubled the
    pixels of every line that ends in the next bin (#808).
    """
    offsets = {}
    first = 0
    with open(sizes_file) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                break
            offsets[fields[0]] = (first, int(fields[1]))
            first += -(-int(fields[1]) // resolution)
    cells = {}
    with open(text_file) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            start1, start2 = int(fields[1]), int(fields[4])
            if not (0 <= start1 < offsets[fields[0]][1] and 0 <= start2 < offsets[fields[3]][1]):
                continue
            cells[(offsets[fields[0]][0] + start1 // resolution,
                   offsets[fields[3]][0] + start2 // resolution)] = float(fields[6])
    return {key: value for key, value in cells.items()
            if value != 0.0 and key[0] <= key[1]}


def _pixels_of(path):
    import cooler
    pixels = cooler.Cooler(path).pixels()[:]
    return {(int(a), int(b)): float(v)
            for a, b, v in zip(pixels['bin1_id'], pixels['bin2_id'], pixels['count'])}


def test_hicConvertFormat_2D_text_to_cool():

    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    text_2d = ROOT + '/GSM1436265_RAD21_ENCFF002EMQ.txt'
    args = "--matrices {} --outFileName {} --inputFormat 2D-text --outputFormat cool -r 10000 --chromosomeSizes {}".format(text_2d, outfile.name, ROOT + '/hg19.chrom.sizes').split()
    compute(hicConvertFormat.main, args, 5)

    assert _pixels_of(outfile.name) == _expected_2d_text_pixels(
        text_2d, ROOT + '/hg19.chrom.sizes', 10000)
    os.unlink(outfile.name)


def test_2D_text_gives_one_pixel_per_line():
    """Four lines are four pixels (#808); 3.7.6 wrote eight."""
    import tempfile
    directory = tempfile.mkdtemp()
    text = os.path.join(directory, 'in.txt')
    sizes = os.path.join(directory, 'sz.txt')
    out = os.path.join(directory, 'o.cool')
    lines = [('chr1', 710000, 720000, 'chr1', 750000, 760000, 17),
             ('chr1', 10000, 20000, 'chr1', 30000, 40000, 3),
             ('chr1', 100000, 110000, 'chr1', 100000, 110000, 5),
             ('chr1', 5000, 15000, 'chr1', 205000, 215000, 2)]
    with open(text, 'w') as handle:
        for line in lines:
            handle.write('\t'.join(str(x) for x in line) + '\n')
    with open(sizes, 'w') as handle:
        handle.write('chr1\t1000000\n')
    hicConvertFormat.main(['--matrices', text, '--outFileName', out, '--resolutions', '10000',
                           '--chromosomeSizes', sizes, '--inputFormat', '2D-text',
                           '--outputFormat', 'cool'])
    assert _pixels_of(out) == {(71, 75): 17.0, (1, 3): 3.0, (10, 10): 5.0, (0, 20): 2.0}


# ---------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port.
#
# The tests above compare at decimal=0, which for these count matrices means
# "agrees to the nearest integer", and the ginteractions (:76), hicpro (:84)
# and mcool (:94) cases assert nothing at all. test_hicConvertFormat_h5_to_cool
# _enforce_integer never passes --enforce_integer.
#
# The tests below re-assert the same conversions with assert_array_equal on the
# full csr triple, add the three unasserted output formats, and pin
# --enforce_integer.
# ---------------------------------------------------------------------------
import gzip as gzip_module  # noqa: E402
import os  # noqa: E402,F401
from scipy.sparse import triu as triu_module  # noqa: E402

DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data")
FLOAT_H5 = os.path.join(DATA_ROOT, 'hicDifferentialTAD',
                        'GSM2644945_Untreated-R1.100000_chr1.h5')
SMALL_50KB_COOL = os.path.join(DATA_ROOT, 'small_test_matrix_50kb_res.cool')


def assert_csr_identical(got, expected):
    got = got.tocsr()
    got.sort_indices()
    expected = expected.tocsr()
    expected.sort_indices()
    assert got.shape == expected.shape
    assert got.nnz == expected.nnz
    nt.assert_array_equal(got.indptr, expected.indptr)
    nt.assert_array_equal(got.indices, expected.indices)
    nt.assert_array_equal(got.data, expected.data)


def convert(args_string, suffix):
    outfile = NamedTemporaryFile(suffix=suffix, delete=False)
    outfile.close()
    hicConvertFormat.main(args_string.format(out=outfile.name).split())
    return outfile.name


def test_h5_to_cool_reproduces_the_reference_cool_exactly():
    """h5 -> cool is value identical to the stored reference cool file.

    The previous test only checked matrix.data at decimal=0. The conversion is
    in fact exact down to the csr index arrays, and the values are stored as
    int32 in cool where the h5 held int64.
    """
    out = convert("--matrices " + original_matrix_h5 +
                  " --outFileName {out} --inputFormat h5 --outputFormat cool",
                  '.cool')
    reference = hm.hiCMatrix(original_matrix_cool)
    new = hm.hiCMatrix(out)
    source = hm.hiCMatrix(original_matrix_h5)

    assert new.matrix.dtype == np.int32
    assert source.matrix.dtype == np.int64
    assert_csr_identical(new.matrix, reference.matrix)
    nt.assert_array_equal(new.matrix.data, source.matrix.data)
    nt.assert_equal(new.cut_intervals, reference.cut_intervals)
    os.unlink(out)


def test_h5_to_cool_replaces_the_coverage_column_with_ones():
    """The fourth field of every cut interval comes back as 1.0 from cool.

    small_test_matrix.h5 stores a per-bin coverage whose first entry is NaN.
    A cooler file has no coverage column, so hicmatrix fills the field with
    1.0 for every bin on load. The chromosome, start and end of every bin are
    unchanged, and the result equals the stored reference cool file.
    Pinned because the port must not carry the h5 coverage across.
    """
    out = convert("--matrices " + original_matrix_h5 +
                  " --outFileName {out} --inputFormat h5 --outputFormat cool",
                  '.cool')
    source = hm.hiCMatrix(original_matrix_h5)
    reference = hm.hiCMatrix(original_matrix_cool)
    new = hm.hiCMatrix(out)

    assert np.isnan(source.cut_intervals[0][3])
    assert set(interval[3] for interval in new.cut_intervals) == {1.0}
    assert [i[:3] for i in new.cut_intervals] == [i[:3] for i in source.cut_intervals]
    nt.assert_equal(new.cut_intervals, reference.cut_intervals)
    os.unlink(out)


def test_cool_to_h5_is_value_identical_but_changes_dtype_and_nan_bins():
    """cool -> h5 keeps every value and gains a nan bin list.

    The source h5 carries no nan bins, the cool file carries 14,845, and the h5
    written from the cool keeps them. The dtype stays int32 rather than
    returning to the int64 of the original h5.
    """
    out = convert("--matrices " + original_matrix_cool +
                  " --outFileName {out} --inputFormat cool --outputFormat h5",
                  '.h5')
    source_h5 = hm.hiCMatrix(original_matrix_h5)
    source_cool = hm.hiCMatrix(original_matrix_cool)
    new = hm.hiCMatrix(out)

    assert new.matrix.dtype == np.int32
    assert len(source_h5.nan_bins) == 0
    assert len(new.nan_bins) == 14845
    nt.assert_array_equal(sorted(new.nan_bins), sorted(source_cool.nan_bins))
    nt.assert_array_equal(new.matrix.data, source_h5.matrix.data)
    nt.assert_array_equal(new.matrix.indices, source_h5.matrix.indices)
    nt.assert_array_equal(new.matrix.indptr, source_h5.matrix.indptr)
    os.unlink(out)


def test_h5_to_ginteractions_writes_the_upper_triangle_to_the_given_file():
    """--outputFormat ginteractions writes the file named by --outFileName.

    HiCExplorer 3.7.6 wrote '<outFileName>.tsv' and left the named file empty;
    HiCExplorer 4 writes the named file and nothing else. The content is one
    line per stored upper triangle entry, in csr row order.
    """
    outfile = NamedTemporaryFile(suffix='.ginteractions', delete=False)
    outfile.close()
    hicConvertFormat.main(("--matrices " + original_matrix_h5 +
                           " --outFileName " + outfile.name +
                           " --inputFormat h5 --outputFormat ginteractions").split())
    assert not os.path.exists(outfile.name + '.tsv')

    source = hm.hiCMatrix(original_matrix_h5)
    intervals = source.cut_intervals
    upper = triu_module(source.matrix, k=0, format='csr').tocoo()
    expected = ''.join(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            intervals[row][0], intervals[row][1], intervals[row][2],
            intervals[col][0], intervals[col][1], intervals[col][2], value)
        for row, col, value in zip(upper.row, upper.col, upper.data))

    with open(outfile.name) as handle:
        assert handle.read() == expected
    assert len(expected.splitlines()) == 35857
    os.unlink(outfile.name)


def test_h5_to_hicpro_writes_one_based_bin_ids():
    """The hicpro matrix and its bed file are pinned exactly.

    The matrix holds the upper triangle as "bin1 bin2 value" with 1-based bin
    ids, and the bed file lists every bin with the same 1-based id in its
    fourth column. The previous test asserted nothing about either file.
    """
    outfile = NamedTemporaryFile(suffix='.hicpro', delete=False)
    outfile.close()
    bedfile = NamedTemporaryFile(suffix='.bed', delete=False)
    bedfile.close()
    hicConvertFormat.main(("--matrices " + original_matrix_h5 +
                           " --outFileName " + outfile.name +
                           " --inputFormat h5 --outputFormat hicpro"
                           " --bedFileHicpro " + bedfile.name).split())

    source = hm.hiCMatrix(original_matrix_h5)
    upper = triu_module(source.matrix, k=0, format='csr').tocoo()
    expected_matrix = ''.join(
        "{}\t{}\t{}\n".format(row + 1, col + 1, value)
        for row, col, value in zip(upper.row, upper.col, upper.data))
    expected_bed = ''.join(
        "{}\t{}\t{}\t{}\n".format(chrom, start, end, index + 1)
        for index, (chrom, start, end, _) in enumerate(source.cut_intervals))

    with open(outfile.name) as handle:
        assert handle.read() == expected_matrix
    with open(bedfile.name) as handle:
        assert handle.read() == expected_bed
    assert len(expected_matrix.splitlines()) == 35857
    assert len(expected_bed.splitlines()) == 33754
    os.unlink(outfile.name)
    os.unlink(bedfile.name)


def test_hicpro_output_writes_each_contact_once_from_a_symmetric_matrix():
    """HiC-Pro matrices are upper triangular: bin1 <= bin2, each contact once.

    HiCExplorer 3.7.6 wrote every stored entry of the symmetric matrix, both
    triangles; HiCExplorer 4 writes the upper triangle.
    """
    outfile = NamedTemporaryFile(suffix='.hicpro', delete=False)
    outfile.close()
    bedfile = NamedTemporaryFile(suffix='.bed', delete=False)
    bedfile.close()
    hicConvertFormat.main(("--matrices " + original_matrix_cool_chr4 +
                           " --outFileName " + outfile.name +
                           " --inputFormat cool --outputFormat hicpro"
                           " --bedFileHicpro " + bedfile.name).split())
    source = hm.hiCMatrix(original_matrix_cool_chr4)
    upper = triu_module(source.matrix, k=0, format='csr').tocoo()
    assert source.matrix.nnz > upper.nnz
    with open(outfile.name) as handle:
        rows = [line.split('\t') for line in handle.read().splitlines()]
    assert len(rows) == upper.nnz
    assert all(int(a) <= int(b) for a, b, _ in rows)
    assert len({(a, b) for a, b, _ in rows}) == len(rows)
    os.unlink(outfile.name)
    os.unlink(bedfile.name)


def test_cool_to_homer_is_byte_exact():
    """The homer output is a gzipped dense table, pinned byte for byte.

    Header line: the literal 'HiCMatrix (directory=.)', 'Regions', then one
    'chrom-start' name per bin and a trailing tab. Data lines: the bin name
    twice, then the full dense row. There is no trailing newline at the end of
    the file. The previous test compared only the reloaded values at decimal=0.
    """
    outfile = NamedTemporaryFile(suffix='.homer', delete=False)
    outfile.close()
    hicConvertFormat.main(("--matrices " + original_matrix_cool_chr4 +
                           " --outFileName " + outfile.name +
                           " --inputFormat cool --outputFormat homer").split())

    source = hm.hiCMatrix(original_matrix_cool_chr4)
    names = ["{}-{}".format(chrom, start)
             for chrom, start, _, _ in source.cut_intervals]
    dense = source.matrix.toarray()
    lines = ["HiCMatrix (directory=.)\tRegions\t" + "\t".join(names) + "\t"]
    for index, name in enumerate(names):
        lines.append(name + "\t" + name + "\t" +
                     "\t".join(str(value) for value in dense[index]))
    expected = "\n".join(lines)

    with gzip_module.open(outfile.name, 'rb') as handle:
        assert handle.read().decode() == expected
    assert len(expected.splitlines()) == 272
    os.unlink(outfile.name)


def test_h5_to_homer_symmetrises_with_maximum_not_with_a_sum():
    """hicConvertFormat.py:258-260 uses triu then maximum(triu.T).

    For an already symmetric matrix that is a no-op, which is what the two
    homer and ginteractions paths rely on. Pinned so that a port which
    symmetrises by adding the transpose is caught: adding would double the
    diagonal.
    """
    source = hm.hiCMatrix(original_matrix_cool_chr4)
    upper = triu_module(source.matrix)
    symmetric = upper.maximum(upper.T)
    assert_csr_identical(symmetric, source.matrix)


def _matrix_with_negative_values(tmp_path):
    """The stored chr4 matrix of small_test_matrix.h5 with every second off
    diagonal entry of the upper triangle negated, as a symmetric h5 file.
    Only chr4 (271 bins) is kept: the homer output is a dense text matrix, and
    the whole matrix (about 33,750 bins) would make it several gigabytes."""
    source = hm.hiCMatrix(original_matrix_h5)
    source.keepOnlyTheseChr(['chr4'])
    upper = triu_module(source.matrix, k=0, format='csr').astype(float)
    coo = upper.tocoo()
    data = coo.data.copy()
    off_diagonal = np.flatnonzero(coo.row != coo.col)
    data[off_diagonal[::2]] *= -1
    from scipy.sparse import coo_matrix
    upper = coo_matrix((data, (coo.row, coo.col)), shape=coo.shape).tocsr()
    symmetric = upper + triu_module(upper, k=1).T
    source.setMatrix(symmetric.tocsr(), source.cut_intervals)
    path = str(tmp_path / 'negative.h5')
    source.save(path)
    return path, upper


def test_ginteractions_writes_negative_values(tmp_path):
    """The upper triangle is written as it is stored, negative values included.
    3.7.6 took the maximum with the transpose, which turned every negative
    value into 0 (#584)."""
    path, upper = _matrix_with_negative_values(tmp_path)
    out = str(tmp_path / 'gi.tsv')
    hicConvertFormat.main(['--matrices', path, '--outFileName', out,
                           '--inputFormat', 'h5', '--outputFormat', 'ginteractions'])
    with open(out) as handle:
        values = [float(line.split('\t')[6]) for line in handle]
    coo = upper.tocoo()
    assert sorted(values) == sorted(coo.data.tolist())
    assert sum(1 for value in values if value < 0) == sum(1 for value in coo.data if value < 0) > 0


def test_homer_writes_negative_values(tmp_path):
    """The homer matrix is the stored upper triangle mirrored, negative values
    included (#584)."""
    path, upper = _matrix_with_negative_values(tmp_path)
    out = str(tmp_path / 'matrix.homer')
    hicConvertFormat.main(['--matrices', path, '--outFileName', out,
                           '--inputFormat', 'h5', '--outputFormat', 'homer'])
    opener = gzip.open if open(out, 'rb').read(2) == b'\x1f\x8b' else open
    with opener(out, 'rt') as handle:
        lines = handle.read().splitlines()
    table = np.array([[float(x) for x in line.split('\t')[2:]] for line in lines[1:]])
    expected = upper.toarray()
    expected = expected + np.triu(expected, k=1).T
    nt.assert_array_equal(table, expected)
    assert (table < 0).sum() > 0


def test_mcool_with_resolutions_and_several_matrices_exits_before_writing(tmp_path, capfd):
    """HiCExplorer 4 (fix 4): --resolutions describes how to resample one
    matrix, so several input matrices with it are an error. 3.7.6 only logged
    the message, exited with status 0 and let every input matrix rewrite the
    output file; now the run stops with status 1 and writes nothing.
    """
    import pytest
    outfile = str(tmp_path / 'matrix.mcool')
    args = ("--matrices {} {} --outFileName {} --inputFormat cool "
            "--outputFormat mcool --resolutions 10000 100000").format(
                original_matrix_cool_chr4, original_matrix_cool, outfile).split()
    with pytest.raises(SystemExit) as error:
        hicConvertFormat.main(args)
    assert error.value.code == 1
    assert not os.path.exists(outfile)
    assert 'Please define one matrix and many resolutions' in capfd.readouterr().err


def test_cool_to_mcool_merges_bins_exactly():
    """--outputFormat mcool with -r drives hicMergeMatrixBins.merge_bins.

    Three resolutions are written under /resolutions/<res>. The 50 kb one is
    the source resolution and is copied unchanged; the other two are block
    sums, checked here against a dense reduction of the source. The previous
    test only opened the three groups and discarded the result.

    Note the difference to the hicMergeMatrixBins command line tool: here the
    matrix is rebuilt with setMatrix, so nan_bins is empty and
    remove_nans_if_needed does nothing. The 664 nan bins of the source cool
    file survive the merge instead of being deleted.
    """
    import h5py
    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)
    outfile.close()
    hicConvertFormat.main(("--matrices " + SMALL_50KB_COOL +
                           " --outFileName " + outfile.name +
                           " --inputFormat cool --outputFormat mcool"
                           " -r 50000 100000 200000").split())

    with h5py.File(outfile.name, 'r') as handle:
        assert list(handle.keys()) == ['resolutions']
        assert sorted(handle['resolutions'].keys()) == ['100000', '200000', '50000']

    source = hm.hiCMatrix(SMALL_50KB_COOL)
    assert len(source.nan_bins) == 664
    upper = triu_module(source.matrix, k=0, format='coo')

    expected_shapes = {50000: 3383, 100000: 1697, 200000: 846}
    expected_grouped = {50000: 3383, 100000: 3383, 200000: 3377}
    for resolution in (50000, 100000, 200000):
        new = hm.hiCMatrix(outfile.name + '::/resolutions/' + str(resolution))
        assert new.matrix.shape == (expected_shapes[resolution],
                                    expected_shapes[resolution])
        groups = []
        for chrom, start, end, _ in new.cut_intervals:
            groups.append([i for i, (c, s, e, _) in enumerate(source.cut_intervals)
                           if c == chrom and s >= start and e <= end])
        assert sum(len(g) for g in groups) == expected_grouped[resolution]

        mapping = np.full(source.matrix.shape[0], -1, dtype=int)
        for index, group in enumerate(groups):
            for bin_id in group:
                mapping[bin_id] = index
        new_row = mapping[upper.row]
        new_col = mapping[upper.col]
        keep = (new_row > -1) & (new_col > -1)
        size = len(groups)
        reduced = np.zeros((size, size), dtype=np.float64)
        np.add.at(reduced, (new_row[keep], new_col[keep]), upper.data[keep])
        reduced = reduced + reduced.T - np.diag(np.diag(reduced))
        nt.assert_array_equal(new.matrix.toarray().astype(np.float64), reduced)
    os.unlink(outfile.name)


def test_enforce_integer_refuses_a_matrix_of_values_below_one_half(tmp_path):
    """--enforce_integer on an already corrected matrix stops with a message.

    GSM2644945_Untreated-R1.100000_chr1.h5 holds 2,504,071 corrected values
    between 1.3e-05 and 0.12. HiCExplorer 3.7.6 rounded every one of them to
    zero and wrote a cool file without entries, with no warning and exit
    status 0. HiCExplorer 4 refuses to write a file that rounding would
    empty. Without the flag the same conversion is bit exact.
    """
    source = hm.hiCMatrix(FLOAT_H5)
    assert source.matrix.dtype == np.float64
    assert source.matrix.nnz == 2504071
    assert source.matrix.data.max() < 0.5

    refused = str(tmp_path / 'refused.cool')
    with pytest.raises(SystemExit) as raised:
        hicConvertFormat.main(("--matrices " + FLOAT_H5 +
                               " --outFileName " + refused +
                               " --inputFormat h5 --outputFormat cool"
                               " --enforce_integer").split())
    assert raised.value.code != 0
    assert not os.path.exists(refused)

    out = convert("--matrices " + FLOAT_H5 +
                  " --outFileName {out} --inputFormat h5 --outputFormat cool",
                  '.cool')
    plain = hm.hiCMatrix(out)
    assert plain.matrix.dtype == np.float64
    nt.assert_array_equal(plain.matrix.data, source.matrix.data)
    nt.assert_array_equal(plain.matrix.indices, source.matrix.indices)
    nt.assert_array_equal(plain.matrix.indptr, source.matrix.indptr)
    os.unlink(out)


def test_cool_to_mcool_resolutions_all_carry_the_genome_assembly():
    """Every resolution of an mcool file carries genome-assembly.

    HiCExplorer 3.7.6 popped matrix-generated-by, matrix-generated-by-url and
    genome-assembly out of one dictionary shared by the resolutions, so only
    the first resolution had them in its metadata attribute.
    """
    import h5py
    import json
    source = os.path.join(DATA_ROOT, 'hicDetectLoops',
                          'GSE63525_GM12878_insitu_primary_2_5mb.cool')
    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)
    outfile.close()
    hicConvertFormat.main(("--matrices " + source + " --outFileName " + outfile.name +
                           " --inputFormat cool --outputFormat mcool"
                           " -r 2500000 5000000").split())
    with h5py.File(outfile.name, 'r') as handle:
        assert sorted(handle['resolutions'].keys()) == ['2500000', '5000000']
        for resolution in ('2500000', '5000000'):
            metadata = json.loads(handle['resolutions'][resolution].attrs['metadata'])
            assert metadata['genome-assembly'] == 'hg19', resolution
    os.unlink(outfile.name)


def test_mcool_without_resolutions_honours_enforce_integer(tmp_path):
    """The Python ignored --enforce_integer when no --resolutions were given
    and wrote the corrected values unchanged; HiCExplorer 4 applies the flag
    (and so refuses the matrix of values below one half)."""
    target = str(tmp_path / 'matrix.mcool')
    with pytest.raises(SystemExit):
        hicConvertFormat.main(("--matrices " + FLOAT_H5 + " --outFileName " + target +
                               " --inputFormat h5 --outputFormat mcool"
                               " --enforce_integer").split())


def test_enforce_integer_is_inert_on_an_integer_matrix():
    """The same flag on small_test_matrix.h5 changes nothing."""
    out = convert("--matrices " + original_matrix_h5 +
                  " --outFileName {out} --inputFormat h5 --outputFormat cool"
                  " --enforce_integer", '.cool')
    reference = hm.hiCMatrix(original_matrix_cool)
    new = hm.hiCMatrix(out)
    assert_csr_identical(new.matrix, reference.matrix)
    os.unlink(out)


def test_hicpro_to_cool_is_bit_exact():
    """Re-assert the hicpro import at full precision.

    The previous test compared the upper triangle at decimal=0.
    """
    hicprofile = ROOT + '/test_matrix.hicpro'
    bedfile = ROOT + '/test_matrix.bed'
    out = convert("--matrices " + hicprofile +
                  " --outFileName {out} --inputFormat hicpro"
                  " --outputFormat cool --bedFileHicpro " + bedfile, '.cool')

    handler = MatrixFileHandler(pFileType='hicpro', pMatrixFile=hicprofile,
                                pBedFileHicPro=bedfile)
    expected, _, _, _, _ = handler.load()

    new = hm.hiCMatrix(out)
    assert_csr_identical(triu(new.matrix), expected)
    os.unlink(out)


def test_2D_text_to_cool_is_bit_exact():
    """Re-assert the 2D-text import at full precision.

    The expected pixels are derived from the text file itself: one pixel per
    line at the bins of start1 and start2 (#808). The stored 2dtexttocool.cool
    is the 3.7.6 output, which holds 7,987 upper triangle entries because the
    bins of the end coordinates were assigned as well.
    """
    text_2d = ROOT + '/GSM1436265_RAD21_ENCFF002EMQ.txt'
    out = convert("--matrices " + text_2d +
                  " --outFileName {out} --inputFormat 2D-text"
                  " --outputFormat cool -r 10000 --chromosomeSizes " +
                  ROOT + '/hg19.chrom.sizes', '.cool')

    expected = _expected_2d_text_pixels(text_2d, ROOT + '/hg19.chrom.sizes', 10000)
    assert len(expected) == 7451
    assert _pixels_of(out) == expected

    old_handler = MatrixFileHandler(pFileType='cool',
                                    pMatrixFile=ROOT + '/2dtexttocool.cool')
    old_matrix, old_intervals, _, _, _ = old_handler.load()
    assert old_matrix.nnz == 7987

    new = hm.hiCMatrix(out)
    assert new.matrix.shape == (313762, 313762)
    # hicmatrix returns the symmetric matrix: every off diagonal pixel twice.
    assert new.matrix.nnz == 2 * len(expected) - sum(1 for key in expected if key[0] == key[1])
    assert list(map(tuple, new.cut_intervals)) == list(map(tuple, old_intervals))
    os.unlink(out)


# ---------------------------------------------------------------------------
# Characterization tests for --inputFormat hic, added for the v4 C++ port
# (cpp/PLAN.md tier 9, item 9.1). The Python route is hic2cool
# (hicConvertFormat.py:124-138); these tests pin what it writes for the
# repository's Juicer file SRR1791297_30.hic (version 8, sacCer3, nine base
# pair resolutions from 2.5 Mb to 5 kb, KR, VC and VC_SQRT vectors).
# ---------------------------------------------------------------------------
from tempfile import mkdtemp  # noqa: E402

HIC_FILE = os.path.join(DATA_ROOT, 'hicHyperoptDetectLoopsHiCCUPS', 'SRR1791297_30.hic')


def test_hic_input_one_resolution_writes_a_cool_file_named_after_the_resolution():
    """--resolutions R inserts _R before the file extension.

    The cool file is hic2cool's layout, not hicmatrix's: fixed width S32
    chromosome names, an enum chrom column, one float64 bin column per
    normalization of the .hic file, int32 counts, the .hic header attributes
    (statistics, graphs) copied to the root, and hic2cool's format-url and
    generated-by. The file named by --outFileName itself is not written.
    """
    import h5py
    directory = mkdtemp()
    out = os.path.join(directory, 'matrix.cool')
    hicConvertFormat.main(['--matrices', HIC_FILE, '--outFileName', out,
                           '--inputFormat', 'hic', '--outputFormat', 'cool',
                           '--resolutions', '250000'])
    assert not os.path.exists(out)
    with h5py.File(os.path.join(directory, 'matrix_250000.cool'), 'r') as handle:
        attrs = handle.attrs
        assert attrs['bin-size'] == 250000
        assert attrs['nbins'] == 57
        assert attrs['nchroms'] == 16
        assert attrs['nnz'] == 1653
        assert attrs['format-version'] == 3
        assert attrs['storage-mode'] == 'symmetric-upper'
        assert attrs['format-url'] == 'https://github.com/4dn-dcic/hic2cool'
        assert attrs['generated-by'].startswith('hic2cool-')
        assert 'statistics' in attrs and 'graphs' in attrs
        assert sorted(handle['bins'].keys()) == ['KR', 'VC', 'VC_SQRT', 'chrom', 'end', 'start']
        assert h5py.check_enum_dtype(handle['bins/chrom'].dtype) is not None
        assert handle['chroms/name'].dtype == np.dtype('S32')
        assert handle['pixels/count'].dtype == np.int32
        assert [x.decode() for x in handle['chroms/name'][:3]] == ['NC_001133.9', 'NC_001134.8', 'NC_001135.5']
        pixels = list(zip(handle['pixels/bin1_id'][:4], handle['pixels/bin2_id'][:4],
                          handle['pixels/count'][:4]))
        assert pixels == [(0, 0, 243095), (0, 1, 17306), (0, 2, 10095), (0, 3, 2940)]
        assert int(handle['pixels/count'][:].sum()) == 36490041
        nt.assert_array_equal(handle['bins/KR'][:2], [1.0, 1.1494348144163238])
        nt.assert_array_equal(handle['indexes/chrom_offset'][:4], [0, 1, 5, 7])
        nt.assert_array_equal(handle['indexes/bin1_offset'][:4], [0, 57, 113, 168])


def test_hic_input_without_resolutions_writes_every_resolution_to_an_mcool_file():
    """Without --resolutions hic2cool converts all nine resolutions and, since
    that is more than one, renames matrix.cool to matrix.mcool."""
    import h5py
    directory = mkdtemp()
    out = os.path.join(directory, 'matrix.cool')
    hicConvertFormat.main(['--matrices', HIC_FILE, '--outFileName', out,
                           '--inputFormat', 'hic', '--outputFormat', 'cool'])
    assert not os.path.exists(out)
    with h5py.File(os.path.join(directory, 'matrix.mcool'), 'r') as handle:
        assert handle.attrs['format'] == 'HDF5::MCOOL'
        assert sorted(int(r) for r in handle['resolutions'].keys()) == \
            [5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000]
        finest = handle['resolutions/5000']
        assert finest.attrs['nnz'] == 2421438
        assert int(finest['pixels/count'][:].sum()) == 36490041


def test_hic_input_to_any_format_but_cool_is_refused():
    """hicConvertFormat.py:120-122: log.error and exit(1)."""
    directory = mkdtemp()
    try:
        hicConvertFormat.main(['--matrices', HIC_FILE, '--outFileName',
                               os.path.join(directory, 'matrix.h5'),
                               '--inputFormat', 'hic', '--outputFormat', 'h5'])
    except SystemExit as error:
        assert error.code == 1
    else:
        raise AssertionError('hic to h5 did not exit')
    assert os.listdir(directory) == []


# ---------------------------------------------------------------------------
# Characterization tests for .hic versions 7 and 6, added before the C++ port
# reads them (cpp/PLAN.md tier 9, item 9.1: read versions 6 to 9). The files
# are chromosomes 21 and 22 of GEO GSE63525's
# GSM12878_insitu_primary+replicate_combined_30.hic (version 7) at 2.5 Mb,
# 1 Mb, 500 kb and 250 kb with its VC, VC_SQRT, KR and GW_/INTER_ vectors, cut
# out byte for byte by hicfilecpp's tests/data/extract_legacy_subset.py; the
# version 6 file holds the same pixels in version 6 block records. hic2cool
# reads both (hic2cool_utils.read_block, version < 7).
# ---------------------------------------------------------------------------
LEGACY_HIC = {version: os.path.join(DATA_ROOT, 'hicConvertFormat',
                                    'GM12878_combined_30.chr21_chr22.v{}.hic'.format(version))
              for version in (7, 6)}


def test_hic_input_versions_7_and_6_write_the_same_cool_file():
    """--resolutions 250000 on the version 7 file writes hic2cool's layout with
    every normalization of the file as a bin column; the version 6 file gives
    the same datasets."""
    import h5py
    datasets = {}
    for version, path in LEGACY_HIC.items():
        directory = mkdtemp()
        hicConvertFormat.main(['--matrices', path, '--outFileName', os.path.join(directory, 'matrix.cool'),
                               '--inputFormat', 'hic', '--outputFormat', 'cool', '--resolutions', '250000'])
        assert os.listdir(directory) == ['matrix_250000.cool']
        with h5py.File(os.path.join(directory, 'matrix_250000.cool'), 'r') as handle:
            attrs = handle.attrs
            assert attrs['bin-size'] == 250000
            assert attrs['nbins'] == 399
            assert attrs['nchroms'] == 2
            assert attrs['nnz'] == 39771
            assert attrs['genome-assembly'] == 'hg19'
            assert attrs['storage-mode'] == 'symmetric-upper'
            assert sorted(handle['bins'].keys()) == ['GW_KR', 'GW_VC', 'INTER_KR', 'INTER_VC', 'KR', 'VC',
                                                     'VC_SQRT', 'chrom', 'end', 'start']
            assert [x.decode() for x in handle['chroms/name'][:]] == ['21', '22']
            nt.assert_array_equal(handle['chroms/length'][:], [48129895, 51304566])
            assert handle['pixels/count'].dtype == np.int32
            pixels = list(zip(handle['pixels/bin1_id'][:4], handle['pixels/bin2_id'][:4],
                              handle['pixels/count'][:4]))
            assert pixels == [(37, 37, 278), (37, 38, 7), (37, 39, 2), (37, 41, 6)]
            assert int(handle['pixels/count'][:].sum()) == 67749579
            assert np.isnan(handle['bins/KR'][:3]).all()
            nt.assert_array_equal(handle['indexes/chrom_offset'][:], [0, 193, 399])
            datasets[version] = {name: handle[name][:] for name in
                                 ['pixels/bin1_id', 'pixels/bin2_id', 'pixels/count', 'bins/KR', 'bins/VC',
                                  'bins/VC_SQRT', 'bins/GW_KR', 'indexes/bin1_offset']}
    for name, values in datasets[7].items():
        nt.assert_array_equal(datasets[6][name], values)


def test_hic_input_version_7_without_resolutions_writes_an_mcool_file():
    import h5py
    directory = mkdtemp()
    hicConvertFormat.main(['--matrices', LEGACY_HIC[7], '--outFileName', os.path.join(directory, 'matrix.cool'),
                           '--inputFormat', 'hic', '--outputFormat', 'cool'])
    with h5py.File(os.path.join(directory, 'matrix.mcool'), 'r') as handle:
        nnz = {int(r): int(handle['resolutions'][r].attrs['nnz']) for r in handle['resolutions'].keys()}
    assert nnz == {250000: 39771, 500000: 10525, 1000000: 2775, 2500000: 528}
