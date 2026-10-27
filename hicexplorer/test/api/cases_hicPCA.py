"""Equivalence cases for hicexplorer.hicPCA (apieq.py).

correlateEigenvectorWithGeneTrack runs on hiCMatrix objects of the
repository's matrices with the repository's gene tracks (dm3_genes.bed.gz,
gzipped, against small_test_matrix_50kb_res; hicPCA/mm9_genes_sorted.bed12
against hicPCA/mm9_reduced_chr1.cool) and eigenvectors drawn from a seeded
generator, as an array and as the list of complex lists main passes.
correlateEigenvectorWithHistonMarkTrack runs on the only bigWig track of the
corpus, bigwig_chrx_2e6_5e6.bw (chrX 2 to 5 Mb), opened with pyBigWig in the
driver, with the 5 kb bins of hicPCA/obsexp_norm.h5 around it; the driver
returns the eigenvectors the function changed in place.

Deviation from 3.7.6 (cpp/python/src/api/hicPCA.cpp): eigenvectors of other
dtypes than float32, float64, complex64 and complex128 raise TypeError naming
HiCExplorer 4.
"""
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicPCA.'
SMALL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
DM3_GENES = os.path.join(TEST_DATA, 'dm3_genes.bed.gz')
MM9 = os.path.join(TEST_DATA, 'hicPCA', 'mm9_reduced_chr1.cool')
MM9_GENES = os.path.join(TEST_DATA, 'hicPCA', 'mm9_genes_sorted.bed12')
OBS_EXP = os.path.join(TEST_DATA, 'hicPCA', 'obsexp_norm.h5')
TRACK = os.path.join(TEST_DATA, 'bigwig_chrx_2e6_5e6.bw')


def load(path):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(path)


def eigenvectors(n, k, seed, kind):
    values = np.random.RandomState(seed).standard_normal((n, k))
    if kind == 'complex':
        return values.astype(np.complex128)
    if kind == 'complex_list':
        return values.astype(np.complex128).tolist()
    if kind == 'list':
        return values.tolist()
    return values.astype(kind)


def gene_track(target, matrix, track, k=2, seed=0, kind='float64', keep=None):
    hic = load(matrix)
    if keep is not None:
        hic.keepOnlyTheseChr(keep)
    vectors = eigenvectors(len(hic.cut_intervals), k, seed, kind)
    return target(hic, vectors, track)


def gene_track_constant(target):
    """An eigenvector that is constant on a chromosome: pearsonr warns and
    returns NaN, which does not flip."""
    hic = load(SMALL)
    vectors = eigenvectors(len(hic.cut_intervals), 2, 4, 'float64')
    first, last = hic.getChrBinRange('chr2L')
    vectors[first:last, 0] = 0.25
    return target(hic, vectors, DM3_GENES)


def gene_track_file(target, content, k=1):
    """A small gene track in the case directory."""
    hic = load(SMALL)
    vectors = eigenvectors(len(hic.cut_intervals), k, 7, 'float64')
    return target(hic, vectors, content)


def histone(target, mark_type='active', seed=0, kind='float64', chromosome='chrX', first=300, last=1100,
            k=2, as_list=False):
    import pyBigWig
    hic = load(OBS_EXP)
    intervals = hic.cut_intervals[first:last]
    start = [x[1] for x in intervals]
    end = [x[2] for x in intervals]
    vectors = eigenvectors(k, len(intervals), seed, kind)
    if as_list:
        vectors = [row.copy() for row in vectors]
    track = pyBigWig.open(TRACK)
    try:
        returned = target(vectors, track, chromosome, start, end, TRACK, mark_type)
    finally:
        track.close()
    return returned, vectors


GENES_SMALL = b"# genes\nchr2L\t100000\t120000\tg1\t0\t+\nchr2L\t5000000\t5100000\tg2\t0\t-\n" \
              b"chrX\t100\t200\tg3\t0\t+\nchrX\t30000000\t30000100\tg4\t0\t+\nchr9\t1\t2\tg5\t0\t+\n" \
              b"chr3R\t27000000\t29000000\tg6\t1\t-\nchr3R\t27900000\t27900010\tg7\t1\t-\n"
GENES_MIXED = b"chr2L\t100000\t120000\tg1\t0\t+\nchrX\t100\t200\tg3\n"


CASES = [
    Case('gene/dm3_float64', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track', args=(SMALL, DM3_GENES)),
    Case('gene/dm3_three_vectors', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, DM3_GENES, 3, 1)),
    Case('gene/dm3_complex', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, DM3_GENES, 2, 2, 'complex')),
    Case('gene/dm3_complex_list', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, DM3_GENES, 2, 3, 'complex_list')),
    Case('gene/dm3_list', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, DM3_GENES, 1, 5, 'list')),
    Case('gene/dm3_float32', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, DM3_GENES, 2, 6, 'float32')),
    Case('gene/dm3_some_chromosomes', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, DM3_GENES, 2, 8, 'float64', ['chr2L', 'chrX'])),
    Case('gene/mm9_bed12', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track', args=(MM9, MM9_GENES, 2, 9)),
    Case('gene/constant_vector', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track_constant'),
    Case('gene/small_track', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track_file',
         args=(File('genes.bed'),), files={'genes.bed': GENES_SMALL}),
    Case('gene/mixed_track', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track_file',
         args=(File('genes.bed'),), files={'genes.bed': GENES_MIXED}),
    Case('gene/missing_track', M + 'correlateEigenvectorWithGeneTrack', driver='gene_track',
         args=(SMALL, os.path.join(TEST_DATA, 'no_such_genes.bed'))),

    Case('histone/active', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone'),
    Case('histone/inactive', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone', args=('inactive',)),
    Case('histone/seeds', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone', args=('active', 3)),
    Case('histone/inactive_seed', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('inactive', 5)),
    Case('histone/complex', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('active', 1, 'complex')),
    Case('histone/complex_inactive', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('inactive', 2, 'complex')),
    Case('histone/float32', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('inactive', 4, 'float32')),
    Case('histone/list_of_vectors', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('active', 6, 'float64', 'chrX', 300, 1100, 3, True)),
    Case('histone/no_signal', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('active', 0, 'float64', 'chrX', 0, 300)),
    Case('histone/other_chromosome', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone',
         args=('active', 0, 'float64', 'chr2L')),
    Case('histone/bad_type', M + 'correlateEigenvectorWithHistonMarkTrack', driver='histone', args=('neutral',)),
]

