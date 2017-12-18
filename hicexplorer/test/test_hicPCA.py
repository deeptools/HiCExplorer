from hicexplorer import hicPCA

from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
import pyBigWig
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

import logging
log = logging.getLogger(__name__)


def are_files_equal(file1, file2):
    equal = True
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                equal = False
                break
    return equal


def are_files_equal_bigwig(file1, file2):
    chrom_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr2RHet', 'chr3RHet', 'chr2LHet', 'chr4', 'chrYHet',
                  'chrU', 'chrX', 'chrXHet', 'chrUextra', 'chrM', 'chr3LHet']
    bw_file1 = pyBigWig.open(file1)
    bw_file2 = pyBigWig.open(file2)

    for chrom in chrom_list:
        try:
            bins_list_file1 = bw_file1.intervals(chrom)
        except Exception:
            log.debug("Chrom not found: {}", chrom)
        try:
            bins_list_file2 = bw_file2.intervals(chrom)
        except Exception:
            log.debug("Chrom not found: {}", chrom)

        nt.assert_equal(bins_list_file1, bins_list_file2)

    return True


def test_pca_bedgraph():
    pca1 = NamedTemporaryFile(suffix='.bedgraph', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bedgraph', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    args = "--matrix {} --outputFileName {} {} -f bedgraph -noe 2".format(matrix, pca1.name, pca2.name).split()
    hicPCA.main(args)

    assert are_files_equal(ROOT + "hicPCA/pca1.bedgraph", pca1.name)
    assert are_files_equal(ROOT + "hicPCA/pca2.bedgraph", pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bigwig():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    args = "--matrix {} --outputFileName {} {} -f bigwig -noe 2".format(matrix, pca1.name, pca2.name).split()
    hicPCA.main(args)

    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1.bw", pca1.name)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2.bw", pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)
