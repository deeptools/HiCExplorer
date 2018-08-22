from hicexplorer import hicPCA

from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
import numpy as np
import pyBigWig
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

import logging
log = logging.getLogger(__name__)


def are_files_equal(file1, file2):
    equal = True
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                # handle the case of flipped values
                split_x = x.split('\t')
                split_y = y.split('\t')
                if split_x[0] == split_y[0] and split_x[1] == split_y[1] and split_x[2] == split_y[2]:
                    # to ignore rounding errors after 2th digit
                    if 0 <= abs(abs(float(split_x[3].strip())) - abs(float(split_y[3].strip()))) <= 0.01:
                        continue
                    else:
                        log.debug('wtf: {}'.format(abs(abs(float(split_x[3].strip())) - abs(float(split_y[3].strip())))))

                        log.debug('split_x {} split_y {}'.format(split_x, split_y))
                equal = False
                break
    return equal


def are_files_equal_bigwig(pFile1, pFile2, pChromosomeList):

    bw_file1 = pyBigWig.open(pFile1)
    bw_file2 = pyBigWig.open(pFile2)

    for chrom in pChromosomeList:
        try:
            bins_list_file1 = bw_file1.intervals(chrom)
        except Exception:
            log.debug("Chrom not found: {}", chrom)
        try:
            bins_list_file2 = bw_file2.intervals(chrom)
        except Exception:
            log.debug("Chrom not found: {}", chrom)
        # sometimes the values are + / - flipped

        if bins_list_file1 is not None and bins_list_file1[0][2] != bins_list_file2[0][2]:
            bins_list_file1 = np.array(bins_list_file1)
            bins_list_file2 = np.array(bins_list_file2)
            bins_list_file1[:][2] *= -1
        if bins_list_file1 is None and bins_list_file2 is None:
            return True
        nt.assert_array_almost_equal(np.absolute(bins_list_file1), np.absolute(bins_list_file2), decimal=2)
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

    chrom_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr2RHet', 'chr3RHet', 'chr2LHet', 'chr4',
                  'chrU', 'chrX', 'chrXHet', 'chr3LHet']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1.bw", pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2.bw", pca2.name, chrom_list)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bedgraph_gene_density():
    pca1 = NamedTemporaryFile(suffix='.bedgraph', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bedgraph', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    chromosomes = 'chrX chrXHet'
    args = "--matrix {} --outputFileName {} {} -f bedgraph -noe 2 --geneTrack {} --chromosomes {}".format(matrix, pca1.name, pca2.name, gene_track, chromosomes).split()
    hicPCA.main(args)

    assert are_files_equal(ROOT + "hicPCA/pca1_gene_track.bedgraph", pca1.name)
    assert are_files_equal(ROOT + "hicPCA/pca2_gene_track.bedgraph", pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bigwig_gene_density():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    chromosomes = 'chrX chrXHet'
    args = "--matrix {} --outputFileName {} {} -f bigwig -noe 2 --geneTrack {} --chromosomes {}".format(matrix, pca1.name, pca2.name, gene_track, chromosomes).split()
    hicPCA.main(args)

    chrom_list = ['chrX', 'chrXHet']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1_gene_track.bw", pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2_gene_track.bw", pca2.name, chrom_list)

    os.unlink(pca1.name)
    os.unlink(pca2.name)
