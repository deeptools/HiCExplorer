from shutil import Error
from hicexplorer import chicExportData


from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import tarfile
import pyBigWig
import numpy.testing as nt
import numpy as np
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
# mpl.use('agg')
import logging
log = logging.getLogger(__name__)
ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")

tolerance = 50
DELTA_DECIMAL = 0


def are_files_equal(file1, file2, delta=1, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            # if x.startswith('File'):
            #     continue
            if i < skip:
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
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
        nt.assert_array_almost_equal(np.absolute(bins_list_file1),
                                     np.absolute(bins_list_file2),
                                     decimal=DELTA_DECIMAL)
    return True


def test_interactions_all():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {}".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', outfile.name, 'all').split()
    chicExportData.main(args)

    file_obj_new = tarfile.open(outfile.name, "r")

    namelist = file_obj_new.getnames()
    assert len(namelist) == 6

    output_folder_new = mkdtemp(prefix="output_")
    output_folder_test_data = mkdtemp(prefix="output_")

    file_obj_new.extractall(output_folder_new)

    file_obj_test_data = tarfile.open(ROOT + 'chicExportData/interaction_all_txt.tar.gz', "r")
    file_obj_test_data.extractall(output_folder_test_data)

    files_new = os.listdir(output_folder_new)
    files_test_data = os.listdir(output_folder_test_data)

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, output_folder_test_data + '/' + file_test_data, skip=1)


# must fail
@pytest.mark.xfail(reason='Test case should fail because of wrong input.')
def test_interaction_bigwig_fail():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {} -oft {}".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', outfile.name, 'all', 'bigwig').split()
    chicExportData.main(args)


@pytest.mark.xfail(reason='Test case should fail because of wrong input.')
def test_interaction_bigwig_fail_range():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {} -oft {} --range {} {}".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', outfile.name, 'all', 'bigwig', 20000, 200000).split()
    chicExportData.main(args)


@pytest.mark.xfail(reason='Test case should fail because of wrong input.')
def test_interaction_bigwig_fail_range_background():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {} -oft {} --range {} {} --backgroundModelFile {}".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', outfile.name, 'all', 'bigwig', 20000, 200000, 'background.txt').split()
    chicExportData.main(args)


@pytest.mark.xfail(reason='Test case should fail because of wrong input.')
def test_significant_one_gene_bigwig_fail():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {} -omn {} -oft {}".format(
        ROOT + 'chicSignificantInteractions/significantInteractions_dual.hdf5', outfile.name, 'geneName', 'Eya1', 'bigwig').split()
    chicExportData.main(args)


@pytest.mark.xfail(reason='Test case should fail because of wrong input.')
def test_target_one_gene_bigwig_fail():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {} -omn {} -oft {}".format(
        ROOT + 'chicSignificantInteractions/targetFile_dual.hdf5', outfile.name, 'geneName', 'Eya1', 'bigwig').split()
    chicExportData.main(args)


def test_interaction_bigwig_all():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {} -oft {} --range {} {} --backgroundModelFile {} --chromosomeSizes {}".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', outfile.name, 'all', 'bigwig', 20000, 200000, ROOT + 'background.txt', ROOT + 'hg19.chrom.sizes').split()
    chicExportData.main(args)
    file_obj_new = tarfile.open(outfile.name, "r")

    namelist = file_obj_new.getnames()
    assert len(namelist) == 12

    output_folder_new = mkdtemp(prefix="output_")
    output_folder_test_data = mkdtemp(prefix="output_")

    file_obj_new.extractall(output_folder_new)

    file_obj_test_data = tarfile.open(ROOT + 'chicExportData/interaction_all_bigwig.tar.gz', "r")
    file_obj_test_data.extractall(output_folder_test_data)
    namelist = file_obj_test_data.getnames()
    assert len(namelist) == 12
    files_new = os.listdir(output_folder_new)
    files_test_data = os.listdir(output_folder_test_data)

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal_bigwig(output_folder_new + '/' + file_new, output_folder_test_data + '/' + file_test_data, ['chr1'])


def test_interaction_bigwig_oneGene():
    output_folder_new = mkdtemp(prefix="output")
    args = "-f {} -o {} -om {} -omn {} -oft {} --range {} {} --backgroundModelFile {} --chromosomeSizes {}".format(
        ROOT + 'chicViewpoint/two_matrices.hdf5', output_folder_new + '/data.tar.gz', 'geneName', 'Eya1', 'bigwig', 20000, 200000, ROOT + 'background.txt', ROOT + 'hg19.chrom.sizes').split()
    chicExportData.main(args)
    files_new = os.listdir(output_folder_new)

    assert len(files_new) == 4
    name_list_oneGene_bigWig = ["background_FL-E13-5_chr1_genes_Eya1_interactions.bigwig",
                                'background_MB-E10-5_chr1_genes_Eya1_interactions.bigwig',
                                'FL-E13-5_chr1_genes_Eya1.bigwig',
                                'MB-E10-5_chr1_genes_Eya1.bigwig']
    files_test_data = name_list_oneGene_bigWig
    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal_bigwig(output_folder_new + '/' + file_new, ROOT + '/chicExportData/' + file_test_data, ['chr1'])


def test_significant_all():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {}".format(
        ROOT + 'chicSignificantInteractions/significantInteractions_dual.hdf5', outfile.name, 'all').split()
    chicExportData.main(args)

    file_obj_new = tarfile.open(outfile.name, "r")

    namelist = file_obj_new.getnames()
    assert len(namelist) == 6

    output_folder_new = mkdtemp(prefix="output_")
    output_folder_test_data = mkdtemp(prefix="output_")

    file_obj_new.extractall(output_folder_new)

    file_obj_test_data = tarfile.open(ROOT + 'chicExportData/significant_all.tar.gz', "r")
    file_obj_test_data.extractall(output_folder_test_data)

    files_new = os.listdir(output_folder_new)
    files_test_data = os.listdir(output_folder_test_data)

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, output_folder_test_data + '/' + file_test_data, skip=1)


def test_significant_oneGene():
    output_folder_new = mkdtemp(prefix="output")
    args = "-f {} -o {} -om {} -omn Eya1".format(
        ROOT + 'chicSignificantInteractions/significantInteractions_dual.hdf5', output_folder_new, 'geneName').split()
    chicExportData.main(args)

    files_new = os.listdir(output_folder_new)
    files_test_data = ['FL-E13-5_chr1_genes_Eya1_significant.txt', 'MB-E10-5_chr1_genes_Eya1_significant.txt']

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, ROOT + '/chicExportData/' + file_test_data, skip=1)


def test_target_all():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {}".format(
        ROOT + 'chicSignificantInteractions/targetFile_dual.hdf5', outfile.name, 'all').split()
    chicExportData.main(args)

    file_obj_new = tarfile.open(outfile.name, "r")

    namelist = file_obj_new.getnames()
    assert len(namelist) == 3

    output_folder_new = mkdtemp(prefix="output_")
    output_folder_test_data = mkdtemp(prefix="output_")

    file_obj_new.extractall(output_folder_new)

    file_obj_test_data = tarfile.open(ROOT + 'chicExportData/target_all.tar.gz', "r")
    file_obj_test_data.extractall(output_folder_test_data)

    files_new = os.listdir(output_folder_new)
    files_test_data = os.listdir(output_folder_test_data)

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, output_folder_test_data + '/' + file_test_data, skip=1)


def test_target_oneGene():
    output_folder_new = mkdtemp(prefix="output")
    args = "-f {} -o {} -om {} -omn Eya1".format(
        ROOT + 'chicSignificantInteractions/targetFile_dual.hdf5', output_folder_new, 'geneName').split()
    chicExportData.main(args)

    files_new = os.listdir(output_folder_new)
    files_test_data = ['FL-E13-5_chr1_MB-E10-5_chr1_genes_Eya1_target.txt']

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, ROOT + '/chicExportData/' + file_test_data, skip=1)


def test_aggregate_all():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {}".format(
        ROOT + 'chicAggregateStatistic/aggregate.hdf5', outfile.name, 'all').split()
    chicExportData.main(args)

    file_obj_new = tarfile.open(outfile.name, "r")

    namelist = file_obj_new.getnames()
    assert len(namelist) == 6

    output_folder_new = mkdtemp(prefix="output_")
    output_folder_test_data = mkdtemp(prefix="output_")

    file_obj_new.extractall(output_folder_new)

    file_obj_test_data = tarfile.open(ROOT + 'chicExportData/aggregate_all.tar.gz', "r")
    file_obj_test_data.extractall(output_folder_test_data)

    files_new = os.listdir(output_folder_new)
    files_test_data = os.listdir(output_folder_test_data)

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, output_folder_test_data + '/' + file_test_data, skip=1)


def test_aggregate_oneGene():
    output_folder_new = mkdtemp(prefix="output")
    args = "-f {} -o {} -om {} -omn Eya1".format(
        ROOT + 'chicAggregateStatistic/aggregate.hdf5', output_folder_new, 'geneName').split()
    chicExportData.main(args)

    files_new = os.listdir(output_folder_new)
    files_test_data = ['FL-E13-5_chr1_MB-E10-5_chr1_FL-E13-5_chr1_genes_Eya1_aggregate.txt',
                       'FL-E13-5_chr1_MB-E10-5_chr1_MB-E10-5_chr1_genes_Eya1_aggregate.txt']
    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, ROOT + '/chicExportData/' + file_test_data, skip=1)


def test_differential_all():
    outfile = NamedTemporaryFile(suffix='.tar.gz', delete=False)
    outfile.close()
    args = "-f {} -o {} -om {}".format(
        ROOT + 'chicDifferentialTest/differential.hdf5', outfile.name, 'all').split()
    chicExportData.main(args)

    file_obj_new = tarfile.open(outfile.name, "r")

    namelist = file_obj_new.getnames()
    assert len(namelist) == 9

    output_folder_new = mkdtemp(prefix="output_")
    output_folder_test_data = mkdtemp(prefix="output_")

    file_obj_new.extractall(output_folder_new)

    file_obj_test_data = tarfile.open(ROOT + 'chicExportData/differential_all.tar.gz', "r")
    file_obj_test_data.extractall(output_folder_test_data)

    files_new = os.listdir(output_folder_new)
    files_test_data = os.listdir(output_folder_test_data)

    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, output_folder_test_data + '/' + file_test_data, skip=1)


def test_differential_one_gene():
    output_folder_new = mkdtemp(prefix="output")
    args = "-f {} -o {} -om {} -omn Eya1".format(
        ROOT + 'chicDifferentialTest/differential.hdf5', output_folder_new, 'geneName').split()
    chicExportData.main(args)

    files_new = os.listdir(output_folder_new)
    files_test_data = ['FL-E13-5_chr1_MB-E10-5_chr1_genes_Eya1_accepted_differential.txt',
                       'FL-E13-5_chr1_MB-E10-5_chr1_genes_Eya1_all_differential.txt',
                       'FL-E13-5_chr1_MB-E10-5_chr1_genes_Eya1_rejected_differential.txt']
    files_new = sorted(files_new)
    files_test_data = sorted(files_test_data)

    for file_new, file_test_data in zip(files_new, files_test_data):
        assert are_files_equal(output_folder_new + '/' + file_new, ROOT + '/chicExportData/' + file_test_data, skip=1)
