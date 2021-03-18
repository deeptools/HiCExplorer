from hicexplorer import hicTrainTADClassifier
from hicexplorer.test.test_compute_function import compute
from hicexplorer.lib_hicTADClassifier import TADClassifier
import numpy.testing as nt
import os
import shutil
from tempfile import mkdtemp
from hicmatrix import HiCMatrix
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

# from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__))),
    "test_data/")
# ROOT = 'test_data/hicTADClassifier/'


def test_hicTrainClassifier():

    test_folder = mkdtemp(prefix="testfolder_hicTrainClassifier")

    args = ['--mode',
            'train_test',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrix_file',
            ROOT + 'unittest_matrix_2.cool',
            '--out_file',
            test_folder + 'train_test.txt',
            '-n',
            'range',
            '-r',
            '10000',
            '--estimators_per_step',
            '10']

    compute(hicTrainClassifier.main, args, 1)
    f = open(test_folder + 'train_test.txt', "r")
    assert f.readline().split()[0] == 'accuracy'

    args = ['--mode',
            'train_new',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrix_file',
            ROOT + 'unittest_matrix_2.cool',
            '--out_file',
            test_folder + 'unittest_classifier_new',
            '-n',
            'range',
            '-r',
            '10000',
            '--estimators_per_step',
            '10']

    compute(hicTrainClassifier.main, args, 1)
    # if estimator is build will be checked by the next run

    args = ['--mode',
            'train_existing',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrix_file',
            ROOT + 'unittest_matrix_2.cool',
            '--out_file',
            test_folder + 'unittest_classifier_new',
            '-n',
            'range',
            '-r',
            '10000',
            '--saved_classifier',
            test_folder + 'unittest_classifier_new.BIN']

    compute(hicTrainClassifier.main, args, 1)

    args = ['--mode',
            'predict_test',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrix_file',
            ROOT + 'unittest_matrix_2.cool',
            '--out_file',
            test_folder + 'predict_test',
            '-n',
            'range',
            '-r',
            '10000',
            '--saved_classifier',
            test_folder + 'unittest_classifier_new.BIN']

    compute(hicTrainClassifier.main, args, 1)
    f = open(test_folder + 'predict_test_results.txt', "r")
    assert f.readline().split()[0] == 'accuracy'

    shutil.rmtree(test_folder)
