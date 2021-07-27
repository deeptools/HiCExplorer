from hicexplorer import hicTrainTADClassifier
from hicexplorer.test.test_compute_function import compute
# from hicexplorer.lib.tadClassifier import TADClassifier
import numpy.testing as nt
import os
import shutil
from tempfile import mkdtemp
import tempfile
from hicmatrix import HiCMatrix
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

# from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__))),
    "test_data/hicTADClassifier/")
# ROOT = 'test_data/hicTADClassifier/'


def test_hicTrainClassifier_train_test():

    test_folder = mkdtemp(prefix="testfolder_hicTrainClassifier")
    print(tempfile.gettempdir())

    args = ['--mode',
            'train_test',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrices',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            test_folder + 'train_test.txt',
            '-n',
            'obs_exp',
            '-r',
            '10000',
            '--estimators_per_step',
            '10']
    # '--chrPrefixProtein', 'remove']
    # hicTrainTADClassifier.main(args)
    compute(hicTrainTADClassifier.main, args, 5)
    f = open(test_folder + 'train_test.txt', "r")
    assert f.readline().split()[0] == 'accuracy'

# def test_hicTrainClassifier_train_new():
    args = ['--mode',
            'train_new',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrices',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            test_folder + 'unittest_classifier_new',
            '-n',
            'range',
            '-r',
            '10000',
            '--estimators_per_step',
            '10']

    compute(hicTrainTADClassifier.main, args, 5)
    # if estimator is build will be checked by the next run

    args = ['--mode',
            'train_existing',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrices',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            test_folder + 'unittest_classifier_new',
            '-n',
            'range',
            '-r',
            '10000',
            '--saved_classifier',
            test_folder + 'unittest_classifier_new.BIN']

    compute(hicTrainTADClassifier.main, args, 5)

    args = ['--mode',
            'predict_test',
            '--domain_file',
            ROOT + 'unittest_domains.bed',
            '--matrices',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            test_folder + 'predict_test',
            '-n',
            'range',
            '-r',
            '10000',
            '--saved_classifier',
            test_folder + 'unittest_classifier_new.BIN']

    compute(hicTrainTADClassifier.main, args, 5)
    f = open(test_folder + 'predict_test_results.txt', "r")
    assert f.readline().split()[0] == 'accuracy'

#     shutil.rmtree(test_folder)
