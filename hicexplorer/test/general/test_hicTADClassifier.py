from hicexplorer import hicTADClassifier
from hicexplorer.test.test_compute_function import compute
from hicexplorer.lib.tadClassifier import TADClassifier
import math
import pandas as pd
import numpy as np
import numpy.testing as nt
import os
import shutil
from tempfile import mkdtemp, NamedTemporaryFile
from hicmatrix import HiCMatrix as hm
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


def test_MP_Domain_Data():

    protein_df = TADClassifier.MP_Domain_Data.readProtein(
        ROOT + 'unittest.regionPeak', None)
    assert protein_df.shape == (28, 6) and protein_df.iloc[0, 0] == 'chr1'

    domain_df = TADClassifier.MP_Domain_Data.read_domain_file(
        ROOT + 'unittest_tads.bed')
    assert domain_df.shape == (11, 9) and domain_df.iloc[0, 0] == 'chr1'

    domain_protein = TADClassifier.MP_Domain_Data.check_domains_against_protein(
        domain_df, protein_df, 10000, threshold=1.0)
    assert len(domain_protein) >= 1 and np.nanmin(
        domain_protein['signal_value']) >= 1.0

    tad_dict = TADClassifier.MP_Domain_Data.build_tad_dictionary(domain_df, [
                                                                 'chr1'])
    assert tad_dict['chr1'][1060000]

    domain_df_2 = TADClassifier.MP_Domain_Data.read_domain_file(
        ROOT + 'unittest_tads_2.bed')
    matrix_2 = TADClassifier.MP_Matrix(
        ROOT + 'unittest_matrix_2.cool', method='range')

    domain_dict_2 = TADClassifier.MP_Domain_Data.build_tad_dictionary(domain_df_2, [
                                                                      'chr1'])
    is_boundary_2 = matrix_2.get_boundary_positions(domain_dict_2)

    assert (is_boundary_2[is_boundary_2].shape[0] == 62 and is_boundary_2.shape == (
        2000,) and is_boundary_2[778] == True and is_boundary_2[779] == False)


def test_MP_Matrix():

    matrix = TADClassifier.MP_Matrix(
        ROOT + 'unittest_matrix.cool', method='range')
    assert (
        matrix.numpy_matrix.shape == (
            500,
            500) and
        matrix.numpy_matrix.max() <= 1 and matrix.numpy_matrix.min() >= 0)

    X_unittest = pd.read_csv(ROOT + 'X_unittest.csv')
    y_unittest = pd.read_csv(ROOT + 'y_unittest.csv')
    X_unittest = X_unittest.to_numpy()
    y_unittest = y_unittest.to_numpy()[:, 0]

    # unselect border cases
    X_unittest = TADClassifier.MP_Matrix.unselect_border_cases(X_unittest, 2)
    y_unittest = TADClassifier.MP_Matrix.unselect_border_cases_list(
        y_unittest, 2)
    assert (X_unittest.shape == (2, 8) and y_unittest.shape == (2,))

    # get positions
    positions = matrix.get_positions()
    assert (positions.shape == (500, 3) and positions[0, 0] == 'chr1')

    # build features
    matrix_2 = hm.hiCMatrix(ROOT + 'unittest_matrix_2.cool')
    features = TADClassifier.MP_Matrix.build_features(
        matrix_2.matrix, 2)
    test_value = 37.03932717672038

    assert (features.shape == (2000,
                               15) and math.isclose(features[10,
                                                             3],
                                                    test_value,
                                                    rel_tol=1e-5) and math.isclose(np.array(matrix_2.getMatrix())[8,
                                                                                                                  8],
                                                                                   features[10,
                                                                                            0],
                                                                                   rel_tol=1e-5))


def test_MP_Classifier():

    X_unittest = pd.read_csv(ROOT + 'X_unittest.csv')
    y_unittest = pd.read_csv(ROOT + 'y_unittest.csv')
    X_impure = X_unittest.copy(deep=True).to_numpy()
    X_unittest = X_unittest.to_numpy()
    y_unittest = y_unittest.to_numpy()[:, 0]

    # impute
    X_impure[2, 2] = np.nan
    X_imputed = TADClassifier.MP_Classifier.impute(X_impure, fill_value=-1.0)
    assert (X_imputed[2, 2] == -1.0 and X_imputed[3, 3] == 43.0)

    # resample
    X, y = TADClassifier.MP_Classifier.resample(
        X_unittest, y_unittest, method='undersample_random', passed_method=None, threads=4)
    assert y[y == True].shape == y[y == False].shape

    X, y = TADClassifier.MP_Classifier.resample(
        X_unittest, y_unittest, method='undersample_cluster_centroids', passed_method=None, threads=4)
    assert y[y == True].shape == y[y == False].shape

    # resample_test
    positions = np.zeros((5, 3))
    positions, X, y = TADClassifier.MP_Classifier.resample_test(
        positions, X_unittest, y_unittest, threads=4)
    assert y[y == True].shape == y[y == False].shape


# def test_get_domains():

#     # get_domains
#     y_unittest = pd.read_csv(
#         ROOT + 'y_unittest.csv',
#         header=None).to_numpy()[
#         :,
#         0].astype(bool)
#     positions_unittest = pd.read_csv(
#         ROOT + 'positions_unittest.csv',
#         header=None).to_numpy()

#     domain_df = TADClassifier.get_domains(positions_unittest, y_unittest)

#     print(domain_df.head)
#     assert domain_df['Chrom'].iloc[0] == 'chr1' and domain_df['Start'].iloc[0] == 20000


def test_hicTADClassifier():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()

    args = ['--matrices',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            outfile.name,
            '-n',
            'range',
            '--saved_classifier',
            ROOT + 'trained_model.bin.BIN',
            '--threads',
            '4'
            ]
    hicTADClassifier.main(args)
    # compute(hicTADClassifier.main, args, 1)
    domain_df = TADClassifier.MP_Domain_Data.read_domain_file(
        outfile.name)
    assert domain_df['Chrom'].iloc[0] == 1


def test_hicTADClassifier_two_matrices():
    #    args = ['--matrices',
    #      ROOT + 'unittest_matrix_2.cool',
    #      '--out_file',
    #      test_folder + 'domains_test.bed',
    #      '-n',
    #      'obs_exp',
    #      '--saved_classifier',
    #      'unittest_classifier.data',
    #      '--unselect_border_cases',
    #      'True',
    #      '--threads',
    #      '4'
    #     ]
    #
    #    compute(main, args, 1)
    #    domain_df = TADClassifier.MP_Domain_Data.read_domain_file(test_folder + 'domains_test.bed')
    #    assert domain_df['Chrom'].iloc[0] = 'chr1'
    outfile_2 = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_2.close()
    outfile_3 = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_3.close()
    args = ['--matrices',
            ROOT + 'gm12878_chr1.cool',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            outfile_2.name,
            outfile_3.name,
            '-n',
            'range',
            '--saved_classifier',
            ROOT + 'trained_model.bin.BIN',
            '--unselect_border_cases',
            # 'True',
            '--threads',
            '4'
            ]
    print(args)
    hicTADClassifier.main(args)
    # compute(hicTADClassifier.main, args, 1)
    domain_df = TADClassifier.MP_Domain_Data.read_domain_file(
        outfile_2.name)
    assert domain_df['Chrom'].iloc[0] == 1
    domain_df_3 = TADClassifier.MP_Domain_Data.read_domain_file(
        outfile_3.name)
    assert domain_df_3['Chrom'].iloc[0] == 1

    # shutil.rmtree(test_folder)


def test_hicTADClassifier_load_model_obs_exp():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()

    args = ['--matrices',
            ROOT + 'gm12878_chr1.cool',
            '--out_file',
            outfile.name,
            '-n',
            'obs_exp',
            '--threads',
            '4'
            ]
    hicTADClassifier.main(args)
    # compute(hicTADClassifier.main, args, 1)
    domain_df = TADClassifier.MP_Domain_Data.read_domain_file(
        outfile.name)
    assert domain_df['Chrom'].iloc[0] == 1


# def test_hicTADClassifier_load_model_range():
#     outfile = NamedTemporaryFile(suffix='.bed', delete=False)
#     outfile.close()

#     args = ['--matrices',
#             ROOT + 'gm12878_chr1.cool',
#             '--out_file',
#             outfile.name,
#             '-n',
#             'range',
#             '--threads',
#             '4'
#             ]
#     hicTADClassifier.main(args)
#     # compute(hicTADClassifier.main, args, 1)
#     domain_df = TADClassifier.MP_Domain_Data.read_domain_file(
#         outfile.name)
#     assert domain_df['Chrom'].iloc[0] == 1
