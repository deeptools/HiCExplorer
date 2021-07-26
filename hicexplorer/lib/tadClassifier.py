# import statements

# classic python libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import Generator
from os import EX_CANTCREAT, walk
from os import listdir
import os
import re
import math
import argparse
import warnings
import logging
import scipy.stats as statsw
from scipy.sparse import csr_matrix, lil_matrix
import pickle
import site
# hicexplorer and pybedtools
from hicexplorer import hicFindTADs
from hicexplorer import hicTransform
from hicmatrix import HiCMatrix as hm
from pybedtools import BedTool
import cooler
from hicexplorer.utilities import obs_exp_matrix
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros

# machine learning libraries
# https://imbalanced-learn.readthedocs.io/en/stable/
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import BaggingClassifier
from sklearn import metrics
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.impute import SimpleImputer
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc

from cleanlab.classification import LearningWithNoisyLabels
from imblearn.under_sampling import *

# PCA
# from sklearn import decomposition


log = logging.getLogger(__name__)


class TADClassifier:
    '''wrapper class for hicTADClassifier program'''

    pretrained_classifier_10000_obsexp = 'pretrained_classifier_10000_obsexp.BIN'
    pretrained_classifier_10000_range = 'pretrained_classifier_10000_range.BIN'

    class MP_Domain_Data:
        '''represents domain and protein information and implements helper functions for their preparation'''

        def __init__(self, domain_file, protein_file=None,
                     threshold=None, leniency=0, resolution=10000, pAddRemoveChrPrexix=None):
            '''read the necessary files and check domains against proteins'''

            # read domain file
            self.domain_df = TADClassifier.MP_Domain_Data.read_domain_file(
                domain_file)

            # read protein file and intersect if necessary
            if(protein_file is not None):
                # if pAddRemoveChrPrexix:

                self.protein_df = TADClassifier.MP_Domain_Data.readProtein(
                    protein_file, pAddRemoveChrPrexix)
                TADClassifier.MP_Domain_Data.apply_binning_and_leniency(
                    self.protein_df, pBinSize=resolution, leniency=leniency)
                self.domain_df = TADClassifier.MP_Domain_Data.check_domains_against_protein(
                    self.domain_df, self.protein_df, resolution=resolution, threshold=threshold)

            # build dictionary for tads from domain file
            self.chromosomes = self.domain_df['Chrom'].unique().tolist()
            self.domain_dict = TADClassifier.MP_Domain_Data.build_tad_dictionary(
                self.domain_df, chromosomes=self.chromosomes)

        def readProtein(pFile, pAddChr):
            '''read in a bed protein file (add pAddChr, if chr are in single digit form)'''

            # read
            protein_df = pd.read_csv(pFile, sep='\t', header=None)[
                [0, 1, 2, 6, 7, 8]]
            if pAddChr is None:
                pass
            elif pAddChr:
                protein_df[0] = 'chr' + protein_df[0].astype(str)
            elif not pAddChr:
                # protein_df[0] = 'chr' + protein_df[0].astype(str)
                protein_df[0] = protein_df[0].str.lstrip('chr')
            log.debug('protein_df {}'.format(protein_df))
            if(protein_df.size < 1):
                raise ValueError('empty protein file passed')

            # use bedtools sorting and apply proper column names
            protein_df_bedtool = BedTool.from_dataframe(protein_df)
            protein_df = protein_df_bedtool.sort().to_dataframe(
                disable_auto_names=True, header=None)
            protein_df.columns = TADClassifier.MP_Domain_Data.get_protein_col_names()
            return protein_df

        def read_domain_file(tad_file):
            '''read in the domain.bed output file of hicFindTads'''

            # read and sort
            domains_by_hicfindtads = pd.read_csv(
                tad_file, sep="\t", header=None)
            domains_by_hicfindtads.columns = TADClassifier.MP_Domain_Data.get_domain_col_names()
            domains_by_hicfindtads = domains_by_hicfindtads.sort_values(
                by=["Chrom", "Start"])

            if(domains_by_hicfindtads.size < 1):
                raise ValueError('empty domain file passed')

            return domains_by_hicfindtads

        def check_domains_against_protein(
                domain_df, protein_df, resolution, threshold=None):
            '''check the given tads in the domain df against the protein peaks in the protein file with a threshold'''

            # build bedtools from dfs
            domain_df['Name'] = domain_df['End']
            domain_df['End'] = domain_df['Start'] + resolution

            tad_bedtool = BedTool.from_dataframe(
                domain_df[["Chrom", "Start", "End", "Name"]])
            protein_bedtool = BedTool.from_dataframe(protein_df)

            # intersect domains and proteins
            domain_protein = tad_bedtool.intersect(
                protein_bedtool, wa=True, wb=True).to_dataframe()
            domain_protein.columns = TADClassifier.MP_Domain_Data.get_domain_protein_col_names()

            # aggregate multiple peaks at one domain
            domain_protein = domain_protein.groupby(['Chrom', 'Start', 'End', 'Name']).agg(
                signal_value=pd.NamedAgg(column='signalValue', aggfunc=max),
                p_value=pd.NamedAgg(column='pValue', aggfunc=max),
                q_value=pd.NamedAgg(column='qValue', aggfunc=max)
            )

            # select only domains with certain threshold if given
            if threshold is not None:
                mask = domain_protein['signal_value'] >= threshold
                domain_protein = domain_protein[mask]

            # build output
            domain_protein.reset_index(inplace=True)
            domain_protein['End'] = domain_protein['Name']

            # debug information; TODO: move
            log.debug('protein peaks: {}'.format(len(protein_df)))
            log.debug('TADs: {}'.format(len(domain_df)))
            log.debug('Matched TADs: {}'.format(len(domain_protein)))

            return domain_protein

        def apply_binning_and_leniency(pDataFrame, pBinSize=10000, leniency=0):
            """bin the given protein file and make peaks wider if necessary"""

            pDataFrame_out = pDataFrame.copy()
            pDataFrame_out['Start'] = (
                pDataFrame['Start'] / pBinSize).astype(int) * pBinSize
            pDataFrame_out['End'] = (
                (pDataFrame['End'] / pBinSize).astype(int) + 1) * pBinSize

            if (leniency > 0):
                pDataFrame_out['Start'] = np.maximum(
                    0, pDataFrame_out['Start'] - int(pBinSize * leniency))
                pDataFrame_out['End'] = pDataFrame_out['End'] + \
                    int(pBinSize * leniency)

            pDataFrame_out.drop_duplicates()
            bedtools_data = BedTool.from_dataframe(pDataFrame_out)
            bedtools_data = bedtools_data.merge()
            bedtools_data = bedtools_data.sort()

            return bedtools_data.to_dataframe()

        def get_protein_col_names():
            '''get column names'''

            return ['Chrom', 'Start', 'End', 'signalValue', 'pValue', 'qValue']

        def get_domain_col_names():
            '''get column names'''

            return ["Chrom", "Start", "End", "Name", "Score",
                    "Strand", "ThickStart", "ThickEnd", "ItemRGB"]

        def get_domain_protein_col_names():
            '''get column names'''

            return ["Chrom", "Start", "End", "Name", "Chrom_p",
                    "Start_p", "End_p", 'signalValue', 'pValue', 'qValue']

        def build_tad_dictionary(domains, chromosomes):
            '''get a dictionary based on domain file for multiple chromosomes'''

            # for given chromosome names, build dictionary from dataframe
            # this is done for performance reasons only
            domain_dicts = {}
            for chromosome in chromosomes:

                # add current chromosomes, then add all positions
                chr_domains = domains[domains['Chrom'] == chromosome]
                b = chr_domains["Start"]
                chr_dict = dict((int(position), True) for position in b)
                domain_dicts[str(chromosome)] = chr_dict

            return domain_dicts

    class MP_Matrix:
        '''acts as a wrapper class for a HiCMatrix Object and implements helper functions for matrix data preparation'''
        # currently supports matrices of one chromosome

        def __init__(self, matrix_file, method=None, range_max=None, pChromosome=None, pThreads=None):

            # use input directly, if already normalized
            if method == 'obs_exp':
                hic_ma, hic_ma_np = TADClassifier.MP_Matrix.read_matrix_file(
                    matrix_file, pChromosome)
                hic_ma = TADClassifier.MP_Matrix.obs_exp_normalization(hic_ma, pThreads=pThreads)
                # hic_ma_np = np.array(hic_ma.getMatrix())
                # hic_ma_np = hic_ma.matrix

            # perform range normalization, with given max value or infer from
            # matrix
            elif method == 'range':
                hic_ma, hic_ma_np = TADClassifier.MP_Matrix.read_matrix_file(
                    matrix_file, pChromosome)
                range_min = 0.0
                o_min = 0.0
                o_max = 1.0

                if range_max is None:
                    # range_max = np.nanmax(hic_ma_np)
                    range_max = hic_ma_np.max()

                # elif (range_max < np.nanmax(hic_ma_np)):
                elif (range_max < hic_ma_np.max()):

                    raise ValueError('range maximum too low for input matrix')

                hic_ma_np = TADClassifier.MP_Matrix.range_normalization(
                    hic_ma_np, range_min, range_max, o_min, o_max)

            else:
                raise NotImplementedError

            self.numpy_matrix = hic_ma_np
            self.hic_matrix = hic_ma
            self.positions = self.get_positions()
            self.range_max = range_max
            self.resolution = None

        def get_positions(self):
            '''get positions for matrix'''

            # get start and end position for every bin in the matrix
            indices = np.arange(0, self.numpy_matrix.shape[0])
            # print('indices {}'.format(indices))
            # print('self.hic_matrix.cut_intervals {}'.format(self.hic_matrix.cut_intervals))
            vec_bin_pos = np.vectorize(self.hic_matrix.getBinPos)
            pos = vec_bin_pos(indices)

            # return it as an ordered array
            return np.transpose(np.array(pos)[0:3, :], (1, 0))

        def get_resolution(self):
            """resolution of matrix"""

            if(self.resolution is None):
                self.resolution = self.hic_matrix.getBinSize()

            return self.resolution

        def get_boundary_positions(self, domain_dicts):
            '''return a list of booleans, denoting if a boundary exists at the current index'''

            # build necessary function and df to check if there is a boundary
            # at position

            is_boundary = np.full(self.positions.shape[0], False, dtype=bool)
            # log.debug(domain_dicts)
            # check for all indices
            for k in range(self.positions.shape[0]):
                try:
                    is_boundary[k] = int(self.positions[k, 1]) in domain_dicts[self.positions[k, 0]]
                except KeyError:
                    is_boundary[k] = False
            return is_boundary

        def get_features(self, distance, use_gradient=False):
            '''build features for self'''

            # run build features for self
            features = TADClassifier.MP_Matrix.build_features(
                self.numpy_matrix, distance)
            if features is None:
                return None
            # experimental: use gradient matrix too
            if(use_gradient):
                x_gradient, y_gradient = self.build_gradient_features(distance)
                features = np.concatenate(
                    (features, x_gradient, y_gradient), axis=1)

            return features

        def build_features(numpy_matrix, distance):
            '''select features from input matrix'''

            # build necessary structures
            m_indices = np.arange(
                distance, numpy_matrix.shape[0] - distance, 1)
            m_list = []

            # get triangle at position index
            def get_triangle(index):
                start = index - distance
                end = index + distance + 1

                # log.debug('triu submatrix: {}'.format(numpy_matrix[start:end, start:end]))
                # log.debug('start {} end {}'.format(start, end))
                triangle = np.triu(numpy_matrix[start:end, start:end].todense(), k=0)
                triangle = triangle.astype(float)
                triangle[np.tril_indices(triangle.shape[0], -1)] = np.NINF
                mask = triangle != float('-inf')
                flattened_triangle = np.ndarray.flatten(triangle[mask])
                return flattened_triangle

            # run for all positions
            # can refactor for vectorization, but resampling and fitting takes
            # way longer anyway
            for i in m_indices:
                m_list.append(get_triangle(i))

            # build output
            # at boundaries: use filler, those will be unselected later anyways
            if len(m_list) == 0:
                return None
            features = np.stack(m_list, axis=1)
            sentinel1 = np.full((features.shape[0], distance), np.nan)
            sentinel2 = np.full((features.shape[0], distance), np.nan)
            features = np.concatenate([sentinel1, features, sentinel2], axis=1)

            return np.transpose(features, (1, 0))

        def unselect_border_cases(X, distance):
            '''unselect cases at border of matrix'''

            return X[(distance - 1):(X.shape[0] - distance), :]

        def unselect_border_cases_list(y, distance):
            '''unselect cases at border of matrix'''

            return y[(distance - 1):y.shape[0] - distance]

        def write_input_to_file(positions, X, y, out_file):
            '''save positions,features and is_boundary to file'''

            output = np.concatenate((positions, X, y[:, np.newaxis]), axis=1)
            np.savetxt(out_file, output, delimiter=";")

        def read_matrix_file(matrix_file, pChromosome):
            ''''reads a given cool file and returns its hiCMatrix and numpy representation'''
            # log.debug('matrix_file {}'.format(matrix_file))
            # log.debug('pChromosome {}'.format(pChromosome))
            # check if instance of string or file and load appropriate
            if isinstance(matrix_file, str):
                if pChromosome is not None:
                    hic_ma = hm.hiCMatrix(matrix_file, pChrnameList=[pChromosome])
                else:
                    hic_ma = hm.hiCMatrix(matrix_file)

                # log.debug('hic_ma: {}'.format(hic_ma.matrix))

            else:
                hic_ma = matrix_file

            # hic_ma_np = np.array(hic_ma.getMatrix())
            hic_ma_np = hic_ma.matrix

            return (hic_ma, hic_ma_np)

        def range_normalization(n, n_min, n_max, o_min, o_max):
            '''apply range normalization'''

            return o_min + (n - n_min) * (o_max - o_min) / (n_max - n_min)

        def obs_exp_normalization(hic_ma, pThreads=None):
            '''apply obs_exp normalization'''
            log.debug('obs/exp matrix computation...')

            trasf_matrix = lil_matrix(hic_ma.matrix.shape)

            # from hicTransformTADs
            def _obs_exp(pSubmatrix, pThreads=None):
                obs_exp_matrix_ = obs_exp_matrix(pSubmatrix, pThreads=pThreads, pDistance=100)
                obs_exp_matrix_ = convertNansToZeros(
                    csr_matrix(obs_exp_matrix_))
                obs_exp_matrix_ = convertInfsToZeros(
                    csr_matrix(obs_exp_matrix_))
                # if len(obs_exp_matrix_.data) == 0:
                # return np.array([[]])
                return obs_exp_matrix_  # .todense()

            for chrname in hic_ma.getChrNames():
                chr_range = hic_ma.getChrBinRange(chrname)
                submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
                submatrix.astype(float)
                obs_exp = _obs_exp(submatrix, pThreads)
                if obs_exp.nnz != 0:
                    trasf_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = lil_matrix(obs_exp)

            hic_ma.setMatrix(
                trasf_matrix.tocsr(),
                cut_intervals=hic_ma.cut_intervals)
            log.debug('obs/exp matrix computation... DONE')

            return hic_ma

    class MP_Classifier:
        '''acts as a wrapper for a classifier'''

        def __init__(self, resolution=10000,
                     normalization_method='range',
                     resampling_method='undersample_random',
                     alternative_resampling_method=None,
                     alternative_classifier=None,
                     impute_value=-1.0,
                     estimators_per_step=10,
                     use_cleanlab=False,
                     threads=4,
                     use_gradient=False,
                     distance=15):

            # set some attribute
            self.resolution = resolution
            self.normalization_method = normalization_method
            self.resampling_method = resampling_method
            self.alternative_resampling_method = alternative_resampling_method
            self.impute_value = impute_value
            self.estimators_per_step = estimators_per_step
            self.use_cleanlab = use_cleanlab
            self.threads = threads
            self.use_gradient = use_gradient
            self.distance = distance
            self.is_default_classifier = alternative_classifier is None

            if(alternative_classifier is None):
                self.classifier = BaggingClassifier(
                    base_estimator=AdaBoostClassifier(),
                    n_estimators=0,
                    warm_start=True,
                    n_jobs=threads)

            else:
                self.classifier = alternative_classifier

            if (use_cleanlab):
                estimators = self.classifier.get_params(
                )['n_estimators'] + self.estimators_per_step
                self.classifier.set_params(n_estimators=estimators)
                self.classifier = LearningWithNoisyLabels(
                    clf=self.classifier, n_jobs=self.threads)

        def single_run(self, positions, X, y, test_size_n=0.2):
            '''train and test the classifier on the given dataset'''

            # impute missing values
            X = TADClassifier.MP_Classifier.impute(X, self.impute_value)

            # division into training set and test set
            pos_train, pos_test, X_train, X_test, y_train, y_test = train_test_split(
                positions, X, y, test_size=test_size_n, random_state=0)

            # enact imbalance strategy
            X_train, y_train = TADClassifier.MP_Classifier.resample(
                X_train, y_train, self.resampling_method, self.alternative_resampling_method, threads=self.threads)

            # set estimators to real wanted value
            estimators = self.classifier.get_params(
            )['n_estimators'] + self.estimators_per_step
            self.classifier.set_params(n_estimators=estimators)

            # train model
            # if input contains only one class, an exception is thrown
            log.debug('fitting')
            if(np.unique(y).shape[0] >= 2):
                self.classifier.fit(X, y)

            else:
                warnings.warn(
                    'input does not contain boundaries or no non-boundaries, fitting aborted')

            # test model
            pos_test, X_test, y_test = TADClassifier.MP_Classifier.resample_test(
                pos_test, X_test, y_test, threads=self.threads)
            y_pred = self.classifier.predict(X_test)

            return (y_test, y_pred)

        def incremental_fit(self, X, y, resample=True):
            '''train existing model incrementally on a new dataset, by adding estimators'''

            # impute missing values
            X = TADClassifier.MP_Classifier.impute(X, self.impute_value)

            # enact oversample imbalance strategy
            if(resample):
                X, y = TADClassifier.MP_Classifier.resample(
                    X, y, self.resampling_method, self.alternative_resampling_method, threads=self.threads)

            # introduce additional estimators to accomodate for new dataset
            if(self.use_cleanlab):
                estimators = self.classifier.clf.get_params(
                )['n_estimators'] + self.estimators_per_step
                self.classifier.clf.set_params(n_estimators=estimators)

            else:
                estimators = self.classifier.get_params(
                )['n_estimators'] + self.estimators_per_step
                self.classifier.set_params(n_estimators=estimators)

            # train model
            # if input contains only one class, an exception is thrown
            log.debug('fitting')
            self.classifier.fit(X, y)

        def predict(self, positions, X):
            '''predict with existing model'''

            if X is None or len(X) == 0:
                return None, None
            # impute missing values
            X = TADClassifier.MP_Classifier.impute(X, self.impute_value)

            # test model
            log.debug('predicting')
            y_pred = self.classifier.predict(X)
            return positions, y_pred

        def predict_test(self, positions, X, y_test):

            # impute missing values
            X = TADClassifier.MP_Classifier.impute(X, self.impute_value)

            # enact resampling test strategy
            positions, X, y_test = TADClassifier.MP_Classifier.resample_test(
                positions, X, y_test, threads=self.threads)

            # test model
            log.debug('predicting')
            if(not y_test[y_test == False].shape[0] <= 0 and not y_test[y_test == True].shape[0] <= 0):
                y_pred = self.classifier.predict(X)
                return positions, X, y_test, y_pred
            else:
                raise ValueError(
                    'input does not contain boundaries and non-boundaries')

        def impute(X, fill_value=-1.0):
            '''impute missing nan values by setting them to fill_value'''

            # call Simple Imputer, which will fill all non real numbers with
            # another value
            log.debug('imputing with value: ' + str(fill_value))
            imp = SimpleImputer(
                missing_values=np.nan,
                strategy='constant',
                fill_value=-1.0)
            X = imp.fit_transform(X)

            return X

        def resample(X, y, method='undersample_random',
                     passed_method=None, threads=4):
            '''enact resample method chosen with method'''

            # resample method: KMEANs
            if(method == 'undersample_cluster_centroids'):
                log.debug('resampling with: ' + method)
                cc = ClusterCentroids(random_state=42, n_jobs=threads)
                return cc.fit_resample(X, y)

            # resample method: random
            elif(method == 'undersample_random'):
                log.debug('resampling with: ' + method)
                rus = RandomUnderSampler(random_state=42)
                return rus.fit_resample(X, y)

            # use users method
            elif(method == 'passed_method' and passed_method is not None):
                log.debug('resampling with: ' + method)
                return passed_method.fit_resample(X, y)

            # just return, if none is chosen
            else:
                return X, y

        def resample_test(positions, X, y, threads=4):
            '''call imblearn random undersample for prediction'''

            # resample test set, for correct accuracy score
            if (np.unique(y).shape[0] > 1):
                pos_X = np.concatenate((positions, X), axis=1)
                rus = RandomUnderSampler(random_state=42)
                pos_X, y = rus.fit_resample(pos_X, y)
                positions = pos_X[:, 0:3]
                X = pos_X[:, 3:(pos_X.shape[1])]

            return positions, X, y

        def print_results(self, y_test, y_pred, out_file=None):
            '''print accuracy, confusion matrix and classification report to logger'''

            # important: perform a balanced test predict using y_test!!!
            acc = "accuracy score: " + \
                str(accuracy_score(y_test, y_pred)) + '\n\n'
            classification = "classification report:\n" + \
                classification_report(y_test, y_pred) + '\n'
            conf = "confusion matrix:\n" + \
                str(confusion_matrix(y_test, y_pred))

            log.debug(acc)
            log.debug(classification)
            log.debug(conf)

            if(out_file):

                f = open(out_file, "w")
                f.writelines(acc)
                f.writelines(classification)
                f.writelines(conf)
                f.close()

        def get_feature_importance(self, out_file):
            '''print feature importances'''

            if(self.is_default_classifier and not self.use_cleanlab):
                # only use on Bagging Classifier
                feature_importances = np.mean(
                    [est.feature_importances_ for est in self.classifier.estimators_], axis=0)
                ind = np.arange(0, feature_importances.shape[0])

                barlist = plt.bar(ind, feature_importances)

                d = self.distance * 2
                off = 0
                barlist[off].set_color('red')

                while(d > 0):
                    off = d + off
                    barlist[off].set_color('red')
                    d = d - 1

                plt.savefig(out_file)

        def get_roc(self, X_test, y_test, out_file):
            '''get receiver operating characteristic'''

            if(self.is_default_classifier and not self.use_cleanlab):
                # from
                # https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html
                fpr = dict()
                tpr = dict()
                roc_auc = dict()

                y_score = self.classifier.decision_function(X_test)
                y_test = y_test.astype(int)

                # Compute micro-average ROC curve and ROC area
                fpr["micro"], tpr["micro"], _ = roc_curve(
                    y_test.ravel(), y_score.ravel())
                roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

                plt.figure()
                lw = 2

                plt.plot(fpr["micro"], tpr["micro"], color='aqua',
                         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc["micro"])

                plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Receiver operating characteristic')
                plt.legend(loc="lower right")
                plt.savefig(out_file)

    # TAD classifier program

    def __init__(self, mode,
                 out_file,
                 normalization_method='range',
                 saved_classifier=None,
                 unselect_border_cases=True,
                 threads=4,
                 threshold=None,
                 leniency=0,
                 resolution=10000,
                 distance=15,
                 impute_value=-1.0,
                 resampling_method='undersample_random',
                 alternative_resampling_method=None,
                 alternative_classifier=None,
                 use_cleanlab=False,
                 estimators_per_step=50,
                 concatenate_before_resample=False,
                 pAddRemoveChrPrexix=None
                 ):

        self.mode = mode
        self.out_file = out_file
        self.threshold = threshold
        self.leniency = leniency
        self.unselect_border_cases = unselect_border_cases
        self.concatenate_before_resample = concatenate_before_resample
        self.addRemoveChrPrexix = pAddRemoveChrPrexix

        if(mode == 'predict' or mode == 'train_existing' or mode == 'predict_test'):
            if(saved_classifier is not None):
                try:
                    self.classifier = pickle.load(open(saved_classifier, 'rb'))
                except Exception as exp:
                    log.error('Tried to load ML model but failed! {}'.format(str(exp)))
                    exit(1)

            else:
                model_location = site.getsitepackages()[0] + '/hicexplorer/trained_models/'

                if(normalization_method == 'obs_exp'):
                    if resolution == 10000:
                        model_location += '10kb_model_cleanlab_proteins_obs_exp.BIN'
                    elif resolution == 25000:
                        model_location += '25kb_model_cleanlab_proteins_obs_exp.BIN'
                    elif resolution == 50000:
                        model_location += '50kb_model_cleanlab_proteins_obs_exp.BIN'
                    elif resolution == 100000:
                        model_location += '100kb_model_cleanlab_proteins_obs_exp.BIN'
                    else:
                        log.error('No trained model for this resolution! Please train a model on your own with hicTrainTADClassifier!')
                        exit(1)

                    log.debug('using obs/exp default classifier')
                else:
                    if resolution == 10000:
                        model_location += '10kb_model_cleanlab_proteins_range.BIN'
                    elif resolution == 25000:
                        model_location += '25kb_model_cleanlab_proteins_range.BIN'
                    elif resolution == 50000:
                        model_location += '50kb_model_cleanlab_proteins_range.BIN'
                    elif resolution == 100000:
                        model_location += '100kb_model_cleanlab_proteins_range.BIN'
                    else:
                        log.error('No trained model for this resolution! Please train a model on your own with hicTrainTADClassifier!')
                        exit(1)

                try:
                    self.classifier = pickle.load(
                        open(model_location, 'rb'))
                except Exception as exp:
                    log.error('Tried to load ML model but failed! {}'.format(str(exp)))
                    exit(1)
        if(mode == 'train_new' or mode == 'train_test'):
            if(alternative_classifier is not None and not isinstance(alternative_classifier, sklearn.base.BaseEstimator)):
                raise ValueError('the passed classifier is not valid')

            if(not resampling_method == 'passed_method' and alternative_resampling_method is not None):
                warnings.warn(
                    'default resampling method chosen, but custom method passed; using passed method')
                resampling_method = 'passed_method'

            if(alternative_resampling_method is not None and not isinstance(alternative_resampling_method, imblearn.base.BaseCleaningSampler)):
                raise ValueError('the passed resampling method is not valid')

            self.classifier = TADClassifier.MP_Classifier(
                resolution=resolution,
                normalization_method=normalization_method,
                resampling_method=resampling_method,
                alternative_resampling_method=alternative_resampling_method,
                alternative_classifier=alternative_classifier,
                impute_value=impute_value,
                estimators_per_step=estimators_per_step,
                use_cleanlab=use_cleanlab,
                use_gradient=False,
                distance=distance,
                threads=threads)

    def prepare_train(self, matrix_file, domain_file, protein_file, pChromosome, ):
        '''prepare matrix and its derivatives for the run'''

        log.debug('preparing domain data')
        prep = TADClassifier.MP_Domain_Data(
            domain_file,
            protein_file,
            resolution=self.classifier.resolution,
            threshold=self.threshold,
            leniency=self.leniency,
            pAddRemoveChrPrexix=self.addRemoveChrPrexix)
        domain_dict = prep.domain_dict

        log.debug('loading matrix')
        # ingest matrix
        matrix = TADClassifier.MP_Matrix(matrix_file,
                                         method=self.classifier.normalization_method, pChromosome=pChromosome)

        # build inputs for classifier
        log.debug('build features')
        is_boundary = matrix.get_boundary_positions(domain_dict)
        features = matrix.get_features(
            self.classifier.distance,
            self.classifier.use_gradient)
        if features is None:
            return matrix, None, features, is_boundary
        positions = matrix.positions

        if(not self.classifier.resolution == matrix.get_resolution()):
            warnings.warn(
                'training matrix with resolution {} on classifier with resolution {}'.format(
                    self.classifier.resolution,
                    matrix.get_resolution()))

        # get rid of cases at the border of matrix
        if(self.unselect_border_cases):
            features = TADClassifier.MP_Matrix.unselect_border_cases(
                features, self.classifier.distance)
            positions = TADClassifier.MP_Matrix.unselect_border_cases(
                positions, self.classifier.distance)
            is_boundary = TADClassifier.MP_Matrix.unselect_border_cases_list(
                is_boundary, self.classifier.distance)

        return matrix, positions, features, is_boundary

    def prepare_predict(self, matrix_file, pChromosome):
        '''prepare matrix and its derivatives for the run'''

        log.debug('loading matrix')
        # ingest matrix
        matrix = TADClassifier.MP_Matrix(matrix_file,
                                         method=self.classifier.normalization_method, pChromosome=pChromosome)

        # build inputs for classifier
        log.debug('build features')
        features = matrix.get_features(
            self.classifier.distance,
            self.classifier.use_gradient)
        positions = matrix.positions

        # get rid of cases at the border of matrix
        if(self.unselect_border_cases):
            features = TADClassifier.MP_Matrix.unselect_border_cases(
                features, self.classifier.distance)
            positions = TADClassifier.MP_Matrix.unselect_border_cases(
                positions, self.classifier.distance)

        return matrix, positions, features

    def train_test(self, matrix_list, domain_list, protein_list, pChromosome):
        '''perform train_test program mode'''

        # prepare
        matrix, positions, features, is_boundary = self.prepare_train(
            matrix_list[0], domain_list[0], protein_list[0], pChromosome=pChromosome[0][0])

        out_file_s = self.out_file.split('.')

        if(out_file_s[-1] == 'txt'):
            out_file = self.out_file
        else:
            out_file = self.out_file + '.txt'

        # run and print classifier output
        try:
            y_test, y_pred = self.classifier.single_run(
                positions, features, is_boundary)
            self.classifier.print_results(
                y_test, y_pred, out_file=out_file)
        except ValueError:
            raise ValueError(
                'training or test set does not contain any boundaries')

    def multi_train(self, matrix_list, domain_list,
                    protein_list, nm_conc_features=1, pChromosome=None):
        '''perform train and on multiple matrices in list, resample separatly, but train once'''

        # note, that the parameter nm_conc_features can be used to resample a
        # few matrices together
        nm_conc_features = nm_conc_features - 1
        i = 0
        conc_features = None
        conc_is_boundary = None
        resampled_X = None
        resampled_y = None

        for j, data in enumerate(zip(matrix_list, domain_list, protein_list, pChromosome)):
            for chromosome in pChromosome[j]:
                matrix, positions, features, is_boundary = self.prepare_train(
                    data[0], data[1], data[2], pChromosome=chromosome)

                if features is None:
                    continue
                matrix.numpy_matrix = None
                matrix.hic_matrix = None
                matrix = None

                if(conc_features is None):
                    conc_features = features
                    conc_is_boundary = is_boundary

                else:
                    conc_features = np.concatenate(
                        (conc_features, features), axis=0)
                    conc_is_boundary = np.concatenate(
                        (conc_is_boundary, is_boundary), axis=0)

                if(i == nm_conc_features):
                    resampled_X, resampled_y = self.incremental_resample(
                        resampled_X, resampled_y, conc_features, conc_is_boundary)
                    i = 0
                    conc_features = None
                    conc_is_boundary = None

                else:
                    i = i + 1

                features = None
                is_boundary = None

        if (conc_features is not None):
            resampled_X, resampled_y = self.incremental_resample(
                resampled_X, resampled_y, conc_features, conc_is_boundary)
            conc_features = None
            conc_is_boundary = None

        try:
            self.classifier.incremental_fit(
                resampled_X, resampled_y, resample=False)
            self.persist(self.out_file)

        except ValueError:
            warnings.warn(
                'training matrix does not contain any boundaries, training skipped')

    def multi_train_concatenate_before_resample(
            self, matrix_list, domain_list, protein_list, pChromosome, nm_conc_features=1000000):
        '''perform train on multiple matrice, concatenate the features of each matrices and then resample and train'''

        # use nm_conc_features to reduce the number of matrices, that are
        # concatenated
        nm_conc_features = nm_conc_features - 1
        i = 0
        conc_features = None
        conc_is_boundary = None

        for j, data in enumerate(zip(matrix_list, domain_list, protein_list)):
            for chromosome in pChromosome[j]:
                matrix, positions, features, is_boundary = self.prepare_train(
                    data[0], data[1], data[2], pChromosome=chromosome)

                matrix.numpy_matrix = None
                matrix.hic_matrix = None
                matrix = None

                if(conc_features is None):
                    conc_features = features
                    conc_is_boundary = is_boundary

                else:
                    conc_features = np.concatenate(
                        (conc_features, features), axis=0)
                    conc_is_boundary = np.concatenate(
                        (conc_is_boundary, is_boundary), axis=0)

                if(i == nm_conc_features):
                    try:
                        self.classifier.incremental_fit(
                            conc_features, conc_is_boundary, resample=True)
                    except ValueError:
                        warnings.warn(
                            'training matrix does not contain any boundaries, training skipped')

                    i = 0
                    conc_features = None
                    conc_is_boundary = None

                else:
                    i = i + 1

                features = None
                is_boundary = None

        if (conc_features is not None):
            try:
                self.classifier.incremental_fit(
                    conc_features, conc_is_boundary, resample=True)
            except ValueError:
                warnings.warn(
                    'training matrix does not contain any boundaries, training skipped')

            conc_features = None
            conc_is_boundary = None

        self.persist(self.out_file)

    def multi_test_predict(self, matrix_list, domain_list, protein_list, pChromosome):
        '''predict a balanced dataset for a meaningful classification report'''

        conc_test_boundary = None
        conc_is_boundary = None
        conc_positions = None

        # ROC
        conc_features = None

        for i, data in enumerate(zip(matrix_list, domain_list, protein_list)):

            for chromosome in pChromosome[i]:
                matrix, positions, features, is_boundary = self.prepare_train(
                    data[0], data[1], data[2], pChromosome=chromosome)

                matrix.numpy_matrix = None
                matrix.hic_matrix = None
                matrix = None

                try:
                    positions, features, y_test, y_pred = self.classifier.predict_test(
                        positions, features, y_test=is_boundary)

                    if(conc_test_boundary is None):
                        conc_test_boundary = y_test
                        conc_is_boundary = y_pred
                        conc_positions = positions

                        # ROC
                        conc_features = features

                    else:
                        conc_test_boundary = np.concatenate(
                            (conc_test_boundary, y_test), axis=0)
                        conc_is_boundary = np.concatenate(
                            (conc_is_boundary, y_pred), axis=0)
                        conc_positions = np.concatenate(
                            (conc_positions, positions), axis=0)

                        # ROC
                        conc_features = np.concatenate(
                            (conc_features, features), axis=0)

                except BaseException:
                    log.debug('one of the inputs did not contain boundaries')

                is_boundary = None
                positions = None
                features = None

        out_file_results = self.out_file + '_results.txt'
        out_file_fi = self.out_file + '_feature_importance.png'
        out_file_roc = self.out_file + '_roc.png'

        self.classifier.get_feature_importance(out_file_fi)
        self.classifier.print_results(
            conc_test_boundary,
            conc_is_boundary,
            out_file=out_file_results)
        self.classifier.get_roc(
            conc_features,
            conc_test_boundary,
            out_file_roc)
        # TADClassifier.print_to_bed(TADClassifier.get_domains(conc_positions,conc_is_boundary),self.out_file)

    def incremental_resample(self, X, y, X_i, y_i):
        '''resample without fitting'''

        try:
            if(X is None):
                X, y = TADClassifier.MP_Classifier.resample(
                    X_i, y_i, method=self.classifier.resampling_method, passed_method=self.classifier.alternative_resampling_method, threads=self.classifier.threads)

            else:
                X_i, y_i = TADClassifier.MP_Classifier.resample(
                    X_i, y_i, method=self.classifier.resampling_method, passed_method=self.classifier.alternative_resampling_method, threads=self.classifier.threads)
                X = np.concatenate((X, X_i), axis=0)
                y = np.concatenate((y, y_i), axis=0)

            return X, y

        except ValueError:
            warnings.warn(
                'training matrix does not contain boundaries, training skipped')

            return X, y

    def persist(self, out_file):

        out_file_s = out_file.split('.')

        if(not out_file_s[-1] == 'BIN'):
            out_file = out_file + '.BIN'

        pickle.dump(self.classifier, open(out_file, 'wb'))

    def run_hicTrainClassifier(self, matrix_list, domain_list, protein_list, pChromosomes):
        '''run hicTrainClassifier with specified mode on file'''

        if(isinstance(matrix_list, str)):
            matrix_list = [matrix_list]

        if(isinstance(domain_list, str)):
            domain_list = [domain_list]

        if(len(domain_list) == 1):
            domain_list = domain_list * len(matrix_list)

        if(protein_list is None):
            protein_list = [None] * len(matrix_list)

        elif(isinstance(protein_list, str)):
            protein_list = [protein_list]

        if(len(protein_list) == 1):
            protein_list = protein_list * len(matrix_list)

        if(not (len(matrix_list) == len(domain_list) and len(protein_list) == len(domain_list))):
            raise ValueError(
                'please pass domain (,optional protein) and matrix lists of same length or pass a single domain (and optional protein) file')

        chromosome_list = []
        if pChromosomes is not None:
            chromosome_list = pChromosomes
        else:
            for matrix in matrix_list:
                cooler_obj = cooler.Cooler(matrix)
                chromosome_list.append(cooler_obj.chromnames)

        if(self.mode == 'train_new' or self.mode == 'train_existing'):

            if(self.concatenate_before_resample):
                self.multi_train_concatenate_before_resample(
                    matrix_list, domain_list, protein_list, pChromosome=chromosome_list)

            else:
                self.multi_train(matrix_list, domain_list, protein_list, pChromosome=chromosome_list)

        elif(self.mode == 'train_test'):
            self.train_test(matrix_list, domain_list, protein_list, pChromosome=chromosome_list)

        elif(self.mode == 'predict_test'):
            self.multi_test_predict(matrix_list, domain_list, protein_list, pChromosome=chromosome_list)

    def run_hicTADClassifier(self, matrix_list, pChromosomes):
        '''predict a dataset'''

        if(isinstance(matrix_list, str)):
            matrix_list = [matrix_list]

        chromosome_list = []
        if pChromosomes is None:
            for matrix in matrix_list:
                cooler_obj = cooler.Cooler(matrix)
                chromosome_list.append(cooler_obj.chromnames)
        else:
            chromosome_list = pChromosomes
        # i = 0
        conc_is_boundary = None
        conc_positions = None
        for i, data in enumerate(matrix_list):
            for chromosome in chromosome_list[i]:
                matrix, positions, features = self.prepare_predict(data, pChromosome=chromosome)

                matrix.numpy_matrix = None
                matrix.hic_matrix = None
                matrix = None

                positions, is_boundary = self.classifier.predict(
                    positions, features)
                if positions is None or len(positions) == 0:
                    continue
                if conc_is_boundary is None:
                    conc_positions = positions
                    conc_is_boundary = is_boundary
                else:
                    conc_positions = np.concatenate((conc_positions, positions), axis=0)

                    conc_is_boundary = np.concatenate((conc_is_boundary, is_boundary), axis=0)

                features = None
                is_boundary = None
                positions = None
                # i = i + 1

            # matrix_name = ".".join(os.path.basename(matrix_list[i]).split(".")[:-1])
            out_file_i = self.out_file[i]
            TADClassifier.print_to_bed(TADClassifier.get_domains(
                conc_positions, conc_is_boundary), out_file_i)

    def print_to_bed(domain_df, path):
        '''print domain file'''

        domain_df.to_csv(path, sep='\t', header=None, index=False)

    def filter_domains(domains):
        '''filter out to small domains'''

        min_tad_size = 50000 - 1
        domain_df = pd.DataFrame({'Chrom': domains[:,
                                                   0],
                                  'Start': domains[:,
                                                   1],
                                  'End': domains[:,
                                                 2]})

        domain_df['Start'] = domain_df.Start.apply(np.int64)
        domain_df['End'] = domain_df.End.apply(np.int64)
        domain_df = domain_df.sort_values(by=["Chrom", "Start"])
        domain_df['Test'] = domain_df['Start'].shift(1, fill_value=0)
        domain_df['Chrom_Test'] = domain_df['Chrom'].shift(1, fill_value=0)
        domain_df['Chrom_Test_Bool'] = domain_df['Chrom'] != domain_df['Chrom_Test']
        domain_df['Test_Bool'] = domain_df['Test'] + min_tad_size <= domain_df['Start']
        domain_df['Test_Bool'] = np.logical_or(domain_df['Test_Bool'], domain_df['Chrom_Test_Bool'])
        domain_df = domain_df[domain_df['Test_Bool']]

        return domain_df[['Chrom', 'Start', 'End']].to_numpy()

    def get_domains(positions, y):
        '''return dataframe of predicted TADs'''

        pos_mask = y[:] == True
        # print('pos_mask {}'.format(pos_mask))
        domains = positions[pos_mask]
        # print('domains {}'.format(domains))

        domains = TADClassifier.filter_domains(domains)
        # print('domains {}'.format(domains))

        name = np.arange(1, domains.shape[0] + 1)
        score = np.full(domains.shape[0], -1)
        strand = np.full(domains.shape[0], 0)
        item_rgb = np.tile([0, 1], domains.shape[0])
        item_rgb = item_rgb[0:domains.shape[0]]
        # print('Chrom {}'.format(domains[:, 0]))

        domain_df = pd.DataFrame({'Chrom': domains[:,
                                                   0],
                                  'Start': domains[:,
                                                   1],
                                  'End': domains[:,
                                                 1],
                                  'Name': name,
                                  'Score': score,
                                  'Strand': strand,
                                  'ThickStart': domains[:,
                                                        1],
                                  'ThickEnd': domains[:,
                                                      1],
                                  'ItemRGB': item_rgb})
        # print('domain_df 1 {}'.format(domain_df))

        def rgb_helper(i):
            if(i == 0):
                return '31,120,180'
            else:
                return '51,160,44'

        def name_helper(i):
            return 'ID_0.03_' + str(i)

        def strand_helper(s):
            return '.'

        domain_df["ItemRGB"] = domain_df["ItemRGB"].map(rgb_helper)
        domain_df["Name"] = domain_df["Name"].map(name_helper)
        domain_df["Strand"] = domain_df["Strand"].map(strand_helper)

        domain_df['Start'] = domain_df.Start.apply(np.int64)
        domain_df['End'] = domain_df.End.apply(np.int64)

        domain_df = domain_df.sort_values(by=["Chrom", "Start"])

        domain_df['Start'] = domain_df['Start'].shift(1, fill_value=0)
        domain_df['ThickStart'] = domain_df['ThickStart'].shift(
            1, fill_value=0)
        # print('domain_df 2{}'.format(domain_df))

        # domain_df = domain_df[1:]
        # print('domain_df 3{}'.format(domain_df))

        domain_df = domain_df[domain_df['Start'] < domain_df['End']]

        return domain_df
