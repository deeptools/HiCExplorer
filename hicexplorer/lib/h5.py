from .matrixFile import MatrixFile
import tables
from scipy.sparse import csr_matrix, triu
import numpy as np
from hicexplorer.utilities import toString
from builtins import super

from past.builtins import zip
from os import unlink
import os
import logging
log = logging.getLogger(__name__)


class H5(MatrixFile, object):

    def __init__(self, pMatrixFile):
        super().__init__(pMatrixFile)

    def load(self):
        """
        Loads a matrix stored in h5 format
        :param matrix_filename:
        :return: matrix, cut_intervals, nan_bins, distance_counts, correction_factors
        """
        log.debug('Load in h5 format')

        with tables.open_file(self.matrixFileName) as f:
            parts = {}
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()

            matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                                shape=parts['shape'])
            # matrix = hiCMatrix.fillLowerTriangle(matrix)
            # get intervals
            intvals = {}
            for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
                if toString(interval_part) == toString('chr_list'):
                    chrom_list = getattr(f.root.intervals, interval_part).read()
                    intvals[interval_part] = toString(chrom_list)
                else:
                    intvals[interval_part] = getattr(f.root.intervals, interval_part).read()

            cut_intervals = zip(intvals['chr_list'], intvals['start_list'], intvals['end_list'], intvals['extra_list'])
            assert len(cut_intervals) == matrix.shape[0], \
                "Error loading matrix. Length of bin intervals ({}) is different than the " \
                "size of the matrix ({})".format(len(cut_intervals), matrix.shape[0])

            # get nan_bins
            if hasattr(f.root, 'nan_bins'):
                nan_bins = f.root.nan_bins.read()
            else:
                nan_bins = np.array([])

            # get correction factors
            if hasattr(f.root, 'correction_factors'):
                correction_factors = f.root.correction_factors.read()
                assert len(correction_factors) == matrix.shape[0], \
                    "Error loading matrix. Length of correction factors does not" \
                    "match size of matrix"
            else:
                correction_factors = None

            # get correction factors
            if hasattr(f.root, 'distance_counts'):
                distance_counts = f.root.correction_factors.read()
            else:
                distance_counts = None

            return matrix, cut_intervals, nan_bins, distance_counts, correction_factors

    def save(self, filename, pSymmetric=True, pApplyCorrection=None):
        """
        Saves a matrix using hdf5 format
        :param filename:
        :return: None
        """
        log.debug('Save in h5 format')

        # self.restoreMaskedBins()
        if not filename.endswith(".h5"):
            filename += ".h5"

        # if the file name already exists
        # try to find a new suitable name
        if os.path.isfile(filename):
            log.warning("*WARNING* File already exists {}\n "
                        "Overwriting ...\n".format(filename))

            unlink(filename)
        if self.nan_bins is None:
            self.nan_bins = np.array([])
        elif not isinstance(self.nan_bins, np.ndarray):
            self.nan_bins = np.array(self.nan_bins)

        # save only the upper triangle of the
        if pSymmetric:
            # symmetric matrix
            matrix = triu(self.matrix, k=0, format='csr')
        else:
            matrix = self.matrix
        filters = tables.Filters(complevel=5, complib='blosc')
        with tables.open_file(filename, mode="w", title="HiCExplorer matrix") as h5file:
            matrix_group = h5file.create_group("/", "matrix", )
            # save the parts of the csr matrix
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                arr = np.array(getattr(matrix, matrix_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(matrix_group, matrix_part, atom,
                                          shape=arr.shape,
                                          filters=filters)
                ds[:] = arr

            # save the matrix intervals
            intervals_group = h5file.create_group("/", "intervals", )
            chr_list, start_list, end_list, extra_list = zip(*self.cut_intervals)
            for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
                arr = np.array(eval(interval_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(intervals_group, interval_part, atom,
                                          shape=arr.shape,
                                          filters=filters)
                ds[:] = arr

            # save nan bins
            if len(self.nan_bins):
                atom = tables.Atom.from_dtype(self.nan_bins.dtype)
                ds = h5file.create_carray(h5file.root, 'nan_bins', atom,
                                          shape=self.nan_bins.shape,
                                          filters=filters)
                ds[:] = self.nan_bins

            # save corrections factors
            if self.correction_factors is not None and len(self.correction_factors):
                self.correction_factors = np.array(self.correction_factors)
                atom = tables.Atom.from_dtype(self.correction_factors.dtype)
                ds = h5file.create_carray(h5file.root, 'correction_factors', atom,
                                          shape=self.correction_factors.shape,
                                          filters=filters)
                ds[:] = np.array(self.correction_factors)

            # save distance counts
            if self.distance_counts is not None and len(self.distance_counts):
                atom = tables.Atom.from_dtype(self.distance_counts.dtype)
                ds = h5file.create_carray(h5file.root, 'distance_counts', atom,
                                          shape=self.distance_counts.shape,
                                          filters=filters)
                ds[:] = np.array(self.distance_counts)
