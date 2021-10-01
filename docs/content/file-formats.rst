.. _file-formats:

HiCExplorer file formats
========================

HiCExplorer has a native interaction file format, h5, and several native formats for capture Hi-C analysis.


Interaction matrix format h5
----------------------------

The h5 format is implemented as a hdf5 container and can be accessed either via HiCMatrix, the hdf5 API or a graphical user interface like hdfview.
The data in the hdf5 format is read and written via the PyTables interface.

The format is structured as follows:

.. code-block:: bash

   - intervals
      - chr_list (CArray)
      - end_list (CArray)
      - extra_list (CArray)
      - start_list (CArray)
   - matrix
      - data (CArray)
      - indices (CArray)
      - indptr (CArray)
      - shape (tupel)
    - distance counts (CArray) [optional]
    - nan_bins (CArray) [optional]
    - correction_factors (CArray) [optional]

The group 'matrix' contains the data which is necessary to create a scipy sparse csr matrix: 

.. code-block:: python

    import tables
    with tables.open_file(self.matrixFileName, 'r') as f:
            parts = {}
            try:
                for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                    parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
            except Exception as e:
                log.info('No h5 file. Please check parameters concerning the file type!')
                e
            matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                                shape=parts['shape'])

The group 'intervals' contains the information of the genomic regions associated to a bin in the matrix via its index. To retrieve the data run:

.. code-block:: python

    with tables.open_file(self.matrixFileName, 'r') as f:
        intvals = {}
        for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
            if toString(interval_part) == toString('chr_list'):
                chrom_list = getattr(f.root.intervals, interval_part).read()
                intvals[interval_part] = toString(chrom_list)
            else:
                intvals[interval_part] = getattr(f.root.intervals, interval_part).read()

        cut_intervals = list(zip(intvals['chr_list'], intvals['start_list'], intvals['end_list'], intvals['extra_list']))

To write the format, the above structure needs to be created via PyTables. For example:


.. code-block:: python

    with tables.open_file(filename, mode="w", title="HiCExplorer matrix") as h5file:
        matrix_group = h5file.create_group("/", "matrix", )
        for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                arr = np.array(getattr(matrix, matrix_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(matrix_group, matrix_part, atom,
                                          shape=arr.shape,
                                          filters=filters)

The matrix object is a scipy `csr_matrix`. 

Please refer to HiCMatrix for a reference implementation: https://github.com/deeptools/HiCMatrix/blob/master/hicmatrix/lib/h5.py


To open a cool / h5 file with the HiCMatrix library, use the following code:

.. code-block:: python

   import hicmatrix.HiCMatrix as hm

   hic_ma = hm.hiCMatrix("/path/to/matrix")

   # csr_matrix
   hic_ma.matrix

   # the corresponding genomic regions as a list
   hic_ma.cut_intervals


To write a h5 use the following code. Please consider, once the datatype of a matrix is specified, for example by reading a h5 matrix or by first time writing a h5 matrix, it cannot be changed to a cool file anymore (except by manipulating the internal matrixFilehandler object). 

.. code-block:: python

   import hicmatrix.HiCMatrix as hm
   from scipy.sparse import csr_matrix


   # create a HiCMatrix object
   hic_ma = hm.hiCMatrix()

   # set important data structures

   hic_ma.matrix = csr_matrix()
   hic_ma.cut_intervals = list[(chr1, start1, end1, 1), (chr1, start2, end2, 1), ..., (chr1, startN, endN, 1)]

   # to store a h5 matrix
   hic_ma.save('/path/to/storage/matrix.h5')
   # to store a cool matrix
   hic_ma.save('/path/to/storage/matrix.cool')

Alternatively, the matrixFileHandler object can be accessed directly:


.. code-block:: python

   from hicmatrix.lib import MatrixFileHandler


   # Load a matrix via the MatrixFileHandler class

   matrixFileHandlerInput = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=matrix,
                                                           pCorrectionFactorTable=args.correction_name,
                                                           pCorrectionOperator=correction_operator,
                                                           pChrnameList=chromosomes_to_load,
                                                           pEnforceInteger=args.enforce_integer,
                                                           pApplyCorrectionCoolerLoad=applyCorrectionCoolerLoad)

   # Load data
   _matrix, cut_intervals, nan_bins, \
      distance_counts, correction_factors = matrixFileHandlerInput.load()

   # create matrixFileHandler output object
   matrixFileHandlerOutput = MatrixFileHandler(pFileType=args.outputFormat, pEnforceInteger=args.enforce_integer, pFileWasH5=format_was_h5, pHic2CoolVersion=hic2CoolVersion, pHiCInfo=cool_metadata)

   # set the variables
   matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                             correction_factors, distance_counts)

   matrixFileHandlerOutput.save(args.outFileName, pSymmetric=True, pApplyCorrection=applyCorrection)



Capture Hi-C HDF containers
---------------------------

The capture Hi-C data analysis creates for the scripts `chicViewpoint`, `chicSignificantInteractions`, `chicAggregateStatistic` and `chicDifferentialTest` individual HDF containers to store the processed data.

chicViewpoint
~~~~~~~~~~~~~

.. code-block:: bash

   - averageContactBin (int)
   - fixateRange (int)
   - range (int, int)
   - resolution (int)
   - type='interactions' (string)
   - matrix 1 
    - chromosome 1
       - gene name 1
           - chromosome (String)
           - end_list (array)
           - gene (String)
           - interaction_data_list (array)
           - pvalue (array)
           - raw (array)
           - reference_point_end (int)
           - reference_point_start (int)
           - relative_position_list (array)
           - start_list (array)
           - sum_of_interactions (float)
           - xfold (array)
       - gene name 2
           - ...
           - ...
        - ...
    - chromosome 2
       - gene 1
          - ...
          - ...
       - ...
    - ...
    - ...
    - genes
      - gene name 1 (link to matrix 1 / chromosome 1 / gene name 1
      - ...
   - matrix 1
     - chromosome 1
        - gene name 1
            - chromosome (String)
            - end_list (array)
            - gene (String)
            - interaction_data_list (array)
            - pvalue (array)
            - raw (array)
            - reference_point_end (int)
            - reference_point_start (int)
            - relative_position_list (array)
            - start_list (array)
            - sum_of_interactions (float)
            - xfold (array)
        - gene name 2
            - ...
            - ...
        - ...
     - chromosome 2
        - gene 1
           - ...
           - ...
     - ...
     - ...

    - genes
      - gene name 1 (link to matrix 2 / chromosome 1 / gene name 1
      - ...


chicSignificantInteractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

`chicSignificantInteractions` creates two files: a target file and a file containing the significant interactions:


Depending on the combination mode (single / dual) the structure is slightly different for the target file:

.. code-block:: bash

   - combinationMode = 'single' (String)
   - fixateRange (int)
   - mode_preselection (String)
   - mode_preselection_value (float)
   - peakInteractionsThreshold (float)
   - pvalue (float)
   - range (int, int)
   - truncateZeroPvalues (Boolean)
   - type='target' (String)
   - matrix 1 
    - chromosome 1
       - gene name 1
           - chromosome (String)
           - end_list (array)
           - reference_point_end (int)
           - reference_point_start (int)
           - start_list (array)
       - gene name 2
           - ...
           - ...
        - ...
    - chromosome 2
       - gene 1
          - ...
          - ...
    - ...
    - ...
    - genes
      - gene name 1 (link to matrix 1 / chromosome 1 / gene name 1
      - ...
   - matrix 1
     - chromosome 1
        - gene name 1
           - chromosome (String)
           - end_list (array)
           - reference_point_end (int)
           - reference_point_start (int)
           - start_list (array)
        - gene name 2
            - ...
            - ...
         - ...
     - chromosome 2
        - gene 1
           - ...
           - ...
     - ...
     - ...

    - genes
      - gene name 1 (link to matrix 2 / chromosome 1 / gene name 1
      - ...

Combination mode 'dual' combines the target regions of one viewpoint region of two matrices.

.. code-block:: bash

   - combinationMode = 'dual' (String)
   - fixateRange (int)
   - mode_preselection (String)
   - mode_preselection_value (float)
   - peakInteractionsThreshold (float)
   - pvalue (float)
   - range (int, int)
   - truncateZeroPvalues (Boolean)
   - type='target' (String)
   - matrix 1
      - matrix 2
         - chromosome 1
            - gene name 1
               - chromosome (String)
               - end_list (array)
               - reference_point_end (int)
               - reference_point_start (int)
               - start_list (array)
            - gene name 2
               - ...
               - ...
            - ...
         - chromosome 2
            - gene 1
               - ...
               - ...
         - genes
            - gene name 1 (link to matrix 1 / matrix 2 / chromosome 1 / gene name 1
            - ...
         - matrix 3
            - chromosome 1
               - gene name 1 
                  - ...
            - ...
            - genes
               - gene name 1 (link to matrix 1 / matrix 3 / chromosome 1 / gene name 1
               - ...
   - matrix 2
      - matrix 3
         - chromosome 1
            - gene name 1
               - ...
            - ...
         - ...
      - ...
         - genes
            - gene name 1 (link to matrix 2 / matrix 3 / chromosome 1 / gene name 1
            - ...
   - ...



Significant interactions file in the 'single' and 'dual' mode don't have a difference in their structure:


.. code-block:: bash

   - combinationMode  = 'single' (String) / 'dual' (String)
   - fixateRange (int)
   - mode_preselection (String)
   - mode_preselection_value (float)
   - peakInteractionsThreshold (float)
   - pvalue (float)
   - range (int, int)
   - truncateZeroPvalues (Boolean)
   - type='significant' (String)
   - matrix 1 
    - chromosome 1
       - gene name 1
            - chromosome (String)
            - end_list (array)
            - gene (String)
            - interaction_data_list (array)
            - pvalue (array)
            - raw (array)
            - reference_point_end (int)
            - reference_point_start (int)
            - relative_position_list (array)
            - start_list (array)
            - sum_of_interactions (float)
            - xfold (array)
       - gene name 2
           - ...
           - ...
        . ---
    - chromosome 2
       - gene 1
          - ...
          - ...
    - ...
    - ...
    - genes
      - gene name 1 (link to matrix 1 / chromosome 1 / gene name 1
      - ...
   - matrix 1
     - chromosome 1
        - gene name 1
            - chromosome (String)
            - end_list (array)
            - gene (String)
            - interaction_data_list (array)
            - pvalue (array)
            - raw (array)
            - reference_point_end (int)
            - reference_point_start (int)
            - relative_position_list (array)
            - start_list (array)
            - sum_of_interactions (float)
            - xfold (array)
        - gene name 2
            - ...
            - ...
         . ---
     - chromosome 2
        - gene 1
           - ...
           - ...
     - ...
     - ...

    - genes
      - gene name 1 (link to matrix 2 / chromosome 1 / gene name 1
      - ...



chicAggregateStatistic
~~~~~~~~~~~~~~~~~~~~~~


.. code-block:: bash

   - type='aggregate' (String)
   - matrix_1_matrix_2 
      - matrix 1
         - chromosome 1
            - gene name 1
                  - chromosome (String)
                  - end_list (array)
                  - gene_name (String)
                  - interaction_data_list (array)
                  - raw_target_list (array)
                  - relative_distance_list (array)
                  - start_list (array)
                  - sum_of_interactions (float)
            - gene name 2
           - ...
           - ...
         - genes
            - gene name 1 (link to matrix_1_matrix_2 / matrix 1 / chromosome 1 / gene name 1
            - ...
      - matrix 2
         - chromosome 1
               - gene name 1
                     - chromosome (String)
                     - end_list (array)
                     - gene_name (String)
                     - interaction_data_list (array)
                     - raw_target_list (array)
                     - relative_distance_list (array)
                     - start_list (array)
                     - sum_of_interactions (float)
               - gene name 2
            - ...
            - ...
         - genes
            - gene name 1 (link to matrix_1_matrix_2/ matrix 2 / chromosome 1 / gene name 1
            - ...
   - matrix_2_matrix_2
      - matrix 2
         - ...
      - matrix 3
         - ...



chicDifferentialTest
~~~~~~~~~~~~~~~~~~~~


.. code-block:: bash

   - type='differential' (String)
   - alpha (float)
   - test: fisher / chi2 (String)
   - matrix_1
      - matrix 2
         - chromosome 1
            - gene name 1
                  - accepted
                     - chromosome (String)
                     - end_list (array)
                     - gene (String)
                     - interaction_data_list (array)
                     - pvalue (array)
                     - raw_target_list_1 (array)
                     - raw_target_list_2 (array)
                     - relative_distance_list (array)
                     - start_list (array)
                     - sum_of_interactions_1 (float)
                     - sum_of_interactions_2 (float)
                  - all
                     - chromosome (String)
                     - end_list (array)
                     - gene (String)
                     - interaction_data_list (array)
                     - pvalue (array)
                     - raw_target_list_1 (array)
                     - raw_target_list_2 (array)
                     - relative_distance_list (array)
                     - start_list (array)
                     - sum_of_interactions_1 (float)
                     - sum_of_interactions_2 (float)
                  - rejected 
                     - chromosome (String)
                     - end_list (array)
                     - gene (String)
                     - interaction_data_list (array)
                     - pvalue (array)
                     - raw_target_list_1 (array)
                     - raw_target_list_2 (array)
                     - relative_distance_list (array)
                     - start_list (array)
                     - sum_of_interactions_1 (float)
                     - sum_of_interactions_2 (float)
            - gene name 2
                  - accepted
                  - all
                  - rejected
            - ...
         - genes
            - gene name 1 (link to matrix_1 / matrix 2 / chromosome 1 / gene name 1
            - ...
      - matrix 3
         - chromosome 1
               - gene name 1
                     - accepted
                     - all
                     - rejected
               - gene name 2
            - ...
            - ...
         - genes
            - gene name 1 (link to matrix_1 / matrix 3 / chromosome 1 / gene name 1
            - ...
      - ...
   - matrix 2
      - matrix 3
         - ...

