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
           - reference_point_end (int)
           - reference_point_start (int)
           - start_list (array)
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

Significant interactions file:


.. code-block:: bash

   - combinationMode (String)
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

chicSignificantInteractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~