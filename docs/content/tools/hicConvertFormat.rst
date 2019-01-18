.. _hicConvertFormat:

hicConvertFormat
================

.. argparse::
   :ref: hicexplorer.hicConvertFormat.parse_arguments
   :prog: hicConvertFormat

Background
^^^^^^^^^^

To reproduce analysis and to compare and use different Hi-C analysis software the exchange of the interaction matrix is crucial.
However, the most Hi-C software is supporting only their own data format which makes the exchange difficult. HiCExplorer supports a range of 
interaction matrices, either to import or to export. 

Import:
    - hic
    - cool
    - h5
    - homer
    - HicPro

Export:
    - cool
    - mcool
    - h5
    - homer
    - ginteractions

With HiCExplorer version 2.2 hicConvertFormat and hicAdjustMatrix replace hicExport from HiCExplorer 2.1 and older versions.


Usage example
^^^^^^^^^^^^^

hic2cool
""""""""

HiCExplorer uses the library hic2cool_  to convert **.hic** interaction matrix files to the cool format. Usually .hic files 
have the three correction factors **KR**, **VC** or **VC_SQRT**; however these can not be applied nativly by HiCExplorer tools because 
HiCExplorer expects the correction values to be stored in the column **weight**.
To work with corrected data the correction factors need to applied separatly, see section cool to cool.

.. _hic2cool: https://github.com/4dn-dcic/hic2cool

The following example will convert a hic file which contains the resolution of 1000 to a cool file with 10kb resolution. The desired 
resolution needs to be existing in the hic file. If no resolution parameter is defined a mcool file with all available resolutions is created.

.. code:: bash

    $ hicConvertFormat -m matrix.hic --inputFormat hic --outputFormat cool -o matrix.cool --resolutions 10000

It is only possible to convert from hic to cool format, no other formats are supported.

cool to cool
""""""""""""

The cool file format is developed and maintained by the Mirny_ lab and allows to access interaction matrices in a easy to use data format.

.. _Mirny: https://github.com/mirnylab/cooler


Cool data format allows to use the following options:

- correction_name: In the case correction factors are not stored in 'weight' the correct column name can be defined with this parameter and the resulting matrix will store the values in 'weight'.
- correction_division: Correction factors can be applied by a multiplication or a division. The default behaviour is to use the multiplication, in the case the correction factors are inverted, set this parameter.
- store_applied_correction: Set this parameter if correction factors should be applied on the data and should be written back to colum 'counts' in the corrected form and not as raw. Default: not set.
- chromosomes: Define a list of chromosomes which should be included in the output matrix. All chromosomes which are not defined are not part of the new matrix. This parameter can speed up the processing especiallly if only one chromosome is used.
- enforce_integer: Raw interaction data is stored as integers, after the correction is applied the data is a float. Set a this parameter to enforce integer values in the new matrix.
- load_raw_values: Set this parameter if the interaction data should not be loaded with the correction factors.

Example usage
-------------

.. code:: bash

    $ hicConvertFormat -m matrix.cool --inputFormat cool --outputFormat cool -o matrix.cool --correction_name KR

Homer
"""""

Homer_ is a software for 'motif discovery and next generation sequencing analysis' and supports Hi-C. HiCExplorer is able to read and write the used 
interaction matrix of Homer. Homer stores the interaction matrix in a simple text file as a dense matrix, to write 
large matrices in Homer format needs a lot of space and can take a few ours to days. 

.. _Homer: http://homer.ucsd.edu/homer/index.html

Example usage
-------------

.. code:: bash

    $ hicConvertFormat -m matrix.homer --inputFormat homer --outputFormat cool -o matrix.cool

Hic-Pro
"""""""

HiC-Pro_ file format needs an additional bed file as input:

Example usage
-------------

.. code:: bash

    $ hicConvertFormat -m matrix.hicpro --bedFileHicpro hicpro.bed --inputFormat hicpro --outputFormat cool -o matrix.cool

.. _HiC-Pro: https://github.com/nservant/HiC-Pro

Create a mcool file
"""""""""""""""""""

With HiCExplorer it is possible to create a multiple cool (mcool) file. These mcool files can be used e.g. with HiGlass_.

.. _HiGlass: http://higlass.io/

To create a mcool file use as input either one matrix in one of the supported read formats and define the desired resolutions or define
multiple input matrices. In the second case the matrices should all have different resolutions.

Example usage
-------------

The resolutions need to be a multiple of the input matrix i.e. matrix with 10kb than 20kb and 30kb are possible but not 35kb.

.. code:: bash

    $ hicConvertFormat -m matrix.cool --inputFormat cool --outputFormat mcool
       -o multi_matrix.mcool --resolutions 20000 40000 70000 120000 500000

.. code:: bash

    $ hicConvertFormat -m matrix10kb.cool matrix20kb.cool matrix30kb.cool 
        --inputFormat cool --outputFormat mcool -o multi_matrix.mcool

The mcool matrix contains the individual matrices as follows:


.. code:: 

    multi_matrix.mcool::/resolutions/10000
    multi_matrix.mcool::/resolutions/40000
    multi_matrix.mcool::/resolutions/70000
    multi_matrix.mcool::/resolutions/120000
    multi_matrix.mcool::/resolutions/500000


