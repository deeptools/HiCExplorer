.. _hicInfo:

hicInfo
=======

.. argparse::
   :ref: hicexplorer.hicInfo.parse_arguments
   :prog: hicInfo


An example output looks like this:

.. code-block:: INI

   File:	hic_matrix.h5
   WARNING: bin size is not homogeneous. Median 434
   Size:	339,684
   Sum:	243,779,407
   Bin length:	434
   Chromosomes:	2L, 2R, 3L, 3R, 4, X
   Non-zero elements:	229,122,848
   Minimum:	0
   Maximum:	2116
   NaN bins:	0

With HiCExplorer version 3 we support metadata in cooler files. The default behavior is to use the metadata:

.. code:: bash

    $ hicInfo -m matrix.cool

To use the old method (and default for h5) please add the parameter:

.. code:: bash

    $ hicInfo -m matrix.cool --no_metadata