.. _findRestSite:

findRestSites
=============

.. argparse::
   :ref: hicexplorer.findRestSite.parse_arguments
   :prog: findRestSite

Further usage
^^^^^^^^^^^^^

If multiple restriction enzymes are used in one experiment, ``findRestSite`` can be used to find restirction sites individually per enzyme. Afterwards, all output bed files should be combined. However, it should be noted that the QC report will not be correct for this specific usage.
