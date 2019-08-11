.. _hicQuickQC:

hicQuickQC
==========

Background
^^^^^^^^^^

This tool considers the first 1,000,000 reads (or user defined number) of the mapped bam files to get a quality estimate of the Hi-C data.

Description
^^^^^^^^^^^

.. argparse::
   :ref: hicexplorer.hicQuickQC.parse_arguments
   :prog: hicQuickQC

For more information see `encode <https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf>`__ .