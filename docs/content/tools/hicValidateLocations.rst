.. _hicValidateLocations:

hicValidateLocations
=====================

hicValidateLoops is a tool to compare the detect loops from hicDetectLoops (or from any other software as long as the data format is followed, see below) 
with known peak protein locations to validate if the computed loops do have the expected anchor points. Loops in mammals are usually bound by CTCF or Cohesin, 
therefore it is important to know if the detect loops have protein peaks at their X and Y position.

.. figure:: ../../images/loops_bonev_cavalli.png

    Loops in Hi-C, graphic from Bonev & Cavalli, Nature Reviews Genetics 2016


Data format
===========

The data format of hicDetectLoops output is:

chr_x start_x end_x chr_y start_y end_y p-value

As ``--protein`` the input of narrowPeak or broadPeak files are accepted. However, as long as the ``--protein`` input file contains chromosome, start and end in the first three columns, it should work.

.. argparse::
   :ref: hicexplorer.hicValidateLocations.parse_arguments
   :prog: hicValidateLocations

