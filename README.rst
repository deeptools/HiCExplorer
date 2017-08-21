.. image:: https://travis-ci.org/maxplanck-ie/HiCExplorer.svg?branch=master
   :target: https://travis-ci.org/maxplanck-ie/HiCExplorer
.. image:: https://readthedocs.org/projects/hicexplorer/badge/?version=latest
   :target: http://hicexplorer.readthedocs.io/?badge=latest
.. image:: https://anaconda.org/bioconda/hicexplorer/badges/installer/conda.svg
   :target: https://anaconda.org/bioconda/hicexplorer
.. image:: https://quay.io/repository/biocontainers/hicexplorer/status
   :target: https://quay.io/repository/biocontainers/hicexplorer
.. image:: https://badge.fury.io/py/HiCExplorer.svg
       :target: https://badge.fury.io/py/HiCExplorer

HiCExplorer
===========

Set of programs to process, analyze and visualize Hi-C data
-----------------------------------------------------------

Sequencing techniques that probe the 3D organization of the genome generate large amounts of data whose processing,
analysis and visualization is challenging. Here, we present HiCExplorer, a set of tools for the analysis and
visualization of chromosome conformation data. HiCExplorer facilitates the creation of contact matrices, correction
of contacts, TAD detection, merging, reordering or chromosomes, conversion from different formats and detection of
long-range contacts. Moreover, it allows the visualization of multiple contact matrices along with other types of
data like genes, compartments, ChIP-seq coverage tracks (and in general any type of genomic scores) and long range contacts.


Citation:
^^^^^^^^^

Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening, Kin Chung Lam, Bianca Habermann, Asifa Akhtar, Thomas Manke
**"High-resolution TADs reveal DNA sequences underlying genome organization in flies"**. bioRxiv 115063; doi: https://doi.org/10.1101/115063

.. image:: ./docs/images/hicex2.png

Installation
^^^^^^^^^^^^

HiCExplorer is available for:

-  Command line usage (via pip/anaconda/github)
-  Integration into Galaxy servers (via toolshed/API/web-browser)

There are many easy ways to install HiCExplorer. Details can be found
`here <https://hicexplorer.readthedocs.io/en/latest/content/installation.html>`__


Install with conda
++++++++++++++++++

The easiest way to install HiCExplorer is using `BioConda <http://bioconda.github.io/>`_
::

   $ conda install hicexplorer -c bioconda -c conda-forge



Install with pip
++++++++++++++++
::

   $ pip install HiCExplorer

Install by cloning this repository
++++++++++++++++++++++++++++++++++

You can install any one of the HiCExplorer branches on command line
(linux/mac) by cloning this git repository :

::

    $ git clone https://github.com/maxplanck-ie/HiCExplorer.git
    $ cd HiCExplorer
    $ python setup.py install

If you don't have root permission, you can set a specific folder using the ``--prefix`` option

::

	$ python setup.py install --prefix /User/Tools/hicexplorer


Documentation:
^^^^^^^^^^^^^^
Please visit our complete documentation `Here <http://hicexplorer.readthedocs.org/>`_
