HiCExplorer
===========

HiCExplorer
-----------

***Set of programms to process, analyze and visualize Hi-C data.***

HiCexplorer addresses the common tasks of Hi-C analysis from processing
to visualization.

Usage of HiCExplorer
--------------------

Tools
~~~~~

The following is the list of tools available in HiCExplorer. Click on
each tool name to get more information on what it does, and how to use
it.

+----------------------------------------+----------------------------------------------------------------------------------------------+
| TOOL                                   | DESCRIPTION                                                                                  |
+========================================+==============================================================================================+
| `***hicBuildMatrix*** <>`__            | Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads               |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicCorrectMatrix*** <>`__          | Uses iterative correction to remove biases from a Hi-C matrix                                |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicFindEnrichedContacts*** <>`__   | Identifies enriched Hi-C contacts                                                            |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicCorrelate*** <>`__              | Computes and visualises the correlation of Hi-C matrices                                     |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicFindTADs*** <>`__               | Identifies Topologically Associating Domains (TADs)                                          |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicMergeMatrixBins*** <>`__        | Merges consecutives bins on a Hi-C matrix to reduce resolution                               |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicPlotMatrix*** <>`__             | Plots a Hi-C matrix as a heatmap                                                             |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicPlotTADs*** <>`__               | Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)   |
+----------------------------------------+----------------------------------------------------------------------------------------------+
| `***hicSumMatrices*** <>`__            | Adds Hi-C matrices of the same size                                                          |
+----------------------------------------+----------------------------------------------------------------------------------------------+

Analysis workflow
~~~~~~~~~~~~~~~~~

`An Image here <>`__ Showing the analysis workflow using the tools

Examples
~~~~~~~~

-  `Hi-C analysis of mouse ESCs using HiC-Explorer <>`__.
-  `Hi-C analysis of human ESCs using HiC-Explorer <>`__.

--------------

Installation
------------

-  `Quick start <#quick>`__
-  `General Installation <#general>`__
-  `Installation on a Mac <#mac>`__

Quick install
~~~~~~~~~~~~~

The easiest way to install the latest HiCExplorer version is to use
``pip`` with our github repository.

::

    $ pip install git+ssh://git@github.com/maxplanck-ie/HiCExplorer.git --user

General Installation
~~~~~~~~~~~~~~~~~~~~

::

    $ git clone https://github.com/maxplanck-ie/HiCExplorer.git
    $ cd HiCExplorer
    $ python setup.py install

By default, the script will install the python library and executable
codes globally, which means you need to be root or administrator of the
machine to complete the installation. If you need to provide a
nonstandard install prefix, or any other nonstandard options, you can
provide many command line options to the install script.

::

    $ python setup.py --help

For example, to install under a specific location use:

::

    $ python setup.py install --prefix <target directory>

Installation on a MAC
~~~~~~~~~~~~~~~~~~~~~

The easiest way to get numpy and scipy dependencies is to install the
`Anaconda Scientific Python
Distribution <https://www.continuum.io/downloads>`__. After
installation, open a terminal ("Applications" → "Terminal") and follow
the `General Installation <#general>`__

If individual installation of the dependencies is preferred, follow
those steps:

Requirement: Python 2.7 installed

Download the packages and install them using dmg images: +
http://sourceforge.net/projects/numpy/files/NumPy/ +
http://sourceforge.net/projects/scipy/files/scipy/

Then open terminal ("Applications" → "Terminal") and follow the `General
Installation <#general>`__

This tool suite is developed by the `Bioinformatics
Facility <http://www1.ie-freiburg.mpg.de/bioinformaticsfac>`__ at the
`Max Planck Institute for Immunobiology and Epigenetics,
Freiburg <http://www1.ie-freiburg.mpg.de/>`__.
