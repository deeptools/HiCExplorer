Installation
=============

.. contents::
    :local:

Requirements
-------------

* Python 2.7
* numpy >= 1.8.1
* scipy >= 0.14.0
* pysam >= 0.8.3
* matplotlib >= 1.3.1
* bx-python >= 0.7.1
* biopython >= 1.65

The fastet way to obtain **Python 2.7 together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation. All of the requirements for HiCExplorer can be installed in Anaconda with:

.. code:: bash

    $ conda install -c bioconda hicexplorer=0.1

Command line installation using ``pip``
-----------------------------------------

Install HiCExplorer using the following command:
::

	$ pip install hicexplorer

All python requirements should be automatically installed.

If you need to specify a specific path for the installation of the tools, make use of `pip install`'s numerous options:

.. code:: bash

    $ pip install --install-option="--prefix=/MyPath/Tools/hicexplorer" git+https://github.com/maxplanck-ie/HiCExplorer.git


Command line installation without ``pip``
-------------------------------------------

You are highly recommended to use `pip` rather than these more complicated steps.

1. Install the requirements listed above in the "requirements" section. This is done automatically by `pip`.

2. Download source code
::

	$ git clone https://github.com/maxplanck-ie/HiCExplorer.git

3. To install the source code (if you don't have root permission, you can set
a specific folder using the ``--prefix`` option)
::

	$ python setup.py install --prefix /User/Tools/hicexplorer
