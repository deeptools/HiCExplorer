Installation
=============

.. contents::
    :local:

Requirements
-------------

* python >= 3.6
* numpy >= 1.19.*
* scipy >= 1.5.*
* matplotlib-base >= 3.1.*
* ipykernel >= 5.3.0
* pysam >= 0.16.*
* intervaltree >= 3.1.*
* biopython
* pytables >= 3.6.*
* pandas >= 1.1.*
* pybigwig >= 0.3.*
* jinja2 >= 2.11
* unidecode >= 1.1.*
* hicmatrix >= 15
* hic2cool >= 0.8.3
* psutil >= 5.7.*
* pygenometracks >= 3.5
* fit_nbinom >= 1.1
* cooler >= 0.8.10
* krbalancing >= 0.0.5 (Needs the library eigen; openmp is recommended for linux users. No openmp support on macOS.)
* pybedtools >= 0.8.*
* future >= 0.18
* tqdm >= 4.50
* hyperopt >= 0.2.4
* python-graphviz >= 0.14


**Warning:** Python 2.7 support is discontinued. Moreover, the support for pip is discontinued too. 
**Warning:** We strongly recommend to use the conda package manager and will no longer give support on all issues raising with pip.

Command line installation using ``conda``
-----------------------------------------

The fastet way to obtain **Python 3.6 or 3.7 together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation. All of the requirements for HiCExplorer can be installed in Anaconda with:

.. code:: bash

    $ conda install hicexplorer -c bioconda -c conda-forge

We strongly recommended to use conda to install HiCExplorer. 

Command line installation using ``pip``
-----------------------------------------

The installation via pip is discontinued with version 3.0. The reason for this is that we want to provide a 'one-click' installation. However,
with version 3.0 we added the C++ library eigen as dependency and pip does not support non-Python packages. 

For older versions you can still use pip: 
Install HiCExplorer using the following command:
::

	$ pip install hicexplorer

All python requirements should be automatically installed.

If you need to specify a specific path for the installation of the tools, make use of `pip install`'s numerous options:

.. code:: bash

    $ pip install --install-option="--prefix=/MyPath/Tools/hicexplorer" git+https://github.com/deeptools/HiCExplorer.git

**Warning:** It can be that you have to install additional packages via your system package manager to successfully install HiCExplorer via pip.
**Warning:** We strongly recommend to use the conda package manager and will no longer give support on all issues raising with pip.


Command line installation without ``pip``
-------------------------------------------

You are highly recommended to use `pip` rather than these more complicated steps.

1. Install the requirements listed above in the "requirements" section. This is done automatically by `pip`.

2. Download source code
::

	$ git clone https://github.com/deeptools/HiCExplorer.git

or if you want a particular release, choose one from https://github.com/deeptools/HiCExplorer/releases:
::

	$ wget https://github.com/deeptools/HiCExplorer/archive/1.5.12.tar.gz
	$ tar -xzvf

3. To install the source code (if you don't have root permission, you can set
a specific folder using the ``--prefix`` option)
::

	$ python setup.py install --prefix /User/Tools/hicexplorer




Galaxy installation
--------------------

HiCExplorer can be easily integrated into a local `Galaxy <http://galaxyproject.org>`_, the wrappers are provided at the `Galaxy tool shed <https://toolshed.g2.bx.psu.edu/>`_.

Installation with Docker
^^^^^^^^^^^^^^^^^^^^^^^^

The HiCExplorer Galaxy instance is also available as a docker container, for those wishing to use the Galaxy
framework but who also prefer a virtualized solution. This container is quite simple to install:

::

    $ sudo docker pull quay.io/bgruening/galaxy-hicexplorer

To start and otherwise modify this container, please see the instructions on `the docker-galaxy-stable github repository <https://github.com/bgruening/docker-galaxy-stable>`__. Note that you must use `bgruening/galaxy-hicexplorer` in place of `bgruening/galaxy-stable` in the examples, as the HiCExplorer Galaxy container is built on top of the galaxy-stable container.

.. tip:: For support, or feature requests contact: deeptools@googlegroups.com
