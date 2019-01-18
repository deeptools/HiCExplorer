Installation
=============

.. contents::
    :local:

Requirements
-------------

* Python 2.7 or 3.6
* numpy >= 1.15
* scipy >= 1.1
* matplotlib >= 2.2
* pysam >= 0.14
* intervaltree >= 2.1
* biopython >= 1.72
* pytables >= 3.4
* pyBigWig >= 0.3
* future >= 0.17
* cooler >= 0.7.11
* six >= 1.11
* jinja2 >= 2.10
* pandas >= 0.23
* unidecode >= 1.0
* hicmatrix >= 5
* pygenometracks >= 2.1
* psutil >= 5.4.8

Command line installation using ``conda``
-----------------------------------------

The fastet way to obtain **Python 3.6 together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation. All of the requirements for HiCExplorer can be installed in Anaconda with:

.. code:: bash

    $ conda install hicexplorer -c bioconda -c conda-forge

Command line installation using ``pip``
-----------------------------------------

Install HiCExplorer using the following command:
::

	$ pip install hicexplorer

All python requirements should be automatically installed.

If you need to specify a specific path for the installation of the tools, make use of `pip install`'s numerous options:

.. code:: bash

    $ pip install --install-option="--prefix=/MyPath/Tools/hicexplorer" git+https://github.com/deeptools/HiCExplorer.git


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
