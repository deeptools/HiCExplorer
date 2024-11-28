|azure| |rtd| |conda| |docker| |galaxy|

.. |azure| image:: https://dev.azure.com/wolffj/HiCExplorer/_apis/build/status/deeptools.HiCExplorer?branchName=master
	:target: https://dev.azure.com/wolffj/HiCExplorer/_build/latest?definitionId=2&branchName=master
.. |rtd| image:: https://readthedocs.org/projects/hicexplorer/badge/?version=latest
   :target: http://hicexplorer.readthedocs.io/?badge=latest
.. |conda| image:: https://anaconda.org/bioconda/hicexplorer/badges/installer/conda.svg
   :target: https://anaconda.org/bioconda/hicexplorer
.. |galaxy| image:: https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==
   :target: https://usegalaxy.eu/root?tool_id=hicexplorer_hicplotviewpoint
.. |docker| image:: https://quay.io/repository/biocontainers/hicexplorer/status
   :target: https://quay.io/repository/biocontainers/hicexplorer



HiCExplorer
===========

Set of programs to process, analyze and visualize Hi-C, Micro-C and cHi-C data
------------------------------------------------------------------------------

Sequencing techniques that probe the 3D organization of the genome generate large amounts of data whose processing,
analysis and visualization is challenging. Here, we present HiCExplorer, a set of tools for the analysis and
visualization of chromosome conformation data. HiCExplorer facilitates the creation of contact matrices, correction
of contacts, TAD detection, A/B compartments, merging, reordering or chromosomes, conversion from different formats including
`cooler <https://github.com/mirnylab/cooler>`_ and detection of long-range contacts. Moreover, it allows the visualization of
multiple contact matrices along with other types of data like genes, compartments, ChIP-seq coverage tracks (and in general
any type of genomic scores), long range contacts and the visualization of viewpoints.
With version 3.7.6 we introduce the support for Micro-C data for the build of matrices.

Single-cell Hi-C data
---------------------

We provide the scHiCExplorer to create, manipulate, analyse and visualize single-cell Hi-C data in its own software:
The  `scHiCExplorer <https://github.com/joachimwolff/schicexplorer>`_.

Citation:
^^^^^^^^^

Joachim Wolff, Rolf Backofen, Björn Grüning.
**Loop detection using Hi-C data with HiCExplorer**, GigaScience, Volume 11, 2022, giac061, https://doi.org/10.1093/gigascience/giac061

Joachim Wolff, Leily Rabbani, Ralf Gilsbach, Gautier Richard, Thomas Manke, Rolf Backofen, Björn A Grüning.
**Galaxy HiCExplorer 3: a web server for reproducible Hi-C, capture Hi-C and single-cell Hi-C data analysis, quality control and visualization, Nucleic Acids Research**, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W177–W184, https://doi.org/10.1093/nar/gkaa220

Joachim Wolff, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach, Thomas Manke, Rolf Backofen, Fidel Ramírez, Björn A Grüning. 
**"Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and visualization", Nucleic Acids Research**, Volume 46, Issue W1, 2 July 2018, Pages W11–W16, doi: https://doi.org/10.1093/nar/gky504

Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening, Kin Chung Lam, Bianca Habermann, Asifa Akhtar, Thomas Manke.
**"High-resolution TADs reveal DNA sequences underlying genome organization in flies". Nature Communications**, Volume 9, Article number: 189 (2018), doi: https://doi.org/10.1038/s41467-017-02525-w

.. image:: ./docs/images/hicex3.png

Availability
^^^^^^^^^^^^

HiCExplorer is available as a **command line suite of tools** on this very GitHub repository and also on other platforms (detailed in *Installation* below).

A **Galaxy HiCExplorer version** is directly available to users at http://hicexplorer.usegalaxy.eu. Training material is available at the `Galaxy Training Network <http://galaxyproject.github.io/training-material/topics/epigenetics/tutorials/hicexplorer/tutorial.html>`_,
while a Galaxy Tour is available `here <https://hicexplorer.usegalaxy.eu/tours/hixexplorer>`_ for users not familiar with this platform. Galaxy HiCExplorer is also available as a Docker image at the `Docker Galaxy HiCExplorer GitHub repository <https://github.com/deeptools/docker-galaxy-hicexplorer>`_. Finally, this Galaxy version is available on the `Galaxy Tool Shed <https://toolshed.g2.bx.psu.edu/>`_ and on the corresponding `GitHub repository <https://github.com/galaxyproject/tools-iuc>`_.



Installation
^^^^^^^^^^^^

With version 3.0, HiCExplorer is available for Python 3 only, the Python 2 support is discontinued. HiCExplorer can be installed with conda.

-  Anaconda and GitHub for command line usage.
-  Toolshed and Docker image for its integration on Galaxy servers.

There are many easy ways to install HiCExplorer. Details can be found
`here <https://hicexplorer.readthedocs.io/en/latest/content/installation.html>`_.


**We strongly recommended to use conda to install HiCExplorer.**


Command line version
++++++++++++++++++++

Install with conda
__________________

The easiest way to install HiCExplorer is using `BioConda <http://bioconda.github.io/>`_
::

   $ conda install hicexplorer -c bioconda -c conda-forge
   
   
We highly recommend conda environments to separate software from each other. With it, different versions of dependencies do not interfere with each other.

::

   $ conda create --name hicexplorer hicexplorer=3.6 python=3.8 -c bioconda -c conda-forge
   $ conda activate hicexplorer
   
To deactivate the environment use:

::

   $ conda deactivate
   

To learn more about conda and environments, please consider the following `documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#>`_.
   
   

Install by cloning this repository
__________________________________

You can install any one of the HiCExplorer branches on command line
(linux/mac) by cloning this git repository :

::

    $ git clone https://github.com/deeptools/HiCExplorer.git
    $ cd HiCExplorer
    $ python setup.py install

If you don't have root permission, you can set a specific folder using the ``--prefix`` option

::

	$ python setup.py install --prefix /User/Tools/hicexplorer

If you don't use conda, please take care of all dependencies on your own.

Galaxy version
++++++++++++++

Install with Docker
___________________

Installation instructions as a Docker image can be followed at https://github.com/deeptools/docker-galaxy-hicexplorer.


Install with Tool Shed
______________________

Galaxy HiCExplorer is part of the `Galaxy Tool Shed <https://toolshed.g2.bx.psu.edu/>`_ and can be installed from there to any Galaxy server following `this link <https://toolshed.g2.bx.psu.edu/repository/browse_repository?id=f1554978eeb3da8b>`_.


Documentation:
^^^^^^^^^^^^^^
Please visit our complete documentation `Here <http://hicexplorer.readthedocs.org/>`_. This documentation is also available directly within `Galaxy <http://hicexplorer.usegalaxy.eu/>`_.
