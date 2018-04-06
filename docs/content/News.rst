News and Developments
=====================

Release 2.1.1
------------
**27 March 2018**

This release fixes a problem related to python3 in which chromosome names were of bytes type

Release 2.1
-----------
**5 March 2018**

The 2.1 version of HiCExplorer comes with new features and bugfixes.

- Adding the new feature `hicAggregateContacts`: A tool that allows plotting of aggregated Hi-C sub-matrices of a specified list of positions.
- Many improvements to the documentation and the help text. Thanks to @GinaRe and @gtrichard.
- hicPlotMatrix:
    - supports only bigwig files for an additional data track.
    - the argument `--pca` was renamed to `--bigwig`
    - Smoothing the bigwig values to neighboring bins if no data is present there
    - Fixes to a bug concerning a crash of `tight_layout`
    - Adding the possibility to flip the sign of the values of the bigwig track
    - Adding the possibility to scale the values of the bigwig track 
- hicPlotViewpoint: Adds a feature to plot multiple matrices in one image
- cooler file format:
   - supports mcool files
   - applies correction factors if present
   - optionally reads `bin['weight']`
- fixes:
    - a crash in hicPlotTads if `horizontal lines` were used
    - checks if all characters of a title are ASCII. If not they are converted to the closest looking one.
- Updated and fixate version number of the dependencies


Release 2.0
-----------

**December 21, 2017**

This release makes HiCExplorer ready for the future:

* Python 3 support
* `Cooler <https://github.com/mirnylab/cooler>`_ file format support
* A/B comparment analysis
* Improved visualizations
 * bug fixes for --perChr option in hicPlotMatrix
 * eigenvector track with --pca for hicPlotMatrix
 * visualization of interactions around a reference point or region with hicPlotViewpoint
* Higher test coverage
* re-licensing from GPLv2 to GPLv3

Release 1.8.1
--------------

**November 27, 2017**

Bug fix release:

* a fix concerning the handling chimeric alignments in hicBuildMatrix. Thanks to Aleksander Jankowski @ajank
* handling of dangling ends was too strict
* improved help message in hicBuildMatrix

Release 1.8
-----------

**October 25, 2017**

This release is adding new features and fixes many bugs:

 * hicBuildMatrix: Added multicore support, new parameters --threads and --inputBufferSize
 * hicFindTADs:
  * One call instead of two: hicFindTADs TAD_score and hicFindTADs find_TADs merged to hicFindTADs.
  * New multiple correction method supported: False discovery rate. Call it with --correctForMultipleTesting fdr and --threshold 0.05.
 * Update of the tutorial: mES-HiC analysis.
 * Additional test cases and docstrings to improve the software quality
 * Fixed a bug occurring with bigwig files with frequent NaN values which resulted in only NaN averages
 * hicPlotTADs: Support for plotting points
 * Moved galaxy wrappers to https://github.com/galaxyproject/tools-iuc
 * Fixed multiple bugs with saving matrices
 * hicCorrelate: Changes direction of dendograms to left

Release 1.7.2
-------------

**April 3, 2017**

 * Added option to plot bigwig files as a line hicPlotTADs
 * Updated documentation
 * Improved hicPlotMatrix --region output
 * Added compressed matrices. In our tests the compressed matrices are significantly smaller.


**March 28, 2017**

Release 1.7
-----------

**March 28, 2017**

This release adds a quality control module to check the results from hicBuildMatrix. By default, now hicBuildMatrix
generates a HTML page containing the plots from the QC measures. The results from several runs of hicBuildMatrix can
be combined in one page using the new tool hicQC.

Also, this release added a module called hicCompareMatrices that takes two Hi-C matrices and computes
the difference, the ratio or the log2 ratio. The resulting matrix can be plotted with hicPlotMatrix
to visualize the changes.


Preprint introducing HiCExplorer is now online
----------------------------------------------

**March 8, 2017**

Our #biorXiv preprint on DNA sequences behind Fly genome architecture is online!

Read the article here : `<http://biorxiv.org/content/early/2017/03/08/115063>`_

In this article, we introduce HiCExplorer : Our easy to use tool for Hi-C data analysis, also available in `Galaxy <https://galaxyproject.org/>`_.

We also introduce `HiCBrowser <https://github.com/maxplanck-ie/HiCBrowser>`_ : A standalone software to visualize Hi-C along with other genomic datasets.

Based on HiCExplorer and HiCBrowser, we built a useful resource for anyone to browse and download the chromosome
conformation datasets in Human, Mouse and Flies. It's called `the chorogenome navigator <http://chorogenome.ie-freiburg.mpg.de/>`_

Along with these resources, we present an analysis of DNA sequences behind 3D genome of Flies. Using high-resolution
Hi-C analysis, we find a set of DNA motifs that characterize TAD boundaries in Flies and show the importance of these motifs in genome organization.

We hope that these resources and analysis would be useful for the community and welcome any feedback.


HiCExplorer wins best poster prize at VizBi2016
-----------------------------------------------

**March 20, 2016**

We are excited to announce that HiCExplorer has won
the `NVIDIA Award for Best Scientific Poster <https://vizbi.org/blog/2016/02/11/nvidia-award-for-best-scientific-poster/>`_
in VizBi2016, the international conference on visualization of biological data.

`Read more here <https://vizbi.org/blog/2016/03/20/winner-of-nvidia-best-scientific-poster-award-2/>`_

This was our poster :

.. image:: https://vizbi.org/Posters/Images/2016/B12.png
   :scale: 50 %
   :alt: HiCExplorer
   :align: left
