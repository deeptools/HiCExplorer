News and Developments
=====================

Release 3.4
-----------
**28 January 2020**

- Fixing a bug in hicAdjustMatrix: `keep` option had a bug concerning the cutting before the end of a chromosome or the start position was not from the beginning of the chromosome 
- hicCompartmentPolarization was renamed to hicCompartmentalization and got some bug fixes 
- Extending the option on how the observed vs. Expected matrix is computed and adding the parameter `--ligation_factor` to achieve a rescale behaviour of the values as it is implemented in Homer. The same changes are applied to `hicTransform` 
- Improved the documentation 
- Adding an option in hicAverageRegions to select start, end, center or start_end as start index for up/downstream range. 
- hicBuildMatrix: Removed default value of binSize to enable mutually exclusive group error if not one of them is set. Behaviour so far was that the binSize was taken. 
- hicPlotSVL: adding xlegend to plot of SVL ratios to indicate the data points per boxplots are the chromosome ratios 
- hicQuickQC: Removed binSize option of hicQuickQC because it does not matter for QC calculation and adding a sentence to recommend the usage of restriction enzyme and dangling end sequence. Fixing bug issue #464 
- hicNormalize: Adding option in hicNormalize to remove values after the normalization if values are smaller than a given threshold 
- Capture Hi-C modules: Change background model distribution assumption from negative binomial to continuous negative binomial by using Gamma functions as a replacement for the binomial coefficient. Source: https://stats.stackexchange.com/questions/310676/continuous-generalization-of-the-negative-binomial-distribution/311927 
- hicInfo: Implementing feature request #456. The length of chromosomes is now show in the information too 


Release 3.3.1
-------------
**15 November 2019**

- Fixing a bug in the labeling of chicPlotViewpoints if the value range is counted in MB
- Add an option to chicViewpoint to pre-compute a x-fold of p-value over the maximum value of the relative distance


Release 3.3
-----------
**8 October 2019**

- Fixing many bugs:
   - A bug in hicDetectLoops if a sub-matrix was very small
   - A bug in hicPlotMatrix if the region defined by --region was only a chromosome and loops should be plotted too
   - A bug in hicPlotMatrix if a loop region should be plotted and chromosomeOrder argument was used too
   - A bug in hicAggregateContacts (issue #405) if chromosomes were present in the matrix but not in the bed file
   - A bug in hicAdjustMatrix concerning a bed file and consecutive regions, see issue #343
   - A bug in hicAdjustMatrix if a chromosome is present in the matrix but not in the bed file, see issue #397
   - A bug in hicCompartmentsPolarization concerning the arguments 'quantile' and 'outliers' were interpreted as strings but should be integers
   - A bug in hicAdjustMatrix concerning the 'keep' option and how matrices are reordered internally. Thanks @LeilyR

- Added features as requested:
   - hicPCA ignores now masked bins, see issue #342
   - chicPlotViewpoint: 
      - Better legend handling on x-axis
      - Peaks are now display with their fill width
      - Add option `--pValueSignificantLevels` to categorize the p-values in x levels (e.g. 0.001 0.05 0.1)
   - chicViewpoint:
      - adding sorting via viewpoints and not by samples option (--allViewpointsList)
   - Adding an option to hicNormalize to normalize via multiplication and a use defined value (see issues #385, #424)

- Rearrange hicAdjustMatrix to have a better accessibility to its functions from outside of main
- Improving the documentation and fixing grammar / spelling mistakes. Thanks @simonbray
- New script: hicPlotSVL to investigate short range vs long range ratios.


Release 3.2
-----------
** 22 August 2019**

- Adding the new captured Hi-C module. Viewpoint analysis based on a background model, significant interaction detection and differential analysis are provided.
- Adding documentation for captured Hi-C module and a tutorial on how to use it.
- Adding a module to be able to detect quite fast the quality of a Hi-C data set: hicQuickQC.
- Adding a tool to merge loops of different resolutions.
- Improving validation of locations: Presorting is no longer necessary; adding feature to add 'chr' prefix to loop or protein chromosome name
- Change loop detection slightly to improve results and fixed bugs:
   - preselection p-value was ignored and only p-value was used 
   - adding additional test to the peak region test to decrease false discoveries
   - exchanging pThreshold / ln(distance) to remove too low values by a share of the maximum value of the distance. New parameter 'maximumInteractionPercentageThreshold'
- Removal of the folder 'scripts' and its content. These were outdated scripts and will maybe part of regular Hi-C tools in the future.

Release 3.1
-----------
**9 July 2019**

- KR correction improvements: It is now able to process larger data sets like GM12878 primary+replicate on 10kb resolution.
- Adding script for validation of loop locations with protein peak locations
- Adding script hicCompartmentsPolarization: Rearrange the average interaction frequencies using the first PC values to represent the global compartmentalization signal


Release 3.0.2
-------------
**28 June 2019**

- Pinning dependencies to:

   - hicmatrix version 9: API changes in version 10
   - krbalancing version 0.0.4: API changes in version 0.0.5
   - matplotlib version 3.0: Version 3.1 raises 'Not implemented error' for unknown reasons.

- Set fit_nbinom to version 1.1: Version 1.0 Had deprecated function call of scipy > 1.2.
- Small documentation fixes and improvements.


Release 3.0.1
-------------
**5 April 2019**

- Fixes KR balancing correction factors
- Deactivates log.debug


Release 3.0
-----------
**3 April 2019**

- Python 3 only. Python 2.X is no longer supported
- Additional Hi-C interaction matrix correction algorithm 'Knight-Ruiz' as a C++ module for a faster runtime and less memory usage.
- Enriched regions detection tool: 'hicDetectLoops' based on strict candidate selection, 'hicFindEnrichedContacts' was deleted
- Metadata for cooler files is supported: hicBuildMatrix and hicInfo are using it 
- New options for hicPlotMatrix: --loops to visualize computed loops from hicDetectLoops and --bigwigAdditionalVerticalAxis to display a bigwig track on the vertical axis too.


Release 2.2.3
-------------
**22 March 2019**

- This bug fix release patches an issue with cooler files, hicBuildMatrix and the usage of a restriction sequence file instead of fixed bin size.


Release 2.2.2
--------------
**27 February 2019**

- This bug fix release removes reference to hicExport that were forgotten to delete in 2.2. Thanks @BioGeek for this contribution.

Release 2.2.1
-------------
**7 February 2019**

- Muting log output of matplotlib and cooler
- Set version number of hicmatrix to 7
- Optional parameter for hicInfo to write the result to a file instead to the bash

Release 2.2
-----------
**18 January 2019**

This release contains:

- replaced hicExport by hicConvertFormat and hicAdjustMatrix
- extended functionality for hicConvertFormat

   - read support for homer, hicpro, cool, h5
   - write support for h5, homer, cool
   - convert hic to cool
   - creation of mcool matrices

- hicAdjustMatrix

   - remove, keep or mask specified regions from a file, or chromosomes

- hicNormalize

   - normalize matrices to 0 - 1 range or to the read coverage of the lowest given

- hicBuildMatrix

   - support for build mcool

- restructuring the central class HiCMatrix to object oriented model and moved to its own library: `deeptools/HiCMatrix <https://github.com/deeptools/HiCMatrix>`_.

   - Extended read / write support for file formats
   - better (faster, less memory) support for cool format 
   - remove of old, unused code
   - restrict support to h5 and cool matrices, except hicConvertFormat 

- hicFindTADs: Option to run computation per specified chromosomes
- hicPlotTADs: removed code and calls pyGenomeTracks
- hicAverageRegions: Sum up in a given range around defined reference points. Useful to detect changes in TAD structures between different samples. 
- hicPlotAverageRegions: Plots such a average region
- hicTransform: Restructuring the source code, remove of option 'all' because it was generating confusion. Adding option 'exp_obs', exp_obs_norm and exp_obs_lieberman. These three different options use different expectation matrix computations. 
- hicPCA

  - Adding --norm option to compute the expected matrix in the way HOMER is doing it. Useful for drosophila genomes
  - Adding option to write out the intermediate matrices 'obs_exp' and 'pearson' which are necessary in the computation of the PCA


- hicPlotMatrix

  - Add option to clip bigwig values
  - Add option to scale bigwig values


- Removed hicLog2Ration, functionality is covered by hicCompareMatrices
- Extending test cases to cover more source code and be hopefully more stable.
- Many small bugfixes 

Publication
-----------
**13 June 2018**

We are proud to announce our latest publication:

Joachim Wolff, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach, Thomas Manke, Rolf Backofen, Fidel Ramírez, Björn A Grüning. 
"Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and visualization", 
Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W11–W16, doi: https://doi.org/10.1093/nar/gky504

Release 2.1.4
-------------
**25 May 2018**

- cooler file format correction factors are applied as they should be
- parameter '--region' of hicBuildMatrix works with Python 3

Release 2.1.3
-------------
**7 May 2018**

The third bugfix release of version 2.1 corrects an error in hicPlotViewpoint. It adds a feature requested in issue #169 which should have been included in release 2.1 but was accidentally not.

From 2.1 release note:
hicPlotViewpoint: Adds a feature to plot multiple matrices in one image

Release 2.1.2
-------------
**26 April 2018**

The second bug fix release of 2.1 includes:

- documentation improvements
- fixing broken Readthedocs documentation
- Small bug fix concerning hicPlotMatrix and cooler: --chromosomeOrder is now possible with more than one chromosome
- Small fixes concerning updated dependencies: Fixing version number a bit more specific and not that strict in test cases delta values.

Release 2.1.1
-------------
**27 March 2018**

This release fixes a problem related to python3 in which chromosome names were of bytes type

Release 2.1
-----------
**5 March 2018**

The 2.1 version of HiCExplorer comes with new features and bugfixes.

- Adding the new feature `hicAggregateContacts`: A tool that allows plotting of aggregated Hi-C sub-matrices of a specified list of positions.
- Many improvements to the documentation and the help text. Thanks to Gina Renschler and Gautier Richard from the MPI-IE Freiburg, Germany.
- hicPlotMatrix

    - supports only bigwig files for an additional data track.
    - the argument `--pca` was renamed to `--bigwig`
    - Smoothing the bigwig values to neighboring bins if no data is present there
    - Fixes to a bug concerning a crash of `tight_layout`
    - Adding the possibility to flip the sign of the values of the bigwig track
    - Adding the possibility to scale the values of the bigwig track 

- hicPlotViewpoint: Adds a feature to plot multiple matrices in one image
- cooler file format

   - supports mcool files
   - applies correction factors if present
   - optionally reads `bin['weight']`

- fixes

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

 * bug fixes for ``--perChr`` option in hicPlotMatrix
 * eigenvector track with ``--pca`` for hicPlotMatrix
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
