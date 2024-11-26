News and Developments
=====================

Release 3.7.6
-------------
**27 November 2024**

- Add a new hicBuildMatrixMicroC script to build matrices from Micro-C data. It is the same as hicBuildMatrix but without the enforcement of the restriction cut site file, restriction sequence and dangling end since this is not necessary for Micro-C data.
- Update to include support for Python 3.11 and 3.12


Release 3.7.5
-------------
**June 2024**

- Update the version file.

Release 3.7.4
-------------
**24 April 2024**

- Allow chicAggregateStatistic.py to to extract the aggregated data from the views.hdf based on differential.hdf or differential_target.bed. Now the BED may have the target name in the 4th column. In that case, the aggregation is done per target.
- Allow hicCorrectMatrix.py to write filtered out regions to a BED file

Thanks @pavanvidem

Warning: In this version the version file has not been modified so the tools gives 3.7.3 as version.

Release 3.7.3
-------------
**17 November 2023**

- Maintenance update for HiCExplorer to keep up to date with APIs of dependencies
- Add additional of the polarization ratio to the output of hicCompartmentalization. Thanks @contessoto.


Release 3.7.2
-------------
**1 October 2021**

- Fixing a bug in hicHyperoptDetectLoops concerning the validation of the loop locations
- Fixing a bug in hicValidateLocations concerning the computation of TADs
- Adding the -v option for the version for hicInfo
- Adding documentation concerning the native file formats of HiCExplorer
- Fixing the way A/B compartments are computed based on Lieberman-Aiden 2009. The intermediate Pearson matrix is not used anymore.

Release 3.7.1
-------------
**9 August 2021**

The promised functions for the strand information in 3.7 were not part of the release due to confusion in merging multiple branches. Therefore they are now included in the 3.7.1 release:

- hicAggregateContacts: Option to consider the strand orientation (#633)
- hicAverageRegions: Option to consider the strand orientation (#633)

Additionally, a small bug fix in hicInfo concerning the correct sum of the interaction matrix.

Release 3.7
-----------
**27 July 2021**

- A new TAD prediction method: `hicTADClassifier` and a method to train own classifiers `hicTrainTADClassifier`. Thanks @AlbertLidel for the design and implementation.
- New file formats for capture Hi-C modules. `chicViewpoint`, `chicSignificantInteractions`, `chicAggregateStatistic` and `chicDifferentialTest` use now HDF5 based files. A new export script `chicExportData` is provided to retrieve text files from these HDF5 files.
- Implementing a few feature requests:
    - `hicPlotMatrix`: TADs can be visualized with hicPlotMatrix
    - `hicAdjustMatrix` is able to remove intra- or inter-chromosomal contacts (#664, #704)
    - `hicValidateLocations`: An option to validate TADs and to use additional data stored in `cool` matrices
    - `hicPCA`: Adding a function to select which eigenvector should be used for the output (#669)
    - `hicConvertFormat`: Adding the function to export hicpro file format and to import `2D-text` files.
    - `hicFindRestSites`: Support of multiple restriction cut sequences (#659)
    - `hicPlotMatrix`: Option for loop locations spanning more than one bin to define if the start, center or end should be used for plotting (#640)
    - `hicInterIntraTAD`: A new script to compute the ratio of inter and intra TAD. ($404)
    - `hicAggregateContacts`: Option to consider the strand orientation (#633)
    - `hicAverageRegions`: Option to consider the strand orientation (#633)
    - `hicCompareMatrices`: An option to not normalize the matrices before the computation. (#677, #645) Thanks @lldelisle
    - `hicDifferentialTAD`: Adding rank sum statistics to the output (#728, #727). Thanks @dawe
    - `hicPlotDistVsCounts`: Adding a function to plot the counts vs distance in TAD regions. (#696) Thanks @lldelisle
- Bug fixes:
     - `hicCorrectMatrix`: A bug that lead to wrong correction factors for the KR correction for `cool` files (#724)
     - `hicDifferentialTAD`: Solved multicore issue related to skipping data at the start and end of chromosomes (#725, #685)
     - `hicHyperoptDetectLoops`: Added an option to set if the `chr` prefix should be added or removed (#723)
     - `hicPCA`: Solving an issue if the region defined by the gene track is larger the region stored in the interaction matrix (#655, #710, #716, #719)
     - `hicPCA`: Fixing a bug where the masking of bins was automatically applied which lead to differing matrix dimensions for the e.g. the Pearson correlation matrices (#618)
     - `hicBuildMatrix`: Solving a bug if multiple restriction cut sites have the same dangling ends (#720)
     - `hicBuildMatrix`: Solving a bug that the parameter `--removeSelfLigations` was always set to true. Changed parameter name to `--keepSelfLigations` to keep the functionality. If the parameter is not set, the self ligations are removed.
     - `hicBuildMatrix`: If a region is specified, only the restrictionCutSite file information for that region is loaded to save memory (#646)
     - `hicConvertFormat`: Fixing a bug to copy the genome annotation information in the case of a `cool` to `cool` file conversion (#657)
     - `hicCorrelate`: Correcting the range of colors for the heatmap (#585)
     - `hicCompartmentalization`: Fixed index bug (#635, #637) Thanks  @LeilyR
- Updating `hicBuildMatrix` to be able to work with biopython versions > 1.77. Thanks @lldelisle (#691)


Release 3.6
-----------
**10 November 2020**

- hicAggregateContacts, thanks @LeilyR:
    - hicAggragateContact has been updated to be able to handle inter chromosomal contacts as well as inter chromosomal contacts
    - Added scikit-learn to dependencies
- hicBuildMatrix: Fixing another bug concerning scaffolds without any restriction enzyme cut site
- Updated dependencies
- Adding default values to the documentation. Thanks @lldelisle
- hicTransform: Fixing a bug in case one of the intermediate matrices is empty
- Official Python 3.8 support:
   - Manually setting 'fork' as the start method for multiprocessing because the default on macOS was set to 'spawn'


Release 3.5.3
-------------
**14 October 2020**

- Bug fix release: 
   - Reads from scaffolds without any restriction enzyme cut site are considered as 'same fragment'. An appearance of such a read will not lead to a crash anymore
- Minor documentation improvements


Release 3.5.2
-------------
**06 October 2020**

- Bug fix release: 
   - enforcing version 15 of HiCMatrix. Version 14 had a bug concerning the application of the correction factors of cool files. See issue #595
   - Fixing a bug in hicDetectLoops in single-core mode. Thanks @alecw
   - Fixing a bug in hicDifferentialTAD concerning multiple chromosomes in the bed file. See issue #587
   - Updating dependencies to newest versions, except biopython. Forcing here <1.77 because their API change is breaking our source code. See issue #609
   - Fixing #596
- Changing internal structure of the docs. Navigation should be better now.


Release 3.5.1
--------------
**11 August 2020**

- patch hicCorrelate
- hicBuildMatrix: Better help text
- patch for hicPlotMatrix if matrix does not start at 0 and --region is used
- bug fix for remove option in hicAdjustMatrix

Release 3.5
-----------
**10 July 2020**

- Major update for hicDetectLoops: Results are now closer to HiCCUPS, it is faster and needs less memory. 
- hicHyperoptDetectLoops: New tool to compute the best parameter setting if the number of expected loops is given and a protein boundary file (e.g. CTCF on mammals) is provided
- hicHyperoptDetectLoopsHiCCUPS: New tool to compute the best parameter setting for Juicers HiCCUPS. HiCCUPS and all its dependencies are not part of HiCExplorer and must be provided by the user. Number of expected loops and a protein boundary file (e.g. CTCF on mammals) must be provided.
- hicMergeDomains: New tool to compute a merge of TADs computed on different resolutions. Moreover it provides a cleaning of the boundaries with the help of a protein peak file, and the hierarchical dependencies of TADs can be plotted. This tool is the result of the Bachelor thesis from Sarah Domogalla (@SarahMaria27). Thanks for your contribution!
- hicDifferentialTAD: New tool to compute differential TADs between two samples
- Bug fix for hicFindTADs: The format of the intermediate z-score matrices are now depending on the format of the input Hi-C matrix
- Bug fix for chic*-modules: Fixate range is now correct applied.
- chicSignificantInteractions, hicDetectLoops: Option to use a per relative position p-value with a p-value threshold file
- Adding hicCreateThresholdFile: A script to generate a per position p-value threshold for hicDetectLoops and chicSignificantInteractions 
- Bug fix for hicPlotMatrix:
   - multiple bigwig tracks in the vertical setting are now supported
   - correct plot of bigwig if the given matrix does not start at the beginning of the chromosome
   - new parameters 'increaseFigureWidth' and 'increaseFigureHeight' to add more space to the plot if multiple bigwigs are plotted. Adjust this parameter to correct the plot of the heatmap which may be not quadratic anymore.
   - restriction of the loop regions to the user given range. This effects especially SVGs that will now contain less objects as before.
- New feature for hicBuildMatrix:
   - multiple restriction cut sequences, dangling ends and restriction cut sites files are now supported
   - restriction cut sequences, dangling ends and restriction cut sites files are now mandatory parameters. This is now enforced to guarantee a correct detection of self ligations and self circles
- hicPrepareQCreport: New support for multiple dangling ends
- hicQuickQC: restriction cut sequences, dangling ends and restriction cut sites files are now mandatory parameters
- hicFindRestSite: gz support for fasta file
- Add fallback modes to multiple scripts if the parallelization fails.
- hicAggregate: interactions between two bed files by comparing every row in the first bed file with its corresponding row in the second file. (issue #390)
- hicAdjustMatrix:  fix #341, 454
   - fixed --action remove:  it actually remove the part from the matrix previously was masking it
   - the case where the end of the chromosome need to be removed.
- New Azure testing for the general test cases, the trivial ones run on travis

Publication
-----------
**02 July 2020** 

Galaxy HiCExplorer 3: a web server for reproducible Hi-C, capture Hi-C and single-cell Hi-C data analysis, quality control and visualization 
Joachim Wolff, Leily Rabbani, Ralf Gilsbach, Gautier Richard, Thomas Manke, Rolf Backofen, Björn A Grüning
Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W177–W184


Release 3.4.3
-------------
**23 March 2020**

- Fixing the wrong p-value computation in for chicViewpoint. New method is more accurate for floating points.
- Fixing a bug in chicViewpointBackgroundModel and chicQualityControl if an non-existing reference point was used.
- Improving all chic* modules with a capturing of errors in sub-processes. It is now guaranteed that the computation will terminate. Either successful or by error. 
- Add option 'truncateZero' to chicViewpointBackgroundModel: This removes all zero values for the distributions before fitting to fight over dispersion.
- Add option 'decimalPlaces' to chicViewpoint to adjust the decimal places in the output for all floating values. Helpful for really small p-values
- Add option 'truncateZeroPvalues' to chicSignificantInteractions to set all p-values which are 0 to 1 and are therefore ignored.
- Add option 'truncateZeroPvalues' to chicPlotViewpoint to set all p-values which are 0 to 1 and do not disturb the presentation of real p-values

Release 3.4.2
-------------
**7 March 2020**

- This release fixes the wrong name scheme which was used in the chicModules. The most .bed files are now .txt files.
- hicDetectLoops got an inner chromosome parallelization to decrease the compute time.
- hicPlotMatrix got three new parameters: rotationX, rotationY and fontSize to adjust the position and font size of the labels. We hope this can lead in certain cases to a a better readability
- hicPlotMatrix: fixed a bug that occurred if the list of chromosomes was given and the last chromosome appeared as an additional label. 
- Improving and updating the documentation.


Preprint
--------
**6 March 2020**

The preprint of the loop detection algorithm is online via biorXiv: `<https://www.biorxiv.org/content/10.1101/2020.03.05.979096v1>`_



Release 3.4.1
-------------
**3 February 2020**

- This release fixes a bug in chicViewpoint that caused a crash if the data to be averaged is less than the window size.

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
