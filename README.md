======================================================================
HiCExplorer
======================================================================
### Set of programms to process, analyze and visualize Hi-C data.

HiCexplorer addresses the common tasks of Hi-C analysis from processing to visualization.

![gallery](https://raw.githubusercontent.com/maxplanck-ie/HiCExplorer/master/examples/images/hicexplorer.png?token=AEu_1VmdSzz0lipVV1DMKuYgYcIjUb4qks5U6zbwwA%3D%3D)

Usage of HiCExplorer
==============

Examples of usage

```shell
# build matrix from idependently mated read pairs
$ hiCBuildMatrix --samFiles mate1.sam mate2.sam --binSize 10000 --restrictionSequence GATC -outBam hic.bam -o hic_matrix.npz
# this creates two files, a bam file containing only the valid Hi-C read pairs
# and a matrix summarizing the Hi-C contacts at the given resolution.

# correct Hi-C matrix
$ hiCorrectMatrix -m hic_matrix.npz -o hic_corrected.npz

# visualize the corrected matrix
$ hiCPlotMatrix -m hic_corrected.npz -o hic_plot.png
```



There are seven major functions available in MACS serving as sub-commands.

| tool 				| description	|
| ----------------------------- | ---------------------------------- |
| hicBuildMatrix 		| Creates a Hi-C matrix using the aligned BAM files of the Hi-C sequencing reads. 	|
| hicCorrectMatrix 		| Uses iterative correction to remove biases from a Hi-C matrix 	|
| hicFindEnrichedContacts  	| Identifies enriched Hi-C contacts.            	|
| hicCorrelate 			| Computes and visualises the correlation of Hi-C matrices            	|
| hicFindTADs 			| Identifies Topologically Associating Domains (TADs)            	|
| hicMergeMatrixBins		| Merges consecutives bins on a Hi-C matrix to reduce resolution     |
| hicPlotMatrix			| Plots a Hi-C matrix as a heatmap |
| hicPlotTADs			| Plots TADs as a track that can be combined with other tracks (genes, signal, interactions)|
| hicSumMatrices		| Adds Hi-C matrices of the same size|


-------------------------------------------------------------------------------------------------------------------

<a name="installation"/></a>
Installation
---------------

[General Installation](#general)

[Installation on a Mac](#mac)


<a name="general"/></a>
__A second option is to clone the repository:__
	
	$ git clone https://github.com/maxplanck-ie/HiCExplorer.git
	$ cd HiCExplorer
	$ python setup.py install
	
By default, the script will install the python library and executable
codes globally, which means you need to be root or administrator of
the machine to complete the installation. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script.

	$ python setup.py --help

For example, to install under a specific location use:

	$ python setup.py install --prefix <target directory>

<a name="mac"></a>
### Installation on a MAC

The easiest way to get numpy and scipy dependencies is to install the
[Anaconda Scientific Python Distribution][]. After installation, open
a terminal ("Applications" → "Terminal") and follow the [General Installation](#general)
  	   
If individual installation of the dependencies is preferred, follow 
those steps:

Requirement: Python 2.7 installed

Download the packages and install them using dmg images:
- http://sourceforge.net/projects/numpy/files/NumPy/
- http://sourceforge.net/projects/scipy/files/scipy/

Then open terminal ("Applications" → "Terminal")
and follow the [General Installation](#general)


This tool suite is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).

