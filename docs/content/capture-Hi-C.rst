.. _capture-Hi-C:

Captured Hi-C data analysis
===========================

.. contents::
    :local:
.. toctree::
    :maxdepth: 1

    capture-Hi-C


How we use HiCExplorer to analyse cHi-C data
--------------------------------------------

This How-to for is based on the published dataset from Andrey et al 2017(doi:10.1101/gr.213066.116). For the tutorial we use the samples FL-E13.5 and MB-E-10.5. 


Download the raw data
---------------------

Please download the raw data via the following links or via `NCBI GSE84795 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84795>`__ .

+--------------------------------------+---------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+
| Dataset                              | forward                                                                                           | reverse                                                                                           |
+======================================+===================================================================================================+===================================================================================================+
| CC-FL-E135-Wt-Mm-Rep1                | `SRR3950565_1 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565_1.fastq.gz>`__ | `SRR3950565_2 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565_2.fastq.gz>`__ |
+--------------------------------------+---------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+
| CC-FL-E135-Wt-Mm-Rep2                | `SRR3950566_1 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566_1.fastq.gz>`__ | `SRR3950566_2 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566_2.fastq.gz>`__ |
+--------------------------------------+---------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+
| CC-MB-E105-Wt-Mm-Rep1                | `SRR3950559_1 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559_1.fastq.gz>`__ | `SRR3950559_2 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559_2.fastq.gz>`__ |
+--------------------------------------+---------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+
|CC-MB-E105-Wt-Mm-Rep2                 | `SRR3950560_1 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560_1.fastq.gz>`__ | `SRR3950560_2 <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560_2.fastq.gz>`__ |
+--------------------------------------+---------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------+


Mapping
-------

Map the files with a mapper of your choice against the mm9 reference genome, as an example the mapping with bowtie2 is shown.

.. code:: bash

    bowtie2 -x mm9_index --threads 8 -U SRR3950565_1.fastq.gz --reorder | samtools view -Shb - > SRR3950565_1.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950565_2.fastq.gz --reorder | samtools view -Shb - > SRR3950565_2.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950566_1.fastq.gz --reorder | samtools view -Shb - > SRR3950566_1.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950566_2.fastq.gz --reorder | samtools view -Shb - > SRR3950566_2.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950559_1.fastq.gz --reorder | samtools view -Shb - > SRR3950559_1.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950559_2.fastq.gz --reorder | samtools view -Shb - > SRR3950559_2.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950560_1.fastq.gz --reorder | samtools view -Shb - > SRR3950560_1.bam
    bowtie2 -x mm9_index --threads 8 -U SRR3950560_2.fastq.gz --reorder | samtools view -Shb - > SRR3950560_2.bam

    


Create cHi-C matrices
---------------------

To create cHi-C matrix we use HiCExplorer's hicBuildMatrix for each replicate separately and merge the replicate to one matrix later. Like Andrey et al we use a resolution of 1kb and use the restriction enzyme DpnII.

.. code:: bash

    hicBuildMatrix --samFiles SRR3950565_1.bam SRR3950565_2.bam  --binSize 1000 --restrictionSequence GATC --outFileName SRR3950565.cool --QCfolder SRR3950565_QC --threads 6
    hicBuildMatrix --samFiles SRR3950566_1.bam SRR3950566_2.bam  --binSize 1000 --restrictionSequence GATC --outFileName SRR3950566.cool --QCfolder SRR3950566_QC --threads 6
    hicBuildMatrix --samFiles SRR3950559_1.bam SRR3950559_2.bam  --binSize 1000 --restrictionSequence GATC --outFileName SRR3950559.cool --QCfolder SRR3950559_QC --threads 6
    hicBuildMatrix --samFiles SRR3950560_1.bam SRR3950560_2.bam  --binSize 1000 --restrictionSequence GATC --outFileName SRR3950560.cool --QCfolder SRR3950560_QC --threads 6


.. code:: bash

    hicSumMatrix --matrices SRR3950565.cool SRR3950566.cool --outFileName FL-E13-5.cool
    hicSumMatrix --matrices SRR3950559.cool SRR3950560.cool --outFileName MB-E10-5.cool




Creation of reference point file
--------------------------------

Andrey et al writes they used in total 460 reference points, but 24 were removed due to low sequence coverage or non correspondence to a promoter region, leading to 446 in total.

To reproduce this, we need the all reference points which are published in Supplementary Table `S2 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S2.xlsx>`__ and `S8 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S8.xlsx>`__ .  


It is the easiest to create via Excel the reference point file in the following format and store it as a tab separated file:

.. code:: bash

    chr1	4487435	4487435 Sox17

Or just download the prepared file. We will do the Quality control on our own and compare to the results of Andrey.



Quality control
^^^^^^^^^^^^^^^

.. code:: bash

    chicQualityControl -m FL-E13-5.cool MB-E10-5.cool -rp reference_points.bed -s 0.025 -t 20


Background model
^^^^^^^^^^^^^^^^

The background model computes for both samples all viewpoints given by the reference points in a range defined by the parameter `fixateRange`. We recommend to set it to 500kb because real interactions above the range 
are rarely observed and very low interaction numbers like 1 are already considered as significant. With this setting, only the interactions in a range 500kb up- and downstream of the reference point are considered for each viewpoint.
Based on this data, two background models are computed, the first one simply computes the average per relative distance to the reference point, and second, a negative binomial distribution per relative distance to
 the reference point is fitted. This first one is used for filtering in the significant interaction evaluation by an x-fold factor and for plotting. The negative binomial model is more important, it is used to 
 compute per relative distance in each sample a p-value and based on it the decision if an interaction is considered as significant is made.

.. code:: bash

    chicViewpointBackgroundModel -m FL-E13-5.cool MB-E10-5.cool --fixateRange 500000 -t 20 -rp reference_points.bed -o background_model.bed


Viewpoint computation
^^^^^^^^^^^^^^^^^^^^^

In this step the viewpoint for each reference point listed in `reference_points.bed` is computed, the up- and downstream range can be given via `--range upstream downstream` and please use the same `--fixateRange` as in the background model computation.
For each relative distance the x-fold over the average value of this relative distance is computed and each location gets a p-value based on the background negative binomial distribution for this relative distance.
For each viewpoint one viewpoint file is created and stored in the folder given by the parameter `--outputFolder`. The name of each viewpoint file starts with its gene / promoter name, followed by the sample (given by the name of the matrix) and the
exact location. For example the viewpoint `chr1	4487435	4487435 Sox17` from `MB-E10-5.cool` matrix will be called `Sox17_MB-E10-5_chr1_4487435_4487435.bed` and looks as follows:

..code:: bash

    TODO: example viewpoint file stuff


in case the parameter `--writeFileNamesToFile` is set, the viewpoint file names are written to a file which can be used for batch processing in the later steps. Each sample is combined with each other sample for each viewpoint to enable the batch processing
for the differential analysis. Example: matrices `FL-E13-5.cool` and  `MB-E10-5.cool` with the two reference points x and y:

..code:: bash

    x_FL-E13-5_chr1_1000_2000.bed
    x_MB-E10-5_chr1_1000_2000.bed
    y_FL-E13-5_chr1_3000_4000.bed
    y_MB-E10-5_chr1_3000_4000.bed



.. code:: bash

    chicViewpoint -m FL-E13-5.cool MB-E10-5.cool --averageContactBin 5 --range 1000000 1000000 -rp views10.bed -bmf background_model.bed --writeFileNamesToFile TSS_view1.txt --outputFolder TSS_view1 --fixateRange 500000 --threads 5


Significant interactions detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


To detect significant interactions and to prepare a target file for each viewpoint which will be used for the differential analysis, the script `chicSignificantInteractions` is used. It offers two modi: Either the user can specify 
a x-fold value or a loose p-value. The first one considers all interactions with a minimum x-fold over the average background for its relative distribution as a candidate or second, all interactions with a loose p-value or less are considered. 
These are preselcetion steps to be able to detect wider peaks in the same way as sharp ones. All detected candidates are merged to one peak in the case they are direct neighbors and the sum of all interaction values of this neighborhood
 is used to compute a new p-value. The p-value is computed based the relative distance negative binomial distribution of the interaction with the original highest interaction value. All considered peaks are accepted as significant interactions if
 their p-value is small as the threshold `--pvalue`.

To exclude interaction with an interaction value which is too less the parameter `--peakInteractionsThreshold` can be set.



.. code:: bash


Batch mode
~~~~~~~~~~



Aggregate data for differential test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash



Differential test
^^^^^^^^^^^^^^^^^



.. code:: bash

Plotting of Viewpoints
^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash
