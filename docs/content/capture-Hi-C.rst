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

To reproduce this, we need the all reference points which are published in Supplementary Table `S2 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S2.xlsx>`__ and `S8 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S8.xlsx >`__ .  


It is the easiest to create via Excel the reference point file in the following format and store it as a tab separated file:

chr1	4487435	4487435 Sox17

Or just download the prepared file. We will do the Quality control on our own and compare to the results of Andrey.



Quality control
^^^^^^^^^^^^^^^

.. code:: bash

    chicQualityControl -m FL-E13-5.cool MB-E10-5.cool -rp reference_points.bed -s 0.025 -t 20


Background model
^^^^^^^^^^^^^^^^


.. code:: bash

    chicViewpointBackgroundModel -m FL-E13-5.cool MB-E10-5.cool -t 20 -rp reference_points.bed -o background_model.bed


Viewpoint computation
^^^^^^^^^^^^^^^^^^^^^

.. code:: bash

chicViewpoint -m FL-E13-5.cool MB-E10-5.cool --averageContactBin 5 --range 2000000 2000000 -rp views10.bed -bmf background_model.bed --writeFileNamesToFile TSS_view1.txt --outputFolder TSS_view1 --fixateRange 500000 --threads 5


Significant interactions detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash




Aggregate data for differential test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash



Differential test
^^^^^^^^^^^^^^^^^



.. code:: bash

Plotting of Viewpoints
^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash
