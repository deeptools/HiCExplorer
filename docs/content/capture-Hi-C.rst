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

Please download the raw data via the following links or via NCBI GSE84795.

CC-FL-E135-Wt-Mm-Rep1: `Forward reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565_1.fastq.gz>`__ and 
                        `reverse reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/005/SRR3950565/SRR3950565_2.fastq.gz>`__ .

CC-FL-E135-Wt-Mm-Rep2: `Forward reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566_1.fastq.gz>`__ and 
                        `reverse reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/006/SRR3950566/SRR3950566_2.fastq.gz>`__ .  

CC-MB-E105-Wt-Mm-Rep1: `Forward reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559_1.fastq.gz>`__ and 
                        `reverse reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/009/SRR3950559/SRR3950559_2.fastq.gz>`__ .

CC-MB-E105-Wt-Mm-Rep2: `Forward reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560_1.fastq.gz>`__ and 
                        `reverse reads <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/000/SRR3950560/SRR3950560_2.fastq.gz>`__ .


Mapping
-------

Map the files with a mapper of your choice against the mm9 reference genome. 


Create cHi-C matrices
---------------------


Creation of reference point file
--------------------------------

Andrey et al writes they used in total 460 reference points, but 24 were removed due to low sequence coverage or non correspondence to a promoter region, leading to 446 in total.

To reproduce this, we need the all Viewpoint which are published in Supplementary Table `S2 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S2.xlsx>`__ and `S8 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S8.xlsx >`__.  


It is the easiest to create via Excel the reference point file in the following format and store it as a tab separated file:

chr1	4487435	4487435 Sox17

Or just download the prepared file. We will do the Quality control on our own and compare to the results of Andrey.



Quality control
^^^^^^^^^^^^^^^


Background model
^^^^^^^^^^^^^^^^


Viewpoint computation
^^^^^^^^^^^^^^^^^^^^^


Significant interactions detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Aggregate data for differential test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Differential test
^^^^^^^^^^^^^^^^^


Plotting of Viewpoints
^^^^^^^^^^^^^^^^^^^^^^