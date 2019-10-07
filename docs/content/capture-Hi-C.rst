Captured Hi-C data analysis
===========================

How we use HiCExplorer to analyse cHi-C data
--------------------------------------------

This How-to is based on the published dataset from `Andrey et al. 2017 <https://doi.org/10.1101/gr.213066.116>`__. For the tutorial, we use the samples FL-E13.5 and MB-E-10.5. 


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

Map the files with a mapper of your choice against the mm9 reference genome; as an example, the mapping with bowtie2 is shown.

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

To create a cHi-C matrix we use HiCExplorer's hicBuildMatrix for each replicate separately and merge the replicates into a single matrix later. Like Andrey et al., we use a resolution of 1kb and use the restriction enzyme DpnII.

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

Andrey et al. state that they used a total of 460 reference points, but that 24 were removed due to low sequence coverage or non-correspondence to a promoter region, leading to 446 in total.

To reproduce this, we need all reference points which are published in Supplementary Table `S2 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S2.xlsx>`__ and `S8 <https://genome.cshlp.org/content/suppl/2017/01/20/gr.213066.116.DC1/Supplemental_Table_S8.xlsx>`__.

It is simplest to create the reference point file in the following format using Excel and store it as a tab separated file:

.. code:: bash

    chr1	4487435	4487435 Sox17

Otherwise, just `download <https://drive.google.com/open?id=14kNDI1xuEiP-8S5lssTbRm5d3s2eMNLU>`__ the prepared file. We will do the quality control on our own and compare with the results of Andrey et al.



Quality control
^^^^^^^^^^^^^^^

As a first step we compute the quality of each viewpoint by considering the sparsity. As soon as one viewpoint in one sample is less than the user-defined threshold (`--sparsity`), the reference point is no longer considered.

.. code:: bash

    chicQualityControl -m FL-E13-5.cool MB-E10-5.cool -rp reference_points.bed --sparsity 0.025 --threads 20

The quality control creates five files: two plots showing the sparsity structure of the samples and three files containing the accepted reference points, the rejected ones and one file with all viewpoints and their sparsity level per sample.

In our example the plots look like the following:

.. image:: ../images/chic/sparsity.png 

.. image:: ../images/chic/histogram.png


The first plot shows the sparsity per sample for each viewpoint, while the second one shows the sparsity distribution as a histogram. It can be seen quite clearly that only a minority of the samples are really sparse and therefore need to be removed. The red line indicates the chosen sparsity level.

The reference point `Tdap2b` at `chr1 19198995`, which has a sparsity of 0.018 in FL-E13-5 and 0.016 in MB-E10-5, is considered to be of bad quality. To confirm this result we plot the viewpoint:

.. image:: ../images/chic/Tfap2b_FL-E13-5_MB-E10-5_chr1_19198995_19198995.png

The plot shows there are effectively no interactions except with the reference point itself and confirm the point should be removed from the data.

The result of the quality control rejected 71 reference points as too sparse, but surprisingly the viewpoints rejected by Andrey et al. are accepted. An explanation for this could be that we only consider two samples and not all samples used by Andrey, and therefore we missed the bad quality of some viewpoints.

Download the data: `Filtered reference points <https://drive.google.com/open?id=1y3G1wJRBy0aZPQJ504N94jLE4jco2GAT>`__, `Quality control raw data <https://drive.google.com/open?id=1E0Ii-5QdZDco8NkEXb-rMoBCcGjYUfJg>`__ and `rejected reference points <https://drive.google.com/open?id=1LGDIoT7etslvHfNSPajYszaQlSsQegBx>`__.

Background model
^^^^^^^^^^^^^^^^

The background model computes all viewpoints given by the reference points for both samples in a range defined by the parameter `fixateRange`. We recommend setting it to 500kb because real interactions above the range 
are rarely observed and very low interaction numbers such as 1 are already considered to be significant. With this setting, only the interactions in a range 500kb up- and downstream of the reference point are considered for each viewpoint.
Based on this data, two background models are computed; the first one simply computes the average per relative distance to the reference point, and secondly, a negative binomial distribution per relative distance to
the reference point is fitted. This first model is used for filtering in the significant interaction evaluation by an x-fold factor and for plotting. The negative binomial model is more important; it is used to 
compute a p-value per relative distance in each sample, which is used to make the decision if an interaction is considered as significant.

.. code:: bash

    chicViewpointBackgroundModel -m FL-E13-5.cool MB-E10-5.cool --fixateRange 500000 -t 20 -rp reference_points.bed -o background_model.bed

The background model looks like this:

.. code:: bash

    Relative position	size nbinom     prob nbinom     max value       mean value
    -500000             75.895607451213	0.998528939430	2.333333333333	0.000101543771
    -499000	        90.348171762247	0.998725799952	2.750000000000	0.000104681360
    -498000	        78.512621775755	0.998514111424	2.800000000000	0.000106107536
    -497000	        75.706478185610	0.998327784087	3.800000000000	0.000116147819


You can `download <https://drive.google.com/open?id=1zblxEWa513LGwkjBknt83oZugg-uIJdw>`__ the background model.


Viewpoint computation
^^^^^^^^^^^^^^^^^^^^^

In this step the viewpoints for each reference point listed in a `reference_points.bed`-file is computed, using the quality controlled file created by `chicQualityControl`. The up- and downstream range can be given via `--range upstream downstream`. Please use the same value for `--fixateRange` as in the background model computation.
For each relative distance the x-fold over the average value of this relative distance is computed and each location is assigned a p-value based on the background negative binomial distribution for this relative distance.
For each viewpoint one viewpoint file is created and stored in the folder given by the parameter `--outputFolder`. 

.. code:: bash

    chicViewpoint -m FL-E13-5.cool MB-E10-5.cool --averageContactBin 5 --range 1000000 1000000 -rp referencePoints.bed -bmf background_model.bed --writeFileNamesToFile interactionFiles.txt --outputFolder interactionFilesFolder --fixateRange 500000 --threads 20


The name of each viewpoint file starts with its sample name (given by the name of the matrix), the
exact location and the gene / promoter name. For example, the viewpoint `chr1	4487435	4487435 Sox17` from `MB-E10-5.cool` matrix will be called `MB-E10-5_chr1_4487435_4487435_Sox17.bed` and looks like the following:

.. code:: text

    # Interaction file, created with HiCExplorer's chicViewpoint version 3.2
    # MB-E10-5.cool chr1_4487435_4487435    3.49  Mbp       5.49  Mbp       Sox17   Sum of interactions in fixate range: 978.0
    # Chromosome    Start   End     Gene    Sum of interactions     Relative position       Relative Interactions   p-value x-fold  Raw
    #
    chr1    3487000 3488000 Sox17   978.0   -1000000        0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3488000 3489000 Sox17   978.0   -999000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3489000 3490000 Sox17   978.0   -998000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3490000 3491000 Sox17   978.0   -997000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3491000 3492000 Sox17   978.0   -996000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3492000 3493000 Sox17   978.0   -995000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3493000 3494000 Sox17   978.0   -994000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3494000 3495000 Sox17   978.0   -993000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3495000 3496000 Sox17   978.0   -992000 0.000000000000  0.894286365313  0.000000000000  0.000000000000
    chr1    3496000 3497000 Sox17   978.0   -991000 0.000000000000  0.894286365313  0.000000000000  0.000000000000


Each file contains a header with information about the HiCExplorer version used, the sample, the viewpoint and the content of the different columns. 

If the parameter `--writeFileNamesToFile` is set, the viewpoint file names are written to a file which can be used for batch processing in the later steps. Each sample is combined with every other sample for each viewpoint to enable the batch processing
for the differential analysis. Example: matrices `FL-E13-5.cool` and  `MB-E10-5.cool` with the three reference points:

.. code:: bash

    FL-E13-5_chr1_4487435_4487435_Sox17.bed
    MB-E10-5_chr1_4487435_4487435_Sox17.bed
    FL-E13-5_chr1_14300280_14300280_Eya1.bed
    MB-E10-5_chr1_14300280_14300280_Eya1.bed
    FL-E13-5_chr1_19093103_19093103_Tfap2d.bed
    MB-E10-5_chr1_19093103_19093103_Tfap2d.bed




Significant interactions detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


To detect significant interactions and to prepare a target file for each viewpoint which will be used for the differential analysis, the script `chicSignificantInteractions` is used. It offers two modes: either the user can specify 
an x-fold value or a loose p-value. The first one considers all interactions with a minimum x-fold over the average background for its relative distribution as a candidate or secondly, all interactions with a loose p-value or less are considered. 
These are preselection steps to be able to detect wider peaks in the same way as sharp ones. All detected candidates are merged to one peak if they are direct neighbors, and the sum of all interaction values of this neighborhood
is used to compute a new p-value. The p-value is computed based on the relative distance negative binomial distribution of the interaction with the original highest interaction value. All peaks considered are accepted as significant interactions if
their p-value is as small as the threshold `--pvalue`.

To exclude interactions with an interaction value smaller than desired the parameter `--peakInteractionsThreshold` can be set.

In this example we use the reference point Mstn at location chr1 53118507, a loose p-value of 0.1 and p-value of 0.01:

.. code:: bash

    chicSignificantInteractions --interactionFile interactionFilesFolder/FL-E13-5_chr1_53118507_53118507_Mstn.bed -bmf background_model.bed --range 1000000 1000000 --pValue 0.01 --loosePValue 0.1


This creates two files: 

.. code:: bash

    FL-E13-5_chr1_53118507_53118507_Mstn_target.bed
    FL-E13-5_chr1_53118507_53118507_Mstn__significant_interactions.bed

These files are stored in the folders given by the parameters `--targetFolder` and `--outputFolder`.

The significant interaction files looks like the following:

.. code:: bash

    # FL-E13-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1517.0
    #Chromosome	Start	End	Gene	Sum of interactions	Relative position	Relative interactions	p-value	x-fold	Raw target
    chr1	53318000	53321000	Mstn	1517.0	200000	0.00395517468600000040	0.00000145009991170397	8.37043994897500098773	6.00000000000000000000
    chr1	53329000	53334000	Mstn	1517.0	212000	0.01081081081000000166	0.00000000000000188738	22.37661518795599846499	16.39999999999999857891
    chr1	53348000	53350000	Mstn	1517.0	231000	0.00329597890600000004	0.00001463968364323609	7.37204640642899988734	5.00000000000000000000
    chr1	53351000	53359000	Mstn	1517.0	239000	0.01437046802899999941	0.00000000000000099920	34.20972383882499912033	21.80000000000000071054
    chr1	53393000	53401000	Mstn	1517.0	278000	0.01793012524599999977	0.00000000000000044409	48.20518935066399990319	27.19999999999999928946
    chr1	53408000	53420000	Mstn	1517.0	294000	0.02307185234000000418	0.00000000000001743050	68.05162660125500906361	35.00000000000000000000


The target file looks like:

.. code:: bash

    # Significant interactions result file of HiCExplorer's chicSignificantInteractions version 3.2-dev
    # targetFolder/FL-E13-5_chr1_53118507_53118507_Mstn_target.bed
    # Mode: loose p-value with 0.1
    # Used p-value: 0.01
    #
    chr1	53318000	53321000
    chr1	53329000	53334000
    chr1	53348000	53359000
    chr1	53393000	53401000
    chr1	53408000	53420000


Batch mode
~~~~~~~~~~

The batch mode supports the computation of many viewpoints at once and is able to create one target list for the same viewpoint and two (or n) samples. To do the batch computation the 
parameter `--batchMode` needs to be added and the folder of the viewpoint files needs to be defined. In addition, the list of viewpoints created by `chicViewpoint` with `--writeFileNamesToFile` needs to be 
used as input. One target file is created for n consecutive lines and can be defined via the `--computeSampleNumber` parameter. However, for the differential test where the target file is needed, only 
two samples and one target file is supported.

.. code:: bash

    chicSignificantInteractions --interactionFile interactionFiles.txt --interactionFileFolder interactionFilesFolder/  -bmf background_model.bed --range 1000000 1000000 --pValue 0.01 --loosePValue 0.1 --batchMode

The output is: 

- A folder containing all target files, name defined by `--targetFolder`, default value is `targetFolder`
- A folder with all significant interaction files, name defined by `--outputFolder`, default value is `significantFiles`
- A list containing the file names of all target files, name defined by `--targetFileList`, default value is `targetList.txt`
- A list containing the file names of all significant interaction files, name defined by `--writeFileNamesToFile`, default value is `significantFilesBatch.txt`


Aggregate data for differential test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The process to aggregate data is only necessary if the differential test is used. Either two files and one target list are used to generate the files for the differential test
or the batch mode can be used. `chicAggregateStatistic` takes the created viewpoint files from `chicViewpoint` as input and either the target files per two samples created by `chicSignificantInteractions`
or one target file which applies for all viewpoints. 

.. code:: bash

    chicAggregateStatistic --interactionFile interactionFilesFolder/FL-E13-5_chr1_53118507_53118507_Mstn.bed interactionFilesFolder/MB-E10-5_chr1_53118507_53118507_Mstn.bed --targetFile targetFolder/FL-E13-5_MB-E10-5_chr1_53118507_53118507_Mstn_target.bed

It selects the original data based on the target locations and returns one file per sample which is used for the differential test.

Batch mode
~~~~~~~~~~

In the batch mode the interaction file is the file containing the viewpoint file names, the folder needs to be defined by `--interactionFileFolder`, the same applies to the target file and folder.
Two viewpoint files are match with one target file created by `chicSignificantInteractions` or one target file for all viewpoints. Please notice the output files are written to the folder name
defined by `--outputFolder`, the default is `aggregatedFiles` and it is recommended to write the file names for further batch processing with `hicDifferentialTest` to `--writeFileNamesToFile`. All output files
get the suffix defined by `--outFileNameSuffix`, default value is `_aggregate_target.bed`.

.. code:: bash

     chicAggregateStatistic --interactionFile interactionFiles.txt --interactionFileFolder interactionFilesFolder --targetFile targetList.txt --targetFileFolder targetFolder --batchMode



Differential test
^^^^^^^^^^^^^^^^^

The differential test tests the interaction value of the reference point and the interaction value of the target of two samples for a differential expression. To achieve this,
either Fisher's test or the chi-squared test can be used. H0 is defined as 'both locations are equal', meaning the differential expressed targets can be found in the H0 rejected file.

This can be computed per sample: 


.. code:: bash

    chicDifferentialTest --interactionFile aggregatedFiles/FL-E13-5_chr1_53118507_53118507_Mstn__aggregate_target.bed aggregatedFiles/MB-E10-5_chr1_53118507_53118507_Mstn__aggregate_target.bed --alpha 0.05 --statisticTest fisher

Or via batch mode:

.. code:: bash

    chicDifferentialTest --interactionFile aggregatedFilesBatch.txt --interactionFileFolder aggregatedFiles --alpha 0.05 --statisticTest fisher --batchMode --threads 20

In both cases it is important to set the desired alpha value and the output is written to `--outputFolder` (default `differentialResults`). For each sample three files are created:

- H0 rejected targets 
- H0 accepted targets 
- one file containing both

In the batch mode, the file `--rejectedFileNamesToFile` is also written and contains the file names of the rejected files. This can be used for the batch processing mode of `chicPlotViewpoint`.

.. code:: bash

    # Differential analysis result file of HiCExplorer's chicDifferentialTest version 3.2-dev
    # This file contains the p-values computed by fisher test
    # To test the smoothed (float) values they were rounded up to the next integer
    #
    # Alpha level 0.05
    # Degrees of freedom 1
    #
    # FL-E13-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1517.0
    # MB-E10-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1670.0
    #Chromosome	Start	End	Gene	Relative distance	sum of interactions 1	target_1 raw	sum of interactions 2	target_2 raw	p-value
    chr1	53089000	53091000	Mstn	-28000	1517.00000	5.00000	1670.00000	10.40000		0.21800
    chr1	53131000	53133000	Mstn	14000	1517.00000	18.20000	1670.00000	23.60000		0.75900
    chr1	53156000	53158000	Mstn	39000	1517.00000	3.00000	1670.00000	10.80000		0.06117
    chr1	53251000	53254000	Mstn	135000	1517.00000	4.00000	1670.00000	9.60000		0.18614
    chr1	53287000	53291000	Mstn	172000	1517.00000	7.20000	1670.00000	15.00000		0.29506
    chr1	53305000	53309000	Mstn	190000	1517.00000	6.20000	1670.00000	12.40000		0.36952
    chr1	53318000	53321000	Mstn	202000	1517.00000	6.00000	1670.00000	3.20000		0.53309
    chr1	53326000	53334000	Mstn	215000	1517.00000	25.80000	1670.00000	22.60000		0.47374
    chr1	53346000	53359000	Mstn	240000	1517.00000	31.60000	1670.00000	22.20000		0.13464
    chr1	53408000	53421000	Mstn	302000	1517.00000	36.40000	1670.00000	28.20000		0.21290


.. code:: bash

    # Differential analysis result file of HiCExplorer's chicDifferentialTest version 3.2-dev
    # This file contains the p-values computed by fisher test
    # To test the smoothed (float) values they were rounded up to the next integer
    #
    # Alpha level 0.05
    # Degrees of freedom 1
    #
    # FL-E13-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1517.0
    # MB-E10-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1670.0
    #Chromosome	Start	End	Gene	Relative distance	sum of interactions 1	target_1 raw	sum of interactions 2	target_2 raw	p-value
    chr1	53393000	53401000	Mstn	282000	1517.00000	27.20000	1670.00000	6.40000		0.00012

.. code:: bash

    # Differential analysis result file of HiCExplorer's chicDifferentialTest version 3.2-dev
    # This file contains the p-values computed by fisher test
    # To test the smoothed (float) values they were rounded up to the next integer
    #
    # Alpha level 0.05
    # Degrees of freedom 1
    #
    # FL-E13-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1517.0
    # MB-E10-5.cool	chr1_53118507_53118507	52.12  Mbp	54.12  Mbp	Mstn	Sum of interactions in fixate range: 1670.0
    #Chromosome	Start	End	Gene	Relative distance	sum of interactions 1	target_1 raw	sum of interactions 2	target_2 raw	p-value
    chr1	53089000	53091000	Mstn	-28000	1517.00000	5.00000	1670.00000	10.40000		0.21800
    chr1	53131000	53133000	Mstn	14000	1517.00000	18.20000	1670.00000	23.60000		0.75900
    chr1	53156000	53158000	Mstn	39000	1517.00000	3.00000	1670.00000	10.80000		0.06117
    chr1	53251000	53254000	Mstn	135000	1517.00000	4.00000	1670.00000	9.60000		0.18614
    chr1	53287000	53291000	Mstn	172000	1517.00000	7.20000	1670.00000	15.00000		0.29506
    chr1	53305000	53309000	Mstn	190000	1517.00000	6.20000	1670.00000	12.40000		0.36952
    chr1	53318000	53321000	Mstn	202000	1517.00000	6.00000	1670.00000	3.20000		0.53309
    chr1	53326000	53334000	Mstn	215000	1517.00000	25.80000	1670.00000	22.60000		0.47374
    chr1	53346000	53359000	Mstn	240000	1517.00000	31.60000	1670.00000	22.20000		0.13464
    chr1	53393000	53401000	Mstn	282000	1517.00000	27.20000	1670.00000	6.40000		0.00012
    chr1	53408000	53421000	Mstn	302000	1517.00000	36.40000	1670.00000	28.20000		0.21290



Plotting of Viewpoints
^^^^^^^^^^^^^^^^^^^^^^

`chicPlotViewpoint` can plot `n` viewpoints in one plot, add the mean background, show the p-value per relative distance per sample as an additional heatmap bar and highlight significant interactions or differential expressed regions.

One viewpoint:

.. code:: bash

    chicPlotViewpoint --interactionFile interactionFilesFolder/FL-E13-5_chr1_53118507_53118507_Mstn.bed --range 200000 200000 -o single_plot.png

.. image:: ../images/chic/single_plot.png

Two viewpoints, background, differential expression and p-values:

.. code:: bash

    chicPlotViewpoint --interactionFile interactionFilesFolder/FL-E13-5_chr1_53118507_53118507_Mstn.bed interactionFilesFolder/MB-E10-5_chr1_53118507_53118507_Mstn.bed --range 300000 300000 --pValue --differentialTestResult differentialResults/FL-E13-5_MB-E10-5_chr1_53118507_53118507_Mstn__H0_rejected.bed --backgroundModelFile background_model.bed -o differential_background_pvalue.png


.. image:: ../images/chic/differential_background_pvalue.png

Two viewpoints, background, significant interactions and p-values:

.. code:: bash

    chicPlotViewpoint --interactionFile interactionFilesFolder/FL-E13-5_chr1_53118507_53118507_Mstn.bed interactionFilesFolder/MB-E10-5_chr1_53118507_53118507_Mstn.bed --range 300000 300000 --pValue --significantInteractions significantFiles/FL-E13-5_chr1_53118507_53118507_Mstn__significant_interactions.bed significantFiles/MB-E10-5_chr1_53118507_53118507_Mstn__significant_interactions.bed --backgroundModelFile background_model.bed -o significant_background_pvalue.png

.. image:: ../images/chic/significant_background_pvalue.png


Batch mode
~~~~~~~~~~

The batch mode is able to plot files under the same parameter setting for multiple viewpoints. These viewpoints are given by the file created by `chicViewpoint` with `--writeFileNamesToFile` parameter.
The number of consecutive lines which should be plotted to one image can be defined using `--plotSampleNumber`. If the differentially expressed regions should be highlighted, setting this parameter to 2 is recommended.

For all modes the principle of a file containing the file names and a folder containing them applies for the plotting too, and using all cores available is highly recommended.

.. code:: bash

    chicPlotViewpoint --interactionFile interactionFiles.txt --interactionFileFolder interactionFilesFolder/ --range 300000 300000 --pValue --significantInteractions significantFilesBatch.txt --significantInteractionFileFolder significantFiles --backgroundModelFile background_model.bed --outputFolder plots --threads 20 --batchMode