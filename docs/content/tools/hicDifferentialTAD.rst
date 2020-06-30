.. _hicDifferentialTAD:

hicDifferentialTAD
==================

.. argparse::
   :ref: hicexplorer.hicDifferentialTAD.parse_arguments
   :prog: hicDifferentialTAD

hicDifferentialTAD computes with a treatment Hi-C file, a control Hi-C file and a precomputed TAD domains file if the detected TADs are differential between the treatment and the control sample.
The TADs need to be precomputed on the treatment file with _hicFindTADs.

hicDifferentialTAD extract per TAD three regions: the intra-TAD, the left and right inter-TAD region. In the following image, the upper visualization shows a region with the
detected TADs which are indicated by the black lines. The bottom shows as an example which regions are used for the differential test: the intra-TAD region is highlighted in red,
the left inter-TAD in sandy-color and the right inter-TAD in blue. Between two samples a Wilcoxon rank-sum test is applied for a TAD under H0: 'The regions are equal'.
For all three regions of a TAD the rank-sum test is independently applied. The user has the choice with the two parameters 'mode' and 'modeReject' to define if a) all three
regions should be considered ('all'), or only a subset e.g. the intra-TAD or intra-TAD and left inter-TAD should be considered; and b) if all regions need to have lower p-value than the 
user given to reject H0 or if it is enough that at least one region is rejecting H0 to consider the region as differential.

.. image:: ../../images/hicDifferentialTAD_example1.png

Example usage
--------------

.. code:: bash

    $ hicDifferentialTAD -tm GSM2644945_Untreated-R1.100000_chr1.cool -cm GSM2644947_Auxin2days-R1.100000_chr1.cool -td untreated_R1_domains.bed -o differential -p 0.01 -t 4 -mr all


In this example data from Nora et al. "Targeted Degradation of CTCF Decouples Local Insulation of Chromosome Domains from Genomic Compartmentalization", Cell 169.5 (2017): 930-944 is used [GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98671]

.. image:: ../../images/hicDifferentialTAD.png


The output are two BED-similar files: the '_accepted.diff_tad' and '_rejected.diff_tad' file. The difference to a real BED file is a) the usage of a header starting with '#', the first six columns are BED6 standard, however, there are three additional columns with the p-values of each intra-TAD, left-inter-TAD and right-inter-TAD test.
The score value and name is copied from the domains.bed file which is a output of hicFindTADs.

.. code:: bash

    # Created with HiCExplorer's hicDifferentialTAD version 3.5-dev
    # H0 'regions are equal' H0 is rejected for all p-value smaller or equal the user given p-value threshold; i.e. regions in this file are considered as differential.
    # Rejected regions with Wilcoxon rank-sum test to p-value: 0.01  with used mode: all and modeReject: all 
    # Chromosome	start	end	name	score	strand	p-value left-inter-TAD	p-value right-inter-TAD	p-value intra-TAD
    chr1	17100000	18100000	ID_0.01_10	-0.230354	.	0.001796788048357949	5.061883070746839e-16	0.0
    chr1	18100000	19100000	ID_0.01_11	-0.5135365	.	0.00012684365278732503	6.438073769459897e-12	0.0
    chr1	19100000	20600000	ID_0.01_12	-0.4202525	.	0.0007752789863786969	4.48832548055169e-24	0.0
    chr1	20600000	21300000	ID_0.01_13	-0.9980865	.	0.006356145543682152	1.8649307626885386e-06	0.0
    chr1	21300000	22600000	ID_0.01_14	-0.859866	.	0.0014703570719553564	4.71767085913077e-20	0.0
    chr1	22600000	24300000	ID_0.01_15	-0.281294	.	5.271193747664565e-11	5.338337394768526e-16	0.0
    chr1	24300000	25800000	ID_0.01_16	-1.0366175	.	1.8337145545266618e-09	1.4870790186599747e-26	0.0
    chr1	25800000	30700000	ID_0.01_17	-0.1471245	.	2.9643530805137102e-12	1.5436334056034135e-174	0.0
    chr1	30700000	33700000	ID_0.01_18	-1.0134455	.	7.999079859441987e-34	6.409108547496729e-103	0.0
    chr1	33700000	34600000	ID_0.01_19	-1.33809	.	2.565221614708799e-07	2.9731582515337764e-18	0.0
    chr1	34600000	36400000	ID_0.01_20	-0.6325015	.	7.07925502998078e-14	2.1170094214663985e-31	0.0
    chr1	38100000	39000000	ID_0.01_22	-0.911804	.	5.0884800973418685e-08	5.205722443509364e-06	0.0
    chr1	39600000	40900000	ID_0.01_24	-0.483398	.	8.317208409235228e-07	9.70932687854428e-07	0.0


The here used data can be downloaded from out test data: github link