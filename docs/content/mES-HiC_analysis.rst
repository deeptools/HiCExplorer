Hi-C analysis of mouse ESCs using HiC-Explorer
==============================================

The following example shows how we can use HiCExplorer to analyse a
published dataset. Here we are using a HiC dataset from `Marks et. al.
2015 <http://www.genomebiology.com/2015/16/1/149>`__, on mouse ESCs.

**Protocol** Collection of the cells for Hi-C and the Hi-C sample
preparation procedure was performed as previously described
`Lieberman-Aiden et
al. <http://www.sciencemag.org/content/326/5950/289.long>`__, with the
slight modification that *DpnII* was used as restriction enzyme during
initial digestion. Paired-end libraries were prepared according to
Lieberman-Aiden et al. and sequenced on the NextSeq 500 platform using 2
× 75 bp sequencing.

Prepare for analysis
--------------------


Download Raw fastq files
~~~~~~~~~~~~~~~~~~~~~~~~

Fastq files can be downloaded from the EBI archive (or NCBI archive). We will store the files in the directory *original_data*.

.. code:: bash

    mkdir original_data

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/007/SRR1956527/SRR1956527_1.fastq.gz -O original_data/SRR1956527_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/007/SRR1956527/SRR1956527_2.fastq.gz -O original_data/SRR1956527_2.fastq.gz

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1956528/SRR1956528_1.fastq.gz -O original_data/SRR1956528_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1956528/SRR1956528_2.fastq.gz -O original_data/SRR1956528_2.fastq.gz

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/009/SRR1956529/SRR1956529_1.fastq.gz -O original_data/SRR1956529_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/009/SRR1956529/SRR1956529_2.fastq.gz -O original_data/SRR1956529_2.fastq.gz


We start with creating an index for our alignment software for the *GRCm38/mm10* genome. 
As a source we use the mm10 genome from `UCSC <http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/bigZips/>`__

.. code:: bash

    mkdir genome_mm10
    wget http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz -O genome_mm10/chromFa.tar.gz
    tar -xvzf genome_mm10/chromFa.tar.gz
    cat genome_mm10/*.fa > genome_mm10/mm10.fa

We have the mm10 genome stored in one fasta file and can build the index. We tried it sucessfuly with hisat2, bowtie2 and bwa. Run the mapping 
with one of them and do not mix them!

Create an index 
~~~~~~~~~~~~~~~

hisat2
^^^^^^^

.. code:: bash

    hisat2-build -p 8 genome_mm10/mm10.fa hisat2/mm10_index 

You can find more information about `hisat <https://ccb.jhu.edu/software/hisat2/manual.shtml>`__

bowtie2
^^^^^^^^
.. code:: bash

    bowtie2-build genome_mm10/mm10.fa bowtie2/mm10_index --threads 8

You can find more information about `bowtie <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__

bwa
^^^^
.. code:: bash

    bwa index -p bwa/mm10_index genome_mm10/mm10.fa

You can find more information about `bwa <http://bio-bwa.sourceforge.net/bwa.shtml>`__


Mapping the RAW files
~~~~~~~~~~~~~~~~~~~~~

Mates have to be mapped individually to avoid mapper specific heuristics designed
for standard paired-end libraries.

We have used the HiCExplorer sucessfuly with `bwa`, `bowtie2` and `hisat2`. However, it is important to:

 * for either `bowtie2` or `hisat2` use the `--reorder` parameter which tells bowtie2 or hisat2 to output
   the *sam* files in the **exact** same order as in the *.fastq* files.
 * use local mapping, in contrast to end-to-end. A fraction of Hi-C reads are chimeric and will not map end-to-end
   thus, local mapping is important to increase the number of mapped reads.
 * Tune the aligner parameters to penalize deletions and insertions. This is important to avoid aligned reads with
   gaps if they happen to be chimeric.


hisat2
^^^^^^

.. code:: bash

    hisat2 -x hisat2/mm10_index --threads 8 ../original_data/SRR1956527_1.fastq --reorder -S SRR1956527_1.sam 
    hisat2 -x hisat2/mm10_index --threads 8 ../original_data/SRR1956527_2.fastq --reorder -S SRR1956527_2.sam 
    hisat2 -x hisat2/mm10_index --threads 8 ../original_data/SRR1956528_1.fastq --reorder -S SRR1956528_1.sam 
    hisat2 -x hisat2/mm10_index --threads 8 ../original_data/SRR1956528_2.fastq --reorder -S SRR1956528_2.sam 
    hisat2 -x hisat2/mm10_index --threads 8 ../original_data/SRR1956529_1.fastq --reorder -S SRR1956529_1.sam 
    hisat2 -x hisat2/mm10_index --threads 8 ../original_data/SRR1956529_2.fastq --reorder -S SRR1956529_2.sam 
    
bowtie2
^^^^^^^^

.. code:: bash

    bowtie2 -x bowtie2/mm10_index --threads 8 ../original_data/SRR1956527_1.fastq --reorder -S SRR1956527_1.sam 
    bowtie2 -x bowtie2/mm10_index --threads 8 ../original_data/SRR1956527_2.fastq --reorder -S SRR1956527_2.sam 
    bowtie2 -x bowtie2/mm10_index --threads 8 ../original_data/SRR1956528_1.fastq --reorder -S SRR1956528_1.sam 
    bowtie2 -x bowtie2/mm10_index --threads 8 ../original_data/SRR1956528_2.fastq --reorder -S SRR1956528_2.sam 
    bowtie2 -x bowtie2/mm10_index --threads 8 ../original_data/SRR1956529_1.fastq --reorder -S SRR1956529_1.sam 
    bowtie2 -x bowtie2/mm10_index --threads 8 ../original_data/SRR1956529_2.fastq --reorder -S SRR1956529_2.sam 
    

.. code:: bash

    bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 bwa/mm10_index original_data/SRR1956527_1.fastq > SRR1956527_1.sam 
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 bwa/mm10_index original_data/SRR1956527_2.fastq > SRR1956527_2.sam 
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 bwa/mm10_index original_data/SRR1956528_1.fastq > SRR1956528_1.sam 
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 bwa/mm10_index original_data/SRR1956528_2.fastq > SRR1956528_2.sam 
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 bwa/mm10_index original_data/SRR1956529_1.fastq > SRR1956529_1.sam 
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 bwa/mm10_index original_data/SRR1956529_2.fastq > SRR1956529_2.sam 
    


Build, visualize and correct Hi-C matrix
----------------------------------------

Create a Hi-C matrix using the aligned files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
the restriction site bed file [-rs],
restriction sequence [-seq]
#-rs dpnII_positions_GRCm37-sorted.bed

Build Hi-C matrix
^^^^^^^^^^^^^^^^^

:ref:`hicBuildMatrix` builds the matrix of read counts over the bins in the
genome, considering the sites around the given restriction site. We need
to provide the input BAM/SAM files,  binsize [-bs], restriction sequence [-seq],
name of output matrix file
[-o] and the name of output bam file (which contains the accepted
alignments) [-b].

To increase the computation speed, please set the system environment variable `HICEXPLORER_FILE_BUFFER_DIR` to a RAM disk like `/dev/shm`.
Be careful: This will consume a sustainable amount of memory, for 8 threads and an input buffer size of 400000 it is recommended to use a 
system with at least 32 GB of RAM. If you have a system with lower specifications, decrease the inputBufferSize. It is recommended to not use 
than 100000. If the memory is still not enough, use a directory on your local hard drive.

.. code:: bash

    export HICEXPLORER_FILE_BUFFER_DIR=/dev/shm


.. code:: bash

    mkdir hicMatrix
    hicBuildMatrix -s SRR1956527_1.sam SRR1956527_2.sam -bs 1000 -seq GATC -b SRR1956527_ref.bam --outFileName hicMatrix/SRR1956527.h5 --QCfolder hicMatrix/SRR1956527_QC --threads 8 --inputBufferSize 400000
    hicBuildMatrix -s SRR1956528_1.sam SRR1956528_2.sam -bs 1000 -seq GATC -b SRR1956528_ref.bam --outFileName hicMatrix/SRR1956528.h5 --QCfolder hicMatrix/SRR1956528_QC --threads 8 --inputBufferSize 400000
    hicBuildMatrix -s SRR1956529_1.sam SRR1956529_2.sam -bs 1000 -seq GATC -b SRR1956529_ref.bam --outFileName hicMatrix/SRR1956529.h5 --QCfolder hicMatrix/SRR1956529_QC --threads 8 --inputBufferSize 400000
    

.. code:: bash

    mkdir hiCmatrix

    for SRR in SRR1956527 SRR1956528 SRR1956529;
    do hicBuildMatrix \
       -s mapped_files/${SRR}_1.bam mapped_files/${SRR}_2.bam \
       -bs 10000 -seq GATC \
       -b ${SRR}_ref.bam -o hiCmatrix/${SRR}.matrix --QCfolder hiCmatrix/${SRR}_QC  --threads 8 --inputBufferSize 400000 & done

The output bam files show that we have around 34M, 54M and 58M selected
reads for SRR1956527, SRR1956528 & SRR1956529, respectively. Normally
25% of the total reads are selected.

The output matrices have counts for the genomic regions. The extension
of output matrix files is *.h5*.


Merge (sum) matrices from replicates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To increase the depth of reads we merge the counts from these three
replicates.

.. code:: bash

    hicSumMatrices -m hiCmatrix/SRR1956527.matrix.h5 hiCmatrix/SRR1956528.matrix.h5 \
                      hiCmatrix/SRR1956529.matrix.h5 -o hiCmatrix/replicateMerged.matrix.h5

Correct Hi-C Matrix
^^^^^^^^^^^^^^^^^^^

:ref:`hicCorrectMatrix` corrects the matrix counts in an iterative manner.
For correcting the matrix, it's important to remove the unassembled
scaffolds (eg NT\_) and keep only chromosomes, as scaffolds create
problems with matrix correction. Therefore we use the chromosome names
(1-19, X, Y) here. **Important** use 'chr1 chr2 chr3 etc.' if your genome index uses
chromosome names with the 'chr' prefix.

Matrix correction works in two steps: first a histogram containing the sum of contact per bin (row sum) is
produced. This plot needs to be inspected to decide the best threshold for removing bins with lower number of reads. The
second steps removes the low scoring bins and does the correction.

(1-19, X, Y) variant:

.. code:: bash

    hicCorrectMatrix diagnostic_plot \
    --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    -m hiCmatrix/replicateMerged.matrix.npz -o hiCmatrix/diagnostic_plot.png

(chr1-ch19, chrX, chrY) variant:

.. code:: bash

    hicCorrectMatrix diagnostic_plot \
    --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \
    -m hiCmatrix/replicateMerged.matrix.npz -o hiCmatrix/diagnostic_plot.png

The output of the program prints a threshold suggestion that is usually accurate but is better to
revise the histogram plot. See :ref:`example_usage` for an example and for more info.

Next we do the correction using the best filter threshold.

.. code:: bash

    hicCorrectMatrix correct \
    --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    --filterThreshold -1.5 10 \
    -m hiCmatrix/replicateMerged.matrix.npz -o hiCmatrix/replicateMerged.Corrected.npz



.. code:: bash

    hicCorrectMatrix correct \
    --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \
    --filterThreshold -1.5 10 \
    -m hiCmatrix/replicateMerged.matrix.npz -o hiCmatrix/replicateMerged.Corrected.npz



Plot Hi-C matrix
~~~~~~~~~~~~~~~~

A 10kb bin matrix is quite large to plot and is better to reduce the resolution (to know the size
of a Hi-C matrix use the tool :ref:`hicInfo`). For this we use the tool :ref:`hicMergeMatrixBins`

Merge matrix bins for plotting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:ref:`hicMergeMatrixBins` merges the bins into larger bins of given number
(specified by -nb). We will merge 100 bins in the original (uncorrected) matrix and
then correct it. The new bin size is going to be 10.000 bp * 100 = 1.000.000 bp

.. code:: bash

    hicMergeMatrixBins \
    -m hiCmatrix/replicateMerged.matrix.npz -nb 100 \
    -o hiCmatrix/replicateMerged.matrix-100bins.npz

Correct the merged matrix
^^^^^^^^^^^^^^^^^^^^^^^^^

We will now correct the merged matrix before plotting.

.. code:: bash

    hicCorrectMatrix diagnostic_plot \
    --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    -m hiCmatrix/replicateMerged.matrix-100bins.npz -o hiCmatrix/diagnostic_plot_100bins.png

    hicCorrectMatrix correct \
    --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    --filterThreshold 0.9 10 \
    -m hiCmatrix/replicateMerged.matrix-100bins.npz -o hiCmatrix/replicateMerged.Corrected-100bins.npz


.. code:: bash

    hicCorrectMatrix diagnostic_plot \
    --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \
    -m hiCmatrix/replicateMerged.matrix-100bins.npz -o hiCmatrix/diagnostic_plot_100bins.png

    hicCorrectMatrix correct \
    --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \
    --filterThreshold 0.9 10 \
    -m hiCmatrix/replicateMerged.matrix-100bins.npz -o hiCmatrix/replicateMerged.Corrected-100bins.npz



Plot the corrected Hi-C Matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**hicPlotMatrix** can plot the merged matrix. We use options :
**--log1p** to plot the log intensites and **dpi** in increase image
resolution

.. code:: bash

    mkdir plots
    hicPlotMatrix \
    --log1p --dpi 300 \
    -m hiCmatrix/replicateMerged.Corrected-100bins.npz \
    --clearMaskedBins \
    -o plots/replicateMerged_Corrected-100bins_plot.png

.. figure:: ./plots/replicateMerged_Corrected-100bins_plot.png
   :alt: corrected\_100kb\_plot

   corrected\_100kb\_plot

Remove outliers from hic-Matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Outliers can be removed by a cutoff after looking at the diagnostic plot
for :ref:`hicCorrectMatrix` (using **diagnostic\_plot** option). Here we
are using a matrix with 20kb bins (produced by *hicMergeMatrixBins -nb
2*), since 20kb seems to be decent resolution.

Select threshold for outlier removal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following is the disgnostic plot that shows a bimodal distribution. We
should remove the values from both lower and upper end of the
distribution.

.. code:: bash

    hicCorrectMatrix diagnostic_plot -m hiCmatrix/replicateMerged.matrix_20kb.npz -o plots/diagPlot-20kb.png

.. figure:: ./plots/diagPlot-20kb.png
   :alt: diagplot

   diagplot

Correct matrix removing outliers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Looking at the above distribution, we can select the value of -2 (lower
end) and 3 (upper end) to remove. This is given by the **-t** option in
hicCorrectMatrix.

.. code:: bash

    hicCorrectMatrix correct \
    --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    -m hiCmatrix/replicateMerged.matrix_20kb.npz \
    -t -2 3 --perchr -o hiCmatrix/replicateMerged.Corrected_20kb.npz


.. code:: bash

    hicCorrectMatrix correct \
    --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \
    -m hiCmatrix/replicateMerged.matrix_20kb.npz \
    -t -2 3 --perchr -o hiCmatrix/replicateMerged.Corrected_20kb.npz


Plot corrected matrix
^^^^^^^^^^^^^^^^^^^^^

We can now plot the one of the chromosomes (eg. chromosome X) , with the
corrected matrix.

.. code:: bash

    hicPlotMatrix \
    --log1p --dpi 300 \
    -m hiCmatrix/replicateMerged.Corrected_20kb.npz \
    --region X -t "Corrected Hi-C matrix for mESC : chrX" \
    -o plots/replicateMerged_Corrected-20kb_plot-chrX.png

.. figure:: ./plots/replicateMerged_Corrected-20kb_plot-chrX.png
   :alt: correctMatrixPlot

   correctMatrixPlot

Find and plot TADs
------------------

Find TADs
~~~~~~~~~

To call TADs a corrected matrix is needed.
TAD calling works in two steps: First HiCExplorer computes a TAD-separation score based on a z-score matrix for
all bins. Then those bins having a local minimum of the TAD-separation score are evaluated with respect to the
surrounding bins to decide assign a p-value. Then a cutoff is applied to select the bins more likely to be TAD
boundaries.

:ref:`hicFindTADs` tries to identify sensible parameters but those can be change to identify more stringent set of
boundaries.

.. code-block:: bash

    mkdir TADs
    hicFindTADs -m hiCmatrix/replicateMerged.Corrected_20kb.npz \
    --minDepth 40000 --maxDepth 120000 --numberOfProcessors 20 --step 20000 \
    --outPrefix TADs/marks_et-al_TADs_20kb-Bins  --minBoundaryDistance 80000 \ # reduce noise by looking at min 80kb steps
    --pvalue 0.05


As an output we get the boundaries and domains as separate bed files.

Plot TADs
~~~~~~~~~

Build Tracks File
^^^^^^^^^^^^^^^^^

We can plot the TADs for a given chromosomal region. For this we need to
create a tracks file containing the instructions to build the plot. The
:doc:`tools/hicPlotTADs` documentation contains the instructions to build the track file.
A small example of a track file is:

.. code-block:: INI

   [x-axis]

   [hic track]
   file = hic.npz
   title = Hi-C
   colormap = RdYlBu_r
   depth = 1000000
   transform = log1p

   [genes]
   file = genes.bed
   title = genes
   color = darkblue
   width = 5
   type = genes


Plot
^^^^

Here I am plotting the TADs we have found (using 20kb bins) along with
the TADs found by Marks et. al., available as bed file
`here <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1652666>`__
and GRCm37\_genes.bed file (from ensembl).

.. code:: bash

    hicPlotTADs --tracks tracks_toPlot/tracks_2.txt \
    --region X:99974316-101359967 --dpi 300 \
    -out plots/marks_et-al_TADs.png -t "Marks et. al. TADs on X"

.. figure:: ./plots/marks_et-al_TADs.png
   :alt: TADplot

   TADplot

Comparing Marks et. al. and Dixon et. al.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We analysed the mESC Hi-C data from `Dixon et.
al <http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html>`__
using Hi-C explorer, and compared it to Marks et. al. dataset. For this
we mapped the reads using bowtie and prepared 20kb matrices. Following
is the plot showing the TADs on the X chromosomes, at 1.2 MB region
around Xist (the X Inactivation Center).

We have plotted here the Hi-C tracks from both the studies, containing
TADs as triangles, detected by Hi-C explorer, along with the boundaries
as bed files provided with the studies, normalized CTCF signal from
ENCODE, spectrum of Hi-C signal produced by *hicFindTADs*, and a
genes.bed file from ensembl.

.. figure:: ./plots/Marks-Dixon_TADs.png
   :alt: TADplot2

   TADplot2
