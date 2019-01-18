Importing and Exporting HiCExplorer data
========================================

Exporting HiCExplorer output to Bioconductor
--------------------------------------------

It's possible to export Hi-C Matrices produced by HiCExplorer to
`bioconductor <http://bioconductor.org/>`__ in R, which allows us to use
existing bioconductor infrastructure for differential Hi-C analysis. The
tool **hicConvertFormat** allows us to write Hi-C matrices in a format that can
easily be imported in bioconductor as **GInteractions** object. Below
is an example.

.. code:: bash

    ## Assuming HiCExplorer is installed in ~/programs
    hicConvertFormat --matrix ~/programs/HiCExplorer/test/test_data/Li_et_al_2015.h5 \ --inputFormat h5
    -o GInteration_example --outputFormat GInteractions

The output file is in tsv format. It looks like this :

+------+---------+---------+------+---------+---------+--------------+
| V1   | V2      | V3      | V4   | V5      | V6      | V7           |
+======+=========+=========+======+=========+=========+==============+
| X    | 19537   | 20701   | X    | 19537   | 20701   | 1054.47483   |
+------+---------+---------+------+---------+---------+--------------+
| X    | 19537   | 20701   | X    | 20701   | 22321   | 375.86990    |
+------+---------+---------+------+---------+---------+--------------+
| X    | 19537   | 20701   | X    | 22321   | 24083   | 222.53900    |
+------+---------+---------+------+---------+---------+--------------+
| X    | 19537   | 20701   | X    | 24083   | 25983   | 114.26340    |
+------+---------+---------+------+---------+---------+--------------+
| X    | 19537   | 20701   | X    | 25983   | 27619   | 95.87463     |
+------+---------+---------+------+---------+---------+--------------+

This file can now be loaded into R as a **GInteractions** object, as
shown below :

.. code:: r

    ## INSIDE R
    library(GenomicRanges)
    library(InteractionSet)

    hic <- read.delim("GInteraction_example.tsv", header = FALSE)

    # Converting data.frame to GInteraction
    convertToGI <- function(df){
                row.regions <- GRanges(df$V1, IRanges(df$V2,df$V3))# interaction start
                col.regions <- GRanges(df$V4, IRanges(df$V5,df$V6))# interaction end
                gi <- GInteractions(row.regions, col.regions)
                gi$norm.freq <- df$V7 # Interaction frequencies
                return(gi)
    }
                            }
    hic.gi <- convertToGI(hic)

Multiple files can be loaded, and converted to an **InteractionSet**
object. If you have prepared matrices using binning, the intervals in
the matrices must be the same. Therefore it's easy to merge these
matrices together in an InteractionSet object. In case some bins don't
match, we can merge the GInteraction objects based on matching bins, as
follows.

.. code:: r

    # assuming hic.gi is a list of two GInteration objects hic.gi1 and hic.gi2
    # hic.gi <- list(hic.gi1, hic.gi2)

    # Get common regions between the two objects
    combined <- unique(c(hic.gi$hic.gi1, hic.gi$hic.gi2))

    # replace original regions with the common regions
    replaceRegions(hic.gi$hic.gi1) <- regions(combined)
    replaceRegions(hic.gi$hic.gi2) <- regions(combined)

    # Get the matching indexes between the two objects
    matched <- lapply(hic.gi, function(x) {
                match(x, combined)
                })

    # Create a count matrix (for interaction frequencies)
    counts <- matrix(0, ncol = 2, nrow=length(combined)) # counts for unmatched bins set to zero

    # fill in the counts for matched bins
    counts[matched$hic.gi1,1] <- hic.gi$hic.gi1$norm.freq
    counts[matched$hic.gi2,2] <- hic.gi$hic.gi2$norm.freq

    # Finally, create the InteractionSet object
    iset <- InteractionSet(counts, combined)

InteractionSet objects can be used for packages like
`diffHic <https://www.bioconductor.org/packages/release/bioc/html/diffHic.html>`__,
for differential Hi-C analysis.

-  For more information on working with GInteraction and InteractionSet
   objects in bioconductor check out `this
   vignette <https://bioconductor.org/packages/devel/bioc/vignettes/InteractionSet/inst/doc/interactions.html>`__.
