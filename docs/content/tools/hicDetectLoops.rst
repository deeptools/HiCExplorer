.. _hicDetectLoops:

hicDetectLoops
===============

hicDetectLoops can detect enriched interaction regions (peaks / loops) based on a strict candidate selection, negative binomial distributions 
and Wilcoxon rank-sum tests. 

The algorithm was mainly develop on GM12878 cells from Rao 2014 on 10kb and 5kb fixed bin size resolution. 

Example usage
--------------

.. code:: bash

    $ hicDetectLoops -m matrix.cool -o loops.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05


The candidate selection is based on the restriction of the maximum genomic distance, here 2MB. This distance is given by Rao 2014. For each genomic distance 
a continuous negative binomial distribution is computed and only interaction pairs with a threshold less than ``--pValuePreselection`` are accepted. 
In a second step, each candidate is considered compared to its neighborhood. This neighborhood is defined by the ``--windowSize`` parameter in the x and y dimension.
Per neighborhood only one candidate is considered, therefore only the candidate with the highest peak values is accepted. As a last step,
the neighborhood is split into a peak and background region (parameter ``--peakWidth``). The peakWidth can never be larger than the windowSize. However, we recommend 
for 10kb matrices a windowSize of 10 and a peakWidth of 6.

With version 3.5 a major revision of this tool was published. The biggest changes are: 

- The introduction of an observed/expected matrix as one of the first steps. Based on it, the other calculations compute without a distance dependent factor and the results are more accurate.
- The testing of peak vs background region with the donut layout as proposed by HiCCUPS with Wilcoxon rank-sum test. Anderson-Darling test is removed.
- Improving the handling of the parallelization and rewrote of the merging of candidates in one neighborhood. Results in faster execution time and less memory demand.
- Loading only the interactions within the range of `maxLoopDistance`. This is possible with HiCMatrix version 13 and results in faster load time and a reduced memory peak. This improvement is only for cool matrices, h5 matrices do not profit from these changes.


hicDetectLoops has many parameters and it can be quite difficult to search for the best parameter setting for your data. With version 3.5 we introduce therefore two new tools :doc:`hicHyperoptDetectLoops` and :doc:`hicHyperoptDetectLoopsHiCCUPS`.
The first one searches for the optimal parameter setting for your data based on HiCExplorer's hicDetectLoops. However, if you want to compare the results to Juicer HiCCUPS, the second tool provides a parameter search for it. Please note that HiCCUPS and any dependency of it are not provided by HiCExplorer and must be installed by the user on its own.


The output file (´´-o loops.bedgraph``) contains the x and y position of each loop and its corresponding p-value of the statistical test. 

.. code::

    1	120000000	122500000	1	145000000	147500000	0.001


The results can visualized using :doc:`hicPlotMatrix`:


.. code:: bash

    $ hicPlotMatrix -m matrix.cool -o plot.png --log1p --region 1:18000000-22000000 --loops loops.bedgraph


.. image:: ../../images/hicDetectLoops.png


.. argparse::
   :ref: hicexplorer.hicDetectLoops.parse_arguments
   :prog: hicDetectLoops

