# hicHyperoptDetectLoops

!!! warning "Not yet ported to C++"
    `hicHyperoptDetectLoops` exists in the Python HiCExplorer but has not been ported to this C++
    rewrite yet: `cpp/STATUS.md` lists it as tier 8, "not started". There is no `cpp/build/tools/hicHyperoptDetectLoops`
    binary, and consequently no real `--help` output to generate a CLI reference from. This page is a
    placeholder so the tool is not silently missing from the documentation; it will be filled in with a
    real reference once the tool is ported.

Searches for the best [hicDetectLoops](hicDetectLoops.md) parameter setting for a given dataset, using
Bayesian hyperparameter optimization (the Python `hyperopt` library) against a set of known protein peak
locations. Because `hicDetectLoops` has many parameters and finding a good setting by hand is difficult,
the Python HiCExplorer added this tool (and its HiCCUPS counterpart,
[hicHyperoptDetectLoopsHiCCUPS](hicHyperoptDetectLoopsHiCCUPS.md)) in version 3.5.

Until this tool is ported, use the Python HiCExplorer's `hicHyperoptDetectLoops` for parameter search,
and this rewrite's [hicDetectLoops](hicDetectLoops.md) for the loop calling itself once a parameter
setting has been chosen.
