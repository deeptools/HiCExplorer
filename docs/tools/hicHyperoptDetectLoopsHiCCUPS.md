# hicHyperoptDetectLoopsHiCCUPS

!!! warning "Not yet ported to C++"
    `hicHyperoptDetectLoopsHiCCUPS` exists in the Python HiCExplorer but has not been ported to this C++
    rewrite yet: `cpp/STATUS.md` lists it as tier 8, "not started". There is no
    `cpp/build/tools/hicHyperoptDetectLoopsHiCCUPS` binary, and consequently no real `--help` output to
    generate a CLI reference from. This page is a placeholder so the tool is not silently missing from
    the documentation; it will be filled in with a real reference once the tool is ported.

Searches for the best parameter setting for Juicer's HiCCUPS loop caller against a set of known protein
peak locations, the HiCCUPS counterpart of [hicHyperoptDetectLoops](hicHyperoptDetectLoops.md). HiCCUPS
itself, and any of its dependencies, are not provided by HiCExplorer and must be installed separately.

Until this tool is ported, use the Python HiCExplorer's `hicHyperoptDetectLoopsHiCCUPS`.
