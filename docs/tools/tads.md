# TADs

!!! warning "Two tools dropped from this rewrite"
    `hicTADClassifier` and `hicTrainTADClassifier` existed in the Python HiCExplorer's TADs category but were dropped from this C++ rewrite on 2026-09-21 (see `cpp/PLAN.md` tier 8): little real-world usage did not justify an ONNX-incompatible Python shell around the C++ tools. There is no C++ build of either tool and no page for them here; use the Python HiCExplorer if you need machine-learning-based TAD boundary prediction.

| Tool | Description |
|---|---|
| [hicFindTADs](hicFindTADs.md) | Identifies Topologically Associating Domains (TADs). |
| [hicMergeDomains](hicMergeDomains.md) | Merges TAD domains called at different resolutions and their hierarchical relation. |
| [hicDifferentialTAD](hicDifferentialTAD.md) | Identifies differential TADs between two Hi-C matrices. |
| [hicMergeTADbins](hicMergeTADbins.md) | Uses a BED file of domains or TAD boundaries to merge the bin counts of a Hi-C matrix. |
| [hicInterIntraTAD](hicInterIntraTAD.md) | Computes and plots the inter-TAD versus intra-TAD contact ratio. |
