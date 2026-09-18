"""The tools of HiCExplorer and whether the C++ port provides them.

PLANNED_TOOLS lists every tool of cpp/PLAN.md tiers 1 to 8 with its tier. A
tool is available when its executable in the configured tool directory prints
a specification with ``--help-json``; otherwise the catalog says why not.
"""

import os

from .workflow.spec import SpecError, SpecLoader

PLANNED_TOOLS = {
    # tier 1: file layer exercisers
    "hicAdjustMatrix": 1, "hicCompareMatrices": 1, "hicConvertFormat": 1, "hicInfo": 1,
    "hicMergeMatrixBins": 1, "hicSumMatrices": 1,
    # tier 2: interval and text tools
    "hicAverageRegions": 2, "hicCreateThresholdFile": 2, "hicFindRestSite": 2, "hicMergeLoops": 2,
    "hicMergeTADbins": 2, "hicNormalize": 2, "hicValidateLocations": 2,
    # tier 3: float matrix math
    "hicCompartmentalization": 3, "hicCorrectMatrix": 3, "hicDetectLoops": 3, "hicFindTADs": 3,
    "hicInterIntraTAD": 3, "hicPCA": 3, "hicPlotSVL": 3, "hicTransform": 3,
    # tier 4: alignment and matrix construction
    "hicBuildMatrix": 4, "hicBuildMatrixMicroC": 4, "hicQuickQC": 4,
    # tier 5: TAD, loop and differential calling
    "hicAggregateContacts": 5, "hicDifferentialTAD": 5, "hicMergeDomains": 5,
    # tier 6: capture Hi-C suite
    "chicAggregateStatistic": 6, "chicDifferentialTest": 6, "chicExportData": 6,
    "chicPlotViewpoint": 6, "chicQualityControl": 6, "chicSignificantInteractions": 6,
    "chicViewpoint": 6, "chicViewpointBackgroundModel": 6,
    # tier 7: plotting
    "hicCorrelate": 7, "hicPlotAverageRegions": 7, "hicPlotDistVsCounts": 7, "hicPlotMatrix": 7,
    "hicPlotTADs": 7, "hicPlotViewpoint": 7, "hicPrepareQCreport": 7, "hicQC": 7,
    # tier 8: machine learning and hyperparameter search
    "hicHyperoptDetectLoops": 8, "hicHyperoptDetectLoopsHiCCUPS": 8, "hicTADClassifier": 8,
    "hicTrainTADClassifier": 8,
    # tier 9: features without a Python counterpart
    "hicDifferentialAnalysis": 9,
    "hicDetectStripes": 9,
}

TIER_NAMES = {
    1: "file layer", 2: "interval and text tools", 3: "float matrix math",
    4: "alignment and matrix construction", 5: "TAD, loop and differential calling",
    6: "capture Hi-C", 7: "plotting", 8: "machine learning", 9: "new features",
}


class ToolEntry:
    def __init__(self, name, tier, spec=None, reason=""):
        self.name = name
        self.tier = tier
        self.spec = spec
        self.reason = reason

    @property
    def available(self):
        return self.spec is not None

    def __repr__(self):
        return "ToolEntry({!r}, available={})".format(self.name, self.available)


def tool_entries(tools_dir):
    """One ToolEntry per planned tool, sorted by tier and name."""
    loader = SpecLoader(tools_dir or None)
    entries = []
    for name, tier in sorted(PLANNED_TOOLS.items(), key=lambda item: (item[1], item[0])):
        tier_text = "PLAN tier {} ({})".format(tier, TIER_NAMES[tier])
        if not tools_dir:
            entries.append(ToolEntry(name, tier, reason="no C++ tool directory is set (Settings)"))
            continue
        if not os.path.isfile(os.path.join(tools_dir, name)):
            entries.append(ToolEntry(name, tier, reason="not ported to C++ yet, " + tier_text))
            continue
        try:
            entries.append(ToolEntry(name, tier, spec=loader.load(name)))
        except SpecError as exc:
            entries.append(ToolEntry(name, tier, reason="the executable gives no specification: {}".format(exc)))
    return entries, loader
