#!/usr/bin/env python3
"""Python reference for chicSignificantInteractions' C++ only option (PLAN.md
9.7 work item 3).

    py_chicSignificantInteractions_calibrated.py <chicSignificantInteractions arguments>
        [--correctForMultipleTesting {none,fdr,bonferroni}]

Without the option, or with `none`, this is hicexplorer.chicSignificantInteractions
itself. With `fdr` or `bonferroni` it builds what the C++ port must produce out
of the unmodified tool's own functions:

1. The tested p-values. The tool's pipeline runs once in this process with the
   output writers disabled. For every sample viewpoint (the triplet [matrix,
   chromosome, gene] the tool reads) the p-values it would compare against
   --pValue are recorded, the first time the triplet is met and in the order
   the tool meets its candidates:
     - with --xFoldBackground or --loosePValue, the tool's compute_new_p_values
       sets the p-value of every merged candidate whose relative position is
       in the background model; it is called with pValue inf and peak
       threshold -inf, so that it computes them all;
     - without a preselection, filter_by_pvalue compares the stored p-value of
       every position.
2. Those p-values are adjusted across all of them, NaN excluded from the
   count and left NaN (benjamini_hochberg_adjusted and bonferroni_adjusted of
   py_hicDifferentialTAD_calibrated.py, the operations hicx::stats performs).
3. The tool runs again, unmodified except that compute_new_p_values and
   filter_by_pvalue accept a candidate when its adjusted p-value is at most
   the threshold, with the peak condition each of them applies. The p-value
   each one stores stays the tool's own.
4. Every viewpoint group of the significant file gains the dataset
   pvalue_adjusted (float64, gzip 9, aligned with pvalue), and both files the
   root attribute correctForMultipleTesting.
"""
from __future__ import annotations

import copy
import math
import os
import sys
import tempfile

CORRECTIONS = ("none", "fdr", "bonferroni")


def fail(message):
    print("chicSignificantInteractions: error: " + message, file=sys.stderr)
    sys.exit(2)


def split_options(argv):
    correction = "none"
    rest = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token == "--correctForMultipleTesting" or token.startswith("--correctForMultipleTesting="):
            if "=" in token:
                value = token.split("=", 1)[1]
            elif i + 1 < len(argv):
                i += 1
                value = argv[i]
            else:
                fail("argument --correctForMultipleTesting: expected one argument")
            if value not in CORRECTIONS:
                fail("argument --correctForMultipleTesting: invalid choice: '{}' (choose from "
                     "'none', 'fdr', 'bonferroni')".format(value))
            correction = value
        else:
            rest.append(token)
        i += 1
    return correction, rest


def main(argv):
    correction, rest = split_options(argv)
    from hicexplorer import chicSignificantInteractions as tool
    if correction == "none":
        tool.main(rest)
        return

    import h5py
    from hicexplorer.lib import Viewpoint
    from py_hicDifferentialTAD_calibrated import benjamini_hochberg_adjusted, bonferroni_adjusted

    state = {"current": None}
    original_read = Viewpoint.readInteractionFile

    def read(self, pFilePath, pInternalIdentifierTriplet):
        state["current"] = tuple(str(part) for part in pInternalIdentifierTriplet)
        return original_read(self, pFilePath, pInternalIdentifierTriplet)

    Viewpoint.readInteractionFile = read

    original_compute = tool.compute_new_p_values
    original_filter = tool.filter_by_pvalue
    original_call = tool.call_multi_core
    original_write_significant = tool.writeSignificantHDF
    original_write_target = tool.writeTargetHDF

    # 1. the tested p-values
    candidates = []
    seen = set()

    def record_compute(pData, pBackgroundModel, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
        data = copy.deepcopy(pData)
        original_compute(data, pBackgroundModel, math.inf, pMergedLinesDict, -math.inf)
        if state["current"] not in seen:
            seen.add(state["current"])
            for key in data:
                if key in pBackgroundModel:
                    candidates.append(((state["current"], key), data[key][-3]))
        return {}, []

    def record_filter(pData, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
        if state["current"] not in seen:
            seen.add(state["current"])
            for key in pData:
                candidates.append(((state["current"], key), pData[key][-3]))
        return {}, []

    def sequential(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution):
        tool.compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackground,
                                      pFilePath, pResolution, pQueue=None)
        return [], [], [], [], [], []

    tool.compute_new_p_values = record_compute
    tool.filter_by_pvalue = record_filter
    tool.call_multi_core = sequential
    tool.writeSignificantHDF = lambda *args, **kwargs: None
    tool.writeTargetHDF = lambda *args, **kwargs: None
    here = os.getcwd()
    with tempfile.TemporaryDirectory() as scratch:
        # compute_interaction_file appends to errorLog.txt in the working
        # directory; the first pass must not add lines to the real one.
        os.chdir(scratch)
        try:
            tool.main(list(rest))
        finally:
            os.chdir(here)

    # 2. the adjustment
    values = [p for _, p in candidates]
    adjusted_values = benjamini_hochberg_adjusted(values) if correction == "fdr" \
        else bonferroni_adjusted(values)
    adjusted = {key: value for (key, _), value in zip(candidates, adjusted_values)}

    def adjusted_of(key):
        return adjusted.get((state["current"], key), math.nan)

    # 3. the decisions on the adjusted values
    def compute(pData, pBackgroundModel, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
        original_compute(pData, pBackgroundModel, math.inf, pMergedLinesDict, -math.inf)
        accepted = {}
        accepted_lines = []
        for key in pData:
            if key in pBackgroundModel:
                threshold = pPValue if isinstance(pPValue, float) else pPValue[key]
                if adjusted_of(key) <= threshold:
                    if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                        accepted[key] = pData[key]
                        target_content = pMergedLinesDict[key][0][:3]
                        target_content[2] = pMergedLinesDict[key][-1][2]
                        accepted_lines.append(target_content)
        return accepted, accepted_lines

    def filter_(pData, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
        accepted = {}
        accepted_lines = []
        for key in pData:
            threshold = pPValue if isinstance(pPValue, float) else pPValue[key]
            if adjusted_of(key) <= threshold:
                if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                    accepted[key] = pMergedLinesDict[key]
                    accepted_lines.append(pMergedLinesDict[key][:3])
        return accepted, accepted_lines

    # 4. the adjusted values in the output
    def write_significant(pOutFileName, pSignificantDataList, pSignificantKeyList, pViewpointObj,
                          pReferencePointsList, pArgs):
        order = iter([tuple(str(part) for part in key)
                      for key, data in zip(pSignificantKeyList, pSignificantDataList) if len(data) > 0])
        original_writer = pViewpointObj.writeInteractionFileHDF5

        def writer(pInteractionFileGroupH5Object, pFileName, pData, pReferencePoint):
            triplet = next(order)
            name = original_writer(pInteractionFileGroupH5Object, pFileName, pData, pReferencePoint)
            values = [adjusted.get((triplet, key), math.nan) for key in pData[5]]
            pInteractionFileGroupH5Object[name].create_dataset(
                "pvalue_adjusted", data=values, compression="gzip", compression_opts=9)
            return name

        pViewpointObj.writeInteractionFileHDF5 = writer
        try:
            original_write_significant(pOutFileName, pSignificantDataList, pSignificantKeyList,
                                       pViewpointObj, pReferencePointsList, pArgs)
        finally:
            del pViewpointObj.writeInteractionFileHDF5
        with h5py.File(pOutFileName, "a") as handle:
            handle.attrs["correctForMultipleTesting"] = correction

    def write_target(pOutFileName, *args):
        original_write_target(pOutFileName, *args)
        with h5py.File(pOutFileName, "a") as handle:
            handle.attrs["correctForMultipleTesting"] = correction

    tool.compute_new_p_values = compute
    tool.filter_by_pvalue = filter_
    tool.call_multi_core = original_call
    tool.writeSignificantHDF = write_significant
    tool.writeTargetHDF = write_target
    tool.main(list(rest))


if __name__ == "__main__":
    main(sys.argv[1:])
