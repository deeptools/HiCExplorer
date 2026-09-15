#!/usr/bin/env python3
"""Python reference for chicDifferentialTest's C++ only option (PLAN.md 9.7
work item 3).

    py_chicDifferentialTest_calibrated.py <chicDifferentialTest arguments>
        [--correctForMultipleTesting {none,fdr,bonferroni}]

Without the option, or with `none`, this is hicexplorer.chicDifferentialTest
itself. With `fdr` or `bonferroni` the unmodified tool runs twice in this
process, its worker processes replaced by in-process calls in the same order:

1. With the output writer disabled, the p-values of every test
   (fisher_exact_test's or chisquare_test's test_result, NaN for a table scipy
   refuses) are collected over all reference points in computation order.
2. They are adjusted together, NaN excluded from the count and left NaN
   (benjamini_hochberg_adjusted and bonferroni_adjusted of
   py_hicDifferentialTAD_calibrated.py, the operations hicx::stats performs).
3. The tool runs again. The test functions return their own p-values, but a
   location is rejected when its adjusted p-value is at most --alpha, for both
   tests; a refused table stays accepted with 1.0. Every result group gains
   pvalue_adjusted_list, aligned with pvalue_list, and the file the attribute
   correctForMultipleTesting.
"""
from __future__ import annotations

import math
import sys

CORRECTIONS = ("none", "fdr", "bonferroni")


def fail(message):
    print("chicDifferentialTest: error: " + message, file=sys.stderr)
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


class InProcess:
    """multiprocessing.Process running its target when started."""

    def __init__(self, target, kwargs):
        self.target = target
        self.kwargs = kwargs

    def start(self):
        self.target(**self.kwargs)

    def join(self):
        pass

    def terminate(self):
        pass


class ListQueue:
    def __init__(self):
        self.items = []

    def put(self, item):
        self.items.append(item)

    def empty(self):
        return not self.items

    def get(self):
        return self.items.pop(0)


class NoSleep:
    @staticmethod
    def sleep(_seconds):
        pass


def main(argv):
    correction, rest = split_options(argv)
    from hicexplorer import chicDifferentialTest as tool
    if correction == "none":
        tool.main(rest)
        return

    import h5py
    from py_hicDifferentialTAD_calibrated import benjamini_hochberg_adjusted, bonferroni_adjusted

    tool.Process = InProcess
    tool.Queue = ListQueue
    tool.time = NoSleep
    original_tests = {"fisher_exact_test": tool.fisher_exact_test, "chisquare_test": tool.chisquare_test}
    original_write = tool.writeResultHDF

    # 1. the p-values of every test
    collected = []

    def recording(name):
        def test(pDataFile1, pDataFile2, pAlpha):
            result = original_tests[name](pDataFile1, pDataFile2, pAlpha)
            collected.extend(result[0])
            return result
        return test

    tool.fisher_exact_test = recording("fisher_exact_test")
    tool.chisquare_test = recording("chisquare_test")
    tool.writeResultHDF = lambda *args, **kwargs: None
    tool.main(list(rest))

    # 2. the adjustment
    adjusted_values = benjamini_hochberg_adjusted(collected) if correction == "fdr" \
        else bonferroni_adjusted(collected)
    position = {"next": 0}
    adjusted_by_row = {}

    # 3. the decisions on the adjusted values
    def deciding(name):
        def test(pDataFile1, pDataFile2, pAlpha):
            test_result, accepted, rejected = original_tests[name](pDataFile1, pDataFile2, pAlpha)
            start = position["next"]
            position["next"] += len(test_result)
            own = {i: p for i, p in accepted + rejected}
            new_accepted = []
            new_rejected = []
            for i, pvalue in enumerate(test_result):
                adjusted = adjusted_values[start + i]
                adjusted_by_row[id(pDataFile1[i])] = adjusted
                refused = math.isnan(pvalue) and own[i] == 1.0 and [i, 1.0] in accepted
                if refused:
                    new_accepted.append([i, 1.0])
                elif adjusted <= pAlpha:
                    new_rejected.append([i, own[i]])
                else:
                    new_accepted.append([i, own[i]])
            return test_result, new_accepted, new_rejected
        return test

    def write(pOutFileName, pAcceptedData, pRejectedData, pAllResultData, pInputData, pAlpha, pTest):
        original_write(pOutFileName, pAcceptedData, pRejectedData, pAllResultData, pInputData, pAlpha, pTest)
        categories = {'accepted': pAcceptedData, 'rejected': pRejectedData, 'all': pAllResultData}
        with h5py.File(pOutFileName, "a") as handle:
            handle.attrs["correctForMultipleTesting"] = correction
            for i, inputData in enumerate(pInputData):
                group = handle[inputData[0][1]][inputData[1][1]][inputData[0][2]][inputData[0][3]]
                for category in ['accepted', 'rejected', 'all']:
                    rows = categories[category][i]
                    if len(rows) == 0:
                        continue
                    values = [adjusted_by_row[id(row[3])] for row in rows]
                    group[category].create_dataset("pvalue_adjusted_list", data=values,
                                                   compression="gzip", compression_opts=9)

    tool.fisher_exact_test = deciding("fisher_exact_test")
    tool.chisquare_test = deciding("chisquare_test")
    tool.writeResultHDF = write
    tool.main(list(rest))


if __name__ == "__main__":
    main(sys.argv[1:])
