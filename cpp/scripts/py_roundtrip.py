#!/usr/bin/env python3
"""Read a contact matrix with hicmatrix and write it out again.

The Python half of the tier 0 writer validation. It is the reference for
cpp/build/tests/hicx_roundtrip and takes the same path through hicmatrix that
hicConvertFormat takes (hicexplorer/hicConvertFormat.py:230-273): load through
the file handler for the input format, hand the five loader values to the
handler for the output format, save. Nothing above the file layer, in
particular nothing of hiCMatrix, is involved, so the load time field swap of
cpp/PLAN.md 2.7 quirk 1 does not enter.

    py_roundtrip.py <input> <output> [--enforce-integer] [--no-symmetric]
                                     [--no-apply-correction]

Run it with the reference venv interpreter named in cpp/AGENTS_CONTRACT.md and
with PYTHONPATH pointing at the repository root.
"""
from __future__ import annotations

import sys

from hicmatrix.lib import MatrixFileHandler


def main(argv):
    positional = [argument for argument in argv if not argument.startswith("-")]
    flags = {argument for argument in argv if argument.startswith("-")}
    unknown = flags - {"--enforce-integer", "--no-symmetric", "--no-apply-correction"}
    if len(positional) != 2 or unknown:
        print(__doc__, file=sys.stderr)
        return 2
    source, destination = positional
    enforce_integer = "--enforce-integer" in flags
    symmetric = "--no-symmetric" not in flags
    apply_correction = "--no-apply-correction" not in flags

    input_format = "h5" if source.endswith(".h5") else "cool"
    output_format = "h5" if destination.endswith(".h5") else "cool"

    handler_in = MatrixFileHandler(pFileType=input_format, pMatrixFile=source,
                                   pEnforceInteger=enforce_integer)
    loaded = handler_in.load()
    if len(loaded) == 2:
        print(f"failed to load {source}: {loaded[1]}", file=sys.stderr)
        return 1
    matrix, cut_intervals, nan_bins, distance_counts, correction_factors = loaded

    hic2cool_version = None
    cool_metadata = None
    if input_format == "cool":
        hic2cool_version = handler_in.matrixFile.hic2cool_version
        cool_metadata = handler_in.matrixFile.hic_metadata

    handler_out = MatrixFileHandler(pFileType=output_format,
                                    pEnforceInteger=enforce_integer,
                                    pFileWasH5=(input_format == "h5"),
                                    pHic2CoolVersion=hic2cool_version,
                                    pHiCInfo=cool_metadata)
    handler_out.set_matrix_variables(matrix, cut_intervals, nan_bins,
                                     correction_factors, distance_counts)
    handler_out.save(destination, pSymmetric=symmetric,
                     pApplyCorrection=apply_correction)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
