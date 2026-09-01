"""Schema driven comparator for the textual outputs.

Built in schemas:

    plain      byte comparison, E0 only. Used for hicInfo, whose whole output
               including the thousands separators and the float reprs is the
               contract.
    bed        3 to 12 columns, chrom/start/end then free columns
    bedgraph   4 columns, the fourth a float
    bedpe      at least 6 columns
    tsv-header a leading header line, then columns typed by the schema

Line count and line order are always compared exactly. Integer and string
columns are compared exactly; float columns at the requested class.
"""
from __future__ import annotations

from .base import Result, fail, values_agree

SCHEMAS = {
    "plain": None,
    "bed": {"comment": "#", "delimiter": "\t", "header": False,
            "columns": ["str", "int", "int"], "rest": "str"},
    "bedgraph": {"comment": "#", "delimiter": "\t", "header": False,
                 "columns": ["str", "int", "int", "float"], "rest": "float"},
    "bedpe": {"comment": "#", "delimiter": "\t", "header": False,
              "columns": ["str", "int", "int", "str", "int", "int"],
              "rest": "float"},
    "tsv-header": {"comment": "#", "delimiter": "\t", "header": True,
                   "columns": [], "rest": "float"},
}


def _read_bytes(path):
    """The content of the file, decompressed when it is gzipped.

    The homer writer is gzip.open(..., 'wt'), whose header carries the
    modification time, so two runs never produce identical compressed bytes.
    E0 for that format therefore means byte identical *content*, which is also
    what hicexplorer/test/general/test_hicConvertFormat.py asserts.
    """
    with open(path, "rb") as handle:
        raw = handle.read()
    if raw[:2] == b"\x1f\x8b":
        import gzip

        return gzip.decompress(raw)
    return raw


def _read_lines(path, schema):
    raw = _read_bytes(path)
    text = raw.decode("utf-8", errors="replace")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if schema is None:
        return lines
    comment = schema["comment"]
    return [line for line in lines if not line.startswith(comment)]


def _column_type(schema, index):
    if index < len(schema["columns"]):
        return schema["columns"][index]
    return schema["rest"]


def compare(path_a, path_b, cls, opts=None):
    opts = opts or {}
    schema_name = opts.get("schema", "plain")
    if schema_name not in SCHEMAS:
        return fail(f"unknown text schema {schema_name!r}")
    schema = SCHEMAS[schema_name]

    if schema is None or cls == "E0":
        a = _read_bytes(path_a)
        b = _read_bytes(path_b)
        if a == b:
            return Result(True, "E0", {"bytes": len(a)})
        diffs = _byte_diff(a, b)
        return Result(False, None, {"bytes_py": len(a), "bytes_cpp": len(b)}, diffs)

    lines_a = _read_lines(path_a, schema)
    lines_b = _read_lines(path_b, schema)
    if len(lines_a) != len(lines_b):
        return fail(f"line count differs: {len(lines_a)} vs {len(lines_b)}",
                    lines_py=len(lines_a), lines_cpp=len(lines_b))

    diffs = []
    delimiter = schema["delimiter"]
    for index, (line_a, line_b) in enumerate(zip(lines_a, lines_b), start=1):
        if line_a == line_b:
            continue
        fields_a = line_a.split(delimiter)
        fields_b = line_b.split(delimiter)
        if len(fields_a) != len(fields_b):
            diffs.append(f"line {index}: column count {len(fields_a)} vs "
                         f"{len(fields_b)}")
            continue
        for column, (value_a, value_b) in enumerate(zip(fields_a, fields_b)):
            if value_a == value_b:
                continue
            kind = _column_type(schema, column)
            if kind == "float":
                try:
                    if values_agree(float(value_a), float(value_b), cls):
                        continue
                except ValueError:
                    pass
            diffs.append(f"line {index} column {column + 1} ({kind}): "
                         f"{value_a!r} vs {value_b!r}")
        if len(diffs) >= 20:
            break

    if diffs:
        return Result(False, None, {"lines": len(lines_a)}, diffs)
    return Result(True, cls, {"lines": len(lines_a)})


def _byte_diff(a, b):
    """First differing byte plus the surrounding lines, for the report."""
    limit = min(len(a), len(b))
    position = limit
    for i in range(limit):
        if a[i] != b[i]:
            position = i
            break
    line_number = a[:position].count(b"\n") + 1
    lines_a = a.decode("utf-8", "replace").split("\n")
    lines_b = b.decode("utf-8", "replace").split("\n")
    diffs = [f"first difference at byte {position}, line {line_number}"]
    for index in range(max(len(lines_a), len(lines_b))):
        left = lines_a[index] if index < len(lines_a) else "<missing>"
        right = lines_b[index] if index < len(lines_b) else "<missing>"
        if left != right:
            diffs.append(f"line {index + 1} py : {left!r}")
            diffs.append(f"line {index + 1} cpp: {right!r}")
        if len(diffs) >= 21:
            break
    return diffs
