"""Files an external program renders from a source the tool writes.

This exists for hicMergeDomains. Its relation trees are drawn by graphviz's
`dot` binary: the Python writes a DOT source through graphviz.Digraph and runs
`dot -Kdot -T<format> -O <source>`, and the port writes a byte identical source
(pinned by cpp/tests/test_merge_domains.cpp) and runs the same command. The
rendered file is therefore equivalent by construction, and what can still go
wrong is checked here: the file exists under the expected name and carries the
signature of its format, and, where `dot` is reproducible, it is byte
identical.

Measured with graphviz 12.0.0 on the hicMergeDomains trees, 2026-09-13: a png
or svg rendering of the same source is byte identical run to run. A pdf is
not: cairo writes /CreationDate into a compressed object stream and a second
compressed stream changes as well, so normalising the named date field does not
make two renderings equal. A pdf output is declared E7 and only its structure
is checked (comparators.compare calls check_structure for E7).

E0 is a byte comparison. No other class is defined for rendered files.
"""
from __future__ import annotations

from pathlib import Path

from .base import Result, fail

SIGNATURES = {
    "png": (b"\x89PNG\r\n\x1a\n",),
    "pdf": (b"%PDF-",),
    "svg": (b"<?xml", b"<svg"),
    "gif": (b"GIF87a", b"GIF89a"),
    "jpg": (b"\xff\xd8\xff",),
    "jpeg": (b"\xff\xd8\xff",),
    "ps": (b"%!PS",),
    "eps": (b"%!PS",),
}


def _signature_problem(path):
    suffix = Path(path).suffix.lower().lstrip(".")
    expected = SIGNATURES.get(suffix)
    if expected is None:
        return f"{path}: no known signature for the suffix .{suffix}"
    with open(path, "rb") as handle:
        head = handle.read(16)
    if not any(head.startswith(signature) for signature in expected):
        return f"{path}: does not start with the {suffix} signature ({head[:8]!r})"
    return None


def check_structure(path_a, path_b, opts=None):
    problems = [problem for problem in (_signature_problem(path_a),
                                        _signature_problem(path_b)) if problem]
    if problems:
        return Result(False, None, {}, ["rendered file has the wrong format"] + problems)
    return Result(True, "E7", {"bytes_py": Path(path_a).stat().st_size,
                               "bytes_cpp": Path(path_b).stat().st_size})


def compare(path_a, path_b, cls, opts=None):
    structure = check_structure(path_a, path_b, opts)
    if not structure.passed:
        return structure
    if cls != "E0":
        return fail(f"class {cls} is not defined for rendered files; declare E0 "
                    f"where the renderer is reproducible and E7 where it is not")
    with open(path_a, "rb") as handle:
        a = handle.read()
    with open(path_b, "rb") as handle:
        b = handle.read()
    if a == b:
        return Result(True, "E0", {"bytes": len(a)})
    first = next((i for i in range(min(len(a), len(b))) if a[i] != b[i]),
                 min(len(a), len(b)))
    return fail(f"rendered files differ, first at byte {first}",
                bytes_py=len(a), bytes_cpp=len(b))
