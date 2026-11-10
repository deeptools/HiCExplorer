"""Archives of data files (tar or tar.gz), for chicExportData.

The archive bytes are not comparable: tarfile stamps each member with
time.time(), members added from disk carry the owner and modification time of
a temporary file, and gzip writes its own header. What the user gets is the
member names and the files, so those are compared:

  * the member names, in order; with options {"ordered": false} as sets, for
    the bigWig archive chicExportData fills from a temporary directory in
    os.walk order;
  * every member: .bigwig and .bw through comparators/bigwig.py at the
    declared class, anything else byte for byte (the text exports are E0).
"""
from __future__ import annotations

import tarfile
import tempfile
from pathlib import Path

from . import bigwig
from .base import Result, fail


def _members(path):
    with tarfile.open(path, "r:*") as archive:
        return [member.name for member in archive.getmembers() if member.isfile()]


def compare(path_a, path_b, cls, opts=None):
    opts = opts or {}
    try:
        names_a = _members(path_a)
        names_b = _members(path_b)
    except (tarfile.TarError, OSError) as error:
        return fail(f"not a readable tar archive: {error}")
    ordered = opts.get("ordered", True)
    if (names_a != names_b) if ordered else (sorted(names_a) != sorted(names_b)):
        return fail("the archives hold different members" + ("" if ordered else " (unordered)"),
                    members_py=names_a, members_cpp=names_b)
    diffs = []
    metrics = {"members": len(names_a), "bigwig_members": 0, "text_members": 0}
    with tempfile.TemporaryDirectory(prefix="tar-members-") as scratch:
        for side, path in (("py", path_a), ("cpp", path_b)):
            with tarfile.open(path, "r:*") as archive:
                for member in archive.getmembers():
                    if member.isfile():
                        target = Path(scratch) / side / member.name
                        target.parent.mkdir(parents=True, exist_ok=True)
                        target.write_bytes(archive.extractfile(member).read())
        for name in sorted(names_a):
            left = Path(scratch) / "py" / name
            right = Path(scratch) / "cpp" / name
            if name.endswith((".bigwig", ".bw")):
                metrics["bigwig_members"] += 1
                result = bigwig.compare(left, right, cls, opts)
                if not result.passed:
                    diffs.extend(f"{name}: {diff}" for diff in result.diffs[:5])
            else:
                metrics["text_members"] += 1
                if left.read_bytes() != right.read_bytes():
                    lines_a = left.read_bytes().split(b"\n")
                    lines_b = right.read_bytes().split(b"\n")
                    first = next((i for i, (x, y) in enumerate(zip(lines_a, lines_b)) if x != y),
                                 min(len(lines_a), len(lines_b)))
                    a = lines_a[first] if first < len(lines_a) else b"<end>"
                    b = lines_b[first] if first < len(lines_b) else b"<end>"
                    diffs.append(f"{name}: line {first + 1} differs: {a[:120]!r} vs {b[:120]!r}")
    if diffs:
        return Result(False, None, metrics, diffs[:20])
    return Result(True, cls, metrics)
