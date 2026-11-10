"""Archives of figures (tar or tar.gz), for chicPlotViewpoint's plots.tar.gz.

The archive bytes are not comparable: tarfile stamps every member with
time.time() and gzip writes its own timestamp. What the user gets is the
list of member names and the images, so those are compared: the same member
names in the same order, and every member through comparators/image.py at
the declared class (E6 or E0). The metrics record the largest RMS and the
per member results.
"""
from __future__ import annotations

import tarfile
import tempfile
from pathlib import Path

from . import image
from .base import Result, fail


def _members(path):
    with tarfile.open(path, "r:*") as archive:
        return [member.name for member in archive.getmembers() if member.isfile()]


def check_structure(path_a, path_b, opts=None):
    try:
        _members(path_a)
        _members(path_b)
    except (tarfile.TarError, OSError) as error:
        return fail(f"not a readable tar archive: {error}")
    return Result(True, "E7", {})


def compare(path_a, path_b, cls, opts=None):
    if cls not in ("E0", "E6"):
        return fail(f"class {cls} is not defined for archives of figures; use E6 or E0")
    try:
        names_a = _members(path_a)
        names_b = _members(path_b)
    except (tarfile.TarError, OSError) as error:
        return fail(f"not a readable tar archive: {error}")
    if names_a != names_b:
        return fail("the archives hold different members",
                    members_py=names_a, members_cpp=names_b)
    worst = 0.0
    per_member = {}
    class_met = "E0"
    with tempfile.TemporaryDirectory(prefix="tar-images-") as scratch:
        for side, path in (("py", path_a), ("cpp", path_b)):
            with tarfile.open(path, "r:*") as archive:
                for member in archive.getmembers():
                    if member.isfile():
                        target = Path(scratch) / side / member.name
                        target.parent.mkdir(parents=True, exist_ok=True)
                        target.write_bytes(archive.extractfile(member).read())
        for name in names_a:
            result = image.compare(str(Path(scratch) / "py" / name),
                                   str(Path(scratch) / "cpp" / name), cls, opts)
            per_member[name] = result.to_json()
            if not result.passed:
                return Result(False, None, {"members": per_member},
                              [f"{name}: {message}" for message in result.to_json().get("diffs", [])])
            rms = result.metrics.get("rms", 0.0)
            worst = max(worst, rms or 0.0)
            if result.class_met != "E0":
                class_met = "E6"
    return Result(True, class_met, {"members": len(names_a), "max_rms": worst,
                                    "per_member": per_member})
