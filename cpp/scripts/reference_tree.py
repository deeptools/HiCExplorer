"""The Python reference implementation: upstream HiCExplorer at the fork point.

HiCExplorer 4 ships no Python implementation of its tools. The equivalence
harness compares the C++ tools with upstream HiCExplorer at the commit the
C++ port forked from, REFERENCE_COMMIT, which is the one place that commit is
recorded. The harness extracts the `hicexplorer` package of that commit into a
directory outside the repository, the reference root, and puts that directory
first on PYTHONPATH for the reference processes only:

    <root>/hicexplorer/...      git archive REFERENCE_COMMIT hicexplorer,
                                without hicexplorer/test
    <root>/REFERENCE.json       commit, overlay and file digest of the tree

The root is --reference-root of equiv.py, else $HICX_REFERENCE_ROOT, else
<cache>/reference/<commit>-<version> under the harness cache directory
(reference_cache.default_cache_dir). A missing root is created on demand from
the git repository $HICX_REFERENCE_GIT, else this repository, else a shallow
fetch of UPSTREAM_URL. A root whose REFERENCE.json names another commit or
overlay, or whose files no longer match the recorded digest, is an error.

The interpreter is the conda env of cpp/AGENTS_CONTRACT.md
(~/miniconda3/envs/__hicexplorer@3.7.6), which provides the dependencies and
also has a released hicexplorer 3.7.6 in its site-packages. PYTHONPATH entries
come before site-packages, so the reference root wins; check_imported() makes
every reference entry point fail if `hicexplorer` came from anywhere else
(site-packages or the repository's own thin API package).

Overlay. The tree is the commit's files with one change: hicexplorer/_version.py
carries the version of the package under test (hicexplorer/_version.py of this
repository) instead of 3.7.6. The tools write their version into output
headers, HDF5 attributes and --version, and the C++ tools write the same
label, so both sides must carry it. Before the reference moved out of the
package it ran with exactly this version string. Help texts are taken
unchanged from the commit: the only help-output case (hicInfo.help) is for a
tool whose help text the fork did not change, and tool_specs.py does not
compare help texts.
"""

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REFERENCE_COMMIT = "c2ad8630c48e04f0ce4229d75d4b1766df79317c"
UPSTREAM_URL = "https://github.com/deeptools/HiCExplorer.git"
STAMP = "REFERENCE.json"
TREE_FORMAT = 1

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent


def product_version():
    text = (REPO_ROOT / "hicexplorer" / "_version.py").read_text()
    match = re.search(r"""__version__\s*=\s*['"]([^'"]+)['"]""", text)
    if not match:
        raise RuntimeError("no __version__ in hicexplorer/_version.py")
    return match.group(1)


def overlay():
    """{relative path: content} written over the extracted commit."""
    return {"hicexplorer/_version.py":
            "# Reference tree of cpp/scripts/reference_tree.py: upstream "
            f"{REFERENCE_COMMIT[:8]} labelled\n# with the version under test.\n"
            f"__version__ = '{product_version()}'\n"}


def default_root():
    configured = os.environ.get("HICX_REFERENCE_ROOT")
    if configured:
        return Path(configured).expanduser()
    sys.path.insert(0, str(SCRIPT_DIR))
    import reference_cache  # pylint: disable=C0415
    return reference_cache.default_cache_dir() / "reference" / \
        f"{REFERENCE_COMMIT[:12]}-{product_version()}"


def tree_digest(root):
    """SHA-256 over the relative names and contents of the extracted package."""
    package = Path(root) / "hicexplorer"
    digest = hashlib.sha256()
    for path in sorted(package.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        digest.update(str(path.relative_to(root)).encode())
        digest.update(b"\0")
        with open(path, "rb") as handle:
            digest.update(hashlib.sha256(handle.read()).hexdigest().encode())
    return digest.hexdigest()


def _has_commit(repository):
    if not repository or not Path(repository).exists():
        return False
    return subprocess.run(["git", "-C", str(repository), "cat-file", "-e",
                           REFERENCE_COMMIT + "^{commit}"],
                          capture_output=True, check=False).returncode == 0


def _source_repository(scratch):
    for candidate in (os.environ.get("HICX_REFERENCE_GIT"), str(REPO_ROOT)):
        if _has_commit(candidate):
            return candidate
    fetched = Path(scratch) / "upstream.git"
    subprocess.run(["git", "init", "--quiet", "--bare", str(fetched)], check=True)
    subprocess.run(["git", "-C", str(fetched), "fetch", "--quiet", "--depth", "1",
                    UPSTREAM_URL, REFERENCE_COMMIT], check=True)
    return str(fetched)


def _stamp_document():
    return {"format": TREE_FORMAT, "commit": REFERENCE_COMMIT,
            "overlay": {path: hashlib.sha256(text.encode()).hexdigest()
                        for path, text in overlay().items()}}


def create(root):
    """Extracts the reference tree into root, which must not exist."""
    root = Path(root)
    root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".reference-", dir=str(root.parent)))
    try:
        repository = _source_repository(staging)
        archive = staging / "reference.tar"
        with open(archive, "wb") as handle:
            subprocess.run(["git", "-C", repository, "archive", "--format=tar",
                            REFERENCE_COMMIT, "hicexplorer",
                            ":(exclude)hicexplorer/test"], stdout=handle, check=True)
        tree = staging / "tree"
        tree.mkdir()
        with tarfile.open(archive) as tar:
            tar.extractall(tree, filter="data")
        for relative, text in overlay().items():
            (tree / relative).write_text(text)
        stamp = _stamp_document()
        stamp["digest"] = tree_digest(tree)
        (tree / STAMP).write_text(json.dumps(stamp, indent=1, sort_keys=True) + "\n")
        try:
            os.rename(tree, root)
        except OSError:
            if not (root / STAMP).is_file():   # a concurrent creator won the race
                raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return root


_VERIFIED = {}


def ensure(root=None):
    """The reference root, created when missing and verified once per process."""
    root = Path(root).expanduser().resolve() if root else default_root().resolve()
    if str(root) in _VERIFIED:
        return root
    if not (root / STAMP).is_file():
        if root.exists() and any(root.iterdir()):
            raise RuntimeError(f"reference root {root} exists without {STAMP}")
        if root.exists():
            root.rmdir()
        create(root)
    stamp = json.loads((root / STAMP).read_text())
    expected = _stamp_document()
    for key in ("format", "commit", "overlay"):
        if stamp.get(key) != expected[key]:
            raise RuntimeError(f"reference root {root}: {key} {stamp.get(key)!r} "
                               f"differs from {expected[key]!r}; remove it to recreate")
    if stamp.get("digest") != tree_digest(root):
        raise RuntimeError(f"reference root {root}: files changed since extraction; "
                           "remove it to recreate")
    _VERIFIED[str(root)] = stamp["digest"]
    return root


def digest_of(root):
    """The verified digest of the tree at root (for the reference cache key)."""
    root = ensure(root)
    return _VERIFIED[str(root)]


def environment(root, base=None):
    """base (default os.environ) with the reference root first on PYTHONPATH."""
    env = dict(os.environ if base is None else base)
    inherited = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(root) + (os.pathsep + inherited if inherited else "")
    env["HICX_REFERENCE_ROOT"] = str(root)
    return env


def check_imported():
    """Fails unless the imported `hicexplorer` is the reference tree's.

    Called by every reference entry point (reference_bin/*, py_*.py), which
    the harness runs with HICX_REFERENCE_ROOT set."""
    import hicexplorer  # pylint: disable=C0415
    root = os.environ.get("HICX_REFERENCE_ROOT")
    origin = Path(hicexplorer.__file__).resolve()
    if not root or Path(root).resolve() not in origin.parents:
        raise SystemExit(f"reference run imported hicexplorer from {origin}, "
                         f"not from the reference root {root!r} (cpp/scripts/reference_tree.py)")


def run_tool(tool):
    """Entry point of reference_bin/<tool>: the reference tool's main()."""
    check_imported()
    import importlib  # pylint: disable=C0415
    module = importlib.import_module("hicexplorer." + tool)
    return module.main()


def main(argv=None):
    import argparse  # pylint: disable=C0415
    parser = argparse.ArgumentParser(description="Create or verify the reference tree.")
    parser.add_argument("--reference-root", default=None)
    options = parser.parse_args(argv)
    root = ensure(options.reference_root)
    print(f"{root}\t{REFERENCE_COMMIT}\t{_VERIFIED[str(root)]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
