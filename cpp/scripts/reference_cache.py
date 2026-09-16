"""The cache of Python reference results for equiv.py.

A full regression spends hours of CPU in the Python reference and minutes in
the C++ port, and the reference rarely changes between two runs. The cache
stores what one run of the reference side of a case produced, so that the
next run compares the C++ against it instead of running the Python again.

What a verdict needs from the Python side, and so what an entry holds:

  * the whole working directory the reference ran in (out_py, the class EN
    noise runs, files written relative to the working directory, the RSS
    trace of a large case), which covers every file a comparator, a validator
    or the memory gate reads, plus stdout and stderr;
  * the measurement of the reference run: exit status, wall and CPU time,
    peak RSS, so the memory and time gates evaluate exactly as with a fresh
    run;
  * the noise record of a class EN case (every run's exit status).

The key is the SHA-256 of a JSON document of everything that can change a
Python result (key_document):

  * the cache format version;
  * the case id, tool, Python script and the argument list with the {data}
    and {out} placeholders unexpanded, the noise source and the number of
    reference runs;
  * the SHA-256 of the content of every input the arguments name under
    {data}, directories expanded, and of every file a pyGenomeTracks .ini
    input names;
  * the SHA-256 of the reference code the case runs: its entry script (bin/
    <tool> or the case's py_script) and every module of the repository it
    imports, found by following the import statements; every non-Python file
    of the hicexplorer package; for a case with a noise source or validators,
    whose child processes run further tools, the whole package, noise_runner.py
    and the validators;
  * the interpreter: its version and the versions of the scientific packages,
    and the content of single-file modules the reference imports from outside
    the repository (fit_nbinom);
  * the environment the harness sets for the Python side, with the
    repository path written as <repo>, and the inherited PYTHONPATH.

The Python side runs in a working directory whose path equiv.py pads to a
fixed length (make_workdir), whatever --tmpdir is. A file that embeds that
path (a command line in a log, an absolute output path in an HDF5 attribute)
is restored with the new path written over the old one, byte for byte; a
restore into a path of a different length (a --tmpdir too long to pad) is a
miss.

The key hashes the content of the inputs, not where they are, so another
checkout with the same data and code finds the same entry. Tools that print
an input path (hicInfo, hicPlotSVL, chicQualityControl, hicValidateLocations)
embed the data directory or the repository in their outputs, and an entry
restored as it was stored would then compare the Python output of the old
checkout with the C++ output of the new one. So store records, for the data
root and the repository root, which files embed them, and restore writes the
current roots over the recorded ones in those files, both in one pass so that
a new root containing an old one is not rewritten twice. A text file may
change length; a file holding a NUL byte is rewritten only when the length
stays the same, and otherwise the entry is a miss. Format version 3 added this
record, so an entry stored without it is never restored.

Layout under the cache directory:

  entries/<key[:2]>/<key>/meta.json     key document, measurement, provenance
  entries/<key[:2]>/<key>/workdir/      the reference working directory
  interpreters/<sha of path>.json       the interpreter fingerprint and the
                                        stamps that tell when to take it again
  history.json                          the last measured peaks and times per
                                        case, for the scheduler
"""

import ast
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

CACHE_FORMAT_VERSION = 3

# Distribution names whose versions enter every key.
SCIENTIFIC_PACKAGES = (
    "numpy", "scipy", "pandas", "matplotlib", "cooler", "HiCMatrix", "h5py",
    "pyGenomeTracks", "krbalancing", "hic2cool", "tables", "scikit-learn",
    "statsmodels", "pysam", "pyBigWig", "pybedtools", "intervaltree", "biopython",
    "hicstraw", "numexpr", "blosc", "scikit-image", "networkx", "jinja2", "Pillow",
)
# Modules imported from outside the repository without distribution
# metadata; their source enters the key.
EXTERNAL_MODULES = ("fit_nbinom",)

# Files of the working directory that are not part of an entry.
SKIPPED_SUFFIXES = (".time",)
EMBEDDED_PATH_SCAN_LIMIT = 256 * 1024 * 1024

# Reference artefacts the validators create after the tools ran, stored with
# the entry once the case has finished. Relative to the working directory.
VALIDATOR_ARTEFACTS = ("fitted_distributions.npz", "downstream/python",
                       "downstream/python_run2")


def default_cache_dir():
    configured = os.environ.get("HICX_EQUIV_CACHE")
    if configured:
        return Path(configured)
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return Path(base) / "hicexplorer-equiv"


# Content hashes by (path, size, mtime), for the lifetime of one harness run,
# in which the inputs do not change.
_HASH_MEMO = {}


def clear_hash_memo():
    _HASH_MEMO.clear()


def _sha256_file(path):
    stat = os.stat(path)
    token = (str(path), stat.st_size, stat.st_mtime_ns)
    if token not in _HASH_MEMO:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        _HASH_MEMO[token] = digest.hexdigest()
    return _HASH_MEMO[token]


def hash_path(path):
    """The content hash of a file, or of a directory with its relative names."""
    path = Path(path)
    if path.is_file():
        return _sha256_file(path)
    if path.is_dir():
        digest = hashlib.sha256()
        for item in sorted(p for p in path.rglob("*") if p.is_file()):
            digest.update(str(item.relative_to(path)).encode())
            digest.update(b"\0")
            digest.update(_sha256_file(item).encode())
        return "dir:" + digest.hexdigest()
    return "absent"


# --------------------------------------------------------------------------
# inputs


def _ini_references(path):
    """Existing files a pyGenomeTracks .ini names (file = ..., *_file = ...)."""
    found = []
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return found
    for line in text.splitlines():
        if "=" not in line or line.lstrip().startswith(("#", ";", "[")):
            continue
        value = line.split("=", 1)[1].strip()
        for token in value.split():
            candidate = Path(token) if os.path.isabs(token) else Path(path).parent / token
            if candidate.is_file():
                found.append(candidate)
    return found


def input_hashes(args, data):
    """{argument token: content hash} for every argument naming a {data} path."""
    hashes = {}
    for argument in args:
        if "{data}" not in argument:
            continue
        token = argument[argument.index("{data}"):]
        relative = token.split("::", 1)[0]
        path = Path(relative.replace("{data}", str(data)))
        hashes[token] = hash_path(path)
        if path.suffix == ".ini" and path.is_file():
            for reference in _ini_references(path):
                hashes[f"{token} -> {reference.relative_to(path.parent) if reference.is_relative_to(path.parent) else reference}"] = hash_path(reference)
    return hashes


# --------------------------------------------------------------------------
# reference code


def _module_file(name, roots):
    parts = name.split(".")
    for root in roots:
        base = Path(root).joinpath(*parts)
        for candidate in (base.with_suffix(".py"), base / "__init__.py"):
            if candidate.is_file():
                return candidate
    return None


def _imported_names(path, package):
    try:
        tree = ast.parse(Path(path).read_text(errors="replace"), filename=str(path))
    except (SyntaxError, ValueError):
        return set()
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                names.update(".".join(parts[:i]) for i in range(1, len(parts) + 1))
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.split(".") if package else []
                base = base[:len(base) - (node.level - 1)] if node.level > 1 else base
                module = ".".join(base + ([node.module] if node.module else []))
            else:
                module = node.module or ""
            if module:
                parts = module.split(".")
                names.update(".".join(parts[:i]) for i in range(1, len(parts) + 1))
            for alias in node.names:
                names.add(f"{module}.{alias.name}" if module else alias.name)
    return names


def import_closure(entry_scripts, repo_root):
    """Every file of the repository the entry scripts import, transitively."""
    repo_root = Path(repo_root)
    files, queue = set(), [Path(p) for p in entry_scripts if Path(p).is_file()]
    while queue:
        path = queue.pop()
        if path in files:
            continue
        files.add(path)
        try:
            relative = path.relative_to(repo_root)
        except ValueError:
            relative = None
        package = ""
        if relative is not None and relative.parts and relative.parts[0] == "hicexplorer":
            parts = list(relative.with_suffix("").parts)
            package = ".".join(parts if path.name == "__init__.py" else parts[:-1])
            if path.name == "__init__.py":
                package = ".".join(parts[:-1])
        roots = [repo_root, path.parent]
        for name in _imported_names(path, package):
            target = _module_file(name, roots)
            if target is None:
                continue
            try:
                target.relative_to(repo_root)
            except ValueError:
                continue
            if "test" in target.relative_to(repo_root).parts[:2]:
                continue
            queue.append(target)
    return files


def _package_files(repo_root, python):
    package = Path(repo_root) / "hicexplorer"
    for path in sorted(package.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        if path.relative_to(package).parts[0] == "test":
            continue
        if (path.suffix == ".py") == python:
            yield path


def code_hash(case, repo_root, entry_script):
    repo_root = Path(repo_root)
    scripts = Path(__file__).resolve().parent
    files = set(_package_files(repo_root, python=False))
    if case.get("noise") or case.get("validators"):
        files.update(_package_files(repo_root, python=True))
        files.update((repo_root / "bin").glob("*"))
        files.add(scripts / "noise_runner.py")
        files.update((scripts / "validators").glob("*.py"))
    files.update(import_closure([entry_script], repo_root))
    digest = hashlib.sha256()
    for path in sorted(files):
        if not path.is_file():
            continue
        try:
            name = str(path.relative_to(repo_root))
        except ValueError:
            name = str(path)
        digest.update(name.encode())
        digest.update(b"\0")
        digest.update(_sha256_file(path).encode())
    return digest.hexdigest()


# --------------------------------------------------------------------------
# interpreter

_FINGERPRINT_SCRIPT = r"""
import hashlib, importlib.metadata as m, importlib.util as u, json, sys
packages, modules = {}, {}
for name in sys.argv[1].split(","):
    try:
        packages[name] = m.version(name)
    except m.PackageNotFoundError:
        packages[name] = None
for name in sys.argv[2].split(","):
    spec = u.find_spec(name)
    if spec is None or not spec.origin or not spec.origin.endswith(".py"):
        modules[name] = None
    else:
        with open(spec.origin, "rb") as handle:
            modules[name] = hashlib.sha256(handle.read()).hexdigest()
print(json.dumps({"python": sys.version, "executable_prefix": sys.prefix,
                  "packages": packages, "modules": modules,
                  "path": [p for p in sys.path if p]}))
"""


def _stamps(python, paths, repo_root):
    stamps = {}
    real = os.path.realpath(python)
    for candidate in [real, os.path.join(os.path.dirname(os.path.dirname(python)), "pyvenv.cfg")] + list(paths):
        if repo_root and str(candidate).startswith(str(repo_root)):
            continue
        try:
            stamps[str(candidate)] = os.stat(candidate).st_mtime_ns
        except OSError:
            stamps[str(candidate)] = None
    return stamps


def interpreter_fingerprint(python, env, cache_dir, repo_root):
    """The interpreter version and package versions, taken again only when
    the interpreter or one of its sys.path directories changed."""
    memo = Path(cache_dir) / "interpreters" / (
        hashlib.sha256(f"{python}\0{env.get('PYTHONPATH', '')}".encode()).hexdigest() + ".json")
    if memo.is_file():
        try:
            saved = json.loads(memo.read_text())
            if _stamps(python, saved["fingerprint"]["path"], repo_root) == saved["stamps"]:
                return saved["fingerprint"]
        except (ValueError, KeyError, OSError):
            pass
    completed = subprocess.run(
        [str(python), "-c", _FINGERPRINT_SCRIPT, ",".join(SCIENTIFIC_PACKAGES),
         ",".join(EXTERNAL_MODULES)],
        env=env, capture_output=True, text=True, check=True)
    fingerprint = json.loads(completed.stdout)
    memo.parent.mkdir(parents=True, exist_ok=True)
    memo.write_text(json.dumps({"fingerprint": fingerprint,
                                "stamps": _stamps(python, fingerprint["path"], repo_root)}))
    return fingerprint


# --------------------------------------------------------------------------
# key


def key_document(*, case, runs, data, repo_root, entry_script, fingerprint, environment):
    fingerprint = {key: value for key, value in fingerprint.items()
                   if key not in ("path", "executable_prefix")}
    return {
        "format": CACHE_FORMAT_VERSION,
        "id": case["id"],
        "tool": case.get("tool"),
        "py_script": case.get("py_script"),
        "args": list(case["args"]),
        "noise": case.get("noise"),
        "runs": runs,
        "validators": sorted(entry.get("name", "") for entry in case.get("validators", [])),
        "large": bool(case.get("large")),
        "inputs": input_hashes(case["args"], data),
        "code": code_hash(case, repo_root, entry_script),
        "interpreter": fingerprint,
        "environment": environment,
    }


def key_of(document):
    return hashlib.sha256(json.dumps(document, sort_keys=True).encode()).hexdigest()


# --------------------------------------------------------------------------
# entries


class Cache:
    def __init__(self, directory):
        self.directory = Path(directory)

    def entry_dir(self, key):
        return self.directory / "entries" / key[:2] / key

    def lookup(self, key):
        meta = self.entry_dir(key) / "meta.json"
        if not meta.is_file():
            return None
        try:
            return json.loads(meta.read_text())
        except (ValueError, OSError):
            return None

    def store(self, key, document, workdir, measurement, noise_record, excluded=("out_cpp",),
              roots=None):
        """Copies the reference side of workdir into a new entry, atomically.

        roots: {"data": path, "repo": path}, the roots the run's outputs may
        embed; which files do is recorded for restore."""
        target = self.entry_dir(key)
        staging = Path(tempfile.mkdtemp(prefix=f".{key[:8]}-", dir=self._ensure(target.parent)))
        tree = staging / "workdir"
        workdir = Path(workdir)
        tree.mkdir()
        for item in workdir.iterdir():
            if item.name in excluded:
                continue
            _copy(item, tree / item.name)
        embedding = _files_embedding(tree, str(workdir).encode())
        embedded_roots = {}
        for name, root in sorted((roots or {}).items()):
            files = _files_embedding(tree, str(root).encode())
            if files:
                embedded_roots[name] = {"path": str(root), "files": files}
        meta = {
            "key": key,
            "document": document,
            "workdir": str(workdir),
            "embedded_workdir_files": embedding,
            "embedded_roots": embedded_roots,
            "measurement": dict(measurement),
            "noise": noise_record,
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "harness_host": os.uname().nodename,
            "size_bytes": sum(p.stat().st_size for p in tree.rglob("*") if p.is_file()),
        }
        (staging / "meta.json").write_text(json.dumps(meta, indent=1))
        if target.exists():
            shutil.rmtree(target, ignore_errors=True)
        os.replace(staging, target)
        return meta

    def add_artefacts(self, key, workdir, relative_paths):
        """Adds validator reference artefacts to an existing entry."""
        target = self.entry_dir(key)
        meta_path = target / "meta.json"
        if not meta_path.is_file():
            return []
        added = []
        for relative in relative_paths:
            source = Path(workdir) / relative
            destination = target / "workdir" / relative
            if source.exists() and not destination.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                _copy(source, destination)
                added.append(relative)
        if added:
            meta = json.loads(meta_path.read_text())
            tree = target / "workdir"
            meta["embedded_workdir_files"] = _files_embedding(tree, meta["workdir"].encode())
            meta["size_bytes"] = sum(p.stat().st_size for p in tree.rglob("*") if p.is_file())
            meta.setdefault("artefacts", []).extend(added)
            meta_path.write_text(json.dumps(meta, indent=1))
        return added

    def restore(self, meta, workdir, roots=None):
        """Copies an entry into workdir with the old working directory written
        over by the new one, and the recorded data and repository roots by the
        current ones (roots, as store takes them). False, with nothing
        restored, when a working directory path or a root embedded in a
        binary file would change length."""
        old = meta["workdir"].encode()
        new = str(workdir).encode()
        if meta["embedded_workdir_files"] and len(old) != len(new):
            return False
        tree = self.entry_dir(meta["key"]) / "workdir"
        replacements = {}
        files = set()
        for name, entry in (meta.get("embedded_roots") or {}).items():
            current = (roots or {}).get(name)
            if current is None or str(current) == entry["path"]:
                continue
            replacements[entry["path"].encode()] = str(current).encode()
            files.update(entry["files"])
        for relative in sorted(files):
            source = tree / relative
            if b"\0" in source.read_bytes() and any(len(a) != len(b) for a, b in replacements.items()):
                return False
        workdir = Path(workdir)
        for item in tree.iterdir():
            _copy(item, workdir / item.name)
        for relative in meta["embedded_workdir_files"]:
            path = workdir / relative
            content = path.read_bytes()
            path.write_bytes(content.replace(old, new))
        if replacements:
            # One pass, longest match first: a new root may contain an old one.
            pattern = re.compile(b"|".join(re.escape(a) for a in
                                           sorted(replacements, key=len, reverse=True)))
            for relative in sorted(files):
                path = workdir / relative
                path.write_bytes(pattern.sub(lambda match: replacements[match.group(0)],
                                             path.read_bytes()))
        return True

    def entries(self):
        root = self.directory / "entries"
        if not root.is_dir():
            return
        for meta_path in root.glob("*/*/meta.json"):
            try:
                yield json.loads(meta_path.read_text())
            except (ValueError, OSError):
                continue

    def remove(self, key):
        shutil.rmtree(self.entry_dir(key), ignore_errors=True)

    # history of measurements, for the scheduler

    def history(self):
        path = self.directory / "history.json"
        try:
            return json.loads(path.read_text())
        except (ValueError, OSError):
            return {}

    def update_history(self, records):
        if not records:
            return
        path = self.directory / "history.json"
        self._ensure(self.directory)
        history = self.history()
        history.update(records)
        staging = path.with_suffix(f".tmp{os.getpid()}")
        staging.write_text(json.dumps(history, indent=1, sort_keys=True))
        os.replace(staging, path)

    @staticmethod
    def _ensure(directory):
        Path(directory).mkdir(parents=True, exist_ok=True)
        return str(directory)


def _copy(source, destination):
    source, destination = Path(source), Path(destination)
    if source.is_dir():
        shutil.copytree(source, destination, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns(*[f"*{s}" for s in SKIPPED_SUFFIXES]))
    elif not source.name.endswith(SKIPPED_SUFFIXES):
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _files_embedding(tree, needle):
    found = []
    for path in sorted(p for p in Path(tree).rglob("*") if p.is_file()):
        if path.stat().st_size > EMBEDDED_PATH_SCAN_LIMIT:
            continue
        if needle in path.read_bytes():
            found.append(str(path.relative_to(tree)))
    return found


def directory_size(directory):
    return sum(p.stat().st_size for p in Path(directory).rglob("*") if p.is_file())


if __name__ == "__main__":
    sys.exit("reference_cache.py is a module of equiv.py")
