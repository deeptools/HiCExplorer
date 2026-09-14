"""Tool specifications (``<tool> --help-json``, schema ``hicexplorer-tool-spec`` v1)."""

import json
import os
import shutil
import subprocess

SCHEMA = "hicexplorer-tool-spec"
SCHEMA_VERSION = 1

# Actions that cannot be set from a workflow file.
NON_SETTABLE_ACTIONS = ("help", "version", "help_json")
# Option dests recognised as the thread count when a step does not name one.
THREAD_DESTS = ("threads", "numberOfProcessors", "processes")


class SpecError(Exception):
    """A tool could not be found or its specification could not be loaded."""


class ToolSpec:
    """Parsed ``--help-json`` output of one tool."""

    def __init__(self, data, executable=None):
        if not isinstance(data, dict) or data.get("schema") != SCHEMA:
            raise SpecError("not a {} document".format(SCHEMA))
        if data.get("schema_version") != SCHEMA_VERSION:
            raise SpecError("unsupported schema_version {!r}".format(data.get("schema_version")))
        self.data = data
        self.executable = executable
        self.tool = data.get("tool")
        self.version = str(data.get("version", ""))
        self.description = data.get("description", "")
        sub = data.get("subcommands")
        self.subcommand_dest = sub.get("dest") if sub else None
        self.subcommand_required = bool(sub.get("required", True)) if sub else False
        self.commands = {c["name"]: c for c in sub.get("commands", [])} if sub else {}

    @property
    def has_subcommands(self):
        return bool(self.commands)

    def arguments(self, subcommand=None):
        """All arguments in spec order: top-level groups, then the subcommand's."""
        out = []
        for group in self.data.get("groups") or []:
            out.extend(group.get("arguments") or [])
        if subcommand is not None and subcommand in self.commands:
            for group in self.commands[subcommand].get("groups") or []:
                out.extend(group.get("arguments") or [])
        return out

    def arguments_by_dest(self, subcommand=None):
        return {a["dest"]: a for a in self.arguments(subcommand)}

    def mutually_exclusive(self, subcommand=None):
        groups = list(self.data.get("mutually_exclusive") or [])
        if subcommand is not None and subcommand in self.commands:
            groups.extend(self.commands[subcommand].get("mutually_exclusive") or [])
        return groups

    def thread_dest(self, subcommand=None, requested=None):
        """The dest receiving the step's thread count, or None."""
        by_dest = self.arguments_by_dest(subcommand)
        if requested is not None:
            return requested if requested in by_dest else None
        for dest in THREAD_DESTS:
            arg = by_dest.get(dest)
            if arg is not None and arg.get("action", "store") == "store":
                return dest
        return None


def long_flag(arg):
    """The first long flag of an option (the first flag when none is long)."""
    flags = arg.get("flags") or []
    for flag in flags:
        if flag.startswith("--"):
            return flag
    return flags[0] if flags else None


def resolve_executable(tool, tools_dir=None):
    """Absolute path of a tool: ``tools_dir``, else ``$HICX_CPP_BIN``, else PATH."""
    if os.sep in tool or (os.altsep and os.altsep in tool):
        raise SpecError("tool name {!r} must not contain a path separator".format(tool))
    directory = tools_dir if tools_dir else os.environ.get("HICX_CPP_BIN") or None
    if directory:
        path = os.path.abspath(os.path.join(directory, tool))
        if not (os.path.isfile(path) and os.access(path, os.X_OK)):
            raise SpecError("tool {} not found as an executable in {}".format(tool, directory))
        return path
    found = shutil.which(tool)
    if not found:
        raise SpecError("tool {} not found on PATH (set --tools-dir or HICX_CPP_BIN)".format(tool))
    return os.path.abspath(found)


class SpecLoader:
    """Resolves tools and loads their specifications, cached per instance."""

    def __init__(self, tools_dir=None, timeout=60):
        self.tools_dir = tools_dir
        self.timeout = timeout
        self._cache = {}

    def executable(self, tool):
        return resolve_executable(tool, self.tools_dir)

    def load(self, tool):
        if tool in self._cache:
            cached = self._cache[tool]
            if isinstance(cached, SpecError):
                raise cached
            return cached
        try:
            spec = self._load(tool)
        except SpecError as exc:
            self._cache[tool] = exc
            raise
        self._cache[tool] = spec
        return spec

    def _load(self, tool):
        exe = self.executable(tool)
        try:
            proc = subprocess.run([exe, "--help-json"], stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, timeout=self.timeout)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise SpecError("{} --help-json failed: {}".format(tool, exc))
        if proc.returncode != 0:
            raise SpecError("{} --help-json exited {}: {}".format(
                tool, proc.returncode, proc.stderr.decode(errors="replace").strip()[:200]))
        try:
            data = json.loads(proc.stdout.decode())
        except ValueError as exc:
            raise SpecError("{} --help-json printed invalid JSON: {}".format(tool, exc))
        spec = ToolSpec(data, executable=exe)
        if spec.tool and spec.tool != tool:
            raise SpecError("{} --help-json describes tool {!r}".format(tool, spec.tool))
        return spec
