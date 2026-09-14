"""Workflow file model: parsing, ``${...}`` references and the step DAG.

Workflow file, version 1::

    version: 1
    name: hic-basic
    threads: 8                  # thread budget (default: number of CPUs)
    inputs:                     # relative to the workflow file's directory
      r1: data/R1.bam
    steps:
      - id: build
        tool: hicBuildMatrix
        subcommand: null        # required for tools with subcommands
        threads: 4              # default 1
        threads_option: null    # dest receiving the thread count (auto-detected)
        after: []               # extra ordering dependencies
        args: {samFiles: ["${inputs.r1}", "${inputs.r2}"], outFileName: "${outputs.matrix}"}
        outputs: {matrix: build/matrix.h5}   # relative to the workdir

References inside argument strings: ``${inputs.NAME}``, ``${outputs.NAME}``
(this step's outputs) and ``${steps.ID.outputs.NAME}``; ``$$`` is a literal
``$``. Tools run with the workdir as their current directory, so literal paths
in ``args`` are relative to the workdir; referenced paths are rendered relative
to the workdir when they lie inside it and absolute otherwise.
"""

import heapq
import os
import re

import yaml

FORMAT_VERSION = 1
STATE_DIR = ".hicexplorer-workflow"

_TOP_KEYS = {"version", "name", "description", "threads", "inputs", "steps"}
_STEP_KEYS = {"id", "tool", "subcommand", "description", "threads", "threads_option",
              "after", "args", "outputs"}
_ID_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")
_REF_RE = re.compile(r"\$\$|\$\{([^}]*)\}")


class WorkflowError(Exception):
    """The workflow file cannot be read at all (not YAML, not a mapping)."""


class Reference:
    """One ``${...}`` reference: kind is 'inputs', 'outputs' or 'steps'."""

    __slots__ = ("kind", "step", "name", "text")

    def __init__(self, kind, name, step=None, text=""):
        self.kind, self.name, self.step, self.text = kind, name, step, text


def parse_references(text):
    """References in a string; raises ValueError on a malformed one."""
    refs = []
    for match in _REF_RE.finditer(text):
        if match.group(0) == "$$":
            continue
        body = match.group(1).strip()
        parts = body.split(".")
        if len(parts) == 2 and parts[0] in ("inputs", "outputs"):
            refs.append(Reference(parts[0], parts[1], text=match.group(0)))
        elif len(parts) == 4 and parts[0] == "steps" and parts[2] == "outputs":
            refs.append(Reference("steps", parts[3], step=parts[1], text=match.group(0)))
        else:
            raise ValueError("malformed reference {}".format(match.group(0)))
    return refs


def iter_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from iter_strings(item)


def is_under(path, directory):
    path = os.path.abspath(path)
    directory = os.path.abspath(directory)
    return path == directory or path.startswith(directory.rstrip(os.sep) + os.sep)


class Step:
    def __init__(self, index):
        self.index = index
        self.id = None
        self.tool = None
        self.subcommand = None
        self.description = ""
        self.threads = 1
        self.threads_option = None
        self.after = []
        self.args = {}
        self.outputs = {}      # name -> path relative to the workdir (normalised)
        self.references = []   # list of (dest, Reference)
        self.deps = set()      # step ids this step depends on

    def __repr__(self):
        return "Step({!r})".format(self.id)


class Workflow:
    def __init__(self, path, workdir=None):
        self.path = os.path.abspath(path)
        self.base_dir = os.path.dirname(self.path)
        self.workdir = os.path.abspath(workdir) if workdir else self.base_dir
        self.version = None
        self.name = ""
        self.description = ""
        self.threads = None
        self.inputs = {}       # name -> absolute path
        self.steps = []
        self.errors = []       # structural errors, one message per problem
        # Absolute output paths outside the workdir are accepted (single tool
        # runs from the GUI); workflow files keep outputs under the workdir.
        self.external_outputs = False

    # -- lookup ---------------------------------------------------------
    def step(self, step_id):
        for step in self.steps:
            if step.id == step_id:
                return step
        return None

    @property
    def state_dir(self):
        return os.path.join(self.workdir, STATE_DIR)

    def budget(self, override=None):
        if override is not None:
            return max(1, int(override))
        if isinstance(self.threads, int) and self.threads > 0:
            return self.threads
        return os.cpu_count() or 1

    # -- paths ----------------------------------------------------------
    def arg_path(self, absolute):
        """How a path appears on a command line run from the workdir."""
        if is_under(absolute, self.workdir):
            return os.path.relpath(absolute, self.workdir)
        return absolute

    def output_abspath(self, step, name):
        return os.path.join(self.workdir, step.outputs[name])

    def workdir_path(self, path):
        """Absolute path of a command-line path (relative to the workdir)."""
        return path if os.path.isabs(path) else os.path.normpath(os.path.join(self.workdir, path))

    def producer_of(self, rel):
        """(step id, output name) of the declared output that yields ``rel``.

        ``rel`` is a workdir-relative path. It matches an output exactly, a
        file inside a directory output, or a file named after a prefix output
        (``tads/t`` yields ``tads/t_domains.bed``).
        """
        rel = os.path.normpath(rel)
        for step in self.steps:
            for name, out in step.outputs.items():
                if rel == out:
                    return step.id, name
        for step in self.steps:
            for name, out in step.outputs.items():
                inside = rel.startswith(out.rstrip(os.sep) + os.sep)
                named = (os.path.dirname(rel) == os.path.dirname(out)
                         and os.path.basename(rel).startswith(os.path.basename(out)))
                if inside or named:
                    return step.id, name
        return None

    # -- references -----------------------------------------------------
    def substitute(self, step, text):
        """Replace the references in one string; unknown references raise KeyError."""
        def repl(match):
            if match.group(0) == "$$":
                return "$"
            ref = parse_references(match.group(0))[0]
            if ref.kind == "inputs":
                return self.arg_path(self.inputs[ref.name])
            if ref.kind == "outputs":
                return step.outputs[ref.name]
            other = self.step(ref.step)
            if other is None:
                raise KeyError(ref.text)
            return other.outputs[ref.name]
        return _REF_RE.sub(repl, text)

    def resolve(self, step, value):
        if isinstance(value, str):
            return self.substitute(step, value)
        if isinstance(value, list):
            return [self.resolve(step, v) for v in value]
        return value

    def resolved_args(self, step):
        return {dest: self.resolve(step, value) for dest, value in step.args.items()}

    # -- DAG ------------------------------------------------------------
    def dependents(self, step_id):
        """All steps downstream of ``step_id`` (transitively)."""
        out, frontier = set(), [step_id]
        while frontier:
            current = frontier.pop()
            for step in self.steps:
                if current in step.deps and step.id not in out:
                    out.add(step.id)
                    frontier.append(step.id)
        return out

    def ancestors(self, step_id):
        out, frontier = set(), [step_id]
        while frontier:
            step = self.step(frontier.pop())
            for dep in (step.deps if step else ()):
                if dep not in out:
                    out.add(dep)
                    frontier.append(dep)
        return out

    def topological_order(self):
        """Steps in dependency order, ties broken by file position.

        Raises WorkflowError naming the steps on a cycle.
        """
        by_id = {s.id: s for s in self.steps}
        indegree = {s.id: len([d for d in s.deps if d in by_id]) for s in self.steps}
        heap = [(s.index, s.id) for s in self.steps if indegree[s.id] == 0]
        heapq.heapify(heap)
        order = []
        while heap:
            _, sid = heapq.heappop(heap)
            order.append(by_id[sid])
            for other in self.steps:
                if sid in other.deps:
                    indegree[other.id] -= 1
                    if indegree[other.id] == 0:
                        heapq.heappush(heap, (other.index, other.id))
        if len(order) != len(self.steps):
            stuck = sorted((s for s in self.steps if s not in order), key=lambda s: s.index)
            raise WorkflowError("dependency cycle among steps: " + ", ".join(s.id for s in stuck))
        return order


def load_workflow(path, workdir=None, external_outputs=False):
    """Parse a workflow file. Structural problems are collected in ``errors``.

    With ``external_outputs`` an output may be an absolute path outside the
    workdir; it is declared, hashed and resumed like any other output. The GUI
    uses this for single tool runs, whose results go wherever the user keeps
    them. Workflow files are loaded without it, so a workflow and its outputs
    stay together in one relocatable work directory.
    """
    try:
        with open(path) as handle:
            data = yaml.safe_load(handle)
    except OSError as exc:
        raise WorkflowError("cannot read {}: {}".format(path, exc))
    except yaml.YAMLError as exc:
        raise WorkflowError("{} is not valid YAML: {}".format(path, exc))
    if not isinstance(data, dict):
        raise WorkflowError("{}: the top level must be a mapping".format(path))
    wf = Workflow(path, workdir)
    wf.external_outputs = bool(external_outputs)
    err = wf.errors.append

    for key in sorted(set(data) - _TOP_KEYS):
        err("workflow: unknown key '{}'".format(key))
    wf.version = data.get("version")
    if wf.version != FORMAT_VERSION:
        err("workflow: version must be {} (found {!r})".format(FORMAT_VERSION, wf.version))
    wf.name = str(data.get("name") or os.path.splitext(os.path.basename(path))[0])
    wf.description = str(data.get("description") or "")
    threads = data.get("threads")
    if threads is not None and (not isinstance(threads, int) or isinstance(threads, bool) or threads < 1):
        err("workflow: threads must be a positive integer (found {!r})".format(threads))
    else:
        wf.threads = threads

    inputs = data.get("inputs") or {}
    if not isinstance(inputs, dict):
        err("workflow: inputs must be a mapping of name to path")
        inputs = {}
    for name, value in inputs.items():
        if not isinstance(value, str) or not value:
            err("input {}: must be a path string".format(name))
            continue
        wf.inputs[str(name)] = os.path.normpath(os.path.join(wf.base_dir, os.path.expanduser(value)))

    steps = data.get("steps")
    if not isinstance(steps, list) or not steps:
        err("workflow: steps must be a non-empty list")
        steps = []
    seen = set()
    for index, raw in enumerate(steps):
        step = Step(index)
        label = "step #{}".format(index + 1)
        if not isinstance(raw, dict):
            err("{}: must be a mapping".format(label))
            continue
        sid = raw.get("id")
        if not isinstance(sid, str) or not _ID_RE.match(sid):
            err("{}: id must match {} (found {!r})".format(label, _ID_RE.pattern, sid))
            sid = "#{}".format(index + 1)
        elif sid in seen:
            err("step {}: duplicate id".format(sid))
        seen.add(sid)
        step.id = sid
        label = "step {}".format(sid)
        for key in sorted(set(raw) - _STEP_KEYS):
            err("{}: unknown key '{}'".format(label, key))
        step.tool = raw.get("tool")
        if not isinstance(step.tool, str) or not step.tool:
            err("{}: tool is required".format(label))
            step.tool = None
        sub = raw.get("subcommand")
        if sub is not None and not isinstance(sub, str):
            err("{}: subcommand must be a string".format(label))
            sub = None
        step.subcommand = sub
        step.description = str(raw.get("description") or "")
        threads = raw.get("threads", 1)
        if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
            err("{}: threads must be a positive integer (found {!r})".format(label, threads))
            threads = 1
        step.threads = threads
        step.threads_option = raw.get("threads_option")
        if step.threads_option is not None and not isinstance(step.threads_option, str):
            err("{}: threads_option must be an argument dest".format(label))
            step.threads_option = None
        after = raw.get("after") or []
        if isinstance(after, str):
            after = [after]
        if not isinstance(after, list) or not all(isinstance(a, str) for a in after):
            err("{}: after must be a list of step ids".format(label))
            after = []
        step.after = after
        args = raw.get("args") or {}
        if not isinstance(args, dict):
            err("{}: args must be a mapping of argument dest to value".format(label))
            args = {}
        step.args = {str(k): v for k, v in args.items()}
        outputs = raw.get("outputs") or {}
        if not isinstance(outputs, dict):
            err("{}: outputs must be a mapping of name to path".format(label))
            outputs = {}
        for name, value in outputs.items():
            if not isinstance(value, str) or not value:
                err("{}: output {} must be a path string".format(label, name))
                continue
            absolute = os.path.normpath(os.path.join(wf.workdir, value))
            if wf.external_outputs and os.path.isabs(value) and not is_under(absolute, wf.workdir):
                step.outputs[str(name)] = absolute
                continue
            if not is_under(absolute, wf.workdir) or absolute == wf.workdir:
                err("{}: output {} ({}) is not under the workdir".format(label, name, value))
                continue
            rel = os.path.relpath(absolute, wf.workdir)
            if rel.split(os.sep)[0] == STATE_DIR:
                err("{}: output {} lies in the engine's state directory".format(label, name))
                continue
            step.outputs[str(name)] = rel
        wf.steps.append(step)

    # References and dependencies.
    ids = {s.id for s in wf.steps}
    for step in wf.steps:
        label = "step {}".format(step.id)
        for dest, value in step.args.items():
            for text in iter_strings(value):
                try:
                    refs = parse_references(text)
                except ValueError as exc:
                    err("{}: {}: {}".format(label, dest, exc))
                    continue
                for ref in refs:
                    step.references.append((dest, ref))
                    if ref.kind == "inputs" and ref.name not in wf.inputs:
                        err("{}: {}: unknown input in {}".format(label, dest, ref.text))
                    elif ref.kind == "outputs" and ref.name not in step.outputs:
                        err("{}: {}: {} is not declared in this step's outputs".format(label, dest, ref.text))
                    elif ref.kind == "steps":
                        other = wf.step(ref.step)
                        if other is None:
                            err("{}: {}: unknown step in {}".format(label, dest, ref.text))
                        elif other is step:
                            err("{}: {}: {} refers to the step itself; use ${{outputs.{}}}".format(
                                label, dest, ref.text, ref.name))
                        elif ref.name not in other.outputs:
                            err("{}: {}: step {} has no output '{}'".format(label, dest, ref.step, ref.name))
                        else:
                            step.deps.add(ref.step)
        for dep in step.after:
            if dep not in ids:
                err("{}: after: unknown step '{}'".format(label, dep))
            elif dep == step.id:
                err("{}: after: a step cannot run after itself".format(label))
            else:
                step.deps.add(dep)

    owners = {}
    for step in wf.steps:
        for name, rel in step.outputs.items():
            if rel in owners and owners[rel] != (step.id, name):
                err("step {}: output {} ({}) is also output {} of step {}".format(
                    step.id, name, rel, owners[rel][1], owners[rel][0]))
            owners.setdefault(rel, (step.id, name))

    try:
        wf.topological_order()
    except WorkflowError as exc:
        err(str(exc))
    return wf
