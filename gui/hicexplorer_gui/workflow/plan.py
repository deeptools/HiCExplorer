"""Per-step execution plan shared by the runner and the exporters."""

import os

from .cmdline import build_argv, shell_join
from .model import is_under, iter_strings


class StepPlan:
    """Everything needed to run or export one step.

    ``argv`` is the exact command line, ``output_kinds`` maps each declared
    output to 'file', 'directory' or 'prefix' (from the spec argument that
    receives it), ``inputs`` lists (command-line path, kind) of every file the
    step reads, and ``mkdirs`` the workdir-relative directories to create first.
    """

    def __init__(self, workflow, step, spec):
        self.step = step
        self.spec = spec
        self.argv = build_argv(workflow, step, spec)
        self.command = shell_join(self.argv)
        self.args = workflow.resolved_args(step)
        by_dest = spec.arguments_by_dest(step.subcommand)

        self.output_kinds = {name: None for name in step.outputs}
        # Declared outputs by absolute path: relative to the workdir, or
        # absolute outside it for a workflow loaded with external_outputs.
        path_to_name = {os.path.normpath(os.path.join(workflow.workdir, rel)): name
                        for name, rel in step.outputs.items()}
        inputs = {}
        output_paths = []
        for dest, value in self.args.items():
            fileinfo = (by_dest.get(dest) or {}).get("file") or {}
            kind = fileinfo.get("kind") or "file"
            for path in iter_strings(value):
                absolute = workflow.workdir_path(path)
                if fileinfo.get("role") == "output":
                    output_paths.append(absolute)
                    key = os.path.normpath(absolute)
                    if key in path_to_name and self.output_kinds[path_to_name[key]] is None:
                        self.output_kinds[path_to_name[key]] = kind
                elif fileinfo.get("role") == "input":
                    inputs.setdefault(workflow.arg_path(absolute), kind)
        # Referenced workflow inputs and upstream outputs count as inputs
        # whatever the spec says about the argument.
        for _dest, ref in step.references:
            if ref.kind == "inputs" and ref.name in workflow.inputs:
                inputs.setdefault(workflow.arg_path(workflow.inputs[ref.name]), "file")
            elif ref.kind == "steps":
                other = workflow.step(ref.step)
                if other is not None and ref.name in other.outputs:
                    inputs.setdefault(other.outputs[ref.name], None)
        self.inputs = inputs  # path -> kind (None: decided from the upstream plan)

        for name, rel in step.outputs.items():
            output_paths.append(os.path.join(workflow.workdir, rel))
        dirs = []
        for absolute in output_paths:
            parent = os.path.dirname(absolute)
            if is_under(parent, workflow.workdir) and parent != workflow.workdir:
                rel = os.path.relpath(parent, workflow.workdir)
                if rel not in dirs:
                    dirs.append(rel)
        self.mkdirs = dirs

    def output_kind(self, name):
        return self.output_kinds.get(name) or "file"


def make_plans(workflow, loader):
    """StepPlan per step id; specs must load (validate first)."""
    plans = {}
    for step in workflow.steps:
        plans[step.id] = StepPlan(workflow, step, loader.load(step.tool))
    # Resolve the kinds of upstream outputs used as inputs.
    for plan in plans.values():
        for path, kind in list(plan.inputs.items()):
            if kind is None:
                resolved = "file"
                for other in plans.values():
                    for name, rel in other.step.outputs.items():
                        if rel == path:
                            resolved = other.output_kind(name)
                plan.inputs[path] = resolved
    return plans
