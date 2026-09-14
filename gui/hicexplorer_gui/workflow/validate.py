"""Validation of a workflow against the tool specifications."""

import glob
import os

from .model import is_under, iter_strings
from .spec import NON_SETTABLE_ACTIONS, SpecError


class Message:
    __slots__ = ("level", "text")

    def __init__(self, level, text):
        self.level, self.text = level, text

    def __str__(self):
        return "{}: {}".format(self.level, self.text)

    def __repr__(self):
        return "Message({!r}, {!r})".format(self.level, self.text)


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _convert(value, type_name):
    """Convert a scalar to the argument type; raises ValueError."""
    if isinstance(value, (list, dict)) or isinstance(value, bool):
        raise ValueError("expected a {} value, found {!r}".format(type_name or "scalar", value))
    if type_name == "int":
        if _is_int(value):
            return value
        if isinstance(value, float):
            if value.is_integer():
                return int(value)
            raise ValueError("expected an integer, found {!r}".format(value))
        try:
            return int(str(value), 10)
        except ValueError:
            raise ValueError("expected an integer, found {!r}".format(value))
    if type_name == "float":
        if isinstance(value, (int, float)):
            return float(value)
        try:
            return float(str(value))
        except ValueError:
            raise ValueError("expected a number, found {!r}".format(value))
    return value if isinstance(value, str) else str(value) if value is not None else value


def file_exists(path, kind):
    if kind == "directory":
        return os.path.isdir(path)
    if kind == "prefix":
        return bool(glob.glob(glob.escape(path) + "*"))
    return os.path.exists(path)


def validate_workflow(workflow, loader, check_files=True):
    """Return a list of Message (level 'error' or 'info'), in step order."""
    msgs = [Message("error", e) for e in workflow.errors]
    err = lambda text: msgs.append(Message("error", text))
    info = lambda text: msgs.append(Message("info", text))

    for step in workflow.steps:
        label = "step {}".format(step.id)
        if step.tool is None:
            continue
        try:
            spec = loader.load(step.tool)
        except SpecError as exc:
            err("{}: {}".format(label, exc))
            continue

        if spec.has_subcommands:
            if step.subcommand is None:
                if spec.subcommand_required:
                    err("{}: subcommand: required, one of {}".format(label, ", ".join(spec.commands)))
            elif step.subcommand not in spec.commands:
                err("{}: subcommand: unknown '{}', expected one of {}".format(
                    label, step.subcommand, ", ".join(spec.commands)))
                continue
        elif step.subcommand is not None:
            err("{}: subcommand: tool {} has no subcommands".format(label, step.tool))
            continue

        by_dest = spec.arguments_by_dest(step.subcommand)
        thread_dest = spec.thread_dest(step.subcommand, step.threads_option)
        if step.threads_option is not None and thread_dest is None:
            err("{}: threads_option: unknown argument '{}'".format(label, step.threads_option))
        if thread_dest is not None and thread_dest in step.args:
            err("{}: {}: set the step's threads instead; the engine passes it to the tool".format(
                label, thread_dest))
        if thread_dest is None and step.threads > 1:
            info("{}: threads: {} has no thread option; {} threads are reserved but not passed".format(
                label, step.tool, step.threads))

        ancestors = workflow.ancestors(step.id)
        resolved = {}
        for dest, value in step.args.items():
            arg = by_dest.get(dest)
            if arg is None:
                err("{}: {}: unknown argument of {}{}".format(
                    label, dest, step.tool, " " + step.subcommand if step.subcommand else ""))
                continue
            action = arg.get("action") or "store"
            if action in NON_SETTABLE_ACTIONS:
                err("{}: {}: action '{}' cannot be set in a workflow".format(label, dest, action))
                continue
            try:
                value = workflow.resolve(step, value)
            except KeyError:
                continue  # already reported as a structural error
            resolved[dest] = value
            if arg.get("cpp_only"):
                note = arg.get("note")
                info("{}: {}: C++-only option{}".format(label, dest, " ({})".format(note) if note else ""))
            _check_value(err, label, dest, arg, action, value)

        for dest, arg in by_dest.items():
            action = arg.get("action") or "store"
            if arg.get("required") and action not in NON_SETTABLE_ACTIONS and dest != thread_dest:
                if step.args.get(dest) is None:
                    err("{}: {}: required argument missing".format(label, dest))

        for group in spec.mutually_exclusive(step.subcommand):
            dests = group.get("dests") or []
            present = [d for d in dests if d in resolved and _is_set(by_dest[d], resolved[d])]
            if len(present) > 1:
                err("{}: {}: mutually exclusive arguments given together".format(label, ", ".join(present)))
            elif group.get("required") and not present:
                err("{}: {}: one of these arguments is required".format(label, ", ".join(dests)))

        if not check_files:
            continue
        for dest, value in resolved.items():
            fileinfo = by_dest[dest].get("file")
            if not fileinfo:
                continue
            kind = fileinfo.get("kind") or "file"
            for path in iter_strings(value):
                absolute = workflow.workdir_path(path)
                rel = os.path.relpath(absolute, workflow.workdir) if is_under(absolute, workflow.workdir) else None
                if fileinfo.get("role") == "output":
                    if not is_under(absolute, workflow.workdir) or absolute == workflow.workdir:
                        err("{}: {}: output {} is not under the workdir".format(label, dest, path))
                    elif rel not in step.outputs.values():
                        info("{}: {}: output {} is not declared in outputs and is not tracked".format(
                            label, dest, path))
                elif fileinfo.get("role") == "input":
                    owner = workflow.producer_of(rel) if rel is not None else None
                    if owner is not None:
                        if owner[0] == step.id:
                            err("{}: {}: input {} is an output of this step".format(label, dest, path))
                        elif owner[0] not in ancestors:
                            err("{}: {}: input {} is output {} of step {}, which is not upstream; "
                                "use ${{steps.{}.outputs.{}}}".format(
                                    label, dest, path, owner[1], owner[0], owner[0], owner[1]))
                    elif not file_exists(absolute, kind):
                        err("{}: {}: input {} does not exist{}".format(
                            label, dest, path, "" if kind == "file" else " ({})".format(kind)))
    return msgs


def _is_set(arg, value):
    if value is None:
        return False
    if (arg.get("action") or "store") in ("store_true", "store_false"):
        return value != arg.get("default")
    return True


def _check_value(err, label, dest, arg, action, value):
    where = "{}: {}".format(label, dest)
    if value is None:
        if arg.get("required"):
            err("{}: required argument is null".format(where))
        return
    if action in ("store_true", "store_false"):
        if not isinstance(value, bool):
            err("{}: expected true or false, found {!r}".format(where, value))
        return
    nargs = arg.get("nargs")
    if action == "append":
        items = value if isinstance(value, list) else [value]
    elif nargs is None or nargs == "?":
        if isinstance(value, (list, dict)):
            err("{}: expected a single value, found a list".format(where))
            return
        items = [value]
    else:
        if isinstance(value, dict):
            err("{}: expected a list, found a mapping".format(where))
            return
        items = value if isinstance(value, list) else [value]
        if nargs == "+" and not items:
            err("{}: expected at least one value".format(where))
        elif _is_int(nargs) and len(items) != nargs:
            err("{}: expected exactly {} values, found {}".format(where, nargs, len(items)))
    type_name = arg.get("type")
    choices = arg.get("choices")
    for item in items:
        if action == "append" and isinstance(item, list):
            continue
        try:
            converted = _convert(item, type_name if type_name in ("int", "float") else "str")
        except ValueError as exc:
            err("{}: {}".format(where, exc))
            continue
        if choices is not None:
            normalised = []
            for choice in choices:
                try:
                    normalised.append(_convert(choice, type_name if type_name in ("int", "float") else "str"))
                except ValueError:
                    normalised.append(choice)
            if converted not in normalised:
                err("{}: invalid choice {!r}, expected one of {}".format(
                    where, item, ", ".join(str(c) for c in choices)))
