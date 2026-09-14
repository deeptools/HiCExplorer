"""Command line construction from a step and its tool specification."""

import shlex

from .spec import NON_SETTABLE_ACTIONS, long_flag


def format_scalar(value):
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def step_values(workflow, step, spec):
    """Resolved argument values including the injected thread count."""
    values = workflow.resolved_args(step)
    dest = spec.thread_dest(step.subcommand, step.threads_option)
    if dest is not None and dest not in values:
        values[dest] = step.threads
    return values


def argv_for_values(spec, subcommand, values, executable=None):
    """The argv list for resolved values.

    The subcommand comes first, then positionals, then options, each group in
    the spec's argument order. Positionals precede options so that a
    positional cannot be taken as a further value of a variable-length option.
    store_true/store_false options are emitted only when the value differs from
    the default; ``null`` omits an argument.
    """
    argv = [executable or spec.tool]
    if subcommand is not None:
        argv.append(subcommand)
    positionals, options = [], []
    for arg in spec.arguments(subcommand):
        dest = arg["dest"]
        action = arg.get("action") or "store"
        if dest not in values or action in NON_SETTABLE_ACTIONS:
            continue
        value = values[dest]
        if value is None:
            continue
        if arg.get("positional"):
            items = value if isinstance(value, list) else [value]
            positionals.extend(format_scalar(v) for v in items)
            continue
        flag = long_flag(arg)
        if action in ("store_true", "store_false"):
            if value != arg.get("default"):
                options.append(flag)
            continue
        if action == "append":
            for item in (value if isinstance(value, list) else [value]):
                options.extend([flag] + ([format_scalar(v) for v in item]
                                         if isinstance(item, list) else [format_scalar(item)]))
            continue
        if isinstance(value, list):
            options.append(flag)
            options.extend(format_scalar(v) for v in value)
        else:
            options.extend([flag, format_scalar(value)])
    return argv + positionals + options


def build_argv(workflow, step, spec, executable=None):
    return argv_for_values(spec, step.subcommand, step_values(workflow, step, spec),
                           executable or spec.executable)


def shell_join(argv):
    return " ".join(shlex.quote(a) for a in argv)
