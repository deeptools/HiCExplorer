"""Case validators: checks that need more than the two output files.

A comparator answers whether two files agree. A validator answers a question
about a case that the files alone cannot settle, for example whether a fitted
model reaches the likelihood of the reference's fits on the data it was fitted
to, or whether a downstream Python tool calls the same set of regions from it.

A case declares validators in its JSON:

    "validators": [{"name": "chic_background_likelihood",
                    "output": "{out}/background.txt"}]

and every validator module exposes

    validate(context) -> comparators.base.Result

where context is a dict with the keys

    case            the case, as loaded
    data            the test data directory
    output          the declared output name, relative to the work directory
    out_py          the Python run's output directory
    out_cpp         the C++ run's output directory
    noise_dirs      the output directories of the further reference runs
    args_py         the case arguments expanded for out_py
    args_cpp        the case arguments expanded for out_cpp
    workdir         the case work directory, for scratch files
    py_python       the interpreter of the Python reference
    env             the environment the Python reference runs in
    repo_root       the repository root

Validators run after the comparators, outside the timed and memory measured
runs, and a failing validator fails the case exactly as a comparator does.
"""
from __future__ import annotations

from importlib import import_module

_KNOWN = {
    "chic_background_likelihood": "chic_background_model",
    "chic_background_downstream": "chic_background_model",
    "differential_analysis_table": "differential_analysis",
}


def names():
    return sorted(_KNOWN)


def run(name, context):
    from comparators.base import fail  # pylint: disable=C0415

    if name not in _KNOWN:
        return fail(f"no validator named {name!r}; known: {names()}")
    module = import_module(f"validators.{_KNOWN[name]}")
    return getattr(module, name)(context)
