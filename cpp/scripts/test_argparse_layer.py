"""Differential test of the C++ argument layer (cpp/core/src/argparse.cpp)
against Python's argparse.

cpp/tests/argparse_driver.cpp declares the parsers built below; both sides
parse the same command lines and must agree on exit status, stdout and stderr
(the namespace is printed as JSON with sorted keys).

    HICX_ARGPARSE_DRIVER=BUILD/tests/hicx_argparse_driver python -m pytest cpp/scripts/test_argparse_layer.py

Skipped when HICX_ARGPARSE_DRIVER is not set. Run it with the Python of the
reference environment, whose argparse is the one the tools reproduce.
"""

import argparse
import contextlib
import io
import json
import os
import subprocess

import pytest

DRIVER = os.environ.get("HICX_ARGPARSE_DRIVER")


def main_parser():
    parser = argparse.ArgumentParser(prog="argdriver", usage="argdriver [options]", add_help=False,
                                     description="A parser with every feature the tools use.")
    g = parser.add_argument_group("Options")
    g.add_argument("--int", "-i", type=int, default=3)
    g.add_argument("--float", type=float, default=0.5)
    g.add_argument("--name", choices=["a", "b"], default="a")
    g.add_argument("--list", type=int, nargs="+")
    g.add_argument("--star", nargs="*")
    g.add_argument("--pair", type=float, nargs=2)
    g.add_argument("--opt", nargs="?", const="C", default="D")
    g.add_argument("--flag", action="store_true")
    g.add_argument("--no-thing", action="store_false")
    g.add_argument("--append", action="append")
    g.add_argument("--choice-int", type=int, choices=[1, 2])
    g.add_argument("--required-when", default="x")
    m = g.add_mutually_exclusive_group()
    m.add_argument("--left", action="store_true")
    m.add_argument("--right", action="store_true")
    g.add_argument("pos", nargs="?")
    g.add_argument("--help", "-h", action="help")
    g.add_argument("--version", "-v", action="version", version="%(prog)s 1.0")
    return parser


def sub_parser():
    parser = argparse.ArgumentParser(prog="argdriver", usage="argdriver {run,show} ...", add_help=False)
    parser.add_argument_group("Options").add_argument("--verbose", action="store_true")
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", usage="argdriver run [options]", add_help=False)
    r = run.add_argument_group("Run options")
    r.add_argument("--x", type=int, required=True)
    r.add_argument("--mode", choices=["fast", "slow"], default="fast")
    r.add_argument("--help", "-h", action="help")
    show = subparsers.add_parser("show", usage="argdriver show [options]", add_help=False)
    show.add_argument_group("Show options").add_argument("--y", action="store_true")
    return parser


def python_side(mode, argv):
    parser = main_parser() if mode == "main" else sub_parser()
    out, err = io.StringIO(), io.StringIO()
    code = 0
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            ns = parser.parse_args(argv)
            print(json.dumps(vars(ns), sort_keys=True))
        except SystemExit as exit_:
            code = exit_.code or 0
    return code, out.getvalue(), err.getvalue()


def cpp_side(mode, argv):
    env = dict(os.environ, HICX_ARGDRIVER_MODE=mode)
    result = subprocess.run([DRIVER] + argv, capture_output=True, text=True, env=env)
    return result.returncode, result.stdout, result.stderr


MAIN_CASES = [
    [],
    ["--int", "7", "--float", "2.25", "--name", "b"],
    ["-i", "5"], ["-i5"], ["-i=5"], ["--int=8"], ["--in", "4"], ["--fl"],
    ["--n", "a"],                                   # ambiguous: --name, --no-thing
    ["--int", "-3"], ["--pair", "-1", "-2.5"], ["--float", "-.5"],
    ["--int", " 42 "], ["--int", "1_000"], ["--int", "+5"], ["--int", "1.5"], ["--int", "x"],
    ["--float", "1e3"], ["--float", "nan"], ["--float", "-inf"], ["--float", "1_0.5"], ["--float", "0x10"],
    ["--float", "."], ["--float", "5."], ["--float", ".5e-3"],
    ["--name", "c"], ["--choice-int", "2"], ["--choice-int", "3"],
    ["--list", "1", "2", "3"], ["--list"], ["--list", "1", "x"],
    ["--star"], ["--star", "a", "b"], ["--pair", "1"], ["--pair", "1", "2", "3"],
    ["--opt"], ["--opt", "val"], ["--opt", "--flag"],
    ["--flag", "--no-thing"], ["--flag=yes"],
    ["--append", "a", "--append", "b"],
    ["--int", "1", "--int", "2"],
    ["--left", "--right"], ["--left", "--left"],
    ["positional"], ["positional", "extra"], ["--bogus"], ["--bogus", "x", "--also"],
    ["--", "--flag"], ["--int"], ["--int", "--flag"],
    ["--star", "a", "--", "b"],
    ["-h"], ["--help", "--bogus"], ["--bogus", "--help"], ["-v"], ["--vers"],
]

SUB_CASES = [
    [], ["run", "--x", "3"], ["run"], ["run", "--x", "3", "--mode", "slow"], ["run", "--x", "q"],
    ["show", "--y"], ["show"], ["nope"], ["--verbose", "show"], ["run", "-h"],
    ["run", "--x", "1", "--mode", "medium"], ["show", "--bogus"], ["--bogus", "show"],
]


@pytest.mark.skipif(not DRIVER, reason="HICX_ARGPARSE_DRIVER is not set")
@pytest.mark.parametrize("mode,argv", [("main", c) for c in MAIN_CASES] + [("sub", c) for c in SUB_CASES],
                         ids=lambda v: v if isinstance(v, str) else " ".join(v) or "<none>")
def test_layer_matches_argparse(mode, argv):
    cpp, py = cpp_side(mode, argv), python_side(mode, argv)
    if any(a in ("-h", "--help") for a in argv) and py[0] == 0 and not py[2]:
        # The tools print hand-written help texts (set_help), so only the
        # exit status, the empty stderr and the usage line are compared.
        assert (cpp[0], cpp[2], cpp[1].splitlines()[0]) == (py[0], py[2], py[1].splitlines()[0])
    else:
        assert cpp == py
