"""Helpers of the GUI tests (PLAN 10.4 and 10.5).

Environment:
    HICX_CPP_BIN            C++ tool directory (tests needing tools skip without it)
    HICX_REFERENCE_PYTHON   interpreter with the Python HiCExplorer dependencies
                            (argparse namespaces, equiv.py comparisons)
    HICX_PYTHON_MODULE_DIR  directory holding the hicx_matrix module
"""

import json
import os
import subprocess
import sys

import pytest

GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(GUI_DIR)
DATA = os.path.join(REPO, "hicexplorer", "test", "test_data")
EQUIV = os.path.join(REPO, "cpp", "scripts", "equiv.py")
CPP_BIN = os.environ.get("HICX_CPP_BIN")
REFERENCE_PYTHON = os.environ.get("HICX_REFERENCE_PYTHON")

needs_tools = pytest.mark.skipif(not CPP_BIN, reason="HICX_CPP_BIN is not set")
needs_reference = pytest.mark.skipif(not REFERENCE_PYTHON, reason="HICX_REFERENCE_PYTHON is not set")


def python_namespaces(requests, workdir):
    """[result] of argparse_namespace.py for [{"tool", "args"}]."""
    path = os.path.join(str(workdir), "namespace_requests.json")
    with open(path, "w") as handle:
        json.dump(requests, handle)
    env = dict(os.environ, PYTHONPATH=REPO + os.pathsep + os.environ.get("PYTHONPATH", ""))
    proc = subprocess.run([REFERENCE_PYTHON, os.path.join(GUI_DIR, "tests", "argparse_namespace.py"), path],
                          cwd=str(workdir), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError("argparse_namespace.py failed:\n" + proc.stderr[-4000:])
    return json.loads(proc.stdout)


def equiv_compare(fmt, a, b, klass="E0"):
    proc = subprocess.run([REFERENCE_PYTHON, EQUIV, "compare", "--format", fmt, "--class", klass, a, b],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    return proc.returncode, proc.stdout


def import_hicx_matrix():
    module_dir = os.environ.get("HICX_PYTHON_MODULE_DIR")
    if module_dir and module_dir not in sys.path:
        sys.path.insert(0, module_dir)
    return pytest.importorskip("hicx_matrix")
