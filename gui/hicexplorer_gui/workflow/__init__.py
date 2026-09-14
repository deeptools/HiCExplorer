"""Headless workflow engine: a versioned YAML DAG of C++ tool steps.

The engine loads each tool's ``--help-json`` specification, validates the
workflow against it, builds the command lines, runs the steps under a thread
budget with resume and cancellation, and exports the same command lines to a
POSIX shell script or a Snakefile.
"""

from .spec import SpecError, ToolSpec, SpecLoader
from .model import Workflow, Step, WorkflowError, load_workflow
from .validate import validate_workflow
from .cmdline import build_argv
from .runner import Runner, EXIT_OK, EXIT_FAILED, EXIT_USAGE, EXIT_CANCELLED

__all__ = [
    "SpecError", "ToolSpec", "SpecLoader",
    "Workflow", "Step", "WorkflowError", "load_workflow",
    "validate_workflow", "build_argv",
    "Runner", "EXIT_OK", "EXIT_FAILED", "EXIT_USAGE", "EXIT_CANCELLED",
]
