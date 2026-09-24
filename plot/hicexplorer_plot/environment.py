"""The drawing environment.

The figures are drawn by matplotlib, and the check below runs before anything
is drawn: the C++ tools run it (python -m hicexplorer_plot --check TOOL)
before they read any input or write any output. A matplotlib older than the
minimum below, or none at all, is refused.

The minimums repeat the versions of plot/pyproject.toml; the ctest
plot_environment_cli compares the two.

The check looks at the module that is actually imported, not only at the
installed distribution metadata, so a different package earlier on
PYTHONPATH is caught too: a module's __version__ when it has one, otherwise
the version of the installed distribution the imported file belongs to.

HICX_PLOT_ALLOW_UNPINNED=1 draws with an older version anyway and says so
on stderr. A missing package is refused either way.
"""

import importlib
import importlib.metadata
import os
import sys

# module name: (distribution name, minimum version)
PINS = {
    "matplotlib": ("matplotlib", "3.8"),
}
TOOL_PINS = {}

OPT_OUT = "HICX_PLOT_ALLOW_UNPINNED"
# Set by the C++ tools once their check passed, so that the drawing process
# does not repeat the opt-out warning.
CHECKED = "HICX_PLOT_CHECKED"
REFUSAL_EXIT = 3


def _numeric(version):
    parts = []
    for piece in version.split("."):
        digits = ""
        for char in piece:
            if not char.isdigit():
                break
            digits += char
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def requirements(tool):
    required = dict(PINS)
    required.update(TOOL_PINS.get(tool, {}))
    return required


def _imported_version(module_name, distribution):
    """(version, None) or (None, the reason there is no usable version)."""
    try:
        module = importlib.import_module(module_name)
    except ImportError as error:
        return None, "is not installed ({})".format(error)
    version = getattr(module, "__version__", None)
    if version is not None:
        return str(version), None
    try:
        dist = importlib.metadata.distribution(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None, "is imported from {} without installed distribution metadata".format(
            getattr(module, "__file__", module_name))
    imported = getattr(module, "__file__", None)
    installed = dist.locate_file(os.path.join(module_name, "__init__.py"))
    if imported is None or os.path.realpath(imported) != os.path.realpath(str(installed)):
        return None, "is imported from {}, which is not the installed {} {}".format(
            imported, distribution, dist.version)
    return str(dist.version), None


def check(tool, stream=None):
    """True when TOOL may draw; writes the refusal or the warning otherwise."""
    stream = sys.stderr if stream is None else stream
    missing, mismatched = [], []
    for module_name, (distribution, required) in sorted(requirements(tool).items()):
        version, problem = _imported_version(module_name, distribution)
        if problem is not None:
            missing.append("{} {} (required: {} {} or newer)".format(distribution, problem, distribution,
                                                            required))
        elif _numeric(version) < _numeric(required):
            mismatched.append("{} {} is installed, but {} {} or newer is required".format(
                distribution, version, distribution, required))
    if not missing and not mismatched:
        return True
    pins = " and ".join("{} {} or newer".format(d, v) for d, v in sorted(requirements(tool).values()))
    if not missing and os.environ.get(OPT_OUT, "") not in ("", "0"):
        if not os.environ.get(CHECKED):
            stream.write(
                "hicexplorer_plot: warning: {tool} draws with {found}; {opt_out} is set, so the "
                "figures may be wrong, which need {pins}.\n".format(
                    tool=tool, found="; ".join(mismatched), opt_out=OPT_OUT, pins=pins))
        return True
    stream.write(
        "hicexplorer_plot: refusing to draw {tool}: {problems}. The interpreter is {exe}. "
        "The drawing needs {pins}: point HICX_PLOT_PYTHON "
        "at a Python that has them{opt_out}. No figure was written.\n".format(
            tool=tool, problems="; ".join(missing + mismatched), exe=sys.executable, pins=pins,
            opt_out="" if missing else ", or set {}=1 to draw anyway with an "
                                       "older version".format(OPT_OUT)))
    return False
