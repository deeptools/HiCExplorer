"""Locates and runs the C++ tool executables that implement HiCExplorer v4."""
import os
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def _candidate_dirs():
    yield os.path.join(_HERE, '_bin')
    env = os.environ.get('HICX_BIN_DIR')
    if env:
        yield env
    yield os.path.join(os.path.dirname(_HERE), 'cpp', 'build', 'tools')


def find_executable(tool):
    exe = tool + ('.exe' if os.name == 'nt' else '')
    for directory in _candidate_dirs():
        path = os.path.join(directory, exe)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    found = shutil.which(tool)
    if found:
        return found
    raise FileNotFoundError(
        "The C++ executable '{}' was not found in the package, in HICX_BIN_DIR or on PATH.".format(tool))


def _environment(tool):
    env = dict(os.environ)
    if 'HICX_PLOT_PYTHON' not in env:
        env['HICX_PLOT_PYTHON'] = sys.executable
    return env


def run(tool, args=None):
    """Run the C++ executable of `tool` with `args` (default sys.argv[1:]).

    Returns the exit status. Output goes straight to the caller's stdout and stderr."""
    if args is None:
        args = sys.argv[1:]
    elif isinstance(args, str):
        args = args.split()
    args = [str(a) for a in args]
    command = [find_executable(tool)] + args
    env = _environment(tool)
    # A caller that redirected sys.stdout or sys.stderr (capture in tests,
    # notebooks, wrappers) only sees output written through them.
    redirect_out = sys.stdout is not sys.__stdout__
    redirect_err = sys.stderr is not sys.__stderr__
    if not (redirect_out or redirect_err):
        return subprocess.call(command, env=env)
    result = subprocess.run(command, env=env,
                            stdout=subprocess.PIPE if redirect_out else None,
                            stderr=subprocess.PIPE if redirect_err else None)
    if redirect_out:
        sys.stdout.write(result.stdout.decode(errors='replace'))
    if redirect_err:
        sys.stderr.write(result.stderr.decode(errors='replace'))
    return result.returncode


def entry_point(tool):
    """Return the console-script function of `tool`."""
    def main(args=None):
        status = run(tool, args)
        if status != 0:
            sys.exit(status)
    main.__name__ = 'main'
    main.__doc__ = 'Run the C++ implementation of {}.'.format(tool)
    return main
