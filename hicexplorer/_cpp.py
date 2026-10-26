"""Locates and runs the C++ tool executables that implement HiCExplorer v4.

Every tool module hicexplorer/<tool>.py is a thin layer over its executable:
main(args) runs the executable, and parse_arguments() returns an
argparse.ArgumentParser built from the executable's --help-json
specification by build_parser(), so the Python parser always has the options,
defaults and choices of the C++ tool."""
import argparse
import json
import os
import subprocess
import sys
import sysconfig

_HERE = os.path.dirname(os.path.abspath(__file__))


def _candidate_dirs():
    env = os.environ.get('HICX_BIN_DIR')
    if env:
        yield env
    # a source checkout: the executables of cpp/build
    yield os.path.join(os.path.dirname(_HERE), 'cpp', 'build', 'tools')
    # an installed package: the wheel puts the executables into the scripts
    # directory of the environment (cpp/cmake/HicxInstall.cmake)
    yield sysconfig.get_path('scripts')
    yield from os.environ.get('PATH', '').split(os.pathsep)


def _is_native(path):
    """True for an executable that is not a script: a console script of an
    older Python HiCExplorer of the same name starts with '#!'."""
    try:
        with open(path, 'rb') as handle:
            return handle.read(2) != b'#!'
    except OSError:
        return False


def find_executable(tool):
    exe = tool + ('.exe' if os.name == 'nt' else '')
    for directory in _candidate_dirs():
        if not directory:
            continue
        path = os.path.join(directory, exe)
        if os.path.isfile(path) and os.access(path, os.X_OK) and _is_native(path):
            return path
    raise FileNotFoundError(
        "The C++ executable '{}' was not found in HICX_BIN_DIR, the environment's scripts "
        "directory or on PATH.".format(tool))


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
    """Return the console-script function of `tool`: it runs the executable
    and exits with its status when that is not 0."""
    def main(args=None):
        status = run(tool, args)
        if status != 0:
            sys.exit(status)
    main.__name__ = 'main'
    main.__doc__ = 'Run the C++ implementation of {}.'.format(tool)
    return main


# --------------------------------------------------------------------------
# argparse parsers from --help-json

_SPECS = {}
_TYPES = {'str': str, 'int': int, 'float': float}


def tool_spec(tool):
    """The --help-json specification of the executable of `tool` (a dict)."""
    if tool not in _SPECS:
        result = subprocess.run([find_executable(tool), '--help-json'], env=_environment(tool),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            raise RuntimeError('{} --help-json exited {}: {}'.format(
                tool, result.returncode, result.stderr.decode(errors='replace').strip()))
        _SPECS[tool] = json.loads(result.stdout.decode())
    return _SPECS[tool]


class _HelpJsonAction(argparse.Action):
    def __init__(self, option_strings, dest, spec=None, **kwargs):
        self._spec = spec
        super().__init__(option_strings, dest, nargs=0, default=argparse.SUPPRESS, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print(json.dumps(self._spec))
        parser.exit()


def _type_of(arg):
    kind = arg.get('type')
    if kind is None:
        return None
    if kind == 'FileType':
        role = (arg.get('file') or {}).get('role')
        return argparse.FileType('w' if role == 'output' else 'r')
    # Types the C++ validates itself (genomicRegion, writableFile) are strings here.
    return _TYPES.get(kind, str)


def _add_argument(container, arg, root_spec):
    action = arg['action']
    flags = list(arg['flags'])
    kwargs = {'help': arg.get('help')}
    if action == 'help':
        return container.add_argument(*flags, action='help', help=kwargs['help'])
    if action == 'version':
        return container.add_argument(*flags, action='version', help=kwargs['help'],
                                      version='%(prog)s ' + root_spec.get('version', ''))
    if action == 'help_json':
        return container.add_argument(*flags, action=_HelpJsonAction, spec=root_spec,
                                      help=argparse.SUPPRESS)
    if arg.get('positional') or not flags:
        names = [arg['dest']]
    else:
        names = flags
        kwargs['dest'] = arg['dest']
        kwargs['required'] = bool(arg.get('required'))
    kwargs['default'] = arg.get('default')
    if action in ('store_true', 'store_false'):
        kwargs['action'] = action
    else:
        kwargs['action'] = 'append' if action == 'append' else 'store'
        kwargs['type'] = _type_of(arg)
        if arg.get('nargs') is not None:
            kwargs['nargs'] = arg['nargs']
        if arg.get('choices') is not None:
            kwargs['choices'] = arg['choices']
        metavar = arg.get('metavar')
        if metavar is not None:
            kwargs['metavar'] = tuple(metavar) if isinstance(metavar, list) else metavar
    return container.add_argument(*names, **kwargs)


def _populate(parser, spec, root_spec):
    mutex_of = {}
    for group in spec.get('mutually_exclusive') or []:
        for dest in group['dests']:
            mutex_of[dest] = group
    created = {}
    for group_spec in spec.get('groups', []):
        group = parser.add_argument_group(group_spec.get('title'))
        for arg in group_spec.get('arguments', []):
            mutex = mutex_of.get(arg.get('dest'))
            container = group
            if mutex is not None:
                key = id(mutex)
                if key not in created:
                    created[key] = group.add_mutually_exclusive_group(
                        required=bool(mutex.get('required')))
                container = created[key]
            _add_argument(container, arg, root_spec)
    subcommands = spec.get('subcommands')
    if subcommands:
        subparsers = parser.add_subparsers(dest=subcommands.get('dest'))
        subparsers.required = bool(subcommands.get('required'))
        for command in subcommands['commands']:
            sub = subparsers.add_parser(command['name'], add_help=False,
                                        description=command.get('description'),
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
            _populate(sub, command, root_spec)


def build_parser(tool, prog=None):
    """An argparse.ArgumentParser equivalent to the C++ parser of `tool`.

    Built from the executable's --help-json: the options, destinations,
    types, nargs, choices, defaults, required flags, mutually exclusive
    groups and subcommands are the executable's. Parsing with it only checks
    a command line; main() of the tool module runs the tool."""
    spec = tool_spec(tool)
    parser = argparse.ArgumentParser(prog=prog or tool, add_help=False,
                                     description=spec.get('description'),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    _populate(parser, spec, spec)
    return parser
