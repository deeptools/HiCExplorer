"""hicexplorer.list_tools: tool_names, tool_list, parse_arguments and
process_args run in C++ (hicexplorer._hicx.list_tools); they must return what
the Python implementation of HiCExplorer 4 returned before the move. That
implementation (3.7.6 has only parse_arguments, process_args and main, with a
fixed list of tools) is kept here as the oracle and run on the same stubbed
package metadata and tool specifications."""
import argparse
import importlib.metadata as metadata

import pytest

import hicexplorer.list_tools as list_tools
from hicexplorer._version import __version__


# ---------------------------------------------------------------------------
# the oracle: the pure Python functions of hicexplorer/list_tools.py before the move

def oracle_tool_names():
    try:
        distribution = metadata.distribution('HiCExplorer')
    except metadata.PackageNotFoundError:
        return []
    names = {entry.name for entry in distribution.entry_points
             if entry.group == 'console_scripts'}
    for path in distribution.files or []:
        if len(path.parts) > 1 and path.parts[-2] in ('bin', 'Scripts'):
            names.add(path.stem if path.suffix == '.exe' else path.name)
    names.discard('hicexplorer')
    return sorted(name for name in names if name.startswith(('hic', 'chic')))


def oracle_summary(tool):
    from hicexplorer._cpp import tool_spec
    try:
        description = tool_spec(tool).get('description') or ''
    except (OSError, RuntimeError, ValueError):
        return ''
    first = description.strip().split('\n')[0]
    return first.split('. ')[0].rstrip('.')


def oracle_tool_list():
    names = oracle_tool_names()
    if not names:
        return 'The tool list is read from the package metadata; HiCExplorer is not installed.\n'
    width = max(len(name) for name in names) + 1
    return ''.join('{:<{}}{}\n'.format(name, width, oracle_summary(name)) for name in names)


def oracle_parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(list_tools.CITATION + '\nThe following is the list of tools:\n\n' +
                     oracle_tool_list() +
                     '\nFor more information visit: http://hicexplorer.readthedocs.org\n'))
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    return parser


# ---------------------------------------------------------------------------
# stubs

class Entry:
    def __init__(self, name, group):
        self.name = name
        self.group = group


class Distribution:
    def __init__(self, scripts=(), files=None, others=()):
        self.entry_points = ([Entry(n, 'console_scripts') for n in scripts] +
                             [Entry(n, 'gui_scripts') for n in others])
        self.files = None if files is None else [metadata.PackagePath(f) for f in files]


DISTRIBUTIONS = {
    'empty': Distribution(),
    'no_files': Distribution(scripts=['hicPlotMatrix'], files=None),
    'mixed': Distribution(
        scripts=['hicHyperoptDetectLoops', 'chicPlotViewpoint', 'pytest', 'hicexplorer'],
        files=['../../bin/hicInfo', '../../bin/hicexplorer', '../../bin/other', 'bin/hicQC',
               'x/Scripts/hicFindTADs.exe', 'x/Scripts/readme.txt', 'x/Scripts/chicQualityControl.exe',
               'hicexplorer/hicInfo.py', 'bin', 'bin/chicViewpoint', 'hicexplorer/bin/hicNormalize'],
        others=['hicGui']),
    'duplicates': Distribution(scripts=['hicInfo', 'hicInfo'], files=['bin/hicInfo']),
    'unicode': Distribution(scripts=['hicNäme', 'hic'], files=['bin/hicZ']),
}


@pytest.fixture(params=sorted(DISTRIBUTIONS) + ['missing'])
def installed(request, monkeypatch):
    def distribution(name):
        assert name == 'HiCExplorer'
        if request.param == 'missing':
            raise metadata.PackageNotFoundError(name)
        return DISTRIBUTIONS[request.param]
    monkeypatch.setattr(metadata, 'distribution', distribution)
    return request.param


DESCRIPTIONS = [
    'Plots a matrix. It does more.', '\n  First line.\nsecond line', '', None, 'No period',
    'Ends with dots...', 'A. B. C', 'Version 1.2. text', 'Umlaut äö. text',
]


@pytest.fixture
def specs(monkeypatch):
    import hicexplorer._cpp as cpp
    calls = []

    def tool_spec(tool):
        calls.append(tool)
        index = sum(map(ord, tool))
        kind = index % (len(DESCRIPTIONS) + 3)
        if kind == len(DESCRIPTIONS):
            raise OSError('no executable')
        if kind == len(DESCRIPTIONS) + 1:
            raise RuntimeError('exited 1')
        if kind == len(DESCRIPTIONS) + 2:
            raise ValueError('bad json')
        return {'description': DESCRIPTIONS[kind]}
    monkeypatch.setattr(cpp, 'tool_spec', tool_spec)
    return calls


# ---------------------------------------------------------------------------
# tests

def test_tool_names(installed):
    assert list_tools.tool_names() == oracle_tool_names()
    assert isinstance(list_tools.tool_names(), list)


def test_tool_names_not_installed():
    # the checkout is not an installed distribution
    assert list_tools.tool_names() == oracle_tool_names()


def test_tool_list(installed, specs):
    assert list_tools.tool_list() == oracle_tool_list()


def test_tool_list_uses_module_tool_names(monkeypatch, specs):
    monkeypatch.setattr(list_tools, 'tool_names', lambda: ['hicA', 'hicLonger'])
    assert list_tools.tool_list() == oracle_summary_list(['hicA', 'hicLonger'])


def oracle_summary_list(names):
    width = max(len(name) for name in names) + 1
    return ''.join('{:<{}}{}\n'.format(name, width, oracle_summary(name)) for name in names)


def test_parse_arguments(installed, specs):
    mine = list_tools.parse_arguments()
    theirs = oracle_parse_arguments()
    assert type(mine) is type(theirs)
    assert mine.format_help() == theirs.format_help()
    assert mine.format_usage() == theirs.format_usage()
    assert mine.formatter_class is argparse.RawDescriptionHelpFormatter


def test_version(capsys):
    with pytest.raises(SystemExit) as stop:
        list_tools.process_args(['--version'])
    assert stop.value.code == 0
    mine = capsys.readouterr().out
    with pytest.raises(SystemExit) as oracle_stop:
        oracle_parse_arguments().parse_args(['--version'])
    assert oracle_stop.value.code == 0
    assert capsys.readouterr().out == mine
    assert mine.strip().endswith(__version__)


def test_process_args(specs):
    assert vars(list_tools.process_args([])) == {}


def test_process_args_unknown(capsys):
    with pytest.raises(SystemExit) as stop:
        list_tools.process_args(['--nope'])
    assert stop.value.code == 2
    assert 'unrecognized arguments: --nope' in capsys.readouterr().err


def test_process_args_uses_module_parse_arguments(monkeypatch):
    parser = argparse.ArgumentParser()
    parser.add_argument('--x', default=5)
    monkeypatch.setattr(list_tools, 'parse_arguments', lambda: parser)
    assert list_tools.process_args(['--x', '7']).x == '7'
    assert list_tools.process_args([]).x == 5


def test_signatures():
    import inspect
    assert str(inspect.signature(list_tools.tool_names)) == '()'
    assert str(inspect.signature(list_tools.tool_list)) == '()'
    assert str(inspect.signature(list_tools.parse_arguments)) == '()'
    assert str(inspect.signature(list_tools.process_args)) == '(args=None)'
    assert str(inspect.signature(list_tools.main)) == '(args=None)'
