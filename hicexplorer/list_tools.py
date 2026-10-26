"""The `hicexplorer` command: prints the version, the citation and the tools.

The tool list is the commands of the installed HiCExplorer distribution
(package metadata: the executables and console scripts it installed); each
description is the first sentence of the tool's own --help-json description."""
import argparse
import sys
from importlib import metadata

from hicexplorer._version import __version__

CITATION = """
HiCExplorer addresses the common tasks of Hi-C analysis from processing to visualization.
Each tool should be called by its own name as in the following example:

 $ hicPlotMatrix -m hic_matrix.h5 -o plot.pdf

If you find HiCExplorer useful for your research please cite as:

Fidel Ramirez, Vivek Bhardwaj, Jose Villaveces, Laura Arrigoni, Bjoern A Gruening, Kin Chung Lam,
Bianca Habermann, Asifa Akhtar, Thomas Manke.
"High-resolution TADs reveal DNA sequences underlying genome organization in flies".
Nature Communications, Volume 9, Article number: 189 (2018), doi: https://doi.org/10.1038/s41467-017-02525-w

Joachim Wolff, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach,
Thomas Manke, Rolf Backofen, Fidel Ramírez, Björn A Grüning.
"Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and visualization",
Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W11–W16, doi: https://doi.org/10.1093/nar/gky504
"""


def tool_names():
    """The commands of the installed HiCExplorer distribution, sorted: the C++
    executables its RECORD lists in the scripts directory and its console
    scripts."""
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


def _summary(tool):
    from hicexplorer._cpp import tool_spec
    try:
        description = tool_spec(tool).get('description') or ''
    except (OSError, RuntimeError, ValueError):
        return ''
    first = description.strip().split('\n')[0]
    return first.split('. ')[0].rstrip('.')


def tool_list():
    names = tool_names()
    if not names:
        return 'The tool list is read from the package metadata; HiCExplorer is not installed.\n'
    width = max(len(name) for name in names) + 1
    return ''.join('{:<{}}{}\n'.format(name, width, _summary(name)) for name in names)


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(CITATION + '\nThe following is the list of tools:\n\n' + tool_list() +
                     '\nFor more information visit: http://hicexplorer.readthedocs.org\n'))
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    return parser


def process_args(args=None):
    return parse_arguments().parse_args(args)


def main(args=None):
    if args is None and len(sys.argv) == 1:
        args = ["--help"]
    process_args(args)
