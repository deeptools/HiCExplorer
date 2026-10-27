"""The `hicexplorer` command: prints the version, the citation and the tools.

The command is the C++ executable `hicexplorer` (cpp/tools/hicexplorer.cpp),
whose tool list and descriptions are generated at build time from the tools'
own specifications; main() runs it. The other names of this module are kept
for API compatibility and run in C++ (hicexplorer._hicx.list_tools):
tool_list() reads the installed distribution's commands (package metadata:
the executables and console scripts it installed), each described by the first
sentence of the tool's own --help-json description."""
from hicexplorer._version import __version__
from hicexplorer._hicx import list_tools as _impl

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
    return _impl.tool_names()


def tool_list():
    return _impl.tool_list()


def parse_arguments():
    return _impl.parse_arguments(CITATION, __version__)


def process_args(args=None):
    return _impl.process_args(args)


def main(args=None):
    """Run the C++ executable `hicexplorer` with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicexplorer')(args)
