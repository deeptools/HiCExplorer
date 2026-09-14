"""The drawing layer of the HiCExplorer v4 plotting tools (cpp/PLAN.md tier 7,
option (a)).

Each plotting tool is a C++ binary that parses its command line with
hicx::cli (so --help-json and the argparse equality of PLAN 10.1 hold), writes
the tool's data outputs, and writes the data its figure is drawn from as a
JSON file. It then replaces itself with

    python -m hicexplorer_plot TOOL DATA.json --remove-data

which draws the figure with the same matplotlib calls as the Python tool in
hicexplorer/, the oracle, which stays unchanged. The C++ memory is gone by the
time matplotlib starts. One module per tool, named after it, exposes
draw(data).

The interpreter the binaries start is HICX_PLOT_PYTHON, else python3. The
console scripts in hicexplorer_plot.tools run the C++ binary with that
variable set to the interpreter they run in.
"""

MATPLOTLIB_VERSION = "3.8.4"
PYGENOMETRACKS_VERSION = "3.9"
