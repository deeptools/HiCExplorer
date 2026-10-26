"""chicChicagoPlotViewpoint: runs the C++ executable chicChicagoPlotViewpoint (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of chicChicagoPlotViewpoint, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicChicagoPlotViewpoint')


def main(args=None):
    """Run the C++ implementation of chicChicagoPlotViewpoint with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicChicagoPlotViewpoint')(args)
