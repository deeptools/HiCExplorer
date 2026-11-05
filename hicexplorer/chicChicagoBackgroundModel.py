"""chicChicagoBackgroundModel: runs the C++ executable chicChicagoBackgroundModel (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of chicChicagoBackgroundModel, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicChicagoBackgroundModel')


def main(args=None):
    """Run the C++ implementation of chicChicagoBackgroundModel with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicChicagoBackgroundModel')(args)
