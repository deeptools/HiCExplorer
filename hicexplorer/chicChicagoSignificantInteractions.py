"""chicChicagoSignificantInteractions: runs the C++ executable chicChicagoSignificantInteractions (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of chicChicagoSignificantInteractions, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicChicagoSignificantInteractions')


def main(args=None):
    """Run the C++ implementation of chicChicagoSignificantInteractions with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicChicagoSignificantInteractions')(args)
