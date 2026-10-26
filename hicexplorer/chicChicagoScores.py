"""chicChicagoScores: runs the C++ executable chicChicagoScores (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of chicChicagoScores, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicChicagoScores')


def main(args=None):
    """Run the C++ implementation of chicChicagoScores with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicChicagoScores')(args)
