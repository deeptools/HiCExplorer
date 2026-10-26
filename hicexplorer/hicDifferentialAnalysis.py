"""hicDifferentialAnalysis: runs the C++ executable hicDifferentialAnalysis (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of hicDifferentialAnalysis, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicDifferentialAnalysis')


def main(args=None):
    """Run the C++ implementation of hicDifferentialAnalysis with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicDifferentialAnalysis')(args)
