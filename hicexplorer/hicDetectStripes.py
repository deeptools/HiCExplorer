"""hicDetectStripes: runs the C++ executable hicDetectStripes (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of hicDetectStripes, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicDetectStripes')


def main(args=None):
    """Run the C++ implementation of hicDetectStripes with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicDetectStripes')(args)
