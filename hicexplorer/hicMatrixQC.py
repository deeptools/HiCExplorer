"""hicMatrixQC: runs the C++ executable hicMatrixQC (no Python implementation exists)."""


def parse_arguments(args=None):
    """The argparse parser of hicMatrixQC, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicMatrixQC')


def main(args=None):
    """Run the C++ implementation of hicMatrixQC with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicMatrixQC')(args)
