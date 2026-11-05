"""hicCreateThresholdFile: a thin layer over the C++ executable hicCreateThresholdFile.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""


def parse_arguments(args=None):
    """The argparse parser of hicCreateThresholdFile, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicCreateThresholdFile')


def main(args=None):
    """Run the C++ implementation of hicCreateThresholdFile with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicCreateThresholdFile')(args)
