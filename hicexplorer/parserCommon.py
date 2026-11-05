"""hicexplorer.parserCommon: argparse helpers shared by the tools.

writableFile calls hicexplorer._hicx.parserCommon (C++,
cpp/python/src/api/parserCommon.cpp). getParentArgParse and CustomFormatter
only construct argparse objects and compute nothing; they stay as they were
in HiCExplorer 3.7.6 (cpp/scripts/public_api_allowlist.tsv).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import parserCommon as _impl


def getParentArgParse(args=None):
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument('--matrix', '-m',
                        help='path of the  Hi-C matrix.',
                        required=True)

    return parser


def writableFile(string):
    return _impl.writableFile(string)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """
    This class allows to use defaultsHelpFormatter and RawDescription at the same time.

    Usage:

        parser = argparse.ArgumentParser(
            formatter_class=CustomFormatter,
            conflict_handler='resolve',
            usage="tex",
            description="text"

    """
    pass
