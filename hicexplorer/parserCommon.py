import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

import logging
log = logging.getLogger(__name__)


def getParentArgParse(args=None):
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument('--matrix', '-m',
                        help='path of the  Hi-C matrix.',
                        required=True)

    return parser


def writableFile(string):
    try:
        open(string, 'w').close()
    except IOError:
        msg = "{} file can be opened for writting".format(string)
        log.debug(msg)
        raise argparse.ArgumentTypeError(msg)
    return string


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
