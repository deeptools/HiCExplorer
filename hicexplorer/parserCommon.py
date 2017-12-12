import argparse
import logging

logging.basicConfig()
log = logging.getLogger("parserCommon")
log.setLevel(logging.WARN)


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
        log.exception(msg)
        raise argparse.ArgumentTypeError(msg)
    return string
