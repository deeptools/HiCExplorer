import argparse


def getParentArgParse(args=None):
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument('--matrix', '-m',
                        help='path of the  Hi-C matrix.',
                        required=True)

    return parser


def writableFile(string):
    try:
        open(string, 'w').close()
    except:
        msg = "{} file can be opened for writting".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string
