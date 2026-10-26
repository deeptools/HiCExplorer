"""Equivalence cases for hicexplorer.parserCommon (apieq.py).

writableFile is C++; getParentArgParse and CustomFormatter only build argparse
objects and stay Python (cpp/scripts/public_api_allowlist.tsv), and are
compared here through what the objects they build do.

Deviation: open(True, 'w') warns that a bool is used as a file descriptor;
writableFile(True) treats True as descriptor 1 without that RuntimeWarning.
"""
import os
import pathlib

from apieq import Case, File

P = 'hicexplorer.parserCommon.'


def written(target, path, *args):
    """The return value, then what the file holds afterwards."""
    result = target(path, *args)
    with open(path, 'rb') as handle:
        return result, handle.read()


def written_bytes(target, path):
    result = target(os.fsencode(path))
    with open(path, 'rb') as handle:
        return result, handle.read()


def descriptor(target, directory):
    """writableFile on an open descriptor closes it."""
    fd = os.open(os.path.join(directory, 'fd.txt'), os.O_WRONLY | os.O_CREAT)
    result = target(fd)
    try:
        os.fstat(fd)
        still_open = True
    except OSError:
        still_open = False
    return result == fd, still_open


def directory_descriptor(target, directory):
    # a fixed descriptor number, which the message and the log show
    opened = os.open(directory, os.O_RDONLY)
    fd = os.dup2(opened, 213)
    os.close(opened)
    try:
        return target(fd)
    finally:
        os.close(fd)


def parent_parser(target, *argv):
    """parse_args of a parser built on the parent parser."""
    import argparse
    parser = argparse.ArgumentParser(parents=[target()], prog='tool')
    try:
        return vars(parser.parse_args(list(argv)))
    except SystemExit as error:
        return 'SystemExit', error.code


def parent_parser_alone(target, args=None):
    parser = target(args)
    actions = [(a.option_strings, a.dest, a.required, a.help) for a in parser._actions]
    return parser.add_help, actions


def formatter_help(target):
    import argparse
    parser = argparse.ArgumentParser(prog='tool', formatter_class=target,
                                     description="first line\n  indented second line\n")
    parser.add_argument('--value', default=5, help='a value')
    parser.add_argument('--flag', action='store_true', help='a flag')
    return parser.format_help(), [cls.__name__ for cls in target.__mro__]


CASES = [
    Case('writableFile/new', P + 'writableFile', args=(File('out.txt'),), driver='written'),
    Case('writableFile/existing_truncated', P + 'writableFile', args=(File('out.txt'),),
         files={'out.txt': b'content'}, driver='written'),
    Case('writableFile/pathlib', P + 'writableFile', args=lambda: (pathlib.Path('.'),)),
    Case('writableFile/bytes', P + 'writableFile', args=(File('out.txt'),),
         driver='written_bytes', files={'out.txt': b'x'}),
    Case('writableFile/missing_directory', P + 'writableFile', args=(File('no/such/out.txt'),)),
    Case('writableFile/directory', P + 'writableFile', args=(File(),)),
    Case('writableFile/embedded_null', P + 'writableFile', args=('a\0b',)),
    Case('writableFile/not_a_path', P + 'writableFile', args=(3.5,)),
    Case('writableFile/none', P + 'writableFile', args=(None,)),
    Case('writableFile/descriptor', P + 'writableFile', args=(File(),), driver='descriptor'),
    Case('writableFile/directory_descriptor', P + 'writableFile', args=(File(),), driver='directory_descriptor'),
    Case('writableFile/bad_descriptor', P + 'writableFile', args=(987654,)),
    Case('writableFile/negative_descriptor', P + 'writableFile', args=(-3,)),

    Case('getParentArgParse/parse', P + 'getParentArgParse', args=('--matrix', 'm.h5'), driver='parent_parser'),
    Case('getParentArgParse/short', P + 'getParentArgParse', args=('-m', 'm.cool'), driver='parent_parser'),
    Case('getParentArgParse/missing', P + 'getParentArgParse', args=(), driver='parent_parser'),
    Case('getParentArgParse/actions', P + 'getParentArgParse', driver='parent_parser_alone'),
    Case('CustomFormatter/help', P + 'CustomFormatter', driver='formatter_help'),
]
