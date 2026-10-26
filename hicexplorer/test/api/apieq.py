"""Function-level equivalence of the public Python API against HiCExplorer 3.7.6.

Every public function of the hicexplorer package keeps the name and signature
it had in 3.7.6 and calls C++ (cpp/python/API_BINDINGS.md). A case calls one
function, or a short driver that uses it, twice: once here on the package
under test, once on the 3.7.6 reference tree (cpp/scripts/reference_tree.py)
in the reference interpreter. Both results are reduced by canon() to plain
data and compared by compare().

What is compared, per call:
  * the return value, including the types of the values (np.float32 is not
    np.float64, a tuple is not a list, an array's dtype and shape count);
  * the arguments after the call, so in-place changes are compared too;
  * whether the return value is one of the arguments (functions that return
    their argument after changing it);
  * an exception's type and message instead of the return value;
  * the records the call logged on the hicexplorer loggers, at every level.

Floats are compared exactly (bit for bit, NaN equal to NaN, 0.0 and -0.0
distinct) unless the case gives rtol/atol, which must then be justified in
the case's `note`. Strings that contain the case's directory are compared
with the directory replaced by <CASE_DIR>, so the two runs may use different
directories.

This module is imported by both interpreters and needs only the standard
library at import time, numpy and scipy when it canonicalises.

Writing cases: see cases_utilities.py. A cases file is a module named
cases_<anything>.py in this directory with a list CASES of Case objects and
the driver functions they name; it must be importable in both interpreters,
so it imports only numpy, scipy, the standard library and this module at top
level, and never hicexplorer.
"""
import copy
import importlib
import importlib.util
import io
import logging
import math
import os
import shutil
import sys

ROOT_LOGGER = 'hicexplorer'
CASE_DIR = '<CASE_DIR>'
# tags of the float leaves of the canonical form, compared with the tolerances
FLOAT = '<float>'
COMPLEX = '<complex>'


class File(object):
    """An argument that is the path of a file of the case (Case.files), or of
    the case's directory itself with name=None."""

    def __init__(self, name=None):
        self.name = name


class Case(object):
    """One equivalence check.

    id         unique within its cases file, used as the pytest id
    target     dotted path of the public name, e.g.
               'hicexplorer.utilities.toString'
    args       positional arguments, or a callable returning them; File
               markers become paths. Each run gets fresh copies.
    kwargs     keyword arguments (dict, or a callable returning one)
    driver     name of a function of the cases file called as
               driver(target_object, *args, **kwargs) instead of the target
               itself, for classes and multi-step uses; its return value is
               compared.
    files      {name: bytes} written into the case directory before each run
    rtol, atol tolerances for floats, 0 means exact
    compare_args  compare the arguments after the call (default True)
    note       why a tolerance or an exception to the defaults is needed
    """

    def __init__(self, id, target, args=(), kwargs=None, driver=None, files=None,
                 rtol=0.0, atol=0.0, compare_args=True, note=''):
        self.id = id
        self.target = target
        self.args = args
        self.kwargs = kwargs or {}
        self.driver = driver
        self.files = files or {}
        self.rtol = rtol
        self.atol = atol
        self.compare_args = compare_args
        self.note = note

    def materialise(self, directory):
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
        for name, content in self.files.items():
            path = os.path.join(directory, name)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'wb') as handle:
                handle.write(content)

    def fresh_arguments(self, directory):
        args = self.args() if callable(self.args) else copy.deepcopy(self.args)
        kwargs = self.kwargs() if callable(self.kwargs) else copy.deepcopy(self.kwargs)

        def resolve(value):
            if isinstance(value, File):
                return directory if value.name is None else os.path.join(directory, value.name)
            return value
        return [resolve(a) for a in args], {k: resolve(v) for k, v in kwargs.items()}


# ---------------------------------------------------------------------------
# loading cases

def load_cases_module(path):
    name = '_apieq_' + os.path.splitext(os.path.basename(path))[0]
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    ids = [case.id for case in module.CASES]
    duplicates = sorted(set(i for i in ids if ids.count(i) > 1))
    if duplicates:
        raise ValueError('{}: duplicate case ids {}'.format(path, duplicates))
    return module


def cases_files(directory):
    return sorted(os.path.join(directory, name) for name in os.listdir(directory)
                  if name.startswith('cases_') and name.endswith('.py'))


def resolve_target(dotted):
    parts = dotted.split('.')
    for split in range(len(parts) - 1, 0, -1):
        try:
            obj = importlib.import_module('.'.join(parts[:split]))
        except ImportError:
            continue
        for attribute in parts[split:]:
            obj = getattr(obj, attribute)
        return obj
    raise ImportError('cannot resolve ' + dotted)


# ---------------------------------------------------------------------------
# running one case

class _Records(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self, level=logging.DEBUG)
        self.records = []

    def emit(self, record):
        self.records.append((record.name, record.levelname, record.getMessage()))


def run_case(module, case, directory):
    """Runs the case in this interpreter and returns its canonical record."""
    case.materialise(directory)
    args, kwargs = case.fresh_arguments(directory)
    logger = logging.getLogger(ROOT_LOGGER)
    handler = _Records()
    old_level = logger.level
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    record = {}
    try:
        target = resolve_target(case.target)
        if case.driver is not None:
            result = getattr(module, case.driver)(target, *args, **kwargs)
        else:
            result = target(*args, **kwargs)
        record['return'] = canon(result, directory)
        record['identity'] = [i for i, a in enumerate(args) if result is a and not _interned(a)]
    except Exception as error:  # the exception is the result
        record['raise'] = (_qualified_type(error), _replace_dir(str(error), directory))
    finally:
        logger.removeHandler(handler)
        logger.setLevel(old_level)
    if case.compare_args and 'raise' not in record:
        record['args'] = canon(args, directory)
        record['kwargs'] = canon(kwargs, directory)
    record['logs'] = [(n, lvl, _replace_dir(m, directory)) for n, lvl, m in handler.records]
    return record


def _interned(value):
    return value is None or isinstance(value, (bool, int, float, str, bytes))


def _qualified_type(error):
    cls = type(error)
    module = cls.__module__
    return cls.__qualname__ if module == 'builtins' else module + '.' + cls.__qualname__


def _replace_dir(text, directory):
    return text.replace(directory, CASE_DIR) if directory else text


def safe_name(case_id):
    return ''.join(c if c.isalnum() or c in '-_.' else '_' for c in case_id)


def compare_records(expected, actual, case):
    """None when the records of the two runs agree, else what differs."""
    if ('raise' in expected) != ('raise' in actual):
        if 'raise' in expected:
            return 'the reference raised {}: {}, the package returned {}'.format(
                expected['raise'][0], expected['raise'][1], _short(actual.get('return')))
        return 'the package raised {}: {}, the reference returned {}'.format(
            actual['raise'][0], actual['raise'][1], _short(expected.get('return')))
    for key in ('raise', 'return', 'identity', 'args', 'kwargs', 'logs'):
        if key in expected or key in actual:
            message = compare(expected.get(key), actual.get(key), case.rtol, case.atol, key)
            if message:
                return message
    return None


# ---------------------------------------------------------------------------
# canonical form

def _type_name(value):
    cls = type(value)
    module = cls.__module__
    return cls.__qualname__ if module == 'builtins' else module + '.' + cls.__qualname__


def canon(value, directory=None):
    """Plain data describing `value`, types included."""
    try:
        import numpy as np
    except ImportError:  # pragma: no cover
        np = None
    if value is None or isinstance(value, bool):
        return ('const', repr(value))
    if np is not None and isinstance(value, np.generic):
        if isinstance(value, (np.floating, np.complexfloating)):
            return ('npscalar', value.dtype.str, (FLOAT, float(value)) if isinstance(value, np.floating)
                    else (COMPLEX, complex(value)))
        if isinstance(value, np.bytes_):
            return ('npscalar', value.dtype.str, bytes(value))
        if isinstance(value, np.str_):
            return ('npscalar', value.dtype.str, _replace_dir(str(value), directory))
        return ('npscalar', value.dtype.str, value.item())
    if isinstance(value, int):
        return (_type_name(value), int(value))
    if isinstance(value, float):
        return (_type_name(value), (FLOAT, float(value)))
    if isinstance(value, str):
        return (_type_name(value), _replace_dir(str(value), directory))
    if isinstance(value, bytes):
        return (_type_name(value), bytes(value))
    if isinstance(value, tuple) and hasattr(value, '_fields'):
        return ('namedtuple', type(value).__module__, type(value).__name__, tuple(value._fields),
                [canon(v, directory) for v in value])
    if isinstance(value, (list, tuple)):
        return (_type_name(value), [canon(v, directory) for v in value])
    if isinstance(value, dict):
        return (_type_name(value), [(canon(k, directory), canon(v, directory)) for k, v in value.items()])
    if isinstance(value, (set, frozenset)):
        return (_type_name(value), sorted((canon(v, directory) for v in value), key=repr))
    if np is not None and isinstance(value, np.ndarray):
        if value.dtype.kind == 'O':
            items = [canon(v, directory) for v in value.ravel().tolist()]
        else:
            items = value.ravel(order='C').tolist()
            if value.dtype.kind in 'fc':
                items = [(FLOAT, v) if isinstance(v, float) else (COMPLEX, v) for v in items]
            elif value.dtype.kind == 'U':
                items = [_replace_dir(v, directory) for v in items]
        return ('ndarray', _type_name(value), value.dtype.str, tuple(value.shape), items)
    try:
        import scipy.sparse as sparse
    except ImportError:  # pragma: no cover
        sparse = None
    if sparse is not None and sparse.issparse(value):
        if value.format in ('csr', 'csc'):
            parts = (canon(value.indptr), canon(value.indices), canon(value.data))
        elif value.format == 'coo':
            parts = (canon(value.row), canon(value.col), canon(value.data))
        else:
            parts = (canon(value.toarray()),)
        return ('sparse', _type_name(value), value.format, tuple(value.shape), value.dtype.str, parts)
    if isinstance(value, BaseException):
        return ('exception', _qualified_type(value), _replace_dir(str(value), directory))
    if isinstance(value, type):
        return ('type', value.__module__, value.__qualname__)
    if isinstance(value, io.IOBase) or type(value).__name__ == 'GzipFile':
        return ('file', _type_name(value), _replace_dir(str(getattr(value, 'name', '')), directory),
                bool(getattr(value, 'closed', False)))
    return ('object', _type_name(value), _replace_dir(repr(value), directory))


# ---------------------------------------------------------------------------
# comparison

def _float_equal(a, b, rtol, atol):
    if isinstance(a, complex) or isinstance(b, complex):
        a, b = complex(a), complex(b)
        return _float_equal(a.real, b.real, rtol, atol) and _float_equal(a.imag, b.imag, rtol, atol)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if rtol == 0 and atol == 0:
        return a == b and math.copysign(1.0, a) == math.copysign(1.0, b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= atol + rtol * abs(b)


def compare(expected, actual, rtol=0.0, atol=0.0, path='result'):
    """None when equal, else a message naming the first difference."""
    if isinstance(expected, tuple) and len(expected) == 2 and expected[0] in (FLOAT, COMPLEX) \
            and isinstance(actual, tuple) and len(actual) == 2 and actual[0] in (FLOAT, COMPLEX):
        if expected[0] != actual[0]:
            return '{}: {} against {}'.format(path, actual, expected)
        if not _float_equal(expected[1], actual[1], rtol, atol):
            return '{}: {!r} against the reference {!r}'.format(path, actual[1], expected[1])
        return None
    if type(expected) is not type(actual):
        return '{}: {!r} against the reference {!r}'.format(path, _short(actual), _short(expected))
    if isinstance(expected, (list, tuple)):
        if len(expected) != len(actual):
            return '{}: length {} against the reference {} ({!r} against {!r})'.format(
                path, len(actual), len(expected), _short(actual), _short(expected))
        for index, (e, a) in enumerate(zip(expected, actual)):
            message = compare(e, a, rtol, atol, '{}[{}]'.format(path, index))
            if message:
                return message
        return None
    if isinstance(expected, dict):
        if sorted(expected) != sorted(actual):
            return '{}: keys {} against the reference {}'.format(path, sorted(actual), sorted(expected))
        for key in expected:
            message = compare(expected[key], actual[key], rtol, atol, '{}.{}'.format(path, key))
            if message:
                return message
        return None
    if isinstance(expected, float):
        return None if _float_equal(expected, actual, rtol, atol) else \
            '{}: {!r} against the reference {!r}'.format(path, actual, expected)
    if expected != actual:
        return '{}: {!r} against the reference {!r}'.format(path, _short(actual), _short(expected))
    return None


def _short(value, limit=300):
    text = repr(value)
    return text if len(text) <= limit else text[:limit] + '...'
