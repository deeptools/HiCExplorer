# The public Python API on C++: how to convert a module

Owner directive: every public module, class and function of the 3.7.6
package stays importable under the same name with the same signature, and its
body calls C++. No Python, numpy or scipy algorithm code stays behind a public
name; matplotlib drawing may stay Python. Nothing public is deleted.

The work list is `cpp/scripts/public_api_inventory.tsv`. A converted module
is `hicexplorer.utilities` (28 functions) together with `hicexplorer.readBed`
(the class `ReadBed`); copy what they do.

## Layout

```
cpp/core/include/hicx/<topic>.hpp    the computation, no Python headers
cpp/core/src/<topic>.cpp             (listed in cpp/core/CMakeLists.txt)
cpp/tests/test_<topic>.cpp           doctest unit tests (cpp/tests/CMakeLists.txt)
cpp/python/src/api/modules.cmake     the registry: one line per bound module
cpp/python/src/api/<entry>.cpp       void hicx::python::bind_<entry>(py::module_&)
cpp/python/src/api/api.hpp           what a binding includes
cpp/python/src/api/support/          shared helpers, see below
hicexplorer/<module>.py              thin wrappers, same names and signatures
hicexplorer/test/api/cases_<module>.py   equivalence cases against 3.7.6
```

The extension is `hicexplorer._hicx`, one pybind11 module with one submodule
per registry entry: `from hicexplorer._hicx import utilities as _impl`. The
entry is the module path below `hicexplorer` with `.` replaced by `_`
(`hicexplorer.lib.cnb` is `lib_cnb`). The wheel installs it into the
`hicexplorer` package (`cpp/cmake/HicxInstall.cmake`); a build tree has it in
`<build>/python/hicexplorer`, which a source checkout's `hicexplorer` finds
by itself (`HICX_PYTHON_EXT_DIR`, the build tree of `HICX_BIN_DIR`, or
`cpp/build`). Configure with `-DHICX_BUILD_PYTHON=ON`.

## Recipe

1. **Core.** Put the computation into `cpp/core` as plain C++ over vectors and
   strings: reuse an existing kernel where it exists (the tools already port
   many module functions; grep the inventory's last column), otherwise add one.
   Reproduce 3.7.6 including its defects, dtype rules and error conditions;
   write the Python exception the function raised as
   `throw hicx::PythonError("IndexError", message)` (`hicx/python_error.hpp`,
   a dotted name like `"argparse.ArgumentTypeError"` works too). Add doctest
   cases, and where a tool kernel computes the same function, a test that
   the two agree on a real matrix.
2. **Binding.** Add `cpp/python/src/api/<entry>.cpp` with
   `void bind_<entry>(py::module_& m)` and add `<entry>` to `modules.cmake`
   (keep the list sorted; that line is the only shared edit). Define each
   function with the 3.7.6 argument names and defaults (`py::arg("pSubmatrix")`,
   `py::arg("pThreads") = py::none()`), take `py::handle`s, convert, call the
   core, build the 3.7.6 result type. For a class, keep the Python class and
   bind its methods as free functions taking `self`
   (`ReadBed_get_bed_interval(self, bed_line)`); call other methods through
   `self.attr(...)` so subclasses keep working.
3. **Wrapper.** Rewrite `hicexplorer/<module>.py`: module docstring, the
   imports whose side effects 3.7.6 had (warnings filters, `mplt_use('Agg')`),
   `log = logging.getLogger(__name__)`, `from hicexplorer._hicx import <entry> as _impl`,
   then every public function with its original signature and docstring and
   the single statement `return _impl.<name>(<arguments in order>)`. Doctests
   in docstrings keep running (`py.test hicexplorer --doctest-modules`).
4. **Cases.** Add `hicexplorer/test/api/cases_<module>.py` (see
   `cases_utilities.py`): `CASES = [Case(id, 'hicexplorer.<module>.<name>', args=...)]`.
   Cover real data of `hicexplorer/test/test_data`, every dtype the function
   is called with, in-place results, every exception the 3.7.6 code raises,
   and the defects you reproduced. `args` may be a callable for fresh
   mutable inputs; `File('x')` is a path in the case directory, created from
   `files={'x': b'...'}`; `driver='name'` calls `name(target, *args)` from the
   cases file for classes and multi-step uses. A cases file is imported by
   both interpreters: top-level imports only numpy, scipy, the standard
   library and `apieq`, never `hicexplorer`.
5. **Check.**
   `python cpp/scripts/public_api_inventory.py . --output cpp/scripts/public_api_inventory.tsv`
   regenerates the inventory, and `--check` lists what is left; your names
   must now say `C++-backed (hicexplorer._hicx)`. The rule is mechanical:
   every statement of the body is a call (or `return` of a call) of an
   `_impl` function with names and constants as arguments. A matplotlib
   drawing function that stays Python goes into
   `cpp/scripts/public_api_allowlist.tsv` with its reason; a function that
   also computes data does not, its computation moves to C++.

## The shared helpers

`support/convert.hpp`:
`utf8(str)`/`to_str(utf8)` (surrogatepass, so every str survives),
`str_of`, `repr_of`, `type_name`;
`dtype_of`, `dispatch_arithmetic(array, name, []<class T>{...})` for
int32/int64/float32/float64 and `dispatch_numeric` for every bool, integer,
real and complex dtype (anything else raises TypeError naming the dtype);
`to_vector<T>` (C order, any strides), `write_back` (in place, refuses
read-only arrays as numpy does), `from_vector`, `take_elements`,
`numpy_scalar<T>`; `sparse(object, name, "rco")` reads a scipy
csr/csc/coo matrix (`entries()` lists the stored entries in `tocoo()` order),
`set_data`, `apply_eliminate_zeros`; `to_csr_matrix` and `csr_matrix_object`
between scipy and the tools' `hicx::CsrMatrix`; `cut_intervals`,
`cut_intervals_object`, `matrix_data` and `assign_matrix_data` between a
hiCMatrix object and `hicx::MatrixData` (the tools' matrix with its bins,
NaN bins and correction factors), keeping the field types of
`cut_intervals`; `to_int64`, `to_double`, `parse_int` (Python's own `int()`).

`support/errors.hpp`: `raise_error(type, message)`, `raise_os_error(errno, filename)`
(the OSError subclass open() raises), `throw_python_error()` after a failed
C API call. `support/runtime.hpp`: `log(logger, level, message)` for the
3.7.6 log records (the cases compare them), `warn`, `read_file` and
`read_text_lines` (open()'s errors, the locale's codec, universal newlines,
TextIOWrapper's UnicodeDecodeError positions).

Calling the Python runtime for conversions and object protocol is fine
(`int()` parsing, `repr`, `dict(...)`, constructing the namedtuple or scipy
matrix the function returns, `copy.deepcopy` for `pInplace=False`, `open()`
when the function returns a file object). Computing with numpy, scipy or
Python code is not.

## Equivalence tests

`hicexplorer/test/api/` runs every case on the package and on the 3.7.6
reference tree (`cpp/scripts/reference_tree.py`) in the reference interpreter
(`~/miniconda3/envs/__hicexplorer@3.7.6/bin/python3`, or
`HICX_REFERENCE_PYTHON`, `HICX_REFERENCE_ROOT`) and compares return value
with types and dtypes, arguments after the call, return-is-argument identity,
exception type and message, and log records. Floats are exact unless a case
sets `rtol`/`atol` with a `note` saying why. Reference results are cached in
`~/.cache/hicexplorer-equiv/api` keyed by the reference tree and the case
file (`HICX_API_EQUIV_CACHE=refresh|off`). Without a reference interpreter
the tests skip.

```
HICX_BIN_DIR=$PWD/cpp/build/tools PYTHONPATH=$PWD python -m pytest hicexplorer/test/api -q
```

## Deviations

Record anything that cannot be made identical in the cases file's docstring
and in the binding, and add a case pinning the 3.7.6 behaviour you do keep.
For utilities and readBed:

* Error paths of worker processes: `expected_interactions` with `pThreads`
  returns `'Fail: <message>'` plus a traceback; the message is identical, the
  traceback text (Python frames of the worker) is not. 3.7.6 also sleeps one
  second per call there; the C++ does not.
* `ReadBed.get_no_comment_line` recurses in 3.7.6, so about a thousand
  consecutive comment lines raise RecursionError; the C++ reads them.
* Types 3.7.6 accepted only by accident raise TypeError naming HiCExplorer 4:
  float coordinates in `enlarge_bins`, float chromosome sizes in
  `getUserRegion`, integers beyond int64, float16, longdouble and complex
  data in the statistics, non-native byte order, numpy.ndarray subclasses in
  `toString` (returned as plain ndarrays), COO input to the functions that
  call `eliminate_zeros`.
* `check_cooler` on a `file::/group` URI whose group does not exist returns
  False; 3.7.6 raised h5py's KeyError. A cooler missing a required group
  does not warn.
* A UTF-8 decoding error of a text file is located by replaying
  TextIOWrapper's 8192 byte chunks, which is exact for regular files; for a
  locale that is not UTF-8 the file is decoded whole, so the reported
  position may differ.
