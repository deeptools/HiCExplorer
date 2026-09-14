"""Generated tool forms against the specifications and the Python parsers (PLAN 10.4).

For every ported tool (and every subcommand):
- the form renders one field per settable argument of the specification;
- an empty form reports every required argument, a value outside the choices
  and two arguments of a mutually exclusive group;
- a filled form produces a command line, and the Python tool's own
  parse_arguments() parses it into exactly the namespace the form values mean.
"""

import os
import re

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from gui_support import CPP_BIN, REFERENCE_PYTHON, REPO, needs_reference, needs_tools, python_namespaces  # noqa: E402
from hicexplorer_gui.catalog import tool_entries  # noqa: E402
from hicexplorer_gui.forms import ChoiceField, FileField, FlagField, ToolForm, is_list_arg  # noqa: E402
from hicexplorer_gui.workflow.spec import NON_SETTABLE_ACTIONS  # noqa: E402

ENTRIES = tool_entries(CPP_BIN)[0] if CPP_BIN else []
CASES = []
for _entry in ENTRIES:
    if _entry.available:
        for _sub in (list(_entry.spec.commands) if _entry.spec.has_subcommands else [None]):
            CASES.append((_entry.name, _sub))

pytestmark = needs_tools


def _spec(tool):
    return next(e.spec for e in ENTRIES if e.name == tool)


def ported_tools():
    """The tools this revision ports: cpp/tools/<tool>.cpp, plus the
    executables CMake builds from another tool's source under a second name
    (hicQC from hicPrepareQCreport.cpp), as cpp/scripts/tool_specs.py lists
    them."""
    tools_dir = os.path.join(REPO, "cpp", "tools")
    names = {os.path.splitext(n)[0] for n in os.listdir(tools_dir)
             if n.endswith(".cpp") and n.startswith(("hic", "chic"))}
    with open(os.path.join(tools_dir, "CMakeLists.txt")) as handle:
        names.update(re.findall(r"add_executable\(((?:hic|chic)\w+)", handle.read()))
    return sorted(names)


def test_ported_tools_available_and_the_rest_explained():
    available = [e for e in ENTRIES if e.available]
    missing = [e for e in ENTRIES if not e.available]
    assert sorted(e.name for e in available) == ported_tools()
    assert missing and all(e.reason.startswith("not ported to C++ yet, PLAN tier") for e in missing), \
        [(e.name, e.reason) for e in missing]


def settable(spec, sub):
    return [a for a in spec.arguments(sub) if (a.get("action") or "store") not in NON_SETTABLE_ACTIONS]


@pytest.mark.parametrize("tool,sub", CASES, ids=["{}{}".format(t, " " + s if s else "") for t, s in CASES])
def test_form_renders_and_validates(qtbot, tmp_path, tool, sub):
    spec = _spec(tool)
    form = ToolForm(spec)
    qtbot.addWidget(form)
    form.set_subcommand(sub)
    args = settable(spec, sub)
    assert set(form.fields) == {a["dest"] for a in args}
    for arg in args:
        field = form.fields[arg["dest"]]
        action = arg.get("action") or "store"
        if action in ("store_true", "store_false"):
            assert isinstance(field, FlagField)
        elif arg.get("file"):
            assert isinstance(field, FileField)
            assert field.role == arg["file"]["role"]
        elif arg.get("choices") and not is_list_arg(arg):
            assert isinstance(field, ChoiceField)
            offered = [field.combo.itemData(i) for i in range(1, field.combo.count())]
            assert offered == [c for c in arg["choices"] if c is not None]

    errors = form.validate()
    by_dest = {a["dest"]: a for a in args}
    for dest, arg in by_dest.items():
        if arg.get("required"):
            assert errors.get(dest) == "required", (dest, errors)
            assert form.fields[dest].error.isVisibleTo(form)
    for dest in errors:
        assert dest in by_dest and (by_dest[dest].get("required") or any(
            g.get("required") and dest in g["dests"] for g in spec.mutually_exclusive(sub))), (dest, errors)

    for arg in args:
        choices = [c for c in (arg.get("choices") or []) if c is not None]
        if choices and not is_list_arg(arg):
            bad = 999999 if arg.get("type") == "int" else "not-a-choice"
            assert "invalid choice" in form.check_values({arg["dest"]: bad}).get(arg["dest"], "")
    for group in spec.mutually_exclusive(sub):
        dests = group["dests"][:2]
        if len(dests) == 2:
            values = {d: sample_value(by_dest[d], str(tmp_path)) for d in dests}
            errors = form.check_values(values)
            assert all("not allowed together" in errors.get(d, "") for d in dests), errors


# -- filling -------------------------------------------------------------

def sample_value(arg, scratch):
    """A value for an argument that its type, choices and nargs accept."""
    action = arg.get("action") or "store"
    if action in ("store_true", "store_false"):
        return action == "store_true"
    count = {"+": 2, "*": 2, "?": 1}.get(arg.get("nargs"), arg.get("nargs") if isinstance(arg.get("nargs"), int) else 1)
    items = []
    for index in range(count):
        items.append(one_value(arg, scratch, index))
    return items if is_list_arg(arg) else items[0]


def one_value(arg, scratch, index):
    choices = [c for c in (arg.get("choices") or []) if c is not None]
    if choices:
        different = [c for c in choices if c != arg.get("default")]
        return (different or choices)[index % len(different or choices)]
    info = arg.get("file")
    type_name = arg.get("type")
    if info or type_name == "FileType":
        role = (info or {}).get("role", "input")
        kind = (info or {}).get("kind", "file")
        formats = (info or {}).get("formats") or ["txt"]
        name = "{}_{}{}.{}".format(arg["dest"], index, "_out" if role == "output" else "", formats[0])
        path = os.path.join(scratch, name)
        if role == "input":
            if kind == "directory":
                os.makedirs(path, exist_ok=True)
            else:
                with open(path, "w") as handle:
                    handle.write("sample\n")
        return path
    if type_name == "int":
        default = arg.get("default")
        return (default + 1 + index) if isinstance(default, int) and not isinstance(default, bool) else 3 + index
    if type_name == "float":
        return 0.25 + index
    if type_name == "genomicRegion":
        return "chr1:1,000-{},000".format(2 + index)
    return "value{}".format(index)


def fill(spec, sub, scratch):
    """Values for every settable argument the Python parser also has: one
    member of each mutually exclusive group, and no C++-only option."""
    values = {}
    skipped = set()
    for group in spec.mutually_exclusive(sub):
        skipped.update(group["dests"][1:])
    for arg in settable(spec, sub):
        if arg["dest"] in skipped or arg.get("cpp_only"):
            continue
        values[arg["dest"]] = sample_value(arg, scratch)
    return values


def expected_namespace(spec, sub, values):
    """The namespace argparse builds for these values: the given values, the
    defaults of the rest (string defaults go through the type, as argparse
    does), files as their names, and the subcommand."""
    out = {}
    for arg in spec.arguments(sub):
        action = arg.get("action") or "store"
        if action in NON_SETTABLE_ACTIONS or arg.get("cpp_only"):
            continue
        dest = arg["dest"]
        type_name = arg.get("type")
        if dest in values:
            value = values[dest]
            if action in ("store_true", "store_false"):
                out[dest] = value
                continue
        else:
            value = arg.get("default")
            if not (isinstance(value, str) and type_name in ("int", "float", "genomicRegion", "FileType")):
                out[dest] = value
                continue
        out[dest] = typed(value, type_name)
    if sub is not None and spec.subcommand_dest:
        out[spec.subcommand_dest] = sub
    return out


def typed(value, type_name):
    if isinstance(value, list):
        return [typed(v, type_name) for v in value]
    if type_name == "int":
        return int(value)
    if type_name == "float":
        return float(value)
    if type_name == "genomicRegion":
        region = "".join(str(value).split()).translate(str.maketrans("", "", ",;|!{}()")).replace("-", ":")
        return region or None
    if type_name == "FileType":
        return {"__file__": value}
    return value


@needs_reference
def test_filled_forms_parse_into_the_meant_namespace(qtbot, tmp_path):
    requests, expectations = [], []
    for tool, sub in CASES:
        spec = _spec(tool)
        scratch = tmp_path / "{}_{}".format(tool, sub or "main")
        scratch.mkdir()
        form = ToolForm(spec)
        qtbot.addWidget(form)
        form.set_subcommand(sub)
        values = fill(spec, sub, str(scratch))
        form.set_values(values)
        assert form.validate() == {}, (tool, sub, form.validate())
        assert form.values() == values, (tool, sub)
        argv = form.argv()
        requests.append({"tool": tool, "args": argv[1:]})
        expectations.append((tool, sub, argv, expected_namespace(spec, sub, values)))
    results = python_namespaces(requests, tmp_path)
    failures = []
    report = []
    for (tool, sub, argv, expected), result in zip(expectations, results):
        name = tool + (" " + sub if sub else "")
        if not result["ok"]:
            failures.append("{}: Python rejected {}: {}".format(name, argv, result["stderr"].strip()[-300:]))
            continue
        got = result["namespace"]
        if got != expected:
            diff = {k: (expected.get(k, "<missing>"), got.get(k, "<missing>"))
                    for k in set(got) | set(expected) if got.get(k, "<m>") != expected.get(k, "<m>")}
            failures.append("{}: namespace differs (expected, python): {}".format(name, diff))
        report.append("{}: {} arguments set, namespace of {} entries equal".format(name, len(argv) - 1, len(got)))
    print("\n".join(report))
    assert not failures, "\n".join(failures)
