#!/usr/bin/env python
"""Compares every C++ tool's --help-json specification with the argparse parser
of the Python tool it ports (cpp/PLAN.md tier 10, item 10.1).

    python cpp/scripts/tool_specs.py --cpp-bin BUILD/tools [--tool NAME ...] [--json OUT]

Run it with the reference venv interpreter (cpp/AGENTS_CONTRACT.md), with the
repository root on PYTHONPATH so that `hicexplorer.<tool>.parse_arguments()`
imports. Exit status 0 when every difference is an allowed one.

What is compared, per argument (matched by dest): the option strings, the
action, the type, nargs, choices, the default and required; per parser the
mutually exclusive groups; and the subcommands, recursively. Help texts,
metavars and group titles are not compared.

Allowed differences:
- arguments the C++ spec marks `cpp_only` that the Python parser does not
  have (and --help-json itself);
- entries of cpp/scripts/tool_spec_deviations.json, each naming the tool, the
  dest, the fields that differ and the reason. Every entry has a "kind":
  "status" for a deviation already recorded in cpp/STATUS.md, or "finding"
  for a difference first found by this comparison (reported, not fixed).
Everything else is a failure, and so is an allowlist entry that no longer
matches a difference.
"""

import argparse
import importlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEVIATIONS = os.path.join(HERE, "tool_spec_deviations.json")
SCHEMA = "hicexplorer-tool-spec"
COMPARED = ("flags", "action", "type", "nargs", "choices", "default", "required")

ACTION_NAMES = {
    "_StoreAction": "store",
    "_StoreTrueAction": "store_true",
    "_StoreFalseAction": "store_false",
    "_StoreConstAction": "store_const",
    "_AppendAction": "append",
    "_AppendConstAction": "append_const",
    "_CountAction": "count",
    "_HelpAction": "help",
    "_VersionAction": "version",
    "BooleanOptionalAction": "boolean_optional",
}


def ported_tools():
    tools = []
    for name in sorted(os.listdir(os.path.join(REPO, "cpp", "tools"))):
        stem, ext = os.path.splitext(name)
        if ext == ".cpp" and stem.startswith(("hic", "chic")):
            tools.append(stem)
    return tools


def json_value(value):
    """A Python default or choice as the JSON the C++ spec writes."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    return repr(value)


def type_name(action):
    kind = ACTION_NAMES.get(type(action).__name__, type(action).__name__)
    if kind not in ("store", "append"):
        return None
    t = action.type
    if t is None or t is str:
        return "str"
    if t is int:
        return "int"
    if t is float:
        return "float"
    if isinstance(t, argparse.FileType):
        return "FileType"
    return getattr(t, "__name__", type(t).__name__)


def normalized_nargs(kind, nargs):
    if kind in ("store_true", "store_false", "store_const", "help", "version", "append_const",
                "count", "help_json"):
        return None
    return nargs


def python_parser_spec(parser):
    arguments = {}
    subcommands = None
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            subcommands = {
                "dest": None if action.dest == argparse.SUPPRESS else action.dest,
                "required": bool(action.required),
                "commands": {name: python_parser_spec(sub) for name, sub in action.choices.items()},
            }
            continue
        kind = ACTION_NAMES.get(type(action).__name__, type(action).__name__)
        arguments[action.dest] = {
            "flags": list(action.option_strings),
            "action": kind,
            "type": type_name(action),
            "nargs": normalized_nargs(kind, action.nargs),
            "choices": None if action.choices is None else [json_value(c) for c in action.choices],
            "default": json_value(action.default),
            "required": bool(action.required),
        }
    mutex = sorted(sorted(a.dest for a in g._group_actions) for g in parser._mutually_exclusive_groups)
    return {"arguments": arguments, "mutually_exclusive": mutex, "subcommands": subcommands}


def cpp_parser_spec(spec):
    arguments = {}
    cpp_only = set()
    for group in spec["groups"]:
        for arg in group["arguments"]:
            arguments[arg["dest"]] = {
                "flags": list(arg["flags"]),
                "action": arg["action"],
                "type": arg["type"],
                "nargs": normalized_nargs(arg["action"], arg["nargs"]),
                "choices": arg["choices"],
                "default": arg["default"],
                "required": bool(arg["required"]),
            }
            if arg.get("cpp_only") or arg["action"] == "help_json":
                cpp_only.add(arg["dest"])
    mutex = sorted(sorted(g["dests"]) for g in spec.get("mutually_exclusive", []))
    subcommands = None
    if spec.get("subcommands"):
        s = spec["subcommands"]
        subcommands = {"dest": s["dest"], "required": bool(s["required"]),
                       "commands": {c["name"]: cpp_parser_spec(c) for c in s["commands"]}}
    return {"arguments": arguments, "mutually_exclusive": mutex, "subcommands": subcommands,
            "cpp_only": cpp_only}


def same(field, a, b):
    if field == "flags":
        return sorted(a) == sorted(b)
    if isinstance(a, float) or isinstance(b, float):
        return isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) \
            and not isinstance(b, bool) and float(a) == float(b) and type(a) is type(b)
    return a == b


def compare(tool, py, cpp, where=""):
    """Differences as (where, dest, field, python value, C++ value)."""
    out = []
    for dest, p in py["arguments"].items():
        c = cpp["arguments"].get(dest)
        if c is None:
            out.append((where, dest, "missing in C++", p, None))
            continue
        for field in COMPARED:
            if not same(field, p[field], c[field]):
                out.append((where, dest, field, p[field], c[field]))
    for dest in cpp["arguments"]:
        if dest not in py["arguments"] and dest not in cpp["cpp_only"]:
            out.append((where, dest, "only in C++ but not marked cpp_only", None, cpp["arguments"][dest]))
    for dest in sorted(cpp["cpp_only"]):
        if dest in py["arguments"]:
            out.append((where, dest, "marked cpp_only but the Python has it", py["arguments"][dest], None))
    if py["mutually_exclusive"] != cpp["mutually_exclusive"]:
        out.append((where, "", "mutually_exclusive", py["mutually_exclusive"], cpp["mutually_exclusive"]))
    ps, cs = py["subcommands"], cpp["subcommands"]
    if (ps is None) != (cs is None):
        out.append((where, "", "subcommands", ps and sorted(ps["commands"]), cs and sorted(cs["commands"])))
    elif ps is not None:
        if ps["dest"] != cs["dest"] or ps["required"] != cs["required"]:
            out.append((where, "", "subcommand dest/required", [ps["dest"], ps["required"]],
                        [cs["dest"], cs["required"]]))
        if sorted(ps["commands"]) != sorted(cs["commands"]):
            out.append((where, "", "subcommand names", sorted(ps["commands"]), sorted(cs["commands"])))
        for name in ps["commands"]:
            if name in cs["commands"]:
                out.extend(compare(tool, ps["commands"][name], cs["commands"][name], where + name + " "))
    return out


def load_cpp_spec(cpp_bin, tool):
    result = subprocess.run([os.path.join(cpp_bin, tool), "--help-json"], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{tool} --help-json exited {result.returncode}: {result.stderr.strip()[:300]}")
    spec = json.loads(result.stdout)
    if spec.get("schema") != SCHEMA or spec.get("schema_version") != 1 or spec.get("tool") != tool:
        raise RuntimeError(f"{tool} --help-json: unexpected schema header")
    return spec


# Tools whose Python module has no parse_arguments of its own. hicPlotTADs is
# `pygenometracks.plotTracks.main(args)` (hicexplorer/hicPlotTADs.py), so its
# parser is pyGenomeTracks'; bin/hicQC runs hicPrepareQCreport.main, and there
# is no hicexplorer.hicQC. The gui tests read the same list.
PARSER_MODULES = {"hicPlotTADs": "pygenometracks.plotTracks",
                  "hicQC": "hicexplorer.hicPrepareQCreport"}


def load_python_parser(tool):
    module = importlib.import_module(PARSER_MODULES.get(tool, "hicexplorer." + tool))
    return module.parse_arguments()


def allowlist():
    with open(DEVIATIONS) as handle:
        return json.load(handle)


def check_tool(tool, cpp_bin):
    """(unexpected differences, allowed differences, stale allowlist entries)."""
    differences = compare(tool, python_parser_spec(load_python_parser(tool)),
                          cpp_parser_spec(load_cpp_spec(cpp_bin, tool)))
    entries = allowlist().get(tool, [])
    used = [False] * len(entries)
    unexpected, allowed = [], []
    for diff in differences:
        where, dest, field = diff[0].strip(), diff[1], diff[2]
        match = None
        for i, entry in enumerate(entries):
            if entry.get("where", "") == where and entry["dest"] == dest and field in entry["fields"]:
                match = i
                break
        if match is None:
            unexpected.append(diff)
        else:
            used[match] = True
            allowed.append(diff + (entries[match]["kind"], entries[match]["reason"]))
    stale = [e for e, u in zip(entries, used) if not u]
    return unexpected, allowed, stale


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cpp-bin", required=True)
    parser.add_argument("--tool", nargs="*")
    parser.add_argument("--json")
    args = parser.parse_args()
    report = {}
    failed = False
    for tool in args.tool or ported_tools():
        try:
            unexpected, allowed, stale = check_tool(tool, args.cpp_bin)
        except Exception as error:  # noqa: BLE001
            print(f"ERROR {tool}: {error}")
            report[tool] = {"error": str(error)}
            failed = True
            continue
        status = "PASS" if not unexpected and not stale else "FAIL"
        failed |= status == "FAIL"
        print(f"{status} {tool}: {len(allowed)} allowed, {len(unexpected)} unexpected, {len(stale)} stale")
        for d in unexpected:
            print(f"    unexpected {d[0]}{d[1]} {d[2]}: python={d[3]!r} cpp={d[4]!r}")
        for d in allowed:
            print(f"    allowed ({d[5]}) {d[0]}{d[1]} {d[2]}: python={d[3]!r} cpp={d[4]!r}")
        for e in stale:
            print(f"    stale allowlist entry {e}")
        report[tool] = {"unexpected": [list(map(repr, d)) for d in unexpected],
                        "allowed": [list(map(repr, d)) for d in allowed], "stale": stale}
    if args.json:
        with open(args.json, "w") as handle:
            json.dump(report, handle, indent=1)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
