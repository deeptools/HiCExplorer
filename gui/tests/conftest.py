"""Fake tools for the workflow engine tests.

Each fake tool is a Python script written into a temporary directory. It prints
its specification for ``--help-json``; otherwise it parses its command line
with argparse built from that specification (so every command line the engine
builds is also checked for parseability), hashes its input file contents and
non-file parameters, and writes that digest to every output. The thread count
is left out of the digest, as real tools' outputs do not depend on it. With
``FAKE_LOG`` set it appends JSON lines with start and end times.
"""

import json
import os
import sys

import pytest
import yaml

GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if GUI_DIR not in sys.path:
    sys.path.insert(0, GUI_DIR)

FAKE_TOOL_BODY = r'''
import argparse, hashlib, json, os, signal, sys, time
spec = json.loads(SPEC)
argv = sys.argv[1:]
if argv == ["--help-json"]:
    print(json.dumps(spec))
    sys.exit(0)
args = []
for g in spec["groups"]:
    args += g["arguments"]
sub = None
if spec.get("subcommands"):
    sub = argv.pop(0)
    cmds = [c for c in spec["subcommands"]["commands"] if c["name"] == sub]
    if not cmds:
        sys.exit("unknown subcommand %s" % sub)
    for g in cmds[0]["groups"]:
        args += g["arguments"]
parser = argparse.ArgumentParser(prog=spec["tool"], add_help=False)
types = {"int": int, "float": float}
for a in args:
    if a["action"] in ("help", "version", "help_json"):
        continue
    kw = {}
    if a["action"] in ("store_true", "store_false"):
        kw["action"] = a["action"]
    else:
        if a["nargs"] is not None:
            kw["nargs"] = a["nargs"]
        if a["type"] in types:
            kw["type"] = types[a["type"]]
        if a["choices"] is not None:
            kw["choices"] = a["choices"]
        kw["default"] = a["default"]
    if a["positional"]:
        parser.add_argument(a["dest"], **kw)
    else:
        parser.add_argument(*a["flags"], dest=a["dest"], required=a["required"], **kw)
ns = parser.parse_args(argv)
log = os.environ.get("FAKE_LOG")
def note(event):
    if log:
        with open(log, "a") as fh:
            fh.write(json.dumps({"tool": spec["tool"], "event": event, "t": time.time(),
                                 "pid": os.getpid(), "threads": getattr(ns, "threads", None),
                                 "nproc": getattr(ns, "numberOfProcessors", None),
                                 "value": getattr(ns, "value", None), "argv": sys.argv}) + "\n")
if getattr(ns, "ignoreTerm", False):
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
note("start")
h = hashlib.sha256((spec["tool"] + str(sub)).encode())
outputs = []
for a in args:
    if a["action"] in ("help", "version", "help_json"):
        continue
    v = getattr(ns, a["dest"])
    f = a.get("file")
    vals = v if isinstance(v, list) else [v]
    if f and f["role"] == "output":
        outputs += [(x, f["kind"]) for x in vals if x is not None]
    elif f and f["role"] == "input":
        for x in vals:
            if x is None:
                continue
            if os.path.isdir(x):
                for root, dirs, files in os.walk(x):
                    dirs.sort()
                    for name in sorted(files):
                        with open(os.path.join(root, name), "rb") as fh:
                            h.update(fh.read())
            else:
                with open(x, "rb") as fh:
                    h.update(fh.read())
    elif a["dest"] not in ("threads", "numberOfProcessors"):
        h.update(json.dumps([a["dest"], v]).encode())
sys.stdout.write("processing %s\n" % spec["tool"])
sys.stderr.write("log line for %s\n" % spec["tool"])
if getattr(ns, "sleep", None):
    time.sleep(ns.sleep)
if getattr(ns, "fail", False):
    note("fail")
    sys.stderr.write("failing on purpose\n")
    sys.exit(3)
digest = h.hexdigest() + "\n"
for path, kind in outputs:
    if kind == "directory":
        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, "report.txt"), "w") as fh:
            fh.write(digest)
    elif kind == "prefix":
        for suffix in ("_a.txt", "_b.txt"):
            with open(path + suffix, "w") as fh:
                fh.write(suffix + digest)
    else:
        with open(path, "w") as fh:
            fh.write(digest)
note("end")
'''


def arg(dest, flags=None, action="store", type=None, nargs=None, choices=None, default=None,
        required=False, file=None, cpp_only=False, positional=False, note=None):
    if flags is None:
        flags = [] if positional else ["--" + dest]
    if action == "store_true" and default is None:
        default = False
    if action == "store_false" and default is None:
        default = True
    return {"dest": dest, "flags": flags, "positional": positional, "action": action,
            "type": type, "nargs": nargs, "choices": choices, "default": default,
            "required": required, "metavar": None, "help": "help for " + dest,
            "file": file, "cpp_only": cpp_only, "note": note}


def infile(kind="file"):
    return {"role": "input", "formats": ["txt"], "kind": kind}


def outfile(kind="file"):
    return {"role": "output", "formats": ["txt"], "kind": kind}


def make_spec(tool, arguments, mutually_exclusive=(), subcommands=None, version="1.0-fake"):
    return {"schema": "hicexplorer-tool-spec", "schema_version": 1, "tool": tool,
            "version": version, "description": "fake " + tool,
            "groups": [{"title": "Arguments", "arguments":
                        [arg("help", flags=["--help", "-h"], action="help")] + list(arguments)}],
            "mutually_exclusive": [dict(m) for m in mutually_exclusive],
            "subcommands": subcommands}


def common_args():
    return [arg("threads", flags=["--threads", "-t"], type="int", default=1),
            arg("sleep", type="float", default=0.0),
            arg("fail", action="store_true"),
            arg("ignoreTerm", action="store_true")]


def fake_specs():
    specs = []
    specs.append(make_spec("makeA", [
        arg("outFileName", flags=["--outFileName", "-o"], required=True, file=outfile()),
        arg("value", flags=["-v", "--value"], type="int", default=1),
    ] + common_args()))
    specs.append(make_spec("combine", [
        arg("label", positional=True, nargs="?", default="none"),
        arg("inputs", flags=["--inputs", "-i"], nargs="+", required=True, file=infile()),
        arg("outFileName", flags=["--outFileName", "-o"], required=True, file=outfile()),
        arg("mode", choices=["sum", "max"], default="sum"),
        arg("scale", type="float", default=1.0),
        arg("pair", type="float", nargs=2),
        arg("flag", action="store_true"),
        arg("noflag", action="store_false"),
        arg("alpha", type="int"),
        arg("beta", type="int"),
        arg("fast", action="store_true", cpp_only=True, note="SIMD path"),
        arg("value", type="int", default=1),
    ] + common_args(), mutually_exclusive=[{"required": False, "dests": ["alpha", "beta"]}]))
    specs.append(make_spec("subtool", [], subcommands={
        "dest": "command", "required": True, "commands": [
            {"name": "correct", "description": "correct", "mutually_exclusive": [],
             "groups": [{"title": "Required", "arguments": [
                 arg("matrix", flags=["--matrix", "-m"], required=True, file=infile()),
                 arg("outFileName", flags=["--outFileName", "-o"], required=True, file=outfile()),
                 arg("filterThreshold", type="float", nargs=2),
             ] + common_args()}]},
            {"name": "plot", "description": "plot", "mutually_exclusive": [],
             "groups": [{"title": "Required", "arguments": [
                 arg("matrix", flags=["--matrix", "-m"], required=True, file=infile()),
             ]}]},
        ]}))
    specs.append(make_spec("prefixtool", [
        arg("matrix", flags=["--matrix", "-m"], required=True, file=infile()),
        arg("outPrefix", required=True, file=outfile("prefix")),
        arg("numberOfProcessors", flags=["-p", "--numberOfProcessors"], type="int", default=1),
        arg("sleep", type="float", default=0.0),
    ]))
    specs.append(make_spec("dirtool", [
        arg("inputs", flags=["--inputs"], nargs="+", required=True, file=infile()),
        arg("QCfolder", required=True, file=outfile("directory")),
        arg("summary", required=True, file=outfile()),
    ]))
    return specs


@pytest.fixture
def tools_dir(tmp_path):
    directory = tmp_path / "bin"
    directory.mkdir()
    for spec in fake_specs():
        path = directory / spec["tool"]
        path.write_text("#!{}\nSPEC = {!r}\n{}".format(sys.executable, json.dumps(spec), FAKE_TOOL_BODY))
        path.chmod(0o755)
    return str(directory)


@pytest.fixture
def fake_log(tmp_path, monkeypatch):
    path = tmp_path / "fake.log"
    monkeypatch.setenv("FAKE_LOG", str(path))
    return path


def read_log(path):
    if not os.path.exists(str(path)):
        return []
    with open(str(path)) as handle:
        return [json.loads(line) for line in handle if line.strip()]


def write_workflow(directory, data, name="wf.yaml"):
    os.makedirs(str(directory), exist_ok=True)
    path = os.path.join(str(directory), name)
    with open(path, "w") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)
    return path


def chain_workflow(threads=4):
    """a (makeA) -> b (combine a + input) -> d (subtool correct b); c (combine input2) independent."""
    return {
        "version": 1, "name": "chain", "threads": threads,
        "inputs": {"x": "data/x.txt", "y": "data/y.txt"},
        "steps": [
            {"id": "a", "tool": "makeA", "args": {"outFileName": "${outputs.out}", "value": 5},
             "outputs": {"out": "a/out.txt"}},
            {"id": "b", "tool": "combine", "threads": 2,
             "args": {"inputs": ["${steps.a.outputs.out}", "${inputs.x}"], "outFileName": "${outputs.out}",
                      "mode": "max"},
             "outputs": {"out": "b/out.txt"}},
            {"id": "c", "tool": "combine",
             "args": {"inputs": ["${inputs.y}"], "outFileName": "${outputs.out}"},
             "outputs": {"out": "c/out.txt"}},
            {"id": "d", "tool": "subtool", "subcommand": "correct",
             "args": {"matrix": "${steps.b.outputs.out}", "outFileName": "${outputs.out}",
                      "filterThreshold": [-1.5, 5]},
             "outputs": {"out": "d/out.txt"}},
        ],
    }


@pytest.fixture
def chain(tmp_path):
    """Workdir with data files and the chain workflow; returns the workflow path."""
    work = tmp_path / "work"
    (work / "data").mkdir(parents=True)
    (work / "data" / "x.txt").write_text("x-content\n")
    (work / "data" / "y.txt").write_text("y-content\n")
    return write_workflow(work, chain_workflow())
