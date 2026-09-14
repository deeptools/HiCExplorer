"""Parses command lines with the Python tools' own argparse parsers.

    python argparse_namespace.py requests.json

requests.json: [{"tool": "hicInfo", "args": ["-m", "a.cool"]}, ...]. Prints a
JSON list, one entry per request: {"ok": true, "namespace": {...}} or
{"ok": false, "exit": code, "stderr": text}. File objects that FileType opens
are reported as {"__file__": name} and closed. Run it with the reference
interpreter and the repository root on PYTHONPATH, in a scratch directory
(FileType('w') creates its files while parsing).
"""

import contextlib
import importlib
import io
import json
import os
import sys

# A tool that delegates to another package is parsed by that package's parser,
# the one cpp/scripts/tool_specs.py compares the C++ spec against
# (hicPlotTADs: pygenometracks.plotTracks). Every other tool is
# hicexplorer.<tool>.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                                "cpp", "scripts"))
from tool_specs import PARSER_MODULES  # noqa: E402


def plain(value):
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if hasattr(value, "name") and hasattr(value, "close"):
        name = value.name
        if value not in (sys.stdin, sys.stdout, sys.stderr):
            value.close()
        return {"__file__": name}
    return {"__repr__": repr(value)}


def main():
    with open(sys.argv[1]) as handle:
        requests = json.load(handle)
    results = []
    for request in requests:
        module = importlib.import_module(PARSER_MODULES.get(request["tool"], "hicexplorer." + request["tool"]))
        parser = module.parse_arguments()
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
                namespace = parser.parse_args(request["args"])
            results.append({"ok": True, "namespace": {k: plain(v) for k, v in vars(namespace).items()}})
        except SystemExit as exc:
            results.append({"ok": False, "exit": exc.code, "stderr": err.getvalue()})
    json.dump(results, sys.stdout)


if __name__ == "__main__":
    main()
