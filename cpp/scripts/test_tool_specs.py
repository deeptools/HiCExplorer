"""pytest over cpp/scripts/tool_specs.py: one test per ported tool.

    HICX_CPP_BIN=BUILD/tools PYTHONPATH=. $VENV/bin/python -m pytest cpp/scripts/test_tool_specs.py

Skipped when HICX_CPP_BIN is not set.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import tool_specs  # noqa: E402

CPP_BIN = os.environ.get("HICX_CPP_BIN")


@pytest.mark.skipif(not CPP_BIN, reason="HICX_CPP_BIN is not set")
@pytest.mark.parametrize("tool", tool_specs.ported_tools())
def test_spec_matches_argparse(tool):
    unexpected, _allowed, stale = tool_specs.check_tool(tool, CPP_BIN)
    lines = [f"{d[0]}{d[1]} {d[2]}: python={d[3]!r} cpp={d[4]!r}" for d in unexpected]
    lines += [f"stale allowlist entry: {e}" for e in stale]
    assert not lines, "\n".join(lines)


def test_allowlist_entries_are_well_formed():
    for tool, entries in tool_specs.allowlist().items():
        assert tool in tool_specs.ported_tools(), tool
        for entry in entries:
            assert entry["kind"] in ("status", "finding"), entry
            assert entry["fields"] and entry["reason"], entry
