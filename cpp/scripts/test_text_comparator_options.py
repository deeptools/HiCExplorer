"""Tests of the text and interval comparator options used for last-bit and tie-order differences.

    python -m pytest cpp/scripts/test_text_comparator_options.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from comparators import interval, text  # noqa: E402


def write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


PAIRS = "a\t1\t3\nb\t2\t3\nc\t3\t2\nd\t4\t1\ne\t5\t1\n"


def test_tie_column_ignores_the_order_inside_a_run(tmp_path):
    other = "b\t2\t3\na\t1\t3\nc\t3\t2\ne\t5\t1\nd\t4\t1\n"
    result = text.compare(write(tmp_path, "a", PAIRS), write(tmp_path, "b", other), "E0", {"tie_column": -1})
    assert result.passed and result.metrics["reordered_lines"] == 4


def test_tie_column_fails_when_runs_are_reordered(tmp_path):
    other = "c\t3\t2\na\t1\t3\nb\t2\t3\nd\t4\t1\ne\t5\t1\n"
    assert not text.compare(write(tmp_path, "a", PAIRS), write(tmp_path, "b", other), "E0", {"tie_column": -1}).passed


def test_tie_column_fails_for_a_different_line_or_count(tmp_path):
    changed = "a\t1\t3\nb\t2\t3\nc\t3\t2\nd\t4\t1\nx\t5\t1\n"
    assert not text.compare(write(tmp_path, "a", PAIRS), write(tmp_path, "b", changed), "E0", {"tie_column": -1}).passed
    shorter = "a\t1\t3\nb\t2\t3\nc\t3\t2\nd\t4\t1\n"
    assert not text.compare(write(tmp_path, "c", PAIRS), write(tmp_path, "d", shorter), "E0", {"tie_column": -1}).passed


def test_tie_column_fails_when_a_key_differs(tmp_path):
    changed = "a\t1\t3\nb\t2\t3\nc\t3\t2\nd\t4\t2\ne\t5\t1\n"
    assert not text.compare(write(tmp_path, "a", PAIRS), write(tmp_path, "b", changed), "E0", {"tie_column": -1}).passed


TAD = "# header 0.05\nchr1\t100\t200\tID_1\t0.5\t.\t0.123456789012345\n"


def test_e3_with_exact_comments(tmp_path):
    options = {"schema": "tsv-header", "exact_comments": True}
    close = TAD.replace("0.123456789012345", "0.123456789012346")
    assert text.compare(write(tmp_path, "a", TAD), write(tmp_path, "b", close), "E3", options).passed
    far = TAD.replace("0.123456789012345", "0.1234568")
    assert not text.compare(write(tmp_path, "c", TAD), write(tmp_path, "d", far), "E3", options).passed
    moved = TAD.replace("100\t200", "101\t200")
    assert not text.compare(write(tmp_path, "e", TAD), write(tmp_path, "f", moved), "E3", options).passed
    header = TAD.replace("0.05", "0.01")
    assert not text.compare(write(tmp_path, "g", TAD), write(tmp_path, "h", header), "E3", options).passed


def test_interval_columns(tmp_path):
    close = TAD.replace("0.123456789012345", "0.123456789012346")
    result = interval.compare(write(tmp_path, "a", TAD), write(tmp_path, "b", close), "E5", {"columns": 4})
    assert result.passed and result.class_met == "E5"
    moved = TAD.replace("100\t200", "101\t200")
    assert not interval.compare(write(tmp_path, "c", TAD), write(tmp_path, "d", moved), "E5", {"columns": 4}).passed
