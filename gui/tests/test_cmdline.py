import os

from conftest import chain_workflow, write_workflow
from hicexplorer_gui.workflow import SpecLoader, build_argv, load_workflow
from hicexplorer_gui.workflow.model import WorkflowError


def plan_argv(tmp_path, tools_dir, data, step_id="s"):
    path = write_workflow(tmp_path / "work", data)
    wf = load_workflow(path)
    assert wf.errors == []
    loader = SpecLoader(tools_dir)
    step = wf.step(step_id)
    return build_argv(wf, step, loader.load(step.tool))[1:]


def test_order_flags_and_values(tmp_path, tools_dir):
    data = {"version": 1, "inputs": {"x": "data/x.txt"}, "steps": [{
        "id": "s", "tool": "combine", "threads": 3,
        "args": {"value": 2, "noflag": False, "flag": False, "pair": [1.5, -2],
                 "outFileName": "${outputs.out}", "label": "L", "inputs": ["${inputs.x}", "lit.txt"],
                 "scale": 1e-05, "alpha": None, "mode": "sum"},
        "outputs": {"out": "s/out.txt"}}]}
    argv = plan_argv(tmp_path, tools_dir, data)
    # positional first, then options in spec order; long flags; store_false only when it differs;
    # store_true false (== default) omitted; null omitted; threads injected as --threads.
    assert argv == ["L", "--inputs", "data/x.txt", "lit.txt", "--outFileName", "s/out.txt",
                    "--mode", "sum", "--scale", "1e-05", "--pair", "1.5", "-2", "--noflag",
                    "--value", "2", "--threads", "3"]


def test_store_true_emitted_when_true(tmp_path, tools_dir):
    data = {"version": 1, "steps": [{"id": "s", "tool": "makeA",
                                     "args": {"outFileName": "${outputs.out}", "fail": True},
                                     "outputs": {"out": "o.txt"}}]}
    assert plan_argv(tmp_path, tools_dir, data) == ["--outFileName", "o.txt", "--threads", "1", "--fail"]


def test_subcommand_first(tmp_path, tools_dir):
    data = chain_workflow()
    argv = plan_argv(tmp_path, tools_dir, data, "d")
    assert argv == ["correct", "--matrix", "b/out.txt", "--outFileName", "d/out.txt",
                    "--filterThreshold", "-1.5", "5", "--threads", "1"]


def test_thread_option_detection(tmp_path, tools_dir):
    data = {"version": 1, "inputs": {"x": "data/x.txt"}, "steps": [
        {"id": "s", "tool": "prefixtool", "threads": 4,
         "args": {"matrix": "${inputs.x}", "outPrefix": "${outputs.p}"}, "outputs": {"p": "t/t"}},
        {"id": "u", "tool": "makeA", "threads": 2, "threads_option": "value",
         "args": {"outFileName": "${outputs.o}"}, "outputs": {"o": "u.txt"}}]}
    assert plan_argv(tmp_path, tools_dir, data, "s") == [
        "--matrix", "data/x.txt", "--outPrefix", "t/t", "--numberOfProcessors", "4"]
    assert plan_argv(tmp_path, tools_dir, data, "u") == ["--outFileName", "u.txt", "--value", "2"]


def test_paths_relative_to_workdir(tmp_path, tools_dir):
    data = chain_workflow()
    wf_path = write_workflow(tmp_path / "wfdir", data)
    wf = load_workflow(wf_path, workdir=str(tmp_path / "elsewhere"))
    step = wf.step("b")
    argv = build_argv(wf, step, SpecLoader(tools_dir).load("combine"))
    assert argv[0] == os.path.join(tools_dir, "combine")
    # inputs outside the workdir are absolute, outputs relative to the workdir
    assert argv[1:6] == ["--inputs", "a/out.txt", str(tmp_path / "wfdir" / "data" / "x.txt"),
                         "--outFileName", "b/out.txt"]


def test_topological_order_ties_by_file_position(tmp_path):
    data = {"version": 1, "steps": [
        {"id": "late", "tool": "t", "after": ["mid"]},
        {"id": "first", "tool": "t"},
        {"id": "mid", "tool": "t", "after": ["first"]},
        {"id": "other", "tool": "t"},
    ]}
    wf = load_workflow(write_workflow(tmp_path, data))
    assert [s.id for s in wf.topological_order()] == ["first", "mid", "late", "other"]
    assert wf.dependents("first") == {"mid", "late"}


def test_cycle_detection(tmp_path):
    data = {"version": 1, "steps": [
        {"id": "a", "tool": "t", "after": ["c"]},
        {"id": "b", "tool": "t", "after": ["a"]},
        {"id": "c", "tool": "t", "after": ["b"]},
        {"id": "free", "tool": "t"},
    ]}
    wf = load_workflow(write_workflow(tmp_path, data))
    assert "dependency cycle among steps: a, b, c" in wf.errors
    try:
        wf.topological_order()
    except WorkflowError as exc:
        assert "a, b, c" in str(exc)
    else:
        raise AssertionError("cycle not detected")
