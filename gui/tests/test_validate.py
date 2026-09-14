import os
import subprocess
import sys

import pytest

from conftest import GUI_DIR, chain_workflow, write_workflow
from hicexplorer_gui.workflow import SpecLoader, load_workflow, validate_workflow
from hicexplorer_gui.workflow.model import WorkflowError


def errors_for(tmp_path, tools_dir, data, create=("data/x.txt", "data/y.txt")):
    work = tmp_path / "work"
    for rel in create:
        (work / rel).parent.mkdir(parents=True, exist_ok=True)
        (work / rel).write_text(rel)
    path = write_workflow(work, data)
    wf = load_workflow(path)
    msgs = validate_workflow(wf, SpecLoader(tools_dir))
    return [m.text for m in msgs if m.level == "error"], [m.text for m in msgs if m.level == "info"]


def one_step(tool, args, outputs=None, **extra):
    step = {"id": "s", "tool": tool, "args": args, "outputs": outputs or {"out": "s/out.txt"}}
    step.update(extra)
    return {"version": 1, "name": "t", "inputs": {"x": "data/x.txt"}, "steps": [step]}


def combine(**args):
    base = {"inputs": ["${inputs.x}"], "outFileName": "${outputs.out}"}
    base.update(args)
    return one_step("combine", base)


def test_valid_chain_has_no_errors(tmp_path, tools_dir):
    errors, infos = errors_for(tmp_path, tools_dir, chain_workflow())
    assert errors == []


@pytest.mark.parametrize("args, expected", [
    ({"bogus": 1}, "step s: bogus: unknown argument of combine"),
    ({"scale": "abc"}, "step s: scale: expected a number"),
    ({"alpha": 1.5}, "step s: alpha: expected an integer"),
    ({"alpha": "7x"}, "step s: alpha: expected an integer"),
    ({"flag": "yes"}, "step s: flag: expected true or false"),
    ({"mode": "median"}, "step s: mode: invalid choice 'median'"),
    ({"pair": [1, 2, 3]}, "step s: pair: expected exactly 2 values, found 3"),
    ({"inputs": []}, "step s: inputs: expected at least one value"),
    ({"scale": [1.0, 2.0]}, "step s: scale: expected a single value"),
    ({"alpha": 1, "beta": 2}, "step s: alpha, beta: mutually exclusive"),
    ({"help": True}, "step s: help: action 'help' cannot be set"),
    ({"threads": 3}, "step s: threads: set the step's threads instead"),
    ({"inputs": ["data/missing.txt"]}, "step s: inputs: input data/missing.txt does not exist"),
    ({"outFileName": "../outside.txt"}, "step s: outFileName: output ../outside.txt is not under the workdir"),
])
def test_argument_errors(tmp_path, tools_dir, args, expected):
    errors, _ = errors_for(tmp_path, tools_dir, combine(**args))
    assert any(e.startswith(expected) for e in errors), errors


def test_required_argument_missing(tmp_path, tools_dir):
    errors, _ = errors_for(tmp_path, tools_dir, one_step("combine", {"outFileName": "${outputs.out}"}))
    assert "step s: inputs: required argument missing" in errors


def test_required_mutually_exclusive_group(tmp_path, tools_dir, monkeypatch):
    from hicexplorer_gui.workflow import spec as specmod
    original = specmod.ToolSpec.mutually_exclusive
    monkeypatch.setattr(specmod.ToolSpec, "mutually_exclusive",
                        lambda self, sub=None: [dict(g, required=True) for g in original(self, sub)])
    errors, _ = errors_for(tmp_path, tools_dir, combine())
    assert "step s: alpha, beta: one of these arguments is required" in errors


def test_cpp_only_is_info(tmp_path, tools_dir):
    errors, infos = errors_for(tmp_path, tools_dir, combine(fast=True))
    assert errors == []
    assert "step s: fast: C++-only option (SIMD path)" in infos


def test_subcommand_errors(tmp_path, tools_dir):
    base = {"matrix": "${inputs.x}", "outFileName": "${outputs.out}"}
    errors, _ = errors_for(tmp_path, tools_dir, one_step("subtool", base))
    assert any(e.startswith("step s: subcommand: required") for e in errors)
    errors, _ = errors_for(tmp_path, tools_dir, one_step("subtool", base, subcommand="nope"))
    assert any(e.startswith("step s: subcommand: unknown 'nope'") for e in errors)
    errors, _ = errors_for(tmp_path, tools_dir, one_step("combine", {}, subcommand="correct"))
    assert "step s: subcommand: tool combine has no subcommands" in errors
    errors, _ = errors_for(tmp_path, tools_dir, one_step("subtool", dict(base, filterThreshold=[1, 2]),
                                                          subcommand="correct"))
    assert errors == []


def test_unknown_tool(tmp_path, tools_dir):
    errors, _ = errors_for(tmp_path, tools_dir, one_step("noSuchTool", {}))
    assert len(errors) == 1 and errors[0].startswith("step s: tool noSuchTool not found")


def test_structural_errors(tmp_path, tools_dir):
    data = chain_workflow()
    data["steps"][2]["id"] = "a"                                             # duplicate id
    data["steps"][3]["args"]["matrix"] = "${steps.zzz.outputs.out}"           # unknown step
    data["steps"][1]["args"]["outFileName"] = "${outputs.nope}"              # undeclared output
    data["steps"][0]["args"]["value"] = "${inputs.nope}"                      # unknown input
    data["steps"][0]["after"] = ["missing"]
    data["steps"][0]["frobnicate"] = 1
    errors, _ = errors_for(tmp_path, tools_dir, data)
    for expected in ["step a: duplicate id",
                     "step d: matrix: unknown step in ${steps.zzz.outputs.out}",
                     "step b: outFileName: ${outputs.nope} is not declared in this step's outputs",
                     "step a: value: unknown input in ${inputs.nope}",
                     "step a: after: unknown step 'missing'",
                     "step a: unknown key 'frobnicate'"]:
        assert expected in errors, errors


def test_input_from_non_upstream_step(tmp_path, tools_dir):
    data = chain_workflow()
    # c reads b's output by literal path without depending on b
    data["steps"][2]["args"]["inputs"] = ["b/out.txt"]
    errors, _ = errors_for(tmp_path, tools_dir, data)
    assert any(e.startswith("step c: inputs: input b/out.txt is output out of step b, which is not upstream")
               for e in errors), errors
    data["steps"][2]["after"] = ["b"]
    errors, _ = errors_for(tmp_path, tools_dir, data)
    assert errors == []


def test_cycle_is_an_error(tmp_path, tools_dir):
    data = chain_workflow()
    data["steps"][0]["after"] = ["d"]
    errors, _ = errors_for(tmp_path, tools_dir, data)
    assert "dependency cycle among steps: a, b, d" in errors


def test_all_errors_reported_together_by_cli(tmp_path, tools_dir):
    data = combine(bogus=1, mode="median", pair=[1])
    work = tmp_path / "work"
    (work / "data").mkdir(parents=True)
    (work / "data" / "x.txt").write_text("x")
    path = write_workflow(work, data)
    env = dict(os.environ, PYTHONPATH=GUI_DIR)
    proc = subprocess.run([sys.executable, "-m", "hicexplorer_gui.workflow.cli", "validate", path,
                           "--tools-dir", tools_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, env=env)
    assert proc.returncode == 1
    lines = [l for l in proc.stderr.splitlines() if l.startswith("error: ")]
    assert len(lines) == 3 and all(l.startswith("error: step s: ") for l in lines)


def test_cli_usage_errors(tmp_path, tools_dir):
    env = dict(os.environ, PYTHONPATH=GUI_DIR)
    run = lambda *a: subprocess.run([sys.executable, "-m", "hicexplorer_gui.workflow.cli"] + list(a),
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env).returncode
    assert run() == 2
    assert run("validate", str(tmp_path / "missing.yaml")) == 2
    path = write_workflow(tmp_path, chain_workflow())
    assert run("run", path, "--threads", "0") == 2
    assert run("export", path, "--format", "zip", "-o", "x") == 2


def test_external_outputs_only_when_enabled(tmp_path, tools_dir):
    """Workflow files keep outputs under the workdir; single tool runs from
    the GUI (external_outputs) may write any existing, writable directory."""
    work = tmp_path / "work"
    (work / "data").mkdir(parents=True)
    (work / "data" / "x.txt").write_text("x")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    target = str(elsewhere / "out.txt")
    data = one_step("combine", {"inputs": ["${inputs.x}"], "outFileName": "${outputs.out}"},
                    outputs={"out": target})
    path = write_workflow(work, data)
    loader = SpecLoader(tools_dir)

    plain = load_workflow(path)
    assert any("output out ({}) is not under the workdir".format(target) in e for e in plain.errors)

    wf = load_workflow(path, external_outputs=True)
    assert wf.errors == []
    assert wf.steps[0].outputs == {"out": target}
    assert wf.output_abspath(wf.steps[0], "out") == target
    assert [m.text for m in validate_workflow(wf, loader) if m.level == "error"] == []

    data["steps"][0]["outputs"] = {"out": str(tmp_path / "missing" / "out.txt")}
    wf = load_workflow(write_workflow(work, data), external_outputs=True)
    errors = [m.text for m in validate_workflow(wf, loader) if m.level == "error"]
    assert any("does not exist" in e for e in errors), errors


def test_unreadable_workflow(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("- just\n- a list\n")
    with pytest.raises(WorkflowError):
        load_workflow(str(path))
