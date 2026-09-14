import filecmp
import json
import os
import re
import shutil
import subprocess
import sys

import pytest

from conftest import GUI_DIR, write_workflow
from hicexplorer_gui.workflow import EXIT_OK, Runner, SpecLoader, load_workflow
from hicexplorer_gui.workflow.export import export_sh, export_snakemake, marker, rule_name
from hicexplorer_gui.workflow.plan import make_plans


def full_workflow(work):
    (work / "data").mkdir(parents=True)
    (work / "data" / "x.txt").write_text("x-content\n")
    (work / "data" / "y.txt").write_text("y-content\n")
    data = {"version": 1, "name": "export test", "threads": 4,
            "inputs": {"x": "data/x.txt", "y": "data/y.txt"},
            "steps": [
                {"id": "d", "tool": "subtool", "subcommand": "correct",
                 "args": {"matrix": "${steps.b.outputs.out}", "outFileName": "${outputs.out}",
                          "filterThreshold": [-1.5, 5]}, "outputs": {"out": "d/deep/out.txt"}},
                {"id": "a", "tool": "makeA", "args": {"outFileName": "${outputs.out}", "value": 5},
                 "outputs": {"out": "a/out.txt"}},
                {"id": "b", "tool": "combine", "threads": 2,
                 "args": {"label": "it's {braced}", "inputs": ["${steps.a.outputs.out}", "${inputs.x}"],
                          "outFileName": "${outputs.out}"}, "outputs": {"out": "b/out.txt"}},
                {"id": "p", "tool": "prefixtool", "args": {"matrix": "${inputs.y}", "outPrefix": "${outputs.t}"},
                 "outputs": {"t": "tads/t"}},
                {"id": "q", "tool": "dirtool", "after": ["p"],
                 "args": {"inputs": ["tads/t_a.txt", "${steps.d.outputs.out}"], "QCfolder": "${outputs.qc}",
                          "summary": "${outputs.s}"},
                 "outputs": {"qc": "qc/folder", "s": "qc/summary.txt"}},
            ]}
    return write_workflow(work, data)


def tree(root):
    files = {}
    for base, dirs, names in os.walk(root):
        dirs[:] = [d for d in dirs if d != ".hicexplorer-workflow" and d != "data"]
        for name in names:
            if not name.endswith((".sh", ".yaml", "Snakefile")):
                full = os.path.join(base, name)
                files[os.path.relpath(full, root)] = open(full, "rb").read()
    return files


def test_sh_export_reproduces_run(tmp_path, tools_dir):
    work = tmp_path / "work"
    path = full_workflow(work)
    wf = load_workflow(path)
    loader = SpecLoader(tools_dir)
    assert Runner(wf, loader, log=lambda m: None).run() == EXIT_OK
    script = tmp_path / "wf.sh"
    env = dict(os.environ, PYTHONPATH=GUI_DIR)
    subprocess.check_call([sys.executable, "-m", "hicexplorer_gui.workflow.cli", "export", path,
                           "--tools-dir", tools_dir, "--format", "sh", "-o", str(script)], env=env)
    text = script.read_text()
    assert text.startswith("#!/bin/sh\n") and "\nset -eu\n" in text
    plans = make_plans(wf, loader)
    commands = [plans[s.id].command for s in wf.topological_order()]
    positions = [text.index(c) for c in commands]
    assert positions == sorted(positions)

    other = tmp_path / "other"
    shutil.copytree(str(work / "data"), str(other / "data"))
    subprocess.check_call(["sh", str(script), str(other)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    expected = tree(str(work))
    assert sorted(expected) == sorted(["a/out.txt", "b/out.txt", "d/deep/out.txt", "qc/folder/report.txt",
                                       "qc/summary.txt", "tads/t_a.txt", "tads/t_b.txt"])
    assert tree(str(other)) == expected


def test_snakefile_structure(tmp_path, tools_dir):
    path = full_workflow(tmp_path / "work")
    wf = load_workflow(path)
    loader = SpecLoader(tools_dir)
    text = export_snakemake(wf, loader)
    plans = make_plans(wf, loader)
    rules = re.findall(r"^rule (\w+):$", text, flags=re.M)
    assert rules == ["all"] + [rule_name(s.id) for s in wf.topological_order()]
    blocks = dict(zip(rules, re.split(r"^rule \w+:$", text, flags=re.M)[1:]))
    for step in wf.steps:
        block = blocks[rule_name(step.id)]
        shell = json.loads(block.split("    shell:\n")[1].strip().splitlines()[0])
        unescaped = shell.replace("{{", "{").replace("}}", "}")
        assert unescaped.endswith(plans[step.id].command)
        assert "    threads: {}\n".format(step.threads) in block
        assert 'touch("{}")'.format(marker(step.id)) in block
        for dep in step.deps:
            assert '"{}",'.format(marker(dep)) in block
        assert '"{}",'.format(marker(step.id)) in blocks["all"]
    assert "{{braced}}" in blocks[rule_name("b")]
    assert 'directory("qc/folder")' in blocks[rule_name("q")]
    assert '"data/y.txt",' in blocks[rule_name("p")]
    assert "tads/t" not in blocks[rule_name("p")].split("    output:")[1].split("threads")[0]
    assert text.count("workdir: ") == 1


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake is not installed")
def test_snakefile_runs(tmp_path, tools_dir):
    work = tmp_path / "work"
    path = full_workflow(work)
    loader = SpecLoader(tools_dir)
    assert Runner(load_workflow(path), loader, log=lambda m: None).run() == EXIT_OK
    other = tmp_path / "other"
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(export_snakemake(load_workflow(path, workdir=str(other)), loader))
    subprocess.check_call(["snakemake", "-s", str(snakefile), "--cores", "4", "--quiet"])
    assert tree(str(other)) == tree(str(work))


def test_export_sh_matches_function(tmp_path, tools_dir):
    path = full_workflow(tmp_path / "work")
    wf = load_workflow(path)
    text = export_sh(wf, SpecLoader(tools_dir))
    assert "mkdir -p d/deep" in text and "\"it's {braced}\"" not in text
