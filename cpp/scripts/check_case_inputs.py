"""List harness case inputs ({data}/...) that are not tracked in a git revision.

Usage: cpp/scripts/check_case_inputs.py <repo> [<rev>]   (default rev: HEAD)

Contract rule 10 (cpp/AGENTS_CONTRACT.md) requires this check before a commit.

Case files are read from <rev> itself, so the check describes exactly what a
fresh checkout of that revision would see. A {data} path counts as tracked when
it is a tracked file, a tracked directory, or the prefix of tracked files (tools
such as hicFindTADs take an output prefix). A cooler URI 'file::/group' is
checked on 'file'.

An untracked input of a case that expects success fails the check. An untracked
input of a case that declares a non-zero expect_exit is listed separately and
does not fail it: such a case may name a missing file on purpose, for example to
pin an argument error. It is still printed so that a person sees each one,
because an error case can also point at a file that was accidentally left out,
and would then exit for the wrong reason.
"""
import json, re, subprocess, sys

repo = sys.argv[1]; rev = sys.argv[2] if len(sys.argv) > 2 else "HEAD"
git = lambda *a: subprocess.run(["git", "-C", repo, *a], capture_output=True, text=True, check=True).stdout
tracked = git("ls-tree", "-r", "--name-only", rev, "hicexplorer/test/test_data").split()
case_files = [p for p in git("ls-tree", "-r", "--name-only", rev, "cpp/scripts/cases").split() if p.endswith(".json")]

def untracked_refs(obj):
    refs = set(re.findall(r"\{data\}/([^\"'\s,\]]+)", json.dumps(obj)))
    out = []
    for ref in sorted(refs):
        path = "hicexplorer/test/test_data/" + ref.split("::", 1)[0]
        if not any(t == path or t.startswith(path) for t in tracked):
            out.append(ref)
    return out

missing, intentional = [], []
for case_file in case_files:
    data = json.loads(git("show", f"{rev}:{case_file}"))
    cases = data["cases"] if isinstance(data, dict) else data
    for case in cases:
        refs = untracked_refs(case)
        if not refs:
            continue
        expect = case.get("expect_exit", 0)
        bucket = intentional if expect not in (0, None) else missing
        bucket.extend((case["id"], expect, ref) for ref in refs)

for label, rows in (("MISSING (case expects success)", missing),
                    ("absent input in an error case (listed, not failing)", intentional)):
    if rows:
        print(label + ":")
        for cid, expect, ref in rows:
            print(f"    {cid}  [expect_exit {expect}]  {ref}")
print(f"{rev}: {len(case_files)} case files, {len(missing)} missing, {len(intentional)} absent in error cases")
sys.exit(1 if missing else 0)
