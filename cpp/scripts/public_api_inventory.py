"""Inventory of the public Python API of hicexplorer/ (outside tests)."""
import ast, re, sys, subprocess
from pathlib import Path
repo = Path(sys.argv[1]); pkg = repo / "hicexplorer"
cpp_files = [p for p in (repo / "cpp").rglob("*") if p.suffix in (".cpp", ".hpp") and "build" not in p.parts]
cpp_text = {p: p.read_text(errors="replace") for p in cpp_files}

def sig(node):
    a = node.args
    parts = []
    pos = a.posonlyargs + a.args
    defaults = [None] * (len(pos) - len(a.defaults)) + list(a.defaults)
    for arg, d in zip(pos, defaults):
        parts.append(arg.arg + ("=" + ast.unparse(d) if d is not None else ""))
    if a.vararg: parts.append("*" + a.vararg.arg)
    elif a.kwonlyargs: parts.append("*")
    for arg, d in zip(a.kwonlyargs, a.kw_defaults):
        parts.append(arg.arg + ("=" + ast.unparse(d) if d is not None else ""))
    if a.kwarg: parts.append("**" + a.kwarg.arg)
    return "(" + ", ".join(parts) + ")"

def cpp_refs(module_file, name):
    stem = module_file.stem
    hits = set()
    pat_name = re.compile(r"\b" + re.escape(name) + r"\b")
    for p, t in cpp_text.items():
        if pat_name.search(t) and (stem in t or module_file.name in t or len(name) > 12):
            hits.add(str(p.relative_to(repo)))
    return sorted(hits)

rows = []
for path in sorted(pkg.rglob("*.py")):
    rel = path.relative_to(repo)
    if rel.parts[1] == "test":
        continue
    text = path.read_text(); tree = ast.parse(text)
    mod = ".".join(path.relative_to(repo).with_suffix("").parts).replace(".__init__", "")
    tool_module = (repo / "cpp" / "tools" / (path.stem + ".cpp")).is_file()
    rows.append([mod, "module", "", str(len(text.splitlines())), "", "cpp/tools/%s.cpp" % path.stem if tool_module else ""])
    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.ClassDef)) or node.name.startswith("_"):
            continue
        n = node.end_lineno - node.lineno + 1
        src = ast.get_source_segment(text, node)
        if isinstance(node, ast.ClassDef):
            rows.append([mod, "class", node.name, str(n), "pure Python", ";".join(cpp_refs(path, node.name))])
            for sub in node.body:
                if isinstance(sub, ast.FunctionDef) and (not sub.name.startswith("_") or sub.name == "__init__"):
                    rows.append([mod, "method", node.name + "." + sub.name + sig(sub), str(sub.end_lineno - sub.lineno + 1),
                                 "pure Python", ";".join(cpp_refs(path, sub.name))])
            continue
        thin = mod in ("hicexplorer._cpp",) or (n <= 6 and ("build_parser(" in src or "entry_point(" in src))
        kind = "C++-backed (runs/reads the executable)" if thin else "pure Python"
        if mod == "hicexplorer.list_tools":
            kind = "Python (package metadata + --help-json)"
        if mod.startswith(("hicexplorer.hicHyperopt",)):
            kind = "Python driver (hyperopt) around C++ tools"
        refs = "" if thin else ";".join(cpp_refs(path, node.name))
        rows.append([mod, "def", node.name + sig(node), str(n), kind, refs])
print("module\tkind\tname(signature)\tlines\timplementation\tC++ code that names it (candidate backing)")
for r in rows:
    print("\t".join(r))
