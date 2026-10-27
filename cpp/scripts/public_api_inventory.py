"""Inventory of the public Python API of hicexplorer/ (outside tests).

    python public_api_inventory.py REPO                 # the TSV on stdout
    python public_api_inventory.py REPO --output FILE   # the TSV into FILE
    python public_api_inventory.py REPO --check         # remaining work, exit 1 if any

Every public module, class, method and function of the package is a row. The
implementation column says how it is implemented:

  C++-backed (hicexplorer._hicx)   the body only calls into the extension
                                   (cpp/python/API_BINDINGS.md): every
                                   statement is a call, or a return of a call,
                                   of a function of hicexplorer._hicx, with
                                   names and constants as arguments; a class
                                   is C++-backed when all its methods are
  C++-backed (runs/reads the executable)  a tool's main()/parse_arguments()
  Python driver (hyperopt) around C++ tools
  pure Python                      everything else: the work list

--check lists the pure Python rows that cpp/scripts/public_api_allowlist.tsv
does not allow (matplotlib drawing, which may stay Python, each entry with
its reason) and exits 1 while any remain. It also fails on allowlist entries
that no longer name a pure Python row.
"""
import argparse
import ast
import re
import sys
from pathlib import Path

C_EXTENSION = "C++-backed (hicexplorer._hicx)"
PURE = "pure Python"


def sig(node):
    a = node.args
    parts = []
    pos = a.posonlyargs + a.args
    defaults = [None] * (len(pos) - len(a.defaults)) + list(a.defaults)
    for arg, d in zip(pos, defaults):
        parts.append(arg.arg + ("=" + ast.unparse(d) if d is not None else ""))
    if a.vararg:
        parts.append("*" + a.vararg.arg)
    elif a.kwonlyargs:
        parts.append("*")
    for arg, d in zip(a.kwonlyargs, a.kw_defaults):
        parts.append(arg.arg + ("=" + ast.unparse(d) if d is not None else ""))
    if a.kwarg:
        parts.append("**" + a.kwarg.arg)
    return "(" + ", ".join(parts) + ")"


# ---------------------------------------------------------------------------
# C++-backed bodies

def extension_aliases(tree):
    """Names that refer to hicexplorer._hicx or one of its submodules."""
    aliases = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.level == 0:
            if node.module == "hicexplorer._hicx" or (node.module or "").startswith("hicexplorer._hicx."):
                aliases.update(a.asname or a.name for a in node.names)
            elif node.module == "hicexplorer":
                aliases.update(a.asname or a.name for a in node.names if a.name == "_hicx")
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name == "hicexplorer._hicx" or a.name.startswith("hicexplorer._hicx."):
                    aliases.add(a.asname or "hicexplorer")
    return aliases


def _root_name(func):
    while isinstance(func, ast.Attribute):
        func = func.value
    return func.id if isinstance(func, ast.Name) else None


def _extension_call(node, aliases):
    if not isinstance(node, ast.Call) or _root_name(node.func) not in aliases:
        return False
    if isinstance(node.func, ast.Name):
        return False  # a bare alias is a module, not a function
    values = list(node.args) + [k.value for k in node.keywords]
    for value in values:
        if isinstance(value, ast.Starred):
            value = value.value
        if isinstance(value, (ast.Name, ast.Constant)):
            continue
        if _extension_call(value, aliases):
            continue
        return False
    return True


def _body_without_docstring(node):
    body = list(node.body)
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
            and isinstance(body[0].value.value, str):
        body = body[1:]
    return body


def calls_only_extension(node, aliases):
    body = _body_without_docstring(node)
    if not body or not aliases:
        return False
    for statement in body:
        if isinstance(statement, (ast.Return, ast.Expr)) and statement.value is not None:
            if not _extension_call(statement.value, aliases):
                return False
        else:
            return False
    return True


def _returns_self(node):
    body = _body_without_docstring(node)
    return len(body) == 1 and isinstance(body[0], ast.Return) and isinstance(body[0].value, ast.Name) \
        and body[0].value.id == "self"


def class_is_extension_backed(node, aliases):
    for item in _body_without_docstring(node):
        if not isinstance(item, ast.FunctionDef):
            return False
        if not (calls_only_extension(item, aliases) or _returns_self(item)):
            return False
    return True


# ---------------------------------------------------------------------------
# inventory

def inventory(repo):
    pkg = repo / "hicexplorer"
    cpp_files = [p for p in (repo / "cpp").rglob("*") if p.suffix in (".cpp", ".hpp") and "build" not in p.parts]
    cpp_text = {p: p.read_text(errors="replace") for p in cpp_files}

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
        text = path.read_text()
        tree = ast.parse(text)
        aliases = extension_aliases(tree)
        mod = ".".join(path.relative_to(repo).with_suffix("").parts).replace(".__init__", "")
        tool_module = (repo / "cpp" / "tools" / (path.stem + ".cpp")).is_file()
        rows.append([mod, "module", "", str(len(text.splitlines())), "",
                     "cpp/tools/%s.cpp" % path.stem if tool_module else ""])
        for node in tree.body:
            if not isinstance(node, (ast.FunctionDef, ast.ClassDef)) or node.name.startswith("_"):
                continue
            n = node.end_lineno - node.lineno + 1
            src = ast.get_source_segment(text, node)
            if isinstance(node, ast.ClassDef):
                kind = C_EXTENSION if class_is_extension_backed(node, aliases) else PURE
                rows.append([mod, "class", node.name, str(n), kind,
                             "" if kind == C_EXTENSION else ";".join(cpp_refs(path, node.name))])
                for sub in node.body:
                    if isinstance(sub, ast.FunctionDef) and (not sub.name.startswith("_") or sub.name == "__init__"):
                        sub_kind = C_EXTENSION if calls_only_extension(sub, aliases) else PURE
                        rows.append([mod, "method", node.name + "." + sub.name + sig(sub),
                                     str(sub.end_lineno - sub.lineno + 1), sub_kind,
                                     "" if sub_kind == C_EXTENSION else ";".join(cpp_refs(path, sub.name))])
                continue
            thin = mod in ("hicexplorer._cpp",) or (n <= 6 and ("build_parser(" in src or "entry_point(" in src))
            if thin:
                kind = "C++-backed (runs/reads the executable)"
            elif calls_only_extension(node, aliases):
                kind = C_EXTENSION
            else:
                kind = PURE
            if mod.startswith(("hicexplorer.hicHyperopt",)):
                kind = "Python driver (hyperopt) around C++ tools"
            refs = "" if kind != PURE else ";".join(cpp_refs(path, node.name))
            rows.append([mod, "def", node.name + sig(node), str(n), kind, refs])
    return rows


def read_allowlist(path):
    allowed = {}
    if not path.is_file():
        return allowed
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 3 or not fields[2].strip():
            raise SystemExit("{}:{}: expected module<TAB>name<TAB>reason".format(path, number))
        allowed[(fields[0], fields[1])] = fields[2]
    return allowed


def row_key(row):
    return row[0], row[2].split("(")[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repo", type=Path)
    parser.add_argument("--output", type=Path, help="write the TSV here instead of stdout")
    parser.add_argument("--check", action="store_true", help="list the pure Python names left; exit 1 if any")
    parser.add_argument("--allowlist", type=Path, help="default cpp/scripts/public_api_allowlist.tsv")
    args = parser.parse_args()
    rows = inventory(args.repo)
    if args.check:
        allowlist = args.allowlist or args.repo / "cpp" / "scripts" / "public_api_allowlist.tsv"
        allowed = read_allowlist(allowlist)
        pure = [r for r in rows if r[4] == PURE]
        remaining = [r for r in pure if row_key(r) not in allowed]
        stale = sorted(set(allowed) - {row_key(r) for r in pure})
        for r in remaining:
            print("{}\t{}\t{}\t{} lines".format(r[0], r[1], r[2], r[3]))
        for module, name in stale:
            print("allowlist entry {}.{} names no pure Python row".format(module, name), file=sys.stderr)
        backed = sum(1 for r in rows if r[4] == C_EXTENSION)
        print("{} pure Python names remain ({} allowed as drawing code, {} C++-backed through "
              "hicexplorer._hicx)".format(len(remaining), len(pure) - len(remaining), backed), file=sys.stderr)
        return 1 if remaining or stale else 0
    out = open(args.output, "w") if args.output else sys.stdout
    out.write("module\tkind\tname(signature)\tlines\timplementation\tC++ code that names it (candidate backing)\n")
    for r in rows:
        out.write("\t".join(r) + "\n")
    if args.output:
        out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
