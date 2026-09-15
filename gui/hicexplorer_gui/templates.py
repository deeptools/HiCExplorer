"""Workflow templates (cpp/PLAN.md 10.7).

A template file (``templates/<name>.yaml``) holds a version 1 workflow
(``workflow/model.py``) and what the GUI needs to fill it in::

    template:
      name: hic                      # file and default workflow name
      title: Hi-C
      description: ...
      parameters:                    # in the order the form shows them
        - name: r1
          kind: file                 # file, files, directory or value
          description: first mates, name-sorted BAM
        - name: binSize
          kind: value
          default: 100000
          description: bin size in bp
      unavailable:                   # analyses of the template's purpose that
        - tool: differential loops            # the port does not provide: shown
          purpose: loops that differ          # with the catalog's reason (or the
          reason: not implemented, PLAN 9.7   # given one for a tool the catalog
                                              # does not plan), never run
    workflow:
      version: 1
      inputs: {r1: "@{r1}"}
      steps: [...]

``@{name}`` is replaced by the parameter's value: a string that is exactly one
token takes the value with its type (a list is spliced into a list), a token
inside a longer string its text. A parameter without a default must be given.
"""

import copy
import os
import re

import yaml

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
KINDS = ("file", "files", "directory", "value")
_TOKEN = re.compile(r"@\{([A-Za-z_][A-Za-z0-9_]*)\}")


class TemplateError(Exception):
    pass


class Template:
    def __init__(self, path):
        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
        meta = data.get("template") or {}
        if "workflow" not in data or not meta.get("name"):
            raise TemplateError("{} is not a workflow template (template.name and workflow)".format(path))
        self.path = path
        self.name = str(meta["name"])
        self.title = str(meta.get("title") or self.name)
        self.description = str(meta.get("description") or "")
        self.parameters = []
        for item in meta.get("parameters") or []:
            if item.get("kind", "value") not in KINDS:
                raise TemplateError("{}: parameter {} has an unknown kind {!r}".format(
                    path, item.get("name"), item.get("kind")))
            self.parameters.append(dict(item, kind=item.get("kind", "value")))
        self.unavailable = list(meta.get("unavailable") or [])
        self.workflow = data["workflow"]
        names = {p["name"] for p in self.parameters}
        used = set(_TOKEN.findall(yaml.safe_dump(self.workflow)))
        if used - names:
            raise TemplateError("{}: undefined parameters {}".format(path, ", ".join(sorted(used - names))))

    def tools(self):
        return [step["tool"] for step in self.workflow.get("steps") or []]

    def unavailable_steps(self, entries):
        """(tool, purpose, reason) of the declared unported analyses and of
        every step whose tool the configured tool directory does not provide."""
        by_name = {entry.name: entry for entry in entries}
        out = []
        for item in self.unavailable:
            entry = by_name.get(item["tool"])
            if entry is None:
                reason = item.get("reason") or "not a planned tool"
            elif not entry.available:
                reason = entry.reason
            else:
                reason = "provided by the tool directory, but not a step of this template yet"
            out.append((item["tool"], item.get("purpose", ""), reason))
        for step in self.workflow.get("steps") or []:
            entry = by_name.get(step["tool"])
            if entry is not None and not entry.available:
                out.append((step["tool"], step.get("description", step["id"]), entry.reason))
        return out

    def defaults(self):
        return {p["name"]: p["default"] for p in self.parameters if "default" in p}

    def instantiate(self, values=None, name=None):
        """The workflow document with every parameter replaced."""
        values = dict(values or {})
        known = {p["name"] for p in self.parameters}
        unknown = set(values) - known
        if unknown:
            raise TemplateError("unknown parameters: {}".format(", ".join(sorted(unknown))))
        resolved = {}
        for parameter in self.parameters:
            value = values.get(parameter["name"])
            if value is None or value == "" or value == []:
                if "default" not in parameter:
                    raise TemplateError("parameter {} needs a value".format(parameter["name"]))
                value = parameter["default"]
            resolved[parameter["name"]] = value

        def substitute(node):
            if isinstance(node, dict):
                return {key: substitute(item) for key, item in node.items()}
            if isinstance(node, list):
                out = []
                for item in node:
                    match = _TOKEN.fullmatch(item) if isinstance(item, str) else None
                    if match and isinstance(resolved[match.group(1)], list):
                        out.extend(copy.deepcopy(resolved[match.group(1)]))
                    else:
                        out.append(substitute(item))
                return out
            if isinstance(node, str):
                match = _TOKEN.fullmatch(node)
                if match:
                    return copy.deepcopy(resolved[match.group(1)])
                return _TOKEN.sub(lambda m: str(resolved[m.group(1)]), node)
            return node

        workflow = substitute(copy.deepcopy(self.workflow))
        workflow["name"] = name or workflow.get("name") or self.name
        workflow.setdefault("description", self.description)
        return workflow


def load_templates(directory=TEMPLATE_DIR):
    return [Template(os.path.join(directory, name)) for name in sorted(os.listdir(directory))
            if name.endswith((".yaml", ".yml"))]
