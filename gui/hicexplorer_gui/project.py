"""Projects: a directory with workflow files, their work directories and the run history.

Layout::

    <project>/hicexplorer-project.yaml   name and creation date
    <project>/workflows/<name>.yaml      workflows made in the editor
    <project>/work/<name>/               work directory of a workflow (engine state inside)
    <project>/history/<id>/              one entry per run: summary.json, workflow.yaml and
                                         per step run.json, stdout.txt and stderr.txt

Single tool runs use the project directory as their work directory, so their
outputs lie inside the project.
"""

import datetime
import json
import os
import re

import yaml

PROJECT_FILE = "hicexplorer-project.yaml"


class ProjectError(Exception):
    pass


class Project:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        meta = os.path.join(self.path, PROJECT_FILE)
        if not os.path.isfile(meta):
            raise ProjectError("{} is not a HiCExplorer project (no {})".format(self.path, PROJECT_FILE))
        with open(meta) as handle:
            data = yaml.safe_load(handle) or {}
        self.name = str(data.get("name") or os.path.basename(self.path))

    @classmethod
    def create(cls, path, name=None):
        path = os.path.abspath(path)
        os.makedirs(path, exist_ok=True)
        meta = os.path.join(path, PROJECT_FILE)
        if not os.path.exists(meta):
            with open(meta, "w") as handle:
                yaml.safe_dump({"format": "hicexplorer-project", "version": 1,
                                "name": name or os.path.basename(path),
                                "created": datetime.date.today().isoformat()}, handle, sort_keys=False)
        for sub in ("workflows", "work", "history"):
            os.makedirs(os.path.join(path, sub), exist_ok=True)
        return cls(path)

    @property
    def workflows_dir(self):
        return os.path.join(self.path, "workflows")

    @property
    def history_dir(self):
        return os.path.join(self.path, "history")

    def workdir_for(self, workflow_name):
        path = os.path.join(self.path, "work", safe_name(workflow_name))
        os.makedirs(path, exist_ok=True)
        return path

    def workflows(self):
        if not os.path.isdir(self.workflows_dir):
            return []
        return sorted(os.path.join(self.workflows_dir, f) for f in os.listdir(self.workflows_dir)
                      if f.endswith((".yaml", ".yml")))

    def new_history_entry(self, name):
        stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S-%f")
        path = os.path.join(self.history_dir, "{}_{}".format(stamp, safe_name(name)))
        os.makedirs(path)
        return path

    def history(self):
        """Summaries of past runs, newest first."""
        out = []
        if not os.path.isdir(self.history_dir):
            return out
        for entry in sorted(os.listdir(self.history_dir), reverse=True):
            summary = os.path.join(self.history_dir, entry, "summary.json")
            if os.path.isfile(summary):
                try:
                    with open(summary) as handle:
                        data = json.load(handle)
                except (OSError, ValueError):
                    continue
                data["dir"] = os.path.join(self.history_dir, entry)
                out.append(data)
        return out


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "run"
