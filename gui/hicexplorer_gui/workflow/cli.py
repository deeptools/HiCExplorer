"""``hicexplorer-workflow``: run, validate or export a HiCExplorer workflow file.

Exit status: 0 success, 1 a step or the validation failed, 2 usage error,
130 cancelled.
"""

import argparse
import os
import signal
import sys

from .model import WorkflowError, load_workflow
from .runner import EXIT_CANCELLED, EXIT_FAILED, EXIT_OK, EXIT_USAGE, Runner
from .spec import SpecLoader
from .validate import validate_workflow


def _positive_int(text):
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("expected a positive integer, found {!r}".format(text))
    if value < 1:
        raise argparse.ArgumentTypeError("expected a positive integer, found {!r}".format(text))
    return value


def build_parser():
    parser = argparse.ArgumentParser(
        prog="hicexplorer-workflow",
        description="Run, validate or export a HiCExplorer workflow (a versioned YAML DAG of "
                    "C++ tool steps). Steps are resumed from "
                    "<workdir>/.hicexplorer-workflow when their tool, parameters and input "
                    "contents are unchanged.",
        epilog="Exit status: 0 success, 1 a step or the validation failed, 2 usage error, "
               "130 cancelled.")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    def common(p):
        p.add_argument("workflow", metavar="FILE.yaml", help="Workflow file.")
        p.add_argument("--tools-dir", metavar="DIR",
                       help="Directory with the C++ tool executables "
                            "(default: $HICX_CPP_BIN, else PATH).")
        p.add_argument("--workdir", metavar="DIR",
                       help="Working directory for outputs and state "
                            "(default: the workflow file's directory).")

    p_run = sub.add_parser("run", help="Run the workflow.")
    common(p_run)
    p_run.add_argument("--threads", type=_positive_int, metavar="N",
                       help="Thread budget, overriding the file's 'threads'.")
    p_run.add_argument("--force", action="store_true", help="Rerun every step.")

    p_val = sub.add_parser("validate", help="Check the workflow against the tool specifications.")
    common(p_val)

    p_exp = sub.add_parser("export", help="Write the workflow as a shell script or a Snakefile.")
    common(p_exp)
    p_exp.add_argument("--format", choices=["sh", "snakemake"], required=True)
    p_exp.add_argument("-o", "--output", required=True, metavar="OUT", help="Output file.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not os.path.isfile(args.workflow):
        print("hicexplorer-workflow: error: workflow file {} not found".format(args.workflow), file=sys.stderr)
        return EXIT_USAGE
    if args.workdir is not None and os.path.exists(args.workdir) and not os.path.isdir(args.workdir):
        print("hicexplorer-workflow: error: --workdir {} is not a directory".format(args.workdir), file=sys.stderr)
        return EXIT_USAGE
    if args.tools_dir is not None and not os.path.isdir(args.tools_dir):
        print("hicexplorer-workflow: error: --tools-dir {} is not a directory".format(args.tools_dir),
              file=sys.stderr)
        return EXIT_USAGE
    if args.workdir is not None:
        os.makedirs(args.workdir, exist_ok=True)
    try:
        workflow = load_workflow(args.workflow, args.workdir)
    except WorkflowError as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return EXIT_FAILED
    loader = SpecLoader(args.tools_dir)

    if args.command == "run":
        runner = Runner(workflow, loader, threads=args.threads, force=args.force)
        previous = {}

        def handler(signum, _frame):
            runner.cancel()

        for sig in (signal.SIGINT, signal.SIGTERM):
            previous[sig] = signal.signal(sig, handler)
        try:
            return runner.run()
        finally:
            for sig, old in previous.items():
                signal.signal(sig, old)

    messages = validate_workflow(workflow, loader)
    for msg in messages:
        print(msg, file=sys.stderr)
    if any(m.level == "error" for m in messages):
        return EXIT_FAILED
    if args.command == "validate":
        print("workflow {}: valid ({} steps)".format(workflow.name, len(workflow.steps)))
        return EXIT_OK

    from .export import export_sh, export_snakemake
    text = export_sh(workflow, loader) if args.format == "sh" else export_snakemake(workflow, loader)
    with open(args.output, "w") as handle:
        handle.write(text)
    if args.format == "sh":
        os.chmod(args.output, os.stat(args.output).st_mode | 0o111)
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
