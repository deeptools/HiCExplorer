"""Entry point of the HiCExplorer GUI: ``hicexplorer-gui [--project DIR] [--tools-dir DIR]``."""

import argparse
import sys


def build_parser():
    parser = argparse.ArgumentParser(prog="hicexplorer-gui", description="HiCExplorer graphical interface.")
    parser.add_argument("--project", help="open this project directory (created when missing)")
    parser.add_argument("--tools-dir", help="directory of the C++ tool executables (stored in the settings)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    from PySide6 import QtWidgets

    from .main_window import MainWindow
    from .project import PROJECT_FILE
    from .settings import Settings

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    app.setApplicationName("hicexplorer-gui")
    settings = Settings()
    if args.tools_dir:
        settings.tools_dir = args.tools_dir
    window = MainWindow(settings)
    if args.project:
        import os
        if os.path.isfile(os.path.join(args.project, PROJECT_FILE)):
            window.open_project(args.project)
        else:
            window.create_project(args.project)
    window.resize(1280, 720)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
