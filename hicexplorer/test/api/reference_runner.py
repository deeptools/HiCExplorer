"""Runs equivalence cases on the HiCExplorer 3.7.6 reference tree.

    PYTHONPATH=<reference root> <reference python> reference_runner.py request.pkl results.pkl

request.pkl holds {'work': directory, 'reference_root': path,
'files': [(cases file, [case id, ...]), ...]}; results.pkl receives
{(cases file name, case id): record} with the records of apieq.run_case.
Started by reference.py, never by hand in the package's own interpreter:
the hicexplorer it imports must be the reference tree's.
"""
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import apieq  # noqa: E402


def main():
    with open(sys.argv[1], 'rb') as handle:
        request = pickle.load(handle)
    import hicexplorer
    root = os.path.realpath(request['reference_root'])
    if not os.path.realpath(hicexplorer.__file__).startswith(root + os.sep):
        sys.stderr.write('hicexplorer is {}, not the reference tree {}\n'.format(hicexplorer.__file__, root))
        return 2
    results = {}
    for path, ids in request['files']:
        module = apieq.load_cases_module(path)
        stem = os.path.splitext(os.path.basename(path))[0]
        wanted = set(ids)
        for case in module.CASES:
            if case.id in wanted:
                directory = os.path.join(request['work'], stem, apieq.safe_name(case.id))
                results[(os.path.basename(path), case.id)] = apieq.run_case(module, case, directory)
    with open(sys.argv[2], 'wb') as handle:
        pickle.dump(results, handle, protocol=4)
    return 0


if __name__ == '__main__':
    sys.exit(main())
