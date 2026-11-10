"""The reference side of the API equivalence tests: finds the 3.7.6
reference interpreter and tree, runs reference_runner.py there and caches
its results.

Environment:
    HICX_REFERENCE_PYTHON     interpreter with the 3.7.6 dependencies (default
                              $HICX_REFERENCE_ENV/bin/python3 if that is set,
                              else
                              ~/miniconda3/envs/__hicexplorer@3.7.6/bin/python3)
    HICX_REFERENCE_ENV        conda prefix of the pinned reference environment
                              (cpp/reference_env/README.md)
    HICX_REFERENCE_ROOT       the reference tree (default: made by
                              cpp/scripts/reference_tree.py of this checkout)
    HICX_API_EQUIV_CACHE      use (default), refresh or off
    HICX_API_EQUIV_CACHE_DIR  default ~/.cache/hicexplorer-equiv/api
"""
import hashlib
import json
import os
import pickle
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
DEFAULT_PYTHON = os.path.expanduser('~/miniconda3/envs/__hicexplorer@3.7.6/bin/python3')


class Unavailable(Exception):
    pass


def reference_python():
    python = os.environ.get('HICX_REFERENCE_PYTHON')
    if not python:
        prefix = os.environ.get('HICX_REFERENCE_ENV')
        python = os.path.join(os.path.expanduser(prefix), 'bin', 'python3') if prefix else DEFAULT_PYTHON
    if not os.path.isfile(python):
        raise Unavailable('no reference interpreter at {} (set HICX_REFERENCE_PYTHON)'.format(python))
    return python


def reference_root(python):
    root = os.environ.get('HICX_REFERENCE_ROOT')
    if root:
        return root
    script = os.path.join(REPO, 'cpp', 'scripts', 'reference_tree.py')
    if not os.path.isfile(script):
        raise Unavailable('set HICX_REFERENCE_ROOT: {} not found'.format(script))
    out = subprocess.run([python, script], capture_output=True, text=True)
    if out.returncode != 0:
        raise Unavailable('reference_tree.py failed: ' + out.stderr[-2000:])
    return out.stdout.split('\t')[0].strip()


def _cache_key(python, root, path):
    digest = hashlib.sha256()
    for part in (python, os.path.realpath(root)):
        digest.update(part.encode() + b'\0')
    meta = os.path.join(root, 'REFERENCE.json')
    if os.path.isfile(meta):
        with open(meta, 'rb') as handle:
            digest.update(handle.read())
    for source in (os.path.join(HERE, 'apieq.py'), os.path.join(HERE, 'reference_runner.py'), path):
        with open(source, 'rb') as handle:
            digest.update(handle.read())
    return digest.hexdigest()[:20]


def results_for(paths):
    """{(cases file name, case id): record} for every case of `paths`."""
    python = reference_python()
    root = reference_root(python)
    mode = os.environ.get('HICX_API_EQUIV_CACHE', 'use')
    cache_dir = os.environ.get('HICX_API_EQUIV_CACHE_DIR',
                               os.path.expanduser('~/.cache/hicexplorer-equiv/api'))
    results = {}
    missing = []
    cache_files = {}
    for path in paths:
        stem = os.path.splitext(os.path.basename(path))[0]
        cache_file = os.path.join(cache_dir, '{}-{}.pkl'.format(stem, _cache_key(python, root, path)))
        cache_files[path] = cache_file
        if mode == 'use' and os.path.isfile(cache_file):
            with open(cache_file, 'rb') as handle:
                results.update(pickle.load(handle))
        else:
            missing.append(path)
    if missing:
        import apieq
        work = tempfile.mkdtemp(prefix='hicx-api-reference-')
        request = {'work': work, 'reference_root': root,
                   'files': [(p, [c.id for c in apieq.load_cases_module(p).CASES]) for p in missing]}
        request_file = os.path.join(work, 'request.pkl')
        result_file = os.path.join(work, 'results.pkl')
        with open(request_file, 'wb') as handle:
            pickle.dump(request, handle, protocol=4)
        env = dict(os.environ)
        env['PYTHONPATH'] = root
        env['HICX_REFERENCE_ROOT'] = root
        out = subprocess.run([python, os.path.join(HERE, 'reference_runner.py'), request_file, result_file],
                             capture_output=True, text=True, env=env, cwd=work)
        if out.returncode != 0:
            raise RuntimeError('reference run failed ({}):\n{}'.format(out.returncode, out.stderr[-6000:]))
        with open(result_file, 'rb') as handle:
            fresh = pickle.load(handle)
        results.update(fresh)
        if mode != 'off':
            os.makedirs(cache_dir, exist_ok=True)
            for path in missing:
                name = os.path.basename(path)
                part = {key: value for key, value in fresh.items() if key[0] == name}
                temporary = cache_files[path] + '.tmp{}'.format(os.getpid())
                with open(temporary, 'wb') as handle:
                    pickle.dump(part, handle, protocol=4)
                os.replace(temporary, cache_files[path])
        import shutil
        shutil.rmtree(work, ignore_errors=True)
    return results


def describe():
    python = reference_python()
    return json.dumps({'python': python, 'root': reference_root(python)})
