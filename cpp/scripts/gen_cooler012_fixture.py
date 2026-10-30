#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cooler012/{ordinary,int64,empty,atomic}.txt, the oracle of the three
behaviours of cooler 0.12.0/0.12.1 that coolercpp adopts (the rest of
coolercpp keeps reproducing cooler 0.10.2):

  1. chromosomes longer than 2147483647 bp: chroms/length, bins/start and
     bins/end are int64 in the whole file (int32 otherwise),
  2. coolers without pixels: the pixel datasets have no placeholder rows,
  3. coolers are written into a temporary file (mode "w") or a temporary group
     (mode "a") and moved into place when complete; the final file has the
     same groups and datasets as before.

Every case is a list of cooler.create_cooler calls (the `create` lines) run
into one file, followed by a structural dump of that file. The dump is written
in a form that test_cooler012.cpp reproduces from an HDF5 file written by
coolercpp, so that the two texts are compared character by character:

    <python with cooler 0.12.1> cpp/scripts/gen_cooler012_fixture.py

Reference: cooler 0.12.1 (the version is asserted). Format of a case:

  case <name>
  create <file>[::<group>] mode=<w|a> input=<table|ordered|unordered>
         chroms=<name:length,...> binsize=<int> chunks=<n> pixels=<b1:b2:count,...>
  ...
  dump
  G <path> attrs=<name>=<value>|...        (a group; i<int> f<bits> s<text>)
  D <path> dtype=<> shape=<n> chunks=<n|-> maxshape=<n|inf> filters=<id:cd,..;..> values=<...>
  end

  input=table      one DataFrame (create_cooler sorts it)
  input=ordered    an iterator of `chunks` DataFrames, ordered=True
  input=unordered  the same with ordered=False (external merge sort)
  chunks=0         the iterator yields nothing at all
The creation-date and generated-by attributes are left out of the dump.
"""
import os
import sys
import struct
import subprocess
import tempfile

import cooler
import h5py
import numpy as np
import pandas as pd

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests',
                       'fixtures', 'cooler012')

ORDINARY = '0:0:5,0:2:3,1:1:4,2:5:1,3:3:2,4:9:7,10:12:1,13:14:6'
BIG_PIXELS = '0:0:1,0:2:2,1:3:3,5:200:4'

ORDINARY_CASES = [
    # ordinary cases (unchanged by 0.12)
    ('ordinary_table', [('main.cool', 'w', 'table', 'chrA:1000,chrB:500', 100, 1, ORDINARY)]),
    ('ordinary_ordered_chunks', [('main.cool', 'w', 'ordered', 'chrA:1000,chrB:500', 100, 3, ORDINARY)]),
    ('ordinary_unordered_chunks', [('main.cool', 'w', 'unordered', 'chrA:1000,chrB:500', 100, 3,
                                    '13:14:6,4:9:7,0:0:5,10:12:1,2:5:1,0:2:3,3:3:2,1:1:4')]),
    ('ordinary_group_w', [('main.cool::/x/y', 'w', 'table', 'chrA:1000,chrB:500', 100, 1, ORDINARY)]),
    ('mcool_ordinary', [
        ('main.cool::/resolutions/100', 'a', 'ordered', 'chrA:1000,chrB:500', 100, 2, ORDINARY),
        ('main.cool::/resolutions/500', 'a', 'table', 'chrA:1000,chrB:500', 500, 1, '0:0:1,0:2:9,1:2:4'),
        ('main.cool::/resolutions/100', 'a', 'table', 'chrA:1000,chrB:500', 100, 1, '0:1:2,3:4:5'),
    ]),
    ('root_append', [
        ('main.cool', 'a', 'table', 'chrA:1000,chrB:500', 100, 1, ORDINARY),
        ('main.cool', 'a', 'table', 'chrA:1000,chrB:500', 250, 1, '0:0:1,1:3:2'),
        ('main.cool::/other', 'a', 'table', 'chrA:1000,chrB:500', 100, 1, '0:0:1'),
    ]),
]

INT64_CASES = [
    # (1) chromosomes beyond the int32 maximum
    ('big_binsize_1e7', [('main.cool', 'w', 'table', 'chrA:3000000000,chrB:1000000', 10000000, 1, BIG_PIXELS)]),
    ('big_binsize_1e6', [('main.cool', 'w', 'table', 'chrA:3000000000,chrB:1000000', 1000000, 1, BIG_PIXELS)]),
    ('big_second_chrom', [('main.cool', 'w', 'ordered', 'chrA:5000,chrB:2200000000', 1000000000, 2, '0:0:1,0:3:2,2:3:5')]),
    ('big_unordered', [('main.cool', 'w', 'unordered', 'chrA:3000000000', 500000000, 2, '3:5:1,0:1:2,0:0:1')]),
    ('big_mcool', [
        ('main.cool::/resolutions/10000000', 'a', 'table', 'chrA:3000000000,chrB:1000000', 10000000, 1, BIG_PIXELS),
        ('main.cool::/resolutions/20000000', 'a', 'table', 'chrA:3000000000,chrB:1000000', 20000000, 1, '0:0:1,2:100:3'),
    ]),
    ('int32_maximum_exact', [('main.cool', 'w', 'table', 'chrA:2147483647', 1073741824, 1, '0:0:1,0:1:2')]),
    ('int32_maximum_plus_one', [('main.cool', 'w', 'table', 'chrA:2147483648', 1073741824, 1, '0:0:1,0:1:2')]),
    ('mixed_mcool_small_then_big', [
        ('main.cool::/resolutions/1', 'a', 'table', 'chrA:1000', 100, 1, '0:0:1'),
        ('main.cool::/resolutions/2', 'a', 'table', 'chrA:2147483649', 1073741825, 1, '0:0:1'),
    ]),
]

EMPTY_CASES = [
    # (2) coolers without pixels
    ('empty_table', [('main.cool', 'w', 'table', 'chrA:1000,chrB:500', 100, 1, '')]),
    ('empty_ordered_no_chunk', [('main.cool', 'w', 'ordered', 'chrA:1000,chrB:500', 100, 0, '')]),
    ('empty_ordered_empty_chunk', [('main.cool', 'w', 'ordered', 'chrA:1000,chrB:500', 100, 1, '')]),
    ('empty_unordered_no_chunk', [('main.cool', 'w', 'unordered', 'chrA:1000,chrB:500', 100, 0, '')]),
    ('empty_unordered_empty_chunk', [('main.cool', 'w', 'unordered', 'chrA:1000,chrB:500', 100, 1, '')]),
    ('empty_mcool', [
        ('main.cool::/resolutions/100', 'a', 'ordered', 'chrA:1000,chrB:500', 100, 0, ''),
        ('main.cool::/resolutions/500', 'a', 'table', 'chrA:1000,chrB:500', 500, 1, ''),
    ]),
    ('empty_big', [('main.cool', 'w', 'ordered', 'chrA:3000000000', 1000000000, 0, '')]),
]


def parse_pixels(text):
    rows = [tuple(int(v) for v in item.split(':')) for item in text.split(',') if item]
    return pd.DataFrame({
        'bin1_id': np.array([r[0] for r in rows], dtype=np.int64),
        'bin2_id': np.array([r[1] for r in rows], dtype=np.int64),
        'count': np.array([r[2] for r in rows], dtype=np.int64),
    })


def parse_chroms(text):
    items = [item.split(':') for item in text.split(',')]
    return pd.Series({name: int(length) for name, length in items})


def run_create(directory, call):
    uri, mode, kind, chroms, binsize, chunks, pixels = call
    bins = cooler.util.binnify(parse_chroms(chroms), binsize)
    frame = parse_pixels(pixels)
    target = os.path.join(directory, uri)
    if kind == 'table':
        cooler.create_cooler(target, bins, frame, mode=mode)
        return
    edges = [len(frame) * i // chunks for i in range(chunks + 1)] if chunks else []
    parts = [frame.iloc[edges[i]:edges[i + 1]] for i in range(chunks)]
    cooler.create_cooler(target, bins, iter(parts), mode=mode, ordered=(kind == 'ordered'))


def fmt_attr(value):
    if isinstance(value, (bool, np.bool_)):
        return 'b%d' % int(value)
    if isinstance(value, (int, np.integer)):
        return 'i%d' % int(value)
    if isinstance(value, (float, np.floating)):
        return 'f%016x' % struct.unpack('<Q', struct.pack('<d', float(value)))[0]
    if isinstance(value, bytes):
        value = value.decode()
    return 's%s' % value


def dtype_text(ds):
    members = h5py.check_enum_dtype(ds.dtype)
    if members is not None:
        base = ds.dtype
        items = ','.join('%s=%d' % (k, v) for k, v in sorted(members.items(), key=lambda kv: kv[1]))
        return 'enum:%s%d{%s}' % (base.kind, base.itemsize, items)
    return '%s%d' % (ds.dtype.kind, ds.dtype.itemsize)


def values_text(ds):
    data = ds[()] if ds.shape else np.array([])
    n = len(data)
    if ds.dtype.kind == 'S':
        items = [bytes(v).decode() for v in data]
    elif ds.dtype.kind == 'f':
        items = ['%016x' % struct.unpack('<Q', struct.pack('<d', float(v)))[0] for v in data]
    else:
        items = [str(int(v)) for v in data]
    shown = items if n <= 40 else items[:32] + ['...'] + items[-8:]
    text = ','.join(shown)
    if ds.dtype.kind in 'iu' and n:
        text += ' sum=%d' % int(np.sum(data.astype(np.int64)))
    return text


def dump(path):
    lines = []
    with h5py.File(path, 'r', locking=False) as f:
        def visit(name, obj):
            full = '/' + name
            if isinstance(obj, h5py.Group):
                attrs = '|'.join('%s=%s' % (k, fmt_attr(obj.attrs[k]))
                                 for k in sorted(obj.attrs)
                                 if k not in ('creation-date', 'generated-by'))
                lines.append('G %s links=%s attrs=%s' % (full, ','.join(obj.keys()), attrs))
            else:
                plist = obj.id.get_create_plist()
                filters = ';'.join(
                    '%d:%s' % (plist.get_filter(i)[0], ','.join(str(v) for v in plist.get_filter(i)[2]))
                    for i in range(plist.get_nfilters()))
                maxshape = obj.maxshape[0]
                lines.append('D %s dtype=%s shape=%d chunks=%s maxshape=%s filters=%s values=%s' % (
                    full, dtype_text(obj), obj.shape[0],
                    obj.chunks[0] if obj.chunks else '-',
                    'inf' if maxshape is None else maxshape, filters, values_text(obj)))
        root_attrs = '|'.join('%s=%s' % (k, fmt_attr(f.attrs[k])) for k in sorted(f.attrs)
                              if k not in ('creation-date', 'generated-by'))
        lines.append('G / links=%s attrs=%s' % (','.join(f.keys()), root_attrs))
        f.visititems(visit)
    return lines


def listing(directory):
    return ','.join(sorted(os.listdir(directory)))


def group_paths(path):
    """Paths of all groups of the file, by a separate process, which is how a
    file that is being written can be looked at."""
    code = ("import sys, h5py\n"
            "f = h5py.File(sys.argv[1], 'r', locking=False)\n"
            "out = ['/']\n"
            "f.visititems(lambda n, o: out.append('/' + n) if isinstance(o, h5py.Group) else None)\n"
            "print(','.join(out))\n")
    result = subprocess.run([sys.executable, '-c', code, path], capture_output=True, text=True)
    return result.stdout.strip() or 'unreadable'


def read_bytes(path):
    with open(path, 'rb') as f:
        return f.read()


def scenario(name, directory, mode, uri, fail):
    """One create_cooler call of two chunks, the second invalid when `fail`,
    into a directory that holds a cooler already (a plain file for mode w, the
    groups /resolutions/100 and /resolutions/500 for mode a). Records what is
    visible while the call is in progress and when it is over."""
    lines = ['scenario %s' % name]
    bins = cooler.util.binnify(parse_chroms('chrA:1000,chrB:500'), 100)
    main = os.path.join(directory, 'main.cool')
    if mode == 'w':
        cooler.create_cooler(main, bins, parse_pixels('0:0:5,1:2:3'), mode='w')
    else:
        for res in (100, 500):
            cooler.create_cooler(main + '::/resolutions/%d' % res, cooler.util.binnify(
                parse_chroms('chrA:1000,chrB:500'), res), parse_pixels('0:1:1'), mode='a')
    before = read_bytes(main)
    groups_before = group_paths(main)

    def chunks():
        yield parse_pixels('0:0:1')
        if mode == 'w':
            lines.append('mid dir=%s unchanged=%d' % (listing(directory), read_bytes(main) == before))
        else:
            lines.append('mid groups=%s' % group_paths(main))
        yield parse_pixels('1:99:2' if fail else '1:2:2')

    target = main if uri == '' else main + '::' + uri
    error = 0
    try:
        cooler.create_cooler(target, bins, chunks(), mode=mode, ordered=True)
    except Exception:
        error = 1
    lines.append('error=%d dir=%s' % (error, listing(directory)))
    if mode == 'w':
        lines.append('unchanged=%d' % (read_bytes(main) == before))
    lines.append('groups=%s' % group_paths(main))
    lines.append('groups_before=%s' % groups_before)
    lines.append('end')
    return lines


SCENARIOS = [
    ('w_success', 'w', '', 0),
    ('w_failure', 'w', '', 1),
    ('w_group_failure', 'w', '/x/y', 1),
    ('a_new_group_success', 'a', '/resolutions/200', 0),
    ('a_new_group_failure', 'a', '/resolutions/200', 1),
    ('a_replace_group_success', 'a', '/resolutions/100', 0),
    ('a_replace_group_failure', 'a', '/resolutions/100', 1),
    ('a_nested_new_group_failure', 'a', '/other/deep/new', 1),
    ('a_root_success', 'a', '', 0),
    ('a_root_failure', 'a', '', 1),
]


def write_cases(name, cases):
    out = []
    for case, calls in cases:
        out.append('case %s' % case)
        with tempfile.TemporaryDirectory() as directory:
            for call in calls:
                out.append('create %s mode=%s input=%s chroms=%s binsize=%d chunks=%d pixels=%s' % call)
                run_create(directory, call)
            out.append('dump')
            out.extend(dump(os.path.join(directory, calls[0][0].split('::')[0])))
        out.append('end')
    path = os.path.join(OUT_DIR, name)
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %s (%d cases)' % (os.path.normpath(path), len(cases)))


def main():
    assert cooler.__version__ == '0.12.1', cooler.__version__
    os.makedirs(OUT_DIR, exist_ok=True)
    write_cases('ordinary.txt', ORDINARY_CASES)
    write_cases('int64.txt', INT64_CASES)
    write_cases('empty.txt', EMPTY_CASES)
    out = []
    for name, mode, uri, fail in SCENARIOS:
        with tempfile.TemporaryDirectory() as directory:
            out.extend(scenario(name, directory, mode, uri, fail))
    path = os.path.join(OUT_DIR, 'atomic.txt')
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %s (%d scenarios)' % (os.path.normpath(path), len(SCENARIOS)))


if __name__ == '__main__':
    main()
