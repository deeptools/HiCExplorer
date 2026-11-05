"""Equivalence cases for hicexplorer.hicPrepareQCreport (apieq.py).

The tables are built as hicPrepareQCreport.main builds them, from the QC logs
of the repository's test data (parsed here with pandas, in both interpreters),
plus synthetic frames for the edge cases: missing and duplicate columns, text
and float cells, zero denominators, labels that are not strings.

save_html: the driver returns the file it wrote, with the five random
hexadecimal characters of each table id replaced. make_figure_*: the driver
returns the table the function returns and the DataFrames it plotted (pandas'
plot accessor is observed), and whether the figure file exists. The pixels
depend on the matplotlib version of the interpreter and are not compared.
make_sure_path_exists: the driver returns the directory tree afterwards.

Deviations from 3.7.6 (hicexplorer/hicPrepareQCreport.py,
cpp/python/src/api/hicPrepareQCreport.cpp):
  * DataFrames with a single level index and columns; labels str, int or
    float; cells int, float or str (others raise TypeError).
  * Duplicate row labels, or duplicate column labels among the selected ones,
    raise ValueError where pandas returned a cartesian join.
  * save_html takes a path, not a file descriptor.
  * A table that is not a DataFrame raises AttributeError.
"""
import os
import re

from apieq import Case as _Case, File


def Case(*args, **kwargs):
    # a DataFrame's class has another module path in each pandas version; the
    # drivers compare the table after the call instead of the argument
    kwargs.setdefault('compare_args', False)
    return _Case(*args, **kwargs)


HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
Q = 'hicexplorer.hicPrepareQCreport.'


# ---------------------------------------------------------------------------
# tables

def log_path(*parts):
    return os.path.join(TEST_DATA, *parts)


def qc_table(logs, labels=None):
    """hicPrepareQCreport.main from reading the logs to the derived columns."""
    import pandas as pd
    params = dict()
    for path in logs:
        in_log_part = False
        with open(path) as fh:
            lines = fh.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("File"):
                in_log_part = True
            if in_log_part is True:
                if line.strip() == "" or line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) == 1:
                    continue
                if fields[0] not in params:
                    params[fields[0]] = []
                try:
                    params[fields[0]].append(int(fields[1]))
                except ValueError:
                    params[fields[0]].append(fields[1])
    table = pd.DataFrame(params)
    if labels and len(labels) == len(logs):
        table['Labels'] = labels
        table = table.set_index('Labels')
    else:
        table = table.set_index('File')
    if 'Pairs mappable, unique and high quality' not in table.columns:
        table['Pairs mappable, unique and high quality'] = \
            table['Sequenced reads'] - (table['One mate unmapped'] +
                                        table['One mate not unique'] + table['Low mapping quality'])
    if 'same fragment (800 bp)' in table.columns:
        table['same fragment'] = table['same fragment (800 bp)']
    return table


def single():
    return qc_table([log_path('QC', 'QC.log')], ['sample'])


def multi():
    return qc_table([log_path('QC', 'QC.log'), log_path('QC_rc', 'QC.log')], ['a', 'b'])


def multi_restriction():
    return qc_table([log_path('QC_rc_multiple', 'QC.log'), log_path('QC_region', 'QC.log'),
                     log_path('QC_multi_restriction', 'QC.log')], ['m', 'n', 'o'])


def no_restriction():
    return qc_table([log_path('QC_no_restriction_site', 'QC.log')], ['plain'])


def quick():
    return qc_table([log_path('hicQuickQC', 'QC.log')])


def file_labels():
    return qc_table([log_path('QC', 'QC.log'), log_path('QC_rc', 'QC.log')])


def columns_of(*names, extra_columns=None, **extra):
    """A synthetic table: ints for every name in `names` (a different value per
    column and row), plus keyword columns given as lists."""
    import pandas as pd
    data = {}
    for i, name in enumerate(names):
        data[name] = [100 * (i + 1) + 7 * k for k in range(3)]
    data.update(extra)
    data.update(extra_columns or {})
    return pd.DataFrame(data, index=pd.Index(['x', 'y', 'z'], name='Labels'))


USED = ['Hi-C contacts', 'Pairs mappable, unique and high quality', 'Sequenced reads']
UNMAP = ['Hi-C contacts', 'Low mapping quality', 'One mate not unique', 'One mate unmapped',
         'Sequenced reads']
DISCARD = ['One mate not close to rest site', 'dangling end', 'duplicated pairs', 'same fragment',
           'self circle', 'self ligation (removed)', 'Pairs mappable, unique and high quality']
DISTANCE = ['inter chromosomal', 'Intra short range (< 20kb)', 'Intra long range (>= 20kb)',
            'Hi-C contacts']
ORIENT = ['Read pair type: inward pairs', 'Read pair type: outward pairs',
          'Read pair type: left pairs', 'Read pair type: right pairs']


# ---------------------------------------------------------------------------
# drivers

def frame_data(df):
    """A DataFrame as plain data."""
    def kind(dtype):
        name = str(dtype)
        return name if name in ('int64', 'float64') else 'other'
    return {'index': df.index.tolist(), 'index_name': df.index.name,
            'columns': list(df.columns), 'columns_name': df.columns.name,
            'dtypes': [kind(t) for t in df.dtypes],
            'data': [df.iloc[:, c].tolist() for c in range(df.shape[1])]}


def figure(target, table, name, *extra):
    """Calls make_figure_*(table, name, *extra): the table it returns, the
    DataFrames it plotted, and whether the figure was written."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas.plotting._core as core

    seen = []
    original = core.PlotAccessor.__call__

    def call(self, *args, **kwargs):
        seen.append(frame_data(self._parent))
        return original(self, *args, **kwargs)

    core.PlotAccessor.__call__ = call
    try:
        result = target(table, name, *extra)
    finally:
        core.PlotAccessor.__call__ = original
        plt.close('all')
    return {'returned': None if result is None else frame_data(result), 'plotted': seen,
            'written': os.path.getsize(name) > 0, 'table_after': frame_data(table)}


def tables_of(table, plot_dir):
    """The four derived tables of main, computed as the 3.7.6 functions do."""
    import pandas as pd
    del plot_dir
    unmap = table[['Hi-C contacts', 'Low mapping quality', 'One mate not unique',
                   'One mate unmapped']].T / table['Sequenced reads']
    unmap = table[['Hi-C contacts', 'Low mapping quality', 'One mate not unique',
                   'One mate unmapped']].join(unmap.T, rsuffix='_%')
    unmap = unmap[[u'Hi-C contacts', u'Hi-C contacts_%', u'Low mapping quality',
                   u'Low mapping quality_%', u'One mate not unique', u'One mate not unique_%',
                   u'One mate unmapped', u'One mate unmapped_%']]
    names = []
    out = []
    for name in ['One mate not close to rest site', 'dangling end', 'duplicated pairs',
                 'same fragment', 'self circle', 'self ligation (removed)']:
        if name in list(table.columns):
            names.append(name)
            out.append(name)
            out.append(name + ' %')
        else:
            for name_ in list(table.columns):
                if name in name_:
                    names.append(name_)
                    out.append(name_ + ' %')
    prc = table[names].T / table['Pairs mappable, unique and high quality']
    discarded = table[names].join(prc.T, rsuffix=' %')[out]
    cols = ['inter chromosomal', 'Intra short range (< 20kb)', 'Intra long range (>= 20kb)']
    prc = table[cols].T / table['Hi-C contacts']
    distance = table[cols].join(prc.T, rsuffix=' %')
    distance = distance[['inter chromosomal', 'inter chromosomal %', 'Intra short range (< 20kb)',
                         'Intra short range (< 20kb) %', 'Intra long range (>= 20kb)',
                         'Intra long range (>= 20kb) %']]
    cols = [u'Read pair type: inward pairs', u'Read pair type: outward pairs',
            u'Read pair type: left pairs', u'Read pair type: right pairs']
    t = table[cols].T
    prc = t / t.sum(axis=0)
    orientation = table[cols].join(prc.T, rsuffix=' %')
    orientation = orientation[[c for c in
                               [x for name in cols for x in (name, name + ' %')]]]
    del pd
    return unmap, discarded, distance, orientation


def html(target, table, name, replace=None):
    """save_html with the four derived tables of `table` and `table` itself;
    returns the file, ids normalised."""
    unmap, discarded, distance, orientation = tables_of(table, None)
    return write_html(target, name, unmap, discarded, distance, orientation, table)


def write_html(target, name, *tables):
    target(name, *tables)
    with open(name) as handle:
        text = handle.read()
    return re.sub(r'T_[0-9a-f]{5}', 'T_ID', text)


def html_frames(target, name, frames):
    """save_html with explicit frames (five of them)."""
    return write_html(target, name, *frames())


def tree(target, base, *parts, **kwargs):
    """Runs target(path) and returns the tree below base afterwards."""
    path = os.path.join(base, *parts) if parts else base
    if kwargs.get('trailing'):
        path += '/'
    if kwargs.get('as_bytes'):
        path = os.fsencode(path)
    if kwargs.get('as_path'):
        import pathlib
        path = pathlib.Path(path)
    if kwargs.get('prepare'):
        prepare(base, kwargs['prepare'])
    try:
        result = target(path)
    finally:
        if kwargs.get('prepare') == 'readonly':
            os.chmod(os.path.join(base, 'ro'), 0o755)
    listing = []
    for root, dirs, files in os.walk(base):
        for d in dirs:
            listing.append(('dir', os.path.relpath(os.path.join(root, d), base)))
        for f in files:
            listing.append(('file', os.path.relpath(os.path.join(root, f), base)))
    return result, sorted(listing)


def prepare(base, what):
    os.makedirs(base, exist_ok=True)
    if what == 'file':
        with open(os.path.join(base, 'f'), 'w') as handle:
            handle.write('x')
    elif what == 'dir':
        os.makedirs(os.path.join(base, 'd'))
    elif what == 'readonly':
        os.makedirs(os.path.join(base, 'ro'))
        os.chmod(os.path.join(base, 'ro'), 0o555)


def plain(target, *args):
    return target(*args)


CASES = []

# --- make_figure_pairs_used -----------------------------------------------------------
for tag, make in [('single', single), ('multi', multi), ('multi_restriction', multi_restriction),
                  ('file_labels', file_labels), ('no_restriction', no_restriction),
                  ('quick', quick)]:
    CASES += [
        Case('pairs_used/' + tag, Q + 'make_figure_pairs_used', driver='figure',
             args=lambda make=make: (make(), File('f.png'), 72)),
        Case('unmappable/' + tag, Q + 'make_figure_umappable_non_unique_reads', driver='figure',
             args=lambda make=make: (make(), File('f.png'), 72)),
        Case('discarded/' + tag, Q + 'make_figure_pairs_discarded', driver='figure',
             args=lambda make=make: (make(), File('f.png'), 72)),
        Case('distance/' + tag, Q + 'make_figure_distance', driver='figure',
             args=lambda make=make: (make(), File('f.png'))),
        Case('orientation/' + tag, Q + 'make_figure_read_orientation', driver='figure',
             args=lambda make=make: (make(), File('f.png'), 72)),
    ]

CASES += [
    # missing columns: pandas' KeyErrors
    Case('pairs_used/missing_one', Q + 'make_figure_pairs_used', driver='figure',
         args=lambda: (columns_of(*USED[:2]), File('f.png'), 72)),
    Case('pairs_used/missing_all', Q + 'make_figure_pairs_used', driver='figure',
         args=lambda: (columns_of('other'), File('f.png'), 72)),
    Case('pairs_used/missing_two', Q + 'make_figure_pairs_used', driver='figure',
         args=lambda: (columns_of(USED[0]), File('f.png'), 72)),
    Case('unmappable/missing_count', Q + 'make_figure_umappable_non_unique_reads', driver='figure',
         args=lambda: (columns_of(*UNMAP[:3]), File('f.png'), 72)),
    Case('unmappable/missing_denominator', Q + 'make_figure_umappable_non_unique_reads',
         driver='figure', args=lambda: (columns_of(*UNMAP[:4]), File('f.png'), 72)),
    Case('unmappable/missing_all', Q + 'make_figure_umappable_non_unique_reads', driver='figure',
         args=lambda: (columns_of('other'), File('f.png'), 72)),
    Case('discarded/missing_denominator', Q + 'make_figure_pairs_discarded', driver='figure',
         args=lambda: (columns_of(*DISCARD[:6]), File('f.png'), 72)),
    Case('discarded/no_discard_columns', Q + 'make_figure_pairs_discarded', driver='figure',
         args=lambda: (columns_of(DISCARD[6]), File('f.png'), 72)),
    Case('discarded/substring_columns', Q + 'make_figure_pairs_discarded', driver='figure',
         args=lambda: (columns_of('dangling end GATC', 'dangling end AGCT', 'same fragment (800 bp)',
                                  'self circle', DISCARD[6]), File('f.png'), 72)),
    Case('discarded/int_column_label', Q + 'make_figure_pairs_discarded', driver='figure',
         args=lambda: (columns_of('self circle', DISCARD[6], extra_columns={5: [1, 2, 3]}), File('f.png'), 72)),
    Case('discarded/int_column_label_exact_found', Q + 'make_figure_pairs_discarded', driver='figure',
         args=lambda: (columns_of(*DISCARD, extra_columns={5: [1, 2, 3]}), File('f.png'), 72)),
    Case('distance/missing_one', Q + 'make_figure_distance', driver='figure',
         args=lambda: (columns_of(*DISTANCE[1:]), File('f.png'))),
    Case('distance/missing_denominator', Q + 'make_figure_distance', driver='figure',
         args=lambda: (columns_of(*DISTANCE[:3]), File('f.png'))),
    Case('orientation/missing_one', Q + 'make_figure_read_orientation', driver='figure',
         args=lambda: (columns_of(*ORIENT[:3]), File('f.png'), 72)),
    # text, float, zero cells
    Case('pairs_used/text_cell', Q + 'make_figure_pairs_used', driver='figure',
         args=lambda: (columns_of(*USED, **{'Hi-C contacts': ['1', '2', '3']}), File('f.png'), 72)),
    Case('unmappable/text_numerator', Q + 'make_figure_umappable_non_unique_reads', driver='figure',
         args=lambda: (columns_of(*UNMAP, **{'Low mapping quality': ['a', 'b', 'c']}),
                       File('f.png'), 72)),
    Case('unmappable/text_denominator', Q + 'make_figure_umappable_non_unique_reads',
         driver='figure',
         args=lambda: (columns_of(*UNMAP, **{'Sequenced reads': ['a', 'b', 'c']}), File('f.png'), 72)),
    Case('unmappable/float_columns', Q + 'make_figure_umappable_non_unique_reads', driver='figure',
         args=lambda: (columns_of(*UNMAP, **{'Hi-C contacts': [1.5, 2.25, 3.125],
                                             'Sequenced reads': [10.0, 20.5, 30.0]}),
                       File('f.png'), 72)),
    Case('unmappable/zero_denominator', Q + 'make_figure_umappable_non_unique_reads',
         driver='figure',
         args=lambda: (columns_of(*UNMAP, **{'Sequenced reads': [0, 5, 0], 'Hi-C contacts': [0, 3, 1]}),
                       File('f.png'), 72)),
    Case('distance/zero_denominator', Q + 'make_figure_distance', driver='figure',
         args=lambda: (columns_of(*DISTANCE, **{'Hi-C contacts': [0, 0, 5]}), File('f.png'))),
    Case('orientation/float_columns', Q + 'make_figure_read_orientation', driver='figure',
         args=lambda: (columns_of(*ORIENT, **{ORIENT[0]: [1.5, 2.5, 3.5]}), File('f.png'), 72)),
    Case('orientation/zero_sum', Q + 'make_figure_read_orientation', driver='figure',
         args=lambda: (columns_of(*ORIENT, **{n: [0, 0, 1] for n in ORIENT}), File('f.png'), 72)),
    Case('orientation/text_cell', Q + 'make_figure_read_orientation', driver='figure',
         args=lambda: (columns_of(*ORIENT, **{ORIENT[1]: ['a', 'b', 'c']}), File('f.png'), 72)),
    Case('orientation/large_integers', Q + 'make_figure_read_orientation', driver='figure',
         args=lambda: (columns_of(*ORIENT, **{ORIENT[0]: [2 ** 53 + 1, 5, 7],
                                              ORIENT[1]: [3, 2 ** 53 + 3, 9]}), File('f.png'), 72)),
    Case('pairs_used/one_sample', Q + 'make_figure_pairs_used', driver='figure',
         args=lambda: (single(), File('f.png'), 200)),
    # save_html -----------------------------------------------------------------------
    Case('save_html/single', Q + 'save_html', driver='html',
         args=lambda: (single(), File('h.html'))),
    Case('save_html/multi', Q + 'save_html', driver='html',
         args=lambda: (multi(), File('h.html'))),
    Case('save_html/multi_restriction', Q + 'save_html', driver='html',
         args=lambda: (multi_restriction(), File('h.html'))),
    Case('save_html/file_labels', Q + 'save_html', driver='html',
         args=lambda: (file_labels(), File('h.html'))),
]


def frames_basic():
    import pandas as pd
    index = pd.Index(['a', 'b'], name='Labels')
    derived = pd.DataFrame({'c1': [1234567, 1], 'c1 %': [0.5, 0.0001234]}, index=index)
    other = pd.DataFrame({'c2': [1.0, 2.5], 'c2 %': [1.5, 12345.6789]}, index=index)
    third = pd.DataFrame({'c3': [0, 1], 'c3 %': [1, 2]}, index=index)
    fourth = pd.DataFrame({'c4': [-5, 2], 'c4 %': [float('nan'), float('inf')]}, index=index)
    full = pd.DataFrame({'Min rest. site distance': [300, 400], 'Max library insert size': [1000, 1100],
                         'Sequenced reads': [99983, 10], 'text': ['x<b>', 'y&z'],
                         'real': [0.1234567891, 123456.5]}, index=index)
    return derived, other, third, fourth, full


def frames_unnamed():
    import pandas as pd
    derived, other, third, fourth, full = frames_basic()
    for frame in (derived, other, third, fourth, full):
        frame.index.name = None
    return derived, other, third, fourth, full


def frames_named_columns():
    derived, other, third, fourth, full = frames_basic()
    for frame in (derived, other):
        frame.columns.name = 'what'
    return derived, other, third, fourth, full


def frames_int_index():
    import pandas as pd
    derived, other, third, fourth, full = frames_basic()
    for frame in (derived, other, third, fourth, full):
        frame.index = pd.Index([3, 40])
    return derived, other, third, fourth, full


def frames_float_index():
    import pandas as pd
    derived, other, third, fourth, full = frames_basic()
    for frame in (derived, other, third, fourth, full):
        frame.index = pd.Index([0.5, 12345.6789], name=7)
    return derived, other, third, fourth, full


def frames_int_columns():
    derived, other, third, fourth, full = frames_basic()
    derived.columns = [1, 2]
    other.columns = [0.5, 2]
    return derived, other, third, fourth, full


def frames_empty():
    import pandas as pd
    empty = pd.DataFrame({'c': []}, index=pd.Index([], name='Labels'))
    no_columns = pd.DataFrame(index=pd.Index(['a'], name='Labels'))
    return empty, no_columns, empty, no_columns, empty


def frames_unicode():
    import pandas as pd
    index = pd.Index(['äö', 'sample α'], name='Labels ü')
    frame = pd.DataFrame({'cß': [1, 2], 'd': [0.5, 3.5]}, index=index)
    return frame, frame, frame, frame, frame


def frames_text_in_derived():
    import pandas as pd
    derived, other, third, fourth, full = frames_basic()
    derived = derived.astype({'c1': object})
    derived['c1'] = ['a', 'b']
    return derived, other, third, fourth, full


def frames_text_at_end():
    import pandas as pd
    derived, other, third, fourth, full = frames_basic()
    third = third.astype({'c3': object})
    third['c3'] = [1, 'b']
    return derived, other, third, fourth, full


def frames_object_numbers():
    derived, other, third, fourth, full = frames_basic()
    derived = derived.astype(object)
    full = full.astype(object)
    return derived, other, third, fourth, full


def frames_one_row_wide():
    import pandas as pd
    wide = pd.DataFrame([list(range(1, 13))], columns=['c%d' % i for i in range(12)],
                        index=pd.Index(['only'], name='L'))
    return wide, wide, wide, wide, wide


def frames_boundary():
    import pandas as pd
    values = [0, 1, 2, 999, 1000, 1001, 10 ** 12, -1, -1000, 2 ** 62]
    floats = [0.0, 1.0, 1.0000000000000002, 999.9999, 1000.5, 1e16, 1e-7, 0.99999999999, -0.5,
              123456789.125]
    index = pd.Index(['r%d' % i for i in range(10)], name='Labels')
    first = pd.DataFrame({'ints': values, 'floats': floats}, index=index)
    return first, first, first, first, first


for name, make in [('basic', frames_basic), ('unnamed_index', frames_unnamed),
                   ('named_columns', frames_named_columns), ('int_index', frames_int_index),
                   ('float_index', frames_float_index), ('int_columns', frames_int_columns),
                   ('empty', frames_empty), ('unicode', frames_unicode),
                   ('text_in_derived', frames_text_in_derived),
                   ('text_at_end', frames_text_at_end),
                   ('object_numbers', frames_object_numbers),
                   ('one_row_wide', frames_one_row_wide), ('boundary', frames_boundary)]:
    CASES.append(Case('save_html/frames_' + name, Q + 'save_html', driver='html_frames',
                      args=lambda make=make: (File('h.html'), make)))

CASES += [
    Case('save_html/no_directory', Q + 'save_html', driver='html_frames',
         args=lambda: (File('missing/h.html'), frames_basic)),
    Case('save_html/is_directory', Q + 'save_html', driver='html_frames',
         args=lambda: (File(), frames_basic)),
    Case('save_html/bad_type', Q + 'save_html', driver='html_frames',
         args=lambda: (3.5, frames_basic)),
    Case('save_html/overwrite', Q + 'save_html', driver='html_frames',
         args=lambda: (File('h.html'), frames_basic), files={'h.html': b'old content that is longer ' * 50}),
]

# --- make_sure_path_exists -------------------------------------------------------------
CASES += [
    Case('path/new', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', 'b', 'c')),
    Case('path/new_one', Q + 'make_sure_path_exists', driver='tree', args=lambda: (File('base'), 'a')),
    Case('path/existing_dir', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'd'), kwargs={'prepare': 'dir'}),
    Case('path/existing_file', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'f'), kwargs={'prepare': 'file'}),
    Case('path/below_file', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'f', 'x'), kwargs={'prepare': 'file'}),
    Case('path/deep_below_file', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'f', 'x', 'y', 'z'), kwargs={'prepare': 'file'}),
    Case('path/trailing_slash', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', 'b'), kwargs={'trailing': True}),
    Case('path/dot_component', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', '.', 'b')),
    Case('path/trailing_dot', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', 'b', '.')),
    Case('path/dotdot', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', 'b', '..', 'c')),
    Case('path/double_slash', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a//b')),
    Case('path/bytes', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', 'b'), kwargs={'as_bytes': True}),
    Case('path/pathlike', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'a', 'b'), kwargs={'as_path': True}),
    Case('path/readonly_parent', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'ro', 'sub'), kwargs={'prepare': 'readonly'}),
    Case('path/readonly_deep', Q + 'make_sure_path_exists', driver='tree',
         args=lambda: (File('base'), 'ro', 'sub', 'deeper'), kwargs={'prepare': 'readonly'}),
    Case('path/empty', Q + 'make_sure_path_exists', args=('',)),
    Case('path/root', Q + 'make_sure_path_exists', args=('/',)),
    Case('path/none', Q + 'make_sure_path_exists', args=(None,)),
    Case('path/int', Q + 'make_sure_path_exists', args=(5,)),
    Case('path/relative_existing', Q + 'make_sure_path_exists', args=('.',)),
]
