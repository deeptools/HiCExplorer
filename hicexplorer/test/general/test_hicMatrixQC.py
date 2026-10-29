import os
from tempfile import NamedTemporaryFile

from hicexplorer import hicMatrixQC

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "test_data/")


def test_parse_arguments():
    parser = hicMatrixQC.parse_arguments()
    args = parser.parse_args(['--matrices', 'a.cool', '--outFileName', 'qc.html'])
    assert args.matrices == ['a.cool']
    assert args.outFileName == 'qc.html'


def test_report_and_table():
    report = NamedTemporaryFile(suffix='.html', delete=False)
    table = NamedTemporaryFile(suffix='.tsv', delete=False)
    report.close()
    table.close()
    try:
        hicMatrixQC.main(['--matrices', ROOT + 'gm12878_chr1.cool', '--outFileName', report.name,
                          '--tableFileName', table.name])
        with open(report.name) as handle:
            html = handle.read()
        assert 'hicMatrixQC report' in html
        assert '<script' not in html
        with open(table.name) as handle:
            rows = dict(line.rstrip('\n').split('\t')[:2] for line in handle)
        assert int(rows['contacts.total']) == 326878012
        assert float(rows['decay.exponent.genome']) < 0
    finally:
        os.unlink(report.name)
        os.unlink(table.name)
