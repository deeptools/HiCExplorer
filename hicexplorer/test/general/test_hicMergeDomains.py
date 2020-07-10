import os
import sys
from tempfile import NamedTemporaryFile
from tempfile import mkdtemp
from psutil import virtual_memory
import subprocess
import pytest
import logging
log = logging.getLogger(__name__)

import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images

from hicexplorer import hicMergeDomains

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 4
HIGH_MEMORY = 120

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def are_files_equal(file1, file2, delta=None):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    return equal


@pytest.mark.xfail(raises=SystemExit)
def test_main_one_file():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    args = "-d {} -om {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed',
        outfile_domains.name).split()
    hicMergeDomains.main(args)


def test_main_one_file_protein():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)

    args = "-d {} -om {} -p {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed',
        outfile_domains.name, ROOT + 'hicMergeDomains/ctcf_sorted.bed').split()
    hicMergeDomains.main(args)

    are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/one_file', delta=2)


@pytest.mark.xfail(sys.platform == "darwin", reason='Matplotlib plots for reasons a different image size.')
def test_main_two_file_protein():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_ctcf_relation = NamedTemporaryFile(suffix='.txt', delete=True)
    plot_folder = mkdtemp(prefix="plot_relations")

    args = "-d {} {} -om {} -or {} -ot {} -of {} -p {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed', ROOT + 'hicMergeDomains/50kbtad_domains.bed',
        outfile_domains.name, outfile_ctcf_relation.name, plot_folder + '/two_files_plot_ctcf', 'png', ROOT + 'hicMergeDomains/ctcf_sorted.bed').split()
    hicMergeDomains.main(args)
    are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/two_files_ctcf.bed', delta=2)
    are_files_equal(outfile_ctcf_relation.name, ROOT + 'hicMergeDomains/two_files_relation_ctcf.txt', delta=2)

    list_data = list(range(1, 23, 1))
    list_data = ["{}".format(x) for x in list_data]
    list_data.append('X')
    for i in list_data:
        res = compare_images(ROOT + '/hicMergeDomains/ctcf_plot/two_files_plot_ctcf_' + i + '.png', plot_folder + '/two_files_plot_ctcf_' + i + '.png', tol=40)
        assert res is None, res


@pytest.mark.xfail(sys.platform == "darwin", reason='Matplotlib plots for reasons a different image size.')
def test_main_two_file_no_protein():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_ctcf_relation = NamedTemporaryFile(suffix='.txt', delete=True)
    plot_folder = mkdtemp(prefix="plot_relations")

    args = "-d {} {} -om {} -or {} -ot {} -of {}".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed', ROOT + 'hicMergeDomains/50kbtad_domains.bed',
        outfile_domains.name, outfile_ctcf_relation.name, plot_folder + '/two_files_plot', 'png').split()
    hicMergeDomains.main(args)
    are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/two_files.bed', delta=2)
    are_files_equal(outfile_ctcf_relation.name, ROOT + 'hicMergeDomains/two_files_relation.txt', delta=2)

    list_data = list(range(1, 23, 1))
    list_data = ["{}".format(x) for x in list_data]
    list_data.append('X')
    for i in list_data:
        res = compare_images(ROOT + '/hicMergeDomains/no_ctcf_plot/two_files_plot_' + i + '.png', plot_folder + '/two_files_plot_' + i + '.png', tol=40)
        assert res is None, res


@pytest.mark.skipif(sys.platform == "linux", reason='Already tested')
def test_main_two_file_protein_no_image():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_ctcf_relation = NamedTemporaryFile(suffix='.txt', delete=True)
    plot_folder = mkdtemp(prefix="plot_relations")

    args = "-d {} {} -om {} -or {} -ot {} -of {} -p {} ".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed', ROOT + 'hicMergeDomains/50kbtad_domains.bed',
        outfile_domains.name, outfile_ctcf_relation.name, plot_folder + '/two_files_plot_ctcf', 'png', ROOT + 'hicMergeDomains/ctcf_sorted.bed').split()
    hicMergeDomains.main(args)
    are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/two_files_ctcf.bed', delta=2)
    are_files_equal(outfile_ctcf_relation.name, ROOT + 'hicMergeDomains/two_files_relation_ctcf.txt', delta=2)


@pytest.mark.skipif(sys.platform == "linux", reason='Already tested')
def test_main_two_file_no_protein_no_image():
    outfile_domains = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_ctcf_relation = NamedTemporaryFile(suffix='.txt', delete=True)
    plot_folder = mkdtemp(prefix="plot_relations")

    args = "-d {} {} -om {} -or {} -ot {} -of {}".format(
        ROOT + 'hicMergeDomains/10kbtad_domains.bed', ROOT + 'hicMergeDomains/50kbtad_domains.bed',
        outfile_domains.name, outfile_ctcf_relation.name, plot_folder + '/two_files_plot', 'png').split()
    hicMergeDomains.main(args)
    are_files_equal(outfile_domains.name, ROOT + 'hicMergeDomains/two_files.bed', delta=2)
    are_files_equal(outfile_ctcf_relation.name, ROOT + 'hicMergeDomains/two_files_relation.txt', delta=2)
