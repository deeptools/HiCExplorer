from hicexplorer import hicFindRestSite
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicFindRestSite/")


def are_files_equal(file1, file2, delta=1, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            # if x.startswith('File'):
            #     continue
            if i < skip:
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


def test_fasta_gz():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} -o {}".format(ROOT + 'dm3_chrM.fasta.gz', 'AAGCTT', outfile.name).split()
    hicFindRestSite.main(args)

    assert are_files_equal(ROOT + "hindIII_chrM.bed",
                           outfile.name, skip=0)


def test_fasta():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} -o {}".format(ROOT + 'dm3_chrM.fasta', 'AAGCTT', outfile.name).split()
    hicFindRestSite.main(args)

    assert are_files_equal(ROOT + "hindIII_chrM.bed",
                           outfile.name, skip=0)


def test_fasta_gz_two_patterns():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} {} -o {}".format(ROOT + 'dm3_chrM.fasta.gz', 'AAGCTT', 'GATC', outfile.name).split()
    hicFindRestSite.main(args)

    assert are_files_equal(ROOT + "hindIII_DpnII.bed",
                           outfile.name, skip=0)


def test_fasta_two_patterns():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} {} -o {}".format(ROOT + 'dm3_chrM.fasta', 'AAGCTT', 'GATC', outfile.name).split()
    hicFindRestSite.main(args)

    assert are_files_equal(ROOT + "hindIII_DpnII.bed",
                           outfile.name, skip=0)
