import os.path
from tempfile import NamedTemporaryFile, mkdtemp
import shutil

import logging
log = logging.getLogger(__name__)

from hicexplorer import hicQuickQC

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")

sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"


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


def test_main():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --QCfolder {} -seq {} --lines 1000".format(sam_R1, sam_R2,
                                                                qc_folder, 'GATC').split()
    hicQuickQC.main(args)

    print(set(os.listdir(ROOT + "hicQuickQC/")))
    assert are_files_equal(ROOT + "hicQuickQC/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT + "hicQuickQC/")) == set(os.listdir(qc_folder))

    shutil.rmtree(qc_folder)
