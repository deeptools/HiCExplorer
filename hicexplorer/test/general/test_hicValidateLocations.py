from hicexplorer import hicValidateLocations
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicValidateLocations/")


def are_files_equal(file1, file2, delta=1, skip=0):

    lines_file1_dict = {}
    mismatches = 0
    matches = 0
    line_count_file1 = 0
    with open(file1, 'r') as textfile1:
        file_content = textfile1.readlines()

        for i, line in enumerate(file_content):
            if i < skip:
                continue
            lines_file1_dict[line] = True
            line_count_file1 += 1
    with open(file2, 'r') as textfile2:

        file_content = textfile2.readlines()
        for i, line in enumerate(file_content):
            if i < skip:
                continue
            if line in lines_file1_dict:
                matches += 1
            else:
                mismatches += 1
    if mismatches < delta and line_count_file1 - delta <= matches:
        return True
    else:
        return False


def assert_files_identical(pExpected, pObserved, pSkip=0):
    """Line by line equality, in order.

    The original comparison in this file built a set of the lines of the
    reference and asked whether the observed lines were in it. That is blind to
    the row order, and the row order here is not incidental: the tool writes
    out a subset of the loop table in the order `bedtools sort` produced, and
    bedtools sorts on the start coordinate alone with a non-stable sort, so the
    3,730 rows of loops_1.bedgraph that share a start with another row come out
    in an order that is neither the input order nor sorted by any other field.
    A port that sorted stably would have passed the set comparison and written
    a differently ordered file.

    `pSkip` drops the leading comment lines of the statistics files, which
    carry the version and the absolute input paths.
    """
    with open(pExpected) as expected_file:
        expected = expected_file.readlines()[pSkip:]
    with open(pObserved) as observed_file:
        observed = observed_file.readlines()[pSkip:]
    nt.assert_equal(len(expected), len(observed))
    for index, (left, right) in enumerate(zip(expected, observed)):
        assert left == right, 'line {} differs:\n  expected {!r}\n  observed {!r}'.format(
            index + 1 + pSkip, left, right)


def read_statistics(pPath):
    """The four `key: value` lines of a statistics file, as a dictionary."""
    values = {}
    with open(pPath) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            key, _, value = line.partition(': ')
            values[key] = value.strip()
    return values


def test_loop_narrow_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} --chrPrefixLoops {} ".format(ROOT + 'loops_1.bedgraph',
                                                                                                                              ROOT + 'GSM935376_hg19_Gm12878_Smc3.narrowPeak', 'bed',
                                                                                                                              'loops', outfile.name, 10000, 'add').split()
    compute(hicValidateLocations.main, args, 5)
    assert_files_identical(
        ROOT + 'overlap_smc3_matched_locations', outfile.name + '_matched_locations')
    assert_files_identical(ROOT + 'overlap_smc3_statistics',
                           outfile.name + '_statistics', pSkip=3)


def test_loop_broad_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} --chrPrefixProtein {} ".format(ROOT + 'loops_1.bedgraph',
                                                                                                                                ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak', 'bed',
                                                                                                                                'loops', outfile.name, 10000, 'remove').split()
    compute(hicValidateLocations.main, args, 5)

    assert_files_identical(
        ROOT + 'overlap_ctcf_matched_locations', outfile.name + '_matched_locations')
    assert_files_identical(ROOT + 'overlap_ctcf_statistics',
                           outfile.name + '_statistics', pSkip=3)


def test_loop_cool():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} --chrPrefixLoop {} ".format(ROOT + 'loops_1.bedgraph',
                                                                                                                             ROOT + 'GSM1436265_RAD21_ENCFF002EMQ_10kb.cool', 'cool',
                                                                                                                             'loops', outfile.name, 10000, 'add').split()
    compute(hicValidateLocations.main, args, 5)

    assert_files_identical(
        ROOT + 'overlap_rad21_cool_matched_locations', outfile.name + '_matched_locations')
    assert_files_identical(ROOT + 'overlap_rad21_cool_statistics',
                           outfile.name + '_statistics', pSkip=3)


def test_tad():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "--data {} --validationData {} --validationType {} --method {} --outFileName {} -r {} ".format(ROOT + 'untreated_R1_boundaries.bed',
                                                                                                          ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak', 'bed',
                                                                                                          'tad', outfile.name, 20000).split()
    compute(hicValidateLocations.main, args, 5)

    assert_files_identical(ROOT + 'validatedTADs.txt', outfile.name)
    assert_files_identical(ROOT + 'validatedTADs.txt_statistics',
                           outfile.name + '_statistics', pSkip=3)

    # --method tad writes the matched TADs to the bare --outFileName, not to
    # <name>_matched_locations as the option's help text says. Pinned so that a
    # port does not quietly move the file to the documented name.
    assert not os.path.exists(outfile.name + '_matched_locations')


def test_loop_matched_locations_is_a_subsequence_of_the_sorted_input():
    """The output order is the sorted loop table's order, not the input file's.

    Independent of any master: the tool's own sorted table is reconstructed
    here with pybedtools, and the matched rows must appear in it in the same
    order and with the same text. This is the assertion that fails if the
    bedtools sort order is not reproduced, and it does not depend on the
    checked-in reference file staying in step with the bedtools version.
    """
    import pandas as pd
    from pybedtools import BedTool

    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = ("--data {} --validationData {} --validationType bed --method loops "
            "--outFileName {} -r 10000 --chrPrefixLoops add").format(
        ROOT + 'loops_1.bedgraph', ROOT + 'GSM935376_hg19_Gm12878_Smc3.narrowPeak',
        outfile.name).split()
    compute(hicValidateLocations.main, args, 5)

    loops = pd.read_csv(ROOT + 'loops_1.bedgraph', sep='\t', header=None)
    loops[0] = 'chr' + loops[0].astype(str)
    loops[3] = 'chr' + loops[3].astype(str)
    sorted_lines = BedTool.from_dataframe(loops).sort().to_dataframe(
        disable_auto_names=True, header=None).to_csv(
            sep='\t', header=False, index=False).splitlines()

    with open(outfile.name + '_matched_locations') as handle:
        matched = handle.read().splitlines()

    position = 0
    for line in matched:
        while position < len(sorted_lines) and sorted_lines[position] != line:
            position += 1
        assert position < len(sorted_lines), \
            'matched line not found in the sorted table, or out of order: ' + line
        position += 1
    nt.assert_equal(len(matched), 6379)


def test_loop_chr_prefix_protein_add():
    """--chrPrefixProtein add, which no test exercised.

    The narrowPeak file already carries a chr prefix, so adding another one
    gives chromosome names such as chrchr1. The loops get the prefix as well,
    which is what makes the case discriminating: with both options the two
    files no longer share a chromosome name and nothing matches, while a
    version that ignored --chrPrefixProtein would leave the peaks at chr1 and
    match 6,379 loops. Every peak is still counted, and the matched-locations
    file is written and empty.
    """
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = ("--data {} --validationData {} --validationType bed --method loops "
            "--outFileName {} -r 10000 --chrPrefixLoops add "
            "--chrPrefixProtein add").format(
        ROOT + 'loops_1.bedgraph', ROOT + 'GSM935376_hg19_Gm12878_Smc3.narrowPeak',
        outfile.name).split()
    compute(hicValidateLocations.main, args, 5)

    statistics = read_statistics(outfile.name + '_statistics')
    nt.assert_equal(statistics['Protein peaks'], '34474')
    nt.assert_equal(statistics['Matched Loops'], '0')
    nt.assert_equal(statistics['Total Loops'], '11723')
    nt.assert_equal(statistics['Loops match protein'], '0.0')
    nt.assert_equal(os.path.getsize(outfile.name + '_matched_locations'), 0)


def test_loop_chr_prefix_loops_remove():
    """--chrPrefixLoops remove, which no test exercised.

    The loop file has no chr prefix, and `Series.str.lstrip('chr')` strips a
    character set rather than a prefix, so on chromosome names such as 1, 19
    and X it removes nothing. The run must therefore agree with the default in
    every output byte. That is the assertion: a port that implemented lstrip as
    a prefix strip would still pass here, but one that stripped the leading
    character of every name would not, and neither would one that skipped the
    option's re-inference of the column dtype.
    """
    with_option = NamedTemporaryFile(suffix='out', delete=True)
    with_option.close()
    without_option = NamedTemporaryFile(suffix='out', delete=True)
    without_option.close()

    template = ("--data {} --validationData {} --validationType bed --method loops "
                "--outFileName {} -r 10000 --chrPrefixProtein remove")
    compute(hicValidateLocations.main,
            (template + " --chrPrefixLoops remove").format(
                ROOT + 'loops_1.bedgraph',
                ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak',
                with_option.name).split(), 5)
    compute(hicValidateLocations.main,
            template.format(ROOT + 'loops_1.bedgraph',
                            ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak',
                            without_option.name).split(), 5)

    assert_files_identical(without_option.name + '_matched_locations',
                           with_option.name + '_matched_locations')
    nt.assert_equal(read_statistics(with_option.name + '_statistics'),
                    read_statistics(without_option.name + '_statistics'))


def test_no_out_file_name_writes_nothing(capsys):
    """--outFileName is optional and the tool then only prints.

    No test ran without it. The four printed lines are the whole output in that
    case, and the TAD branch spells them differently from the loop branch
    (Matched TADs against Matched Loops), which is pinned here.
    """
    args = ("--data {} --validationData {} --validationType bed --method tad "
            "-r 20000").format(ROOT + 'untreated_R1_boundaries.bed',
                               ROOT + 'GSM733752_hg19_ctcf_GM12878.broadPeak').split()
    compute(hicValidateLocations.main, args, 5)

    captured = capsys.readouterr()
    assert 'Protein peaks: 22394' in captured.out
    assert 'Matched TADs: 99' in captured.out
    assert 'Total TADs: 201' in captured.out
    assert 'TADs match protein: 0.4925373134328358' in captured.out
