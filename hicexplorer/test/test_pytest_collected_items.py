"""
    Testsuite for testing the number of collected items during pytest execution.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import pytest
import re
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data/")


def test_number_of_collected_items(capsys):
    """
    This test will fail if the number of tests is lower than the number of tests when
    running pytest the last time.

    If test run was succesful new number of tests is stored into file for next time.
    """
    # first read file with the number of tests running last time
    path = ROOT + 'number_of_tests.txt'
    with open(path, 'r') as f:
        number_of_tests = int(''.join([x for x in f]))

    # run pytest and capture all tests running this time
    pytest.main(['hicexplorer', '--doctest-modules', '--collect-only'])
    captured = capsys.readouterr()
    output_raw = ''.join(captured.out)
    output = ''.join(output_raw)
    lines = output.split('\n')

    # search for number of collected items with regex
    regex_collected = r"collected (\d+).+"

    for line in lines:
        matches = re.match(regex_collected, line)
        if matches:
            assert int(matches.group(1)) >= number_of_tests
            current_tests = int(matches.group(1))

    # write new number of tests into file
    with open(path, 'w') as outfile:
        outfile.write(str(current_tests))
