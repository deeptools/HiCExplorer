"""Tests that pin behaviour of the Python implementation the C++ tools do not
reproduce.

These characterization tests were written against the Python implementation
(upstream HiCExplorer at c2ad8630, cpp/AGENTS_CONTRACT.md rule 1) before the
tool was ported. The tool modules now run the C++ executables, and each test
below fails against them for the reason of its group, already on
version4-cpp before the Python implementation left the package. They are
strict expected failures: a C++ change that makes one pass shows up as a
failure, so the list stays exact. The comparison with the Python reference is
the equivalence harness's job (cpp/scripts/equiv.py), which runs the reference
tree of cpp/scripts/reference_tree.py.
"""
import pytest

PYTHON_EXCEPTION = (
    "pins an exception or crash of the Python implementation inside the calling process; "
    "the C++ tool reports the error on stderr and exits with status 1 instead, or does not fail")
PYTHON_DEFECT = (
    "pins a defect of the Python implementation (thread splitting, --doTestRun output, "
    "--numberOfProcessors) that the C++ tool does not reproduce")
BIT_EXACT = (
    "requires the bit-identical floating-point result of the Python implementation; the C++ "
    "result differs in the last bits and meets its equivalence class (cpp/STATUS.md)")
CHROMOSOME_SIZES = (
    "hicMergeMatrixBins requires --chromosomeSizes and derives the merged layout from the "
    "chromosome lengths (deliberate deviation of 2026-09-22, cpp/STATUS.md); these tests pin the "
    "Python layout derived from the bins the input holds")

IN_PROCESS_CAPTURE = (
    "records the DOT sources by patching the graphviz Python module inside the calling process; "
    "the C++ tool writes them in its own process, and the same recordings pin the C++ DOT writer "
    "in cpp/tests/test_merge_domains.cpp")

EXPECTED_FAILURES = {
    "test_chicAggregateStatistic.py::test_characterization_failures_and_thread_counts": PYTHON_EXCEPTION,
    "test_chicDifferentialTest.py::test_characterization_calls_and_failures": PYTHON_EXCEPTION,
    "test_chicExportData.py::test_characterization_options_and_failures": PYTHON_EXCEPTION,
    "test_chicSignificantInteractions.py::test_characterization_threshold_file_fails_when_written": PYTHON_EXCEPTION,
    "test_chicSignificantInteractions.py::test_characterization_thread_counts": PYTHON_EXCEPTION,
    "test_hicAdjustMatrix.py::test_maskBadRegions_raises_a_type_error_on_cool": PYTHON_EXCEPTION,
    "test_hicBuildMatrix.py::test_build_matrix_do_test_run_writes_qc_but_no_matrix": PYTHON_DEFECT,
    "test_hicBuildMatrix.py::test_build_matrix_input_buffer_size_crashes_on_an_exact_multiple": PYTHON_EXCEPTION,
    "test_hicCompareMatrices.py::test_log2ratio_h5_is_bit_exact": BIT_EXACT,
    "test_hicCompareMatrices.py::test_log2ratio_noNorm_h5_is_bit_exact": BIT_EXACT,
    "test_hicCompareMatrices.py::test_log2ratio_cool_is_bit_exact": BIT_EXACT,
    "test_hicCompartmentalization.py::test_quantile_one_fails_and_quantile_zero_writes_an_empty_line": PYTHON_EXCEPTION,
    "test_hicCorrectMatrix.py::test_correct_matrix_ICE": BIT_EXACT,
    "test_hicCorrectMatrix.py::test_char_ICE_whole_is_reproducible_and_format_independent": BIT_EXACT,
    "test_hicCorrectMatrix.py::test_char_ICE_inflationCutoff_without_transCutoff_raises": PYTHON_EXCEPTION,
    "test_hicCorrectMatrix.py::test_char_ICE_inflationCutoff_with_transCutoff_raises": PYTHON_EXCEPTION,
    "test_hicCorrectMatrix.py::test_char_ICE_sequencedCountCutoff_rejects_a_cool_input": PYTHON_EXCEPTION,
    "test_hicDifferentialTAD.py::test_threads_small_chromosome_first_drops_last_tad": PYTHON_DEFECT,
    "test_hicDifferentialTAD.py::test_threads_small_chromosome_last_duplicates_rows": PYTHON_DEFECT,
    "test_hicDifferentialTAD.py::test_threads_one_tad_per_thread_loses_last_left_test": PYTHON_DEFECT,
    "test_hicFindTADs.py::test_find_TADs_high_numberOfProcessors_is_a_defect": PYTHON_DEFECT,
    "test_hicFindTADs.py::test_find_TADs_maxDepth_must_exceed_minDepth": PYTHON_EXCEPTION,
    "test_hicFindTADs.py::test_find_TADs_unreadable_matrix_never_reaches_the_suffix_check": PYTHON_EXCEPTION,
    "test_hicInterIntraTAD.py::test_threads_small_chromosome_last_duplicates_rows": PYTHON_DEFECT,
    "test_hicInterIntraTAD.py::test_threads_small_chromosome_first_drops_last_tad": PYTHON_DEFECT,
    "test_hicMergeDomains.py::test_relation_tree_sources_are_pinned[two_files]": IN_PROCESS_CAPTURE,
    "test_hicMergeDomains.py::test_relation_tree_sources_are_pinned[three_files_protein]": IN_PROCESS_CAPTURE,
    "test_hicMergeDomains.py::test_relation_tree_sources_are_pinned[reversed_order]": IN_PROCESS_CAPTURE,
    "test_hicMergeMatrixBins.py::test_correct_matrix": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_running_window_3_is_bit_exact": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_running_window_5_is_bit_exact": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_running_window_keeps_the_resolution_and_the_intervals": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_running_window_of_one_bin_only_drops_the_nan_bins": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_running_window_rejects_an_even_number_of_bins": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_running_window_cool_output_is_bit_exact": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_merge_bins_4_is_bit_exact": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_merge_bins_drops_chromosomes_with_too_few_bins": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_merge_bins_cut_intervals_span_the_removed_nan_bins": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_merge_bins_cool_output_is_bit_exact": CHROMOSOME_SIZES,
    "test_hicMergeMatrixBins.py::test_merge_bins_5_reproduces_the_stored_reference_exactly": CHROMOSOME_SIZES,
    "test_hicMergeTADbins.py::test_domains_outside_the_matrix_raise_a_type_error": PYTHON_EXCEPTION,
    "test_hicTransform.py::test_pearson_matches_numpy_corrcoef_at_full_precision": BIT_EXACT,
    "test_hicTransform.py::test_pearson_perChromosome_matches_numpy_corrcoef_at_full_precision": BIT_EXACT,
    "test_hicTransform.py::test_covariance_matches_numpy_cov_at_full_precision": BIT_EXACT,
    "test_hicTransform.py::test_covariance_perChromosome_matches_numpy_cov_at_full_precision": BIT_EXACT,
    "test_hicTransform.py::test_chromosomes_changes_the_pearson_result_not_just_its_size": BIT_EXACT,
    "test_hicTransform.py::test_unknown_chromosome_is_an_error": PYTHON_EXCEPTION,
}


def pytest_collection_modifyitems(config, items):
    for item in items:
        key = item.nodeid.split("/")[-1]
        reason = EXPECTED_FAILURES.get(key)
        if reason is not None:
            item.add_marker(pytest.mark.xfail(strict=True, reason=reason))
