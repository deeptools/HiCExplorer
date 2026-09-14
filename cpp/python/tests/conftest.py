import support


def pytest_terminal_summary(terminalreporter):
    if not support.SUMMARY:
        return
    terminalreporter.section("hicx_matrix E2 comparisons")
    for oracle, (arrays, values, nonzero, max_rel) in sorted(support.SUMMARY.items()):
        terminalreporter.write_line(
            f"{oracle}: {arrays} arrays, {values} values ({nonzero} non-zero), "
            f"max relative difference {max_rel:.3e}")
