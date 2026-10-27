// The core code behind hicexplorer.hicPrepareQCreport (hicx/qc_html.hpp, hicx/qc_report_api.hpp).
// Their equivalence with HiCExplorer 3.7.6 is tested from Python
// (hicexplorer/test/api); these tests pin the C++ contracts.

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "hicx/python_error.hpp"
#include "hicx/qc_html.hpp"
#include "hicx/qc_report_api.hpp"


TEST_CASE("qc_styler_html: header rows and formatters") {
    hicx::QcDisplayTable table;
    table.index_name = "Labels";
    table.index = {"a"};
    table.columns = {"count", "fraction"};
    hicx::QcStylerCell count;
    count.is_int = true;
    count.integer = 1234567;
    hicx::QcStylerCell fraction;
    fraction.is_float = true;
    fraction.real = 0.25;
    table.cells = {{count}, {fraction}};
    const std::string derived = hicx::qc_styler_html(table, hicx::QcFormatter::Derived);
    CHECK(derived.find(">1,234,567</td>") != std::string::npos);
    CHECK(derived.find(">25.00%</td>") != std::string::npos);
    CHECK(derived.find("class=\"index_name level0\" >Labels</th>") != std::string::npos);
    const std::string plain = hicx::qc_styler_html(table, hicx::QcFormatter::Plain);
    CHECK(plain.find(">1234567</td>") != std::string::npos);
    CHECK(plain.find(">0.250000</td>") != std::string::npos);
    // without an index name there is no second header row
    table.index_name.reset();
    CHECK(hicx::qc_styler_html(table, hicx::QcFormatter::Plain).find("index_name") ==
          std::string::npos);
    // a text cell has no meaning for the derived formatter
    table.cells[0][0].is_int = false;
    table.cells[0][0].text = "x";
    CHECK_THROWS_AS((void)hicx::qc_styler_html(table, hicx::QcFormatter::Derived),
                    hicx::PythonError);
}

TEST_CASE("qcapi: divisions and sums as pandas does them") {
    using hicx::QcStylerCell;
    QcStylerCell two;
    two.is_int = true;
    two.integer = 2;
    QcStylerCell text;
    text.text = "t";
    CHECK(hicx::qcapi::divide_by_scalar({two}, 4.0) == std::vector<double>{0.5});
    CHECK_THROWS_WITH((void)hicx::qcapi::divide_by_scalar({text}, 4.0),
                      "unsupported operand type(s) for /: 'str' and 'float'");
    CHECK_THROWS_WITH((void)hicx::qcapi::divide_columns({two}, {text}),
                      "unsupported operand type(s) for /: 'int' and 'str'");
    const auto total = hicx::qcapi::sum_columns({{two}, {two}});
    REQUIRE(total.size() == 1);
    CHECK(total[0].is_int);
    CHECK(total[0].integer == 4);
}

