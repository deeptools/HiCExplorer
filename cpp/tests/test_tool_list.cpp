// The core code behind hicexplorer.list_tools and the hicexplorer command (hicx/tool_list.hpp).
// Their equivalence with HiCExplorer 3.7.6 is tested from Python
// (hicexplorer/test/api); these tests pin the C++ contracts.

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "hicx/tool_list.hpp"


TEST_CASE("tool_list: names, summaries, alignment") {
    namespace tl = hicx::toollist;
    tl::DistributionFile in_bin{{"..", "bin", "hicInfo"}, "hicInfo", "hicInfo", ""};
    tl::DistributionFile exe{{"x", "Scripts", "hicQC.exe"}, "hicQC.exe", "hicQC", ".exe"};
    tl::DistributionFile other{{"x", "lib", "hicNo"}, "hicNo", "hicNo", ""};
    tl::DistributionFile itself{{"bin", "hicexplorer"}, "hicexplorer", "hicexplorer", ""};
    CHECK(tl::command_names({"chicX", "pytest", "hicexplorer"}, {in_bin, exe, other, itself}) ==
          std::vector<std::string>{"chicX", "hicInfo", "hicQC"});
    CHECK(tl::summary("\n  First one. Second.\nline") == "First one");
    CHECK(tl::summary("Ends...") == "Ends");
    CHECK(tl::summary("") == "");
    CHECK(tl::format_tool_list({{"hicA", "a"}, {"hicLonger", "b"}}) ==
          "hicA      a\nhicLonger b\n");
}
