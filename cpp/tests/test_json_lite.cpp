#include <doctest/doctest.h>

#include "hicx/json_lite.hpp"

using hicx::json::parse;
using hicx::json::Type;

TEST_CASE("cooler attribute coercion") {
    // 'bin-size' of a variable bin cooler is the literal string "null".
    const auto null_value = parse("null");
    REQUIRE(null_value.has_value());
    CHECK(null_value->is_null());

    // 'format-version' is sometimes written as a string.
    const auto version = parse("3");
    REQUIRE(version.has_value());
    CHECK(version->type() == Type::Int);
    CHECK(version->as_int() == 3);

    // Everything that is not JSON stays a string on the Python side.
    CHECK_FALSE(parse("2018-11-08T20:03:12.737522").has_value());
    CHECK_FALSE(parse("unknown").has_value());
    CHECK_FALSE(parse("cooler-0.7.11").has_value());
    CHECK_FALSE(parse("HiCMatrix-16.dev").has_value());
    CHECK_FALSE(parse("").has_value());
}

TEST_CASE("metadata objects") {
    const auto empty = parse("{}");
    REQUIRE(empty.has_value());
    CHECK(empty->is_object());
    CHECK(empty->find("statistics") == nullptr);

    const auto metadata = parse(
        R"({"format": "HDF5::Cooler", "generated-by": "HiCMatrix-16.dev", )"
        R"("statistics": "Pairs considered\t1000\n", "count": 17, "ratio": 0.5})");
    REQUIRE(metadata.has_value());
    REQUIRE(metadata->is_object());
    const auto* statistics = metadata->find("statistics");
    REQUIRE(statistics != nullptr);
    CHECK(statistics->as_string() == "Pairs considered\t1000\n");
    CHECK(metadata->find("count")->as_int() == 17);
    CHECK(metadata->find("ratio")->as_double() == doctest::Approx(0.5));
    CHECK(metadata->find("missing") == nullptr);
}

TEST_CASE("malformed json is rejected instead of partially parsed") {
    CHECK_FALSE(parse("{\"a\": }").has_value());
    CHECK_FALSE(parse("[1, 2").has_value());
    CHECK_FALSE(parse("3 trailing").has_value());
}
