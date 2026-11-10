// Unit tests for tools/quick_qc_impl.hpp, the NamedTemporaryFile replacement
// hicQuickQC uses for the matrix name that ends up in QC.log.

#include <cstdlib>
#include <filesystem>
#include <regex>
#include <set>
#include <string>

#include <doctest/doctest.h>

#include "../tools/quick_qc_impl.hpp"

namespace qq = hicx::quick_qc;

namespace {

// Sets an environment variable for the lifetime of the object and restores
// the previous value afterwards.
class ScopedEnv {
  public:
    ScopedEnv(const char* name, const std::string& value) : name_(name) {
        const char* previous = std::getenv(name);
        if (previous != nullptr) {
            had_ = true;
            previous_ = previous;
        }
        ::setenv(name, value.c_str(), 1);
    }
    ~ScopedEnv() {
        if (had_) {
            ::setenv(name_, previous_.c_str(), 1);
        } else {
            ::unsetenv(name_);
        }
    }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;

  private:
    const char* name_;
    bool had_ = false;
    std::string previous_;
};

std::filesystem::path scratch_directory(const std::string& tag) {
    const auto base = std::filesystem::temp_directory_path() /
                      ("hicx_test_quick_qc_" + tag + "_" + std::to_string(::getpid()));
    std::filesystem::remove_all(base);
    std::filesystem::create_directories(base);
    return base;
}

}  // namespace

TEST_CASE("temporary names have tempfile's shape") {
    std::mt19937_64 generator(12345);
    const std::regex shape("tmp[a-z0-9_]{8}\\.h5");
    std::set<std::string> names;
    for (int i = 0; i < 2000; ++i) {
        const std::string name = qq::random_temp_basename(generator, ".h5");
        CHECK(std::regex_match(name, shape));
        names.insert(name);
    }
    CHECK(names.size() == 2000);
}

TEST_CASE("create_named_temporary_file creates distinct files exclusively") {
    const auto directory = scratch_directory("create");
    std::set<std::string> paths;
    for (int i = 0; i < 50; ++i) {
        const std::string path = qq::create_named_temporary_file(directory.string(), ".h5");
        CHECK(std::filesystem::is_regular_file(path));
        CHECK(std::filesystem::file_size(path) == 0);
        CHECK(std::filesystem::path(path).parent_path() == directory);
        paths.insert(path);
    }
    CHECK(paths.size() == 50);
    CHECK_THROWS_AS((void)qq::create_named_temporary_file(
                        (directory / "missing" / "deeper").string(), ".h5"),
                    std::runtime_error);
    std::filesystem::remove_all(directory);
}

TEST_CASE("python_gettempdir honours TMPDIR, strips a trailing slash and skips unusable entries") {
    const auto directory = scratch_directory("tempdir");
    {
        const ScopedEnv tmpdir("TMPDIR", directory.string() + "/");
        CHECK(qq::python_gettempdir() == directory.string());
    }
    {
        const ScopedEnv tmpdir("TMPDIR", (directory / "does_not_exist").string());
        const ScopedEnv temp("TEMP", directory.string());
        CHECK(qq::python_gettempdir() == directory.string());
    }
    std::filesystem::remove_all(directory);
}
