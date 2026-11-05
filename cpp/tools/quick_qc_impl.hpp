// hicQuickQC's one piece of its own: the temporary matrix name.
//
// hicQuickQC.py:100 creates NamedTemporaryFile(suffix='.h5', delete=False)
// and hands its name to hicBuildMatrix as --outFileName. A test run never
// writes a matrix, and buildMatrixMethods.py:1366 unlinks the file after the QC
// report, but the name is printed into the first line of QC.log and into the
// first column of every *_table.txt, so it is part of the output. The port
// creates the same kind of name in the same directory, exclusively, and
// removes it afterwards.

#ifndef HICX_TOOLS_QUICK_QC_IMPL_HPP
#define HICX_TOOLS_QUICK_QC_IMPL_HPP

#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace hicx::quick_qc {

// tempfile._RandomNameSequence.characters and its name length.
inline constexpr const char* kTempCharacters = "abcdefghijklmnopqrstuvwxyz0123456789_";
inline constexpr std::size_t kTempNameLength = 8;

// tempfile._candidate_tempdir_list: TMPDIR, TEMP and TMP, then the platform
// directories. The first directory in which a file can be created wins, as in
// tempfile._get_default_tempdir, and it is made absolute.
inline std::string python_gettempdir() {
    std::vector<std::string> candidates;
    for (const char* variable : {"TMPDIR", "TEMP", "TMP"}) {
        const char* value = std::getenv(variable);
        if (value != nullptr && *value != '\0') {
            candidates.emplace_back(value);
        }
    }
    for (const char* directory : {"/tmp", "/var/tmp", "/usr/tmp"}) {
        candidates.emplace_back(directory);
    }
    std::error_code ec;
    candidates.push_back(std::filesystem::current_path(ec).string());
    for (const auto& candidate : candidates) {
        std::filesystem::path directory = std::filesystem::absolute(candidate, ec);
        directory = directory.lexically_normal();
        std::string text = directory.string();
        while (text.size() > 1 && text.back() == '/') {
            text.pop_back();
        }
        if (::access(text.c_str(), W_OK | X_OK) == 0) {
            return text;
        }
    }
    throw std::runtime_error("No usable temporary directory found");
}

// "tmp" + eight characters of kTempCharacters + suffix.
inline std::string random_temp_basename(std::mt19937_64& generator,
                                        const std::string& suffix) {
    std::uniform_int_distribution<int> pick(0, 36);
    std::string name = "tmp";
    for (std::size_t i = 0; i < kTempNameLength; ++i) {
        name += kTempCharacters[pick(generator)];
    }
    return name + suffix;
}

// tempfile.NamedTemporaryFile(suffix=suffix, delete=False).name: the file is
// created with O_EXCL, retried on a collision, and left in place.
inline std::string create_named_temporary_file(const std::string& directory,
                                               const std::string& suffix) {
    std::random_device device;
    std::mt19937_64 generator((static_cast<std::uint64_t>(device()) << 32) ^ device() ^
                              static_cast<std::uint64_t>(::getpid()));
    for (int attempt = 0; attempt < 100; ++attempt) {
        const std::string path = directory + "/" + random_temp_basename(generator, suffix);
        const int fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
        if (fd >= 0) {
            ::close(fd);
            return path;
        }
        if (errno != EEXIST) {
            break;
        }
    }
    throw std::runtime_error("could not create a temporary file in " + directory);
}

}  // namespace hicx::quick_qc

#endif  // HICX_TOOLS_QUICK_QC_IMPL_HPP
