// The C++ implementation behind hicexplorer.parserCommon.writableFile of
// HiCExplorer 3.7.6:
//
//     try:
//         open(string, 'w').close()
//     except IOError:
//         raise argparse.ArgumentTypeError("{} file can be opened for writting")
//
// open(path, 'w') is open(2) with O_WRONLY | O_CREAT | O_TRUNC and mode 0666;
// open(fd, 'w') checks the descriptor with fstat and refuses a directory, and
// close() closes it. Header only, no Python headers.

#ifndef HICX_PARSER_COMMON_HPP
#define HICX_PARSER_COMMON_HPP

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <string>

namespace hicx::parser_common {

// open(path, 'w').close(): 0, or the errno of the OSError open() raises.
[[nodiscard]] inline int truncate_for_writing(const std::string& path) {
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
    if (fd < 0) {
        return errno;
    }
    struct stat info {};
    if (::fstat(fd, &info) == 0 && S_ISDIR(info.st_mode)) {
        ::close(fd);
        return EISDIR;
    }
    ::close(fd);
    return 0;
}

// open(fd, 'w').close() for a descriptor fd >= 0: 0 after closing it, or the
// errno of the OSError open() raises (the descriptor then stays open, as
// io.FileIO leaves a descriptor it does not own).
[[nodiscard]] inline int close_for_writing(int fd) {
    struct stat info {};
    if (::fstat(fd, &info) != 0) {
        return errno;
    }
    if (S_ISDIR(info.st_mode)) {
        return EISDIR;
    }
    ::close(fd);
    return 0;
}

// The message of the ArgumentTypeError, `shown` being format(string, '').
[[nodiscard]] inline std::string writable_file_message(const std::string& shown) {
    return shown + " file can be opened for writting";
}

}  // namespace hicx::parser_common

#endif  // HICX_PARSER_COMMON_HPP
