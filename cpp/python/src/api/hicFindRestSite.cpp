// hicexplorer._hicx.hicFindRestSite: the C++ body of
// hicexplorer.hicFindRestSite.find_pattern. The scan, the sort and the BED
// text are hicx::fasta::find_restriction_sites and restriction_sites_bed
// (cpp/core/include/hicx/fasta_reader.hpp), which the hicFindRestSite
// executable uses too.
//
// Deviations from 3.7.6, recorded in hicexplorer/test/api/cases_hicFindRestSite.py:
//  * 3.7.6 writes the unsorted sites into a NamedTemporaryFile(delete=False)
//    and shells out to `sort`; the temporary file stayed behind in the
//    temporary directory after every call. The C++ sorts in memory and leaves
//    no file.
//  * a pattern that is not letters and '.' is matched by std::regex
//    (ECMAScript), not by Python's re: constructs only Python has behave
//    differently, and an invalid pattern raises ValueError rather than
//    re.error.
//  * a pattern that is not a str raises TypeError naming HiCExplorer 4.

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <regex>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/fasta_reader.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicFindRestSite";

// open(name) of the FASTA (or gzip.open, which opens it with open(name, 'rb')):
// the OSError open() raises, before anything is read.
void check_readable(const py::object& name) {
    const std::string path = fs_path(name);
    // open() names os.fspath(name) in its error
    const py::object shown = py::module_::import("os").attr("fspath")(name);
    const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) {
        raise_os_error(errno, shown);
    }
    struct stat info {};
    const bool directory = ::fstat(fd, &info) == 0 && S_ISDIR(info.st_mode);
    ::close(fd);
    if (directory) {
        raise_os_error(EISDIR, shown);
    }
}

}  // namespace

void bind_hicFindRestSite(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.hicFindRestSite.";

    m.def("find_pattern", [](py::handle pPattern, py::handle fasta_file, py::handle out_file) {
        std::vector<std::string> patterns;
        for (py::handle pattern : py::iter(pPattern)) {
            if (!is_str(pattern)) {
                raise_error("TypeError", "find_pattern: HiCExplorer 4 takes the search patterns as str, not " +
                                             type_name(pattern));
            }
            patterns.push_back(utf8(pattern));
        }
        std::string fasta_path;
        if (!patterns.empty()) {
            // fasta_file.name, or fasta_file itself when it has no name
            py::object name;
            if (py::hasattr(fasta_file, "name")) {
                name = fasta_file.attr("name");
            } else {
                name = py::reinterpret_borrow<py::object>(fasta_file);
            }
            check_readable(name);
            fasta_path = fs_path(name);
        }
        std::string bed;
        try {
            bed = hicx::fasta::restriction_sites_bed(hicx::fasta::find_restriction_sites(patterns, fasta_path));
        } catch (const std::regex_error& error) {
            raise_error("ValueError", std::string("invalid search pattern: ") + error.what());
        }
        log(kLogger, "info", "Sorting file ...");
        out_file.attr("write")(to_str(bed));
        out_file.attr("close")();
    }, py::arg("pPattern"), py::arg("fasta_file"), py::arg("out_file"));
}

}  // namespace hicx::python
