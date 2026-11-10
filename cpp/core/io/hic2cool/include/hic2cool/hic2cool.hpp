// hic2cool: converting Juicer .hic files to cool and mcool files and back.
//
// This is the C++ implementation of the hic2cool fork. hic2cool_convert,
// hic2cool_extractnorms and hic2cool_update take the arguments of the Python
// functions of hic2cool 1.0.1 (https://github.com/4dn-dcic/hic2cool, MIT
// licence), print the same messages and write the same files. Beyond the
// Python package they read .hic versions 6 to 9 (the Python package reads 6
// to 8), use every available core unless told otherwise, and are joined by
// cool2hic_convert, the conversion in the other direction.

#ifndef HIC2COOL_HIC2COOL_HPP
#define HIC2COOL_HIC2COOL_HPP

#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "hic2cool/version.hpp"

namespace hic2cool {

// A condition under which the Python package prints a message to standard
// error and calls sys.exit(1) (hic2cool_utils.force_exit). what() is the
// message.
class ExitError : public std::runtime_error {
  public:
    explicit ExitError(const std::string& message) : std::runtime_error(message) {}
};

// Where messages go. The command line tools print to the process's standard
// streams; the Python module forwards to sys.stdout, sys.stderr and input(),
// so that redirecting those in Python captures the messages as it does for
// the Python package. Only the calling thread uses the console.
struct Console {
    // One line of standard output, without the trailing newline.
    std::function<void(const std::string&)> out;
    // One line of standard error, without the trailing newline.
    std::function<void(const std::string&)> err;
    // input(prompt): prints the prompt and returns the line read.
    std::function<std::string(const std::string&)> input;

    static Console standard();
};

// The number of threads nproc = 0 stands for: the CPUs this process may run
// on (its affinity mask, which taskset, cgroups and batch schedulers set).
[[nodiscard]] int available_threads();

// hic2cool_convert(infile, outfile, resolution=0, nproc=0,
// show_warnings=False, silent=False). resolution 0 writes every base pair
// resolution of the file as a multi resolution (.mcool) file, any other value
// that one resolution as a single resolution .cool file. As in the Python
// package the output name is adjusted to .cool or .mcool; the name written is
// returned. nproc 0 uses available_threads(); the output does not depend on
// nproc. storage_mode "symmetric-upper" (hic2cool's layout) stores the upper
// triangle; "square" stores both triangles, cooler's layout for asymmetric
// matrices, for tools that read square coolers only.
std::string hic2cool_convert(const std::string& infile, const std::string& outfile,
                             std::int64_t resolution = 0, int nproc = 0,
                             bool show_warnings = false, bool silent = false,
                             const std::string& storage_mode = "symmetric-upper",
                             const Console& console = Console::standard());

// The multi resolution layout for a chosen list of resolutions (every
// resolution of the file when the list is empty), written to exactly
// `outfile`: no renaming, and a multi resolution file even for one
// resolution. Returns the path written. Otherwise as hic2cool_convert.
std::string hic2cool_convert_mcool(const std::string& infile, const std::string& outfile,
                                   const std::vector<std::int64_t>& resolutions = {}, int nproc = 0,
                                   bool show_warnings = false, bool silent = false,
                                   const Console& console = Console::standard());

// hic2cool_extractnorms(infile, outfile, exclude_mt=False,
// show_warnings=False, silent=False): adds the normalization vectors of a
// .hic file to the bins tables of the cooler groups of outfile whose
// resolution the .hic file has.
void hic2cool_extractnorms(const std::string& infile, const std::string& outfile,
                           bool exclude_mt = false, bool show_warnings = false,
                           bool silent = false, const Console& console = Console::standard());

// hic2cool_update(infile, outfile='', show_warnings=False, silent=False):
// upgrades a cool or mcool file written by an older hic2cool.
void hic2cool_update(const std::string& infile, const std::string& outfile = "",
                     bool show_warnings = false, bool silent = false,
                     const Console& console = Console::standard());

struct Cool2hicOptions {
    // 0: every resolution of the input (all groups of an mcool file); any
    // other value: that resolution, which the input must have.
    std::int64_t resolution = 0;
    // Coarser resolutions to add, binned from the finest resolution written.
    std::vector<std::int64_t> add_resolutions;
    int nproc = 0;
    // .hic version to write: 8 (Juicer tools 1.22) or 9 (Juicer tools 2.x).
    int hic_version = 9;
    // "auto": carry the normalization columns of the bins table (the ones
    // hic2cool writes) into the .hic file; when there are none, compute VC,
    // VC_SQRT, KR and SCALE as Juicer tools pre does. "none": no
    // normalization. Otherwise a comma separated list of the normalizations
    // to compute, out of VC, VC_SQRT, KR and SCALE.
    std::string normalizations = "auto";
    // When set, cooler's multiplicative "weight" column is written as the
    // divisive normalization vector of this name (1 / weight).
    std::optional<std::string> cooler_weight;
    // The genome id of the .hic header; empty takes the cool file's
    // genome-assembly attribute, or "unknown".
    std::string genome;
    // For square coolers (both triangles stored), which triangle becomes the
    // .hic file's: "auto" checks that the cooler is symmetric and fails when
    // it is not; "upper" or "lower" takes that triangle as it is.
    std::string triangle = "auto";
    bool show_warnings = false;
    bool silent = false;
};

// Writes a .hic file (version 8 or 9) from a .cool or .mcool file or a URI
// "file.mcool::/resolutions/10000". Returns the path written.
std::string cool2hic_convert(const std::string& infile, const std::string& outfile,
                             const Cool2hicOptions& options = {},
                             const Console& console = Console::standard());

}  // namespace hic2cool

#endif  // HIC2COOL_HIC2COOL_HPP
