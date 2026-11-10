// What hic2cool's read_header and read_footer take from a .hic file, and the
// messages they print.

#ifndef HIC2COOL_HIC_HEADER_HPP
#define HIC2COOL_HIC_HEADER_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <hicfilecpp/hicfilecpp.hpp>

#include "hic2cool/hic2cool.hpp"

namespace hic2cool::detail {

struct HicChrom {
    std::int32_t index = 0;
    std::string name;
    std::int64_t length = 0;
};

struct HicHeader {
    // read_header's chrs: chromosomes with a name and a non-zero length,
    // "All" included, in file order.
    std::vector<HicChrom> used;
    std::vector<std::int32_t> resolutions;
    std::string genome;
    // The header attributes as a Python dict: first position, last value.
    std::vector<std::pair<std::string, std::string>> metadata;
    // read_footer's NORMS: the normalization types of the normalized
    // expected values, in file order.
    std::vector<std::string> norms;
};

[[nodiscard]] std::string lowercase(std::string text);

// Opens and reads the file as read_header and read_footer do, printing the
// footer's warning. Throws ExitError for a file that is not a .hic file.
HicHeader read_hic_header(const hicfilecpp::HiCFile& hic, const Console& console);

// The magic string check of read_header.
void check_hic_magic(const std::string& path);

}  // namespace hic2cool::detail

#endif  // HIC2COOL_HIC_HEADER_HPP
