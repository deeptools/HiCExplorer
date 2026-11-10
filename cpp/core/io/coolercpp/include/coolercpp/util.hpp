// cooler.util helpers that operate on bin tables.

#ifndef COOLERCPP_UTIL_HPP
#define COOLERCPP_UTIL_HPP

#include <string>
#include <vector>

#include "coolercpp/json.hpp"
#include "coolercpp/region.hpp"
#include "coolercpp/table.hpp"

namespace coolercpp {

// cooler.util.get_binsize: the single width shared by every bin except the
// last one of each chromosome, as a number carrying the dtype of end - start;
// null when the widths differ or no chromosome has more than one bin.
[[nodiscard]] json::Value get_binsize(const Table& bins);

// cooler.util.get_chromsizes: the end of the last bin of each chromosome, in
// the order of those last bins.
[[nodiscard]] ChromSizes get_chromsizes(const Table& bins);

// cooler.util.natsorted.
[[nodiscard]] std::vector<std::string> natsorted(std::vector<std::string> items);

}  // namespace coolercpp

#endif  // COOLERCPP_UTIL_HPP
