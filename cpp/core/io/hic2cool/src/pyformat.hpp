// The Python text the messages of hic2cool are made of: repr of strings,
// lists, dicts and floats, datetime.utcnow().isoformat(), and the pandas
// DataFrame preview extract-norms prints.

#ifndef HIC2COOL_PYFORMAT_HPP
#define HIC2COOL_PYFORMAT_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace hic2cool::detail {

// repr(str)
[[nodiscard]] std::string py_repr(const std::string& text);
// repr(list) of str and of int
[[nodiscard]] std::string py_repr(const std::vector<std::string>& values);
[[nodiscard]] std::string py_repr(const std::vector<std::int32_t>& values);
[[nodiscard]] std::string py_repr(const std::vector<std::int64_t>& values);
// repr(dict) of int -> str
[[nodiscard]] std::string py_repr(const std::vector<std::pair<std::int64_t, std::string>>& values);
// repr(float)
[[nodiscard]] std::string py_repr(double value);

// datetime.utcnow().isoformat()
[[nodiscard]] std::string utcnow_isoformat();

// str(df.head()) for the bins table extract-norms prints: the columns chrom,
// start, end and one float64 column named `norm`, first rows only.
struct BinsPreviewRow {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    double value = 0.0;
};
[[nodiscard]] std::string pandas_bins_head(const std::string& norm, const std::vector<BinsPreviewRow>& rows);

}  // namespace hic2cool::detail

#endif  // HIC2COOL_PYFORMAT_HPP
