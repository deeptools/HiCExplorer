// Internal: port of cooler/core/_rangequery.py (cooler 0.10.2, BSD-3-Clause):
// CSRReader, DirectRangeQuery2D, FillLowerRangeQuery2D and the conversions of
// their result dictionaries into frames, sparse and dense matrices.
//
// cooler materialises one dict of arrays per task and concatenates them. Here
// every task appends to one PixelDict whose buffers are reserved up front for
// the largest possible result (the pixel count of the row spans, doubled when
// the lower triangle is filled). Untouched reserved pages never become
// resident, so the peak RSS is the size of the actual result, without the
// per-task copies and the concatenation.

#ifndef COOLERCPP_DETAIL_RANGEQUERY_HPP
#define COOLERCPP_DETAIL_RANGEQUERY_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "coolercpp/api.hpp"
#include "h5.hpp"

namespace coolercpp::detail {

struct Bbox {
    std::int64_t i0 = 0;
    std::int64_t i1 = 0;
    std::int64_t j0 = 0;
    std::int64_t j1 = 0;
};

// Bin IDs of a query result: int32 storage whenever every coordinate of the
// query fits, int64 otherwise. The numpy dtype the Python would give the
// column is tracked separately in PixelDict.
class IdArray {
  public:
    explicit IdArray(bool narrow = true) : narrow_(narrow) {}

    [[nodiscard]] bool narrow() const noexcept { return narrow_; }
    [[nodiscard]] std::size_t size() const noexcept { return narrow_ ? n_.size() : w_.size(); }
    [[nodiscard]] std::int64_t operator[](std::size_t k) const {
        return narrow_ ? static_cast<std::int64_t>(n_[k]) : w_[k];
    }
    void reserve(std::size_t n) {
        if (narrow_) {
            n_.reserve(n);
        } else {
            w_.reserve(n);
        }
    }
    [[nodiscard]] std::vector<std::int32_t>& narrow_values() noexcept { return n_; }
    [[nodiscard]] std::vector<std::int64_t>& wide_values() noexcept { return w_; }
    void release() {
        std::vector<std::int32_t>().swap(n_);
        std::vector<std::int64_t>().swap(w_);
    }

  private:
    bool narrow_;
    std::vector<std::int32_t> n_;
    std::vector<std::int64_t> w_;
};

// The dict of arrays a query returns.
struct PixelDict {
    IdArray bin1;
    IdArray bin2;
    DType bin1_dtype = DType::Int64;
    DType bin2_dtype = DType::Int64;
    bool has_parts = false;
    Column value;
    std::optional<std::vector<std::int64_t>> index;
};

class CSRReader {
  public:
    CSRReader(const h5::File& file, std::string pixels_group,
              std::vector<std::int64_t> bin1_offsets);

    [[nodiscard]] std::vector<std::pair<std::int64_t, std::int64_t>> get_spans(
        const Bbox& bbox, std::int64_t chunksize) const;
    [[nodiscard]] PixelDict meta(const std::string& field, bool return_index) const;
    [[nodiscard]] DType dtype_of(const std::string& column) const;
    // Pixels in the stored rows of a span: an upper bound of what a task
    // selects from it in ordinary queries.
    [[nodiscard]] std::size_t span_pixels(std::pair<std::int64_t, std::int64_t> span) const;
    [[nodiscard]] bool narrow_ids(const Bbox& bbox) const;

    // One CSRReader call (cooler's __call__), appended to `out`; `transpose`
    // swaps the two ID columns of this part, as FillLowerRangeQuery2D does.
    void read_into(PixelDict& out, const std::string& field, const Bbox& bbox,
                   std::pair<std::int64_t, std::int64_t> span, bool reflect, bool return_index,
                   bool transpose) const;

  private:
    const h5::File& file_;
    std::string pixels_;
    std::vector<std::int64_t> offsets_;
    std::vector<std::pair<std::string, DType>> dtypes_;
};

// One task of a range query engine: a CSRReader call over a row span of a
// bounding box, optionally transposed.
struct Task {
    Bbox bbox;
    std::pair<std::int64_t, std::int64_t> span;
    bool transpose = false;
};

// The task lists of DirectRangeQuery2D and FillLowerRangeQuery2D.
[[nodiscard]] std::vector<Task> direct_tasks(const CSRReader& reader, const Bbox& bbox,
                                             std::int64_t chunksize);
[[nodiscard]] std::vector<Task> fill_lower_tasks(const CSRReader& reader, const Bbox& bbox,
                                                 std::int64_t chunksize);
// Runs tasks in order into one result (reflect: the lower triangle filling
// of FillLowerRangeQuery2D).
[[nodiscard]] PixelDict run_tasks(const CSRReader& reader, const std::string& field,
                                  const Bbox& query, const std::vector<Task>& tasks,
                                  bool reflect, bool return_index);

[[nodiscard]] PixelDict direct_query(const CSRReader& reader, const std::string& field,
                                     const Bbox& bbox, std::int64_t chunksize,
                                     bool return_index);
[[nodiscard]] PixelDict fill_lower_query(const CSRReader& reader, const std::string& field,
                                         const Bbox& bbox, std::int64_t chunksize,
                                         bool return_index);

[[nodiscard]] Table frame_from_dict(PixelDict&& dict, const std::string& field);
[[nodiscard]] SparseMatrix sparse_from_dict(PixelDict&& dict, const Bbox& bbox);
[[nodiscard]] DenseMatrix dense_from_dict(const PixelDict& dict, const Bbox& bbox);

}  // namespace coolercpp::detail

#endif  // COOLERCPP_DETAIL_RANGEQUERY_HPP
