// The capture Hi-C viewpoint machinery, ported from hicexplorer/lib/viewpoint.py.
//
// This is a reusable core component. chicQualityControl,
// chicViewpointBackgroundModel and chicViewpoint are built on it, and so will be
// chicSignificantInteractions, chicAggregateStatistic, chicDifferentialTest and
// chicExportData. Everything here is written against the Python method it
// replaces, and it reproduces that method's behaviour including the parts that
// look wrong; each of those is named where it happens, and none is fixed.
//
// Python raises three kinds of exception out of this code, and the tools branch
// on which one it was (chicQualityControl turns a TypeError or an IndexError
// into a sparsity of -1 and lets anything else abort the run). The three C++
// exception types below carry that distinction.
//
// Units. A reference point, a region and a range are genomic positions; a
// viewpoint's data are indexed by bin. The background model file is keyed by
// genomic distance. Viewpoint.pvalues looks it up by bin distance, which is one
// of the pinned defects: see p_values().
//
// Determinism: nothing here is threaded or holds global state. The tools run
// independent reference points concurrently and combine the results in file
// order.

#ifndef HICX_CHIC_VIEWPOINT_HPP
#define HICX_CHIC_VIEWPOINT_HPP

#include <cstdint>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/hic_matrix.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx::chic {

// ---------------------------------------------------------------------------
// Python exception categories

class TypeError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};
class IndexError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};
class ValueError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// ---------------------------------------------------------------------------
// Text helpers with Python's semantics

// open(path).readlines(): universal newlines, so "\r\n" and a lone "\r" both
// end a line and are delivered as "\n". Each line keeps its newline.
[[nodiscard]] std::vector<std::string> read_lines(const std::string& path);

// str.strip() for the ASCII whitespace Python recognises.
[[nodiscard]] std::string_view strip(std::string_view text);

// int(text) for a decimal string: surrounding whitespace, a sign and single
// underscores between digits are accepted. Throws ValueError otherwise.
[[nodiscard]] std::int64_t python_int(std::string_view text);

// float(text), including 'nan', 'inf' and surrounding whitespace. Throws
// ValueError.
[[nodiscard]] double python_float(std::string_view text);

// '{:.12f}'.format(value) and friends: correctly rounded fixed notation, with
// Python's spelling of 'nan', 'inf' and '-inf'.
[[nodiscard]] std::string format_fixed(double value, int decimals);

// ---------------------------------------------------------------------------
// Reference points

// One line of a reference point file, as strings: Python keeps them as text
// and converts with int() only where it uses them, so a malformed number is a
// ValueError at that point rather than when the file is read.
struct ReferencePoint {
    std::string chromosome;
    std::string start;
    std::string end;
};

struct ReferencePoints {
    std::vector<ReferencePoint> points;
    std::vector<std::string> genes;
};

// Viewpoint.readReferencePointFile(pBedFile, pGene=True). A line of three
// tab separated fields is (chrom, start, start) with the gene in the third;
// four or more are (chrom, start, end) with the gene in the fourth; anything
// shorter is skipped.
[[nodiscard]] ReferencePoints read_reference_points(const std::string& path);

// '_'.join(str(j) for j in referencePoint)
[[nodiscard]] std::string reference_point_string(const ReferencePoint& point);

// ---------------------------------------------------------------------------
// The matrix a cHi-C tool loads

class ViewpointMatrix {
  public:
    // hm.hiCMatrix(path): symmetric, balancing weights applied when present.
    static ViewpointMatrix load(const std::string& path);

    [[nodiscard]] const BinTable& bins() const noexcept { return matrix_.bins(); }
    [[nodiscard]] const CsrMatrix& matrix() const noexcept { return matrix_.matrix(); }

    // getBinSize(): the median start difference, truncated.
    [[nodiscard]] std::int64_t bin_size() const { return matrix_.bin_size(); }

    // getChrNames(), in bin table order.
    [[nodiscard]] std::vector<std::string> chromosome_names() const;
    [[nodiscard]] bool has_chromosome(const std::string& chromosome) const;

    // getRegionBinRange(chrom, start, end): the bins holding start and end, or
    // nullopt when either position is outside the chromosome's bins (the
    // Python logs an IndexError and returns None). An unknown chromosome raises
    // ValueError.
    [[nodiscard]] std::optional<std::pair<std::int64_t, std::int64_t>> region_bin_range(
        const std::string& chromosome, std::int64_t start, std::int64_t end) const;

    // getReferencePointAsMatrixIndices. A None from the region lookup is
    // unpacked by the caller, which is a TypeError.
    [[nodiscard]] std::pair<std::int64_t, std::int64_t> reference_point_indices(
        const ReferencePoint& point) const;

    // getBinPos(index): (chrom, start, end). ValueError past the last bin.
    [[nodiscard]] const CutInterval& bin_position(std::int64_t index) const;

    // self.hicMatrix.matrix[row, column] on the symmetric matrix.
    [[nodiscard]] double value(std::int64_t row, std::int64_t column) const;

  private:
    HiCMatrix matrix_;
};

// Viewpoint.calculateViewpointRange(pViewpoint, pRange).
struct ViewpointRange {
    std::int64_t region_start = 0;
    std::int64_t region_end = 0;
    // _range, as adjusted when the region is clipped to the chromosome.
    std::int64_t upstream = 0;
    std::int64_t downstream = 0;
};

// The region is the reference point widened by the range, clipped to 0 on the
// left and to one base before the chromosome end on the right. The clipped
// downstream range is (chromosome end - reference end) + bin size, which is
// not a multiple of the bin size in general; that is why a reference point
// near a chromosome end can produce a background of a different length than
// its data (see chicViewpoint's adjustViewpointData). An unknown chromosome
// raises ValueError.
[[nodiscard]] ViewpointRange calculate_viewpoint_range(const ViewpointMatrix& matrix,
                                                       const ReferencePoint& point,
                                                       std::int64_t upstream,
                                                       std::int64_t downstream);

// Viewpoint.computeViewpoint: the rows of the reference point's bins summed
// over the region's columns, then the reference point's own bins collapsed
// into one element. `index_before_viewpoint` is the position of that element.
//
// The three slice assignments that do the collapsing follow numpy's rules, so
// a reference point whose start lies after its end, which makes the lengths
// inconsistent, raises ValueError ("could not broadcast") exactly where numpy
// does, and a region lookup that returns None raises TypeError.
struct ComputedViewpoint {
    std::vector<double> data;
    std::int64_t index_before_viewpoint = 0;
};
[[nodiscard]] ComputedViewpoint compute_viewpoint(const ViewpointMatrix& matrix,
                                                  const ReferencePoint& point,
                                                  const std::string& chromosome,
                                                  std::int64_t region_start,
                                                  std::int64_t region_end);

// Viewpoint.smoothInteractionValues(pData, pWindowSize).
//
// Every element is the mean of a window around it. For an odd window the
// window is symmetric. For an even window it holds one element fewer upstream
// than downstream, except at the right border, where viewpoint.py:561 averages
// pData[-(i + half + 1):] and so takes the full half window upstream. Means are
// numpy's: a pairwise sum divided by the count. The float32 overload computes
// them in single precision, as np.mean does for a float32 slice, and returns
// them widened into the float64 result array.
[[nodiscard]] std::vector<double> smooth_interaction_values(std::span<const double> data,
                                                            std::int64_t window);
[[nodiscard]] std::vector<double> smooth_interaction_values(std::span<const float> data,
                                                            std::int64_t window);

// Viewpoint.computeRelativeValues(pData, pDenominator). `if pDenominator:`
// treats 0.0 as absent, and then the sum of the data is used instead.
[[nodiscard]] std::vector<double> compute_relative_values(std::span<const double> data,
                                                          double denominator);

// ---------------------------------------------------------------------------
// The background model file

// Viewpoint.readBackgroundDataFile: {relative genomic position: values}, in
// the dict's insertion order, which is the file order followed by the
// positions added to cover the range.
class BackgroundModel {
  public:
    void set(std::int64_t key, std::vector<double> values);
    [[nodiscard]] bool contains(std::int64_t key) const;
    // KeyError in Python; a ValueError here, since no caller distinguishes it.
    [[nodiscard]] const std::vector<double>& at(std::int64_t key) const;
    [[nodiscard]] std::int64_t min_key() const;
    [[nodiscard]] std::int64_t max_key() const;
    [[nodiscard]] std::vector<std::int64_t> sorted_keys() const;
    [[nodiscard]] const std::vector<std::int64_t>& insertion_order() const noexcept {
        return order_;
    }
    [[nodiscard]] std::size_t size() const noexcept { return order_.size(); }

  private:
    std::map<std::int64_t, std::vector<double>> values_;
    std::vector<std::int64_t> order_;
};

// pMean=False keeps [size, prob, max value]; pMean=True keeps [mean value].
// Positions beyond --fixateRange are clamped for the extension; the model is
// then extended in steps of |abs(first key) - abs(second key)| until it covers
// pRange on both sides, repeating the value at the clamped end.
[[nodiscard]] BackgroundModel read_background_model(const std::string& path,
                                                    std::int64_t range_upstream,
                                                    std::int64_t range_downstream,
                                                    std::int64_t fixate_range, bool mean);

// Viewpoint.interactionBackgroundData(pBackground, pRange) for a mean model:
// the values at the sorted keys within [-upstream, downstream], flattened.
[[nodiscard]] std::vector<double> interaction_background_data(const BackgroundModel& model,
                                                              std::int64_t upstream,
                                                              std::int64_t downstream);

// Viewpoint.pvalues(pBackgroundModel, pDataList, pIndexReferencePoint):
// 1 - cnb.cdf(x, size, prob) = 1 - betainc(size, x + 1, prob) per element, 1.0
// for a zero, and NaN or infinity replaced by 1.0. betainc is
// hicx::scipy::betainc, Boost's ibeta as scipy 1.14 evaluates it, not the
// cephes translation in stats_ops: in the upper tail the two differ by the one
// ulp that decides whether a stored p-value is 0.0.
//
// The distribution is looked up by `i - pIndexReferencePoint`, a distance in
// bins, in a model keyed by genomic distance (chicViewpointBackgroundModel.py
// writes relative_position * bin_size). A missing key falls back to the
// smallest key upstream and the largest downstream. With bins larger than 1 bp
// the only bin distance that is also a key is 0, so every position except the
// reference point itself is tested against the distribution at -fixateRange
// or +fixateRange. Reproduced, and reported as a defect in the Python.
//
// `data_is_float32` marks data that numpy holds as float32 (chicViewpoint's
// adjustViewpointData without smoothing), for which x + 1 is evaluated in
// single precision before betainc widens it.
[[nodiscard]] std::vector<double> p_values(const BackgroundModel& model,
                                           std::span<const double> data,
                                           std::int64_t index_reference_point,
                                           bool data_is_float32 = false);

// ---------------------------------------------------------------------------
// One viewpoint's content in the interaction file

// Viewpoint.createInteractionFileDataHDF5.
struct InteractionFileData {
    std::string chromosome;
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> ends;
    std::string gene;
    double sum_of_interactions = 0.0;
    std::vector<std::int64_t> relative_positions;
    std::vector<double> interaction_data;
    std::vector<double> pvalues;
    std::vector<double> xfold;
    std::vector<double> raw;
};

// The positions are the region's bins up to the reference point, the reference
// point's first bin, then the bins after its last one, zipped with the data,
// so the shorter of the two decides the length. A relative position is
// start - reference start until the first non negative one, and end -
// reference end from then on, which puts the collapsed element at 0 and gives
// it the coordinates of the reference point's first bin only.
[[nodiscard]] InteractionFileData create_interaction_file_data(
    const ViewpointMatrix& matrix, const ReferencePoint& point, const std::string& chromosome,
    std::int64_t region_start, std::int64_t region_end,
    std::span<const double> interaction_data, std::span<const double> raw,
    const std::string& gene, double sum_of_interactions, std::vector<double> pvalues,
    std::vector<double> xfold);

}  // namespace hicx::chic

#endif  // HICX_CHIC_VIEWPOINT_HPP
