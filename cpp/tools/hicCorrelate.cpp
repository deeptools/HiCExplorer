// Port of hicexplorer/hicCorrelate.py (cpp/PLAN.md tier 7, option (a)).
//
// The matrices are read and reduced here: the chromosome selection, the
// inter-chromosomal and diagonal counts removed, the upper triangle cut to
// --range, log1p, the NaN bins of every matrix removed, one vector per matrix
// over the union of stored positions, the pairwise Pearson or Spearman
// correlations, and the complete linkage clustering of the correlation
// matrix. plot/hicexplorer_plot/hicCorrelate.py draws the scatter grid and the
// heatmap with the reference's calls; the dendrogram is drawn from the linkage
// computed here.
//
// What is reproduced, pinned by the harness cases:
//
//  1. The value types of scipy.sparse: integer matrices stay integers (log1p
//     makes them float64), float32 matrices are reduced and correlated in
//     float32, float64 in float64.
//  2. The per matrix vector is (mat + all_mat).data - all_mat.data, so a
//     stored value v comes back as (v + a) - a with the rounding of that
//     arithmetic, where a is the sum over the matrices in command line order.
//  3. The pairs keep the positions whose two values sum to more than 1.
//  4. pearsonr as scipy 1.14.1 computes it, with numpy's pairwise sums: the
//     mean, the scaled norm xmax * ||xm / xmax||, sum(xm / normxm * ym /
//     normym), the clip to [-1, 1], NaN for a constant input and the sign rule
//     for two points. spearmanr: average ranks, then np.corrcoef's covariance
//     and normalisation order. The covariance goes through a BLAS dot product
//     in numpy, whose summation order is not reproduced, so Spearman values
//     agree at E3.
//  5. linkage(method='complete') and the dendrogram leaf order through
//     hicx::cluster (pdist, nn_chain, label).
//
// Memory: the reference holds every reduced matrix and then the vectors. Here
// a matrix is freed as soon as its stored values are extracted, and its
// values are merged into the union at once: both key lists are walked from
// the back into the enlarged arrays, so earlier matrices are never copied.
// The working set is the keys and one float64 per matrix per position, plus
// the matrix being read. The vectors reach the drawing process as a raw .npy
// written straight from that array.
//
// Deliberate deviations:
//
//  * --labels with a different count than --matrices, and a --range whose
//    upper end is below the bin size, make the reference log an error and
//    exit 0 without writing anything. The port exits 1 (cpp/AGENTS_CONTRACT.md
//    rule 7).
//  * A non-finite value in the vectors (a NaN or infinite count outside the
//    NaN bins) is refused with exit 1. The reference masks it with
//    np.ma.masked_invalid, which changes the pair selection in ways that
//    depend on numpy's masked reductions; no input in the corpus has one.
//
// Mixed value types across matrices (float32 with float64) are accumulated in
// the promoted type from the start, where scipy promotes at the first mixed
// addition; the vectors can then differ in the last bits (E3).

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/clustering.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/npz_file.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kDescription =
    "Computes pairwise correlations between Hi-C matrices data. The correlation is computed "
    "taking the values from each pair of matrices and discarding values that are zero in both "
    "matrices.Parameters that strongly affect correlations are bin size of the Hi-C matrices and "
    "the considered range. The smaller the bin size of the matrices, the finer differences you "
    "score. The --range parameter should be selected at a meaningful genomic scale according to, "
    "for example, the mean size of the TADs in the organism you work with.";

const char* const kUsage =
    "usage: hicCorrelate --matrices MATRICES [MATRICES ...] [--zMin ZMIN]\n"
    "                    [--zMax ZMAX] [--colorMap] [--plotNumbers]\n"
    "                    [--method {pearson,spearman}] [--log1p]\n"
    "                    [--labels sample1 sample2 [sample1 sample2 ...]]\n"
    "                    [--range RANGE] --outFileNameHeatmap OUTFILENAMEHEATMAP\n"
    "                    --outFileNameScatter OUTFILENAMESCATTER\n"
    "                    [--chromosomes CHROMOSOMES [CHROMOSOMES ...]] [--help]\n"
    "                    [--version]\n";

const char* const kHelp =
    "\n"
    "Computes pairwise correlations between Hi-C matrices data. The correlation is\n"
    "computed taking the values from each pair of matrices and discarding values\n"
    "that are zero in both matrices.Parameters that strongly affect correlations\n"
    "are bin size of the Hi-C matrices and the considered range. The smaller the\n"
    "bin size of the matrices, the finer differences you score. The --range\n"
    "parameter should be selected at a meaningful genomic scale according to, for\n"
    "example, the mean size of the TADs in the organism you work with.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        Matrices to correlate (usually .h5 but other formats\n"
    "                        are allowed). hicCorrelate is better used on un-\n"
    "                        corrected matrices in order to exclude any changes\n"
    "                        introduced by the correction. (default: None)\n"
    "\n"
    "Heatmap arguments:\n"
    "  Options for generating the correlation heatmap\n"
    "\n"
    "  --zMin ZMIN, -min ZMIN\n"
    "                        Minimum value for the heatmap intensities. If not\n"
    "                        specified the value is set automatically. (default:\n"
    "                        None)\n"
    "  --zMax ZMAX, -max ZMAX\n"
    "                        Maximum value for the heatmap intensities.If not\n"
    "                        specified the value is set automatically. (default:\n"
    "                        None)\n"
    "  --colorMap            Color map to use for the heatmap. Available values can\n"
    "                        be seen here: http://matplotlib.org/examples/color/col\n"
    "                        ormaps_reference.html (Default: jet).\n"
    "  --plotNumbers         If set, then the correlation number is plotted on top\n"
    "                        of the heatmap. (default: False)\n"
    "\n"
    "Optional arguments:\n"
    "  --method {pearson,spearman}\n"
    "                        Correlation method to use (Default: pearson).\n"
    "  --log1p               If set, then the log1p of the matrix values is used.\n"
    "                        This parameter has no effect for Spearman correlations\n"
    "                        but changes the output of Pearson correlation and, for\n"
    "                        the scatter plot, if set, the visualization of the\n"
    "                        values is easier. (default: False)\n"
    "  --labels sample1 sample2 [sample1 sample2 ...], -l sample1 sample2 [sample1 sample2 ...]\n"
    "                        User defined labels instead of default labels from\n"
    "                        file names. Multiple labels have to be separated by\n"
    "                        space, e.g. --labels sample1 sample2 sample3 (default:\n"
    "                        None)\n"
    "  --range RANGE         In bp with the format low_range:high_range, for\n"
    "                        example 1000000:2000000. If --range is given only\n"
    "                        counts within this range are considered. The range\n"
    "                        should be adjusted to the size of interacting domains\n"
    "                        in the genome you are working with. (default: None)\n"
    "  --outFileNameHeatmap OUTFILENAMEHEATMAP, -oh OUTFILENAMEHEATMAP\n"
    "                        File name to save the resulting heatmap plot.\n"
    "                        Supported file formats are given by matplotlib,\n"
    "                        usually these are: png, pdf, ps, eps and svg.\n"
    "                        (default: heatmap.png)\n"
    "  --outFileNameScatter OUTFILENAMESCATTER, -os OUTFILENAMESCATTER\n"
    "                        File name to save the resulting scatter plot.\n"
    "                        Supported file formats are given by matplotlib,\n"
    "                        usually these are: png, pdf, ps, eps and svg.\n"
    "                        (default: scatter.png)\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        List of chromosomes to be included in the correlation.\n"
    "                        (default: None)\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the vectors, correlations and clustering are computed in C++, and the\n"
    "figures are drawn by the hicexplorer_plot drawing layer with the calls of the\n"
    "Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option\n"
    "--plotData FILE writes the data of the figures as JSON to FILE (and the vectors\n"
    "to FILE.npy) instead of drawing them.\n";

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

enum class Kind { Int, F32, F64 };

Kind promote(Kind a, Kind b) {
    return a == b ? a : Kind::F64;
}

// One matrix after the reference's reduction: the stored positions of the
// upper triangle as row * size + column, ascending, with their values.
struct SampleEntries {
    std::int64_t size = 0;
    Kind kind = Kind::F64;
    std::vector<std::int64_t> keys;
    std::vector<double> values;
    std::vector<std::int64_t> nan_bins;
};

bool contains(const std::vector<std::string>& list, const std::string& value) {
    return std::find(list.begin(), list.end(), value) != list.end();
}

SampleEntries load_sample(const std::string& path,
                          const std::optional<std::vector<std::string>>& chromosomes,
                          const std::optional<std::string>& range, bool log1p) {
    const bool single_chromosome_cool =
        hicx::check_cooler(path) && chromosomes.has_value() && chromosomes->size() == 1;
    const hicx::ToolMatrix matrix =
        single_chromosome_cool ? hicx::ToolMatrix::load(path, chromosomes->front())
                               : hicx::ToolMatrix::load(path);
    const std::vector<hicx::CutInterval>& all = matrix.cut_intervals();

    // keepOnlyTheseChr(chromosomes), then filterOutInterChrCounts.
    std::vector<std::int64_t> new_index(all.size(), -1);
    std::vector<hicx::CutInterval> kept;
    if (!single_chromosome_cool && chromosomes.has_value()) {
        for (const std::string& name : *chromosomes) {
            if (std::none_of(matrix.boundaries().begin(), matrix.boundaries().end(),
                             [&](const auto& entry) { return entry.first == name; })) {
                throw PythonError("ValueError: Chromosome name not in matrix. '" + name + "'");
            }
        }
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (contains(*chromosomes, all[i].chrom)) {
                new_index[i] = static_cast<std::int64_t>(kept.size());
                kept.push_back(all[i]);
            }
        }
    } else {
        kept = all;
        for (std::size_t i = 0; i < all.size(); ++i) {
            new_index[i] = static_cast<std::int64_t>(i);
        }
    }
    const bool filter_inter = !single_chromosome_cool;

    SampleEntries out;
    out.size = static_cast<std::int64_t>(kept.size());
    const hicx::CsrMatrix& csr = matrix.matrix();
    out.kind = csr.dtype_kind() == hicx::DType::Integer   ? Kind::Int
               : csr.dtype_kind() == hicx::DType::Float32 ? Kind::F32
                                                          : Kind::F64;
    for (const std::int64_t bin : matrix.nan_bins()) {
        if (bin >= 0 && static_cast<std::size_t>(bin) < new_index.size() &&
            new_index[static_cast<std::size_t>(bin)] >= 0) {
            out.nan_bins.push_back(new_index[static_cast<std::size_t>(bin)]);
        }
    }

    const hicx::BinTable table(kept);
    const std::int64_t bin_size = table.bin_size();
    std::optional<std::pair<std::int64_t, std::int64_t>> distances;  // [min, max] in bins
    std::int64_t max_depth_in_bins = 0;
    if (range.has_value()) {
        const std::string::size_type colon = range->find(':');
        if (colon == std::string::npos || range->find(':', colon + 1) != std::string::npos) {
            throw PythonError("ValueError: --range needs exactly two values separated by ':'");
        }
        std::int64_t min_dist = 0;
        std::int64_t max_dist = 0;
        if (!cli::python_int(range->substr(0, colon), &min_dist) ||
            !cli::python_int(range->substr(colon + 1), &max_dist)) {
            throw PythonError("ValueError: invalid literal for int() with base 10 in --range " +
                              *range);
        }
        if (max_dist < bin_size) {
            std::fprintf(stderr,
                         "Please specify a max range that is larger than bin size (%lld). The "
                         "reference exits 0 here without writing anything; the C++ port exits 1 "
                         "(cpp/AGENTS_CONTRACT.md rule 7).\n",
                         static_cast<long long>(bin_size));
            std::exit(1);
        }
        max_depth_in_bins =
            static_cast<std::int64_t>(static_cast<double>(max_dist) / static_cast<double>(bin_size));
        auto floor_div = [](std::int64_t a, std::int64_t b) {
            std::int64_t q = a / b;
            return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
        };
        distances = std::make_pair(floor_div(min_dist, bin_size), floor_div(max_dist, bin_size));
    }

    const auto& indptr = csr.indptr();
    const auto& indices = csr.indices();
    const auto& data = csr.data();
    std::vector<std::pair<std::int64_t, double>> row_entries;
    for (std::int64_t row = 0; row < csr.rows(); ++row) {
        const std::int64_t new_row = new_index[static_cast<std::size_t>(row)];
        if (new_row < 0) {
            continue;
        }
        row_entries.clear();
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t col = indices[k];
            if (col <= row || data[k] == 0.0) {
                continue;  // diagflat(0) and triu(k=0)
            }
            const std::int64_t new_col = new_index[static_cast<std::size_t>(col)];
            if (new_col < 0) {
                continue;
            }
            if (filter_inter && kept[static_cast<std::size_t>(new_row)].chrom !=
                                    kept[static_cast<std::size_t>(new_col)].chrom) {
                continue;
            }
            const std::int64_t dist = new_col - new_row;
            if (distances.has_value() &&
                (dist >= max_depth_in_bins || dist > distances->second || dist < distances->first)) {
                continue;
            }
            row_entries.emplace_back(new_row * out.size + new_col, data[k]);
        }
        std::sort(row_entries.begin(), row_entries.end());
        for (const auto& [key, value] : row_entries) {
            out.keys.push_back(key);
            out.values.push_back(value);
        }
    }
    if (log1p) {
        if (out.kind == Kind::F32) {
            for (double& value : out.values) {
                value = static_cast<double>(std::log1p(static_cast<float>(value)));
            }
        } else {
            for (double& value : out.values) {
                value = std::log1p(value);
            }
            out.kind = Kind::F64;
        }
    }
    return out;
}

// Merges one matrix into the union: keys ascending, big with `columns` values
// per key, grown to columns + 1. Both lists are walked from the back into the
// enlarged arrays; a destination is never before its source, so nothing is
// overwritten before it is read.
void merge_sample(std::vector<std::int64_t>& keys, std::vector<double>& big, std::size_t columns,
                  std::size_t total_columns, const SampleEntries& sample) {
    std::size_t count = 0;
    {
        std::size_t i = 0;
        std::size_t j = 0;
        while (i < keys.size() || j < sample.keys.size()) {
            if (j == sample.keys.size() || (i < keys.size() && keys[i] < sample.keys[j])) {
                ++i;
            } else if (i == keys.size() || sample.keys[j] < keys[i]) {
                ++j;
            } else {
                ++i;
                ++j;
            }
            ++count;
        }
    }
    const std::size_t width = columns + 1;
    std::ptrdiff_t i = static_cast<std::ptrdiff_t>(keys.size()) - 1;
    std::ptrdiff_t j = static_cast<std::ptrdiff_t>(sample.keys.size()) - 1;
    // Room for every matrix at once, so that the later merges of matrices that
    // store the same positions grow in place instead of reallocating, which
    // would hold the old and the new array together.
    if (keys.capacity() < count) {
        keys.reserve(count);
    }
    if (big.capacity() < count * width) {
        big.reserve(count * total_columns);
    }
    keys.resize(count);
    big.resize(count * width);
    for (std::ptrdiff_t out = static_cast<std::ptrdiff_t>(count) - 1; out >= 0; --out) {
        bool take_old = false;
        bool take_new = false;
        if (i >= 0 && j >= 0) {
            const std::int64_t old_key = keys[static_cast<std::size_t>(i)];
            const std::int64_t new_key = sample.keys[static_cast<std::size_t>(j)];
            take_old = old_key >= new_key;
            take_new = new_key >= old_key;
        } else {
            take_old = i >= 0;
            take_new = j >= 0;
        }
        const std::int64_t key = take_old ? keys[static_cast<std::size_t>(i)]
                                          : sample.keys[static_cast<std::size_t>(j)];
        const auto destination = static_cast<std::size_t>(out) * width;
        for (std::size_t c = columns; c-- > 0;) {
            big[destination + c] =
                take_old ? big[static_cast<std::size_t>(i) * columns + c] : 0.0;
        }
        big[destination + columns] = take_new ? sample.values[static_cast<std::size_t>(j)] : 0.0;
        keys[static_cast<std::size_t>(out)] = key;
        if (take_old) {
            --i;
        }
        if (take_new) {
            --j;
        }
    }
}

template <typename T>
T pairwise(const std::vector<T>& values) {
    return hicx::npy::pairwise_sum(values);
}

template <typename T>
T sign(T value) {
    return value > 0 ? T(1) : (value < 0 ? T(-1) : T(0));
}

// scipy.stats.pearsonr(x, y)[0] in the floating type T.
template <typename T>
double pearson(const std::vector<T>& x, const std::vector<T>& y) {
    const std::size_t n = x.size();
    if (n < 2) {
        throw PythonError("ValueError: `x` and `y` must have length at least 2.");
    }
    if (n == 2) {
        return static_cast<double>(sign<T>(x[1] - x[0]) * sign<T>(y[1] - y[0]));
    }
    const bool const_x = std::all_of(x.begin(), x.end(), [&](T v) { return v == x[0]; });
    const bool const_y = std::all_of(y.begin(), y.end(), [&](T v) { return v == y[0]; });
    auto centred_norm = [n](const std::vector<T>& v, std::vector<T>* centred) {
        const T mean = pairwise(v) / static_cast<T>(n);
        centred->resize(n);
        T largest = 0;
        for (std::size_t i = 0; i < n; ++i) {
            (*centred)[i] = v[i] - mean;
            largest = std::max(largest, std::abs((*centred)[i]));
        }
        std::vector<T> squares(n);
        for (std::size_t i = 0; i < n; ++i) {
            const T scaled = (*centred)[i] / largest;
            squares[i] = scaled * scaled;
        }
        return largest * std::sqrt(pairwise(squares));
    };
    std::vector<T> xm;
    std::vector<T> ym;
    const T normxm = centred_norm(x, &xm);
    const T normym = centred_norm(y, &ym);
    std::vector<T> products(n);
    for (std::size_t i = 0; i < n; ++i) {
        products[i] = xm[i] / normxm * ym[i] / normym;
    }
    T r = pairwise(products);
    r = std::clamp(r, T(-1), T(1));
    if (const_x || const_y) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(r);
}

// scipy.stats.spearmanr(x, y)[0].
template <typename T>
double spearman(const std::vector<T>& x, const std::vector<T>& y) {
    const std::size_t n = x.size();
    if (n <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::all_of(x.begin(), x.end(), [&](T v) { return v == x[0]; }) ||
        std::all_of(y.begin(), y.end(), [&](T v) { return v == y[0]; })) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<double> rx = hicx::stats::rankdata_average(std::vector<double>(x.begin(), x.end()));
    std::vector<double> ry = hicx::stats::rankdata_average(std::vector<double>(y.begin(), y.end()));
    const double mean_x = hicx::npy::pairwise_sum(rx) / static_cast<double>(n);
    const double mean_y = hicx::npy::pairwise_sum(ry) / static_cast<double>(n);
    double c00 = 0.0;
    double c11 = 0.0;
    double c10 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double a = rx[i] - mean_x;
        const double b = ry[i] - mean_y;
        c00 += a * a;
        c11 += b * b;
        c10 += b * a;
    }
    const double inverse = 1.0 / static_cast<double>(n - 1);
    c00 *= inverse;
    c11 *= inverse;
    c10 *= inverse;
    const double r = c10 / std::sqrt(c11) / std::sqrt(c00);
    return std::clamp(r, -1.0, 1.0);
}

// The pair of columns row and col whose sum exceeds 1, in type T.
template <typename T>
std::pair<std::vector<T>, std::vector<T>> pair_vectors(const std::vector<double>& big,
                                                       std::size_t entries, std::size_t num_files,
                                                       std::size_t row, std::size_t col) {
    std::pair<std::vector<T>, std::vector<T>> out;
    for (std::size_t u = 0; u < entries; ++u) {
        const auto a = static_cast<T>(big[u * num_files + row]);
        const auto b = static_cast<T>(big[u * num_files + col]);
        if (a + b > T(1)) {
            out.first.push_back(a);
            out.second.push_back(b);
        }
    }
    return out;
}

// big_mat.T as a .npy file in the reference's dtype, written from the array.
void write_npy(const std::string& path, const std::vector<double>& big, std::size_t entries,
               std::size_t num_files, Kind kind) {
    const char* dtype = kind == Kind::Int ? "<i8" : kind == Kind::F32 ? "<f4" : "<f8";
    const std::string header = hicx::npz::npy_header(
        dtype, {static_cast<std::int64_t>(entries), static_cast<std::int64_t>(num_files)});
    std::FILE* file = std::fopen(path.c_str(), "wb");
    if (file == nullptr) {
        throw std::runtime_error("cannot write " + path);
    }
    bool ok = std::fwrite(header.data(), 1, header.size(), file) == header.size();
    if (kind == Kind::F64) {
        ok = ok && std::fwrite(big.data(), sizeof(double), big.size(), file) == big.size();
    } else {
        constexpr std::size_t chunk = 1 << 16;
        std::vector<char> buffer(chunk * 8);
        for (std::size_t begin = 0; ok && begin < big.size(); begin += chunk) {
            const std::size_t count = std::min(chunk, big.size() - begin);
            std::size_t bytes = 0;
            for (std::size_t k = 0; k < count; ++k) {
                if (kind == Kind::Int) {
                    const auto value = static_cast<std::int64_t>(big[begin + k]);
                    std::memcpy(buffer.data() + bytes, &value, sizeof value);
                    bytes += sizeof value;
                } else {
                    const auto value = static_cast<float>(big[begin + k]);
                    std::memcpy(buffer.data() + bytes, &value, sizeof value);
                    bytes += sizeof value;
                }
            }
            ok = std::fwrite(buffer.data(), 1, bytes, file) == bytes;
        }
    }
    ok = (std::fclose(file) == 0) && ok;
    if (!ok) {
        throw std::runtime_error("cannot write " + path);
    }
}

std::string basename_of(const std::string& path) {
    const std::string::size_type slash = path.rfind('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("hicCorrelate", kDescription);
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .nargs("+")
        .required()
        .input({"h5", "cool"})
        .help("Matrices to correlate (usually .h5 but other formats are allowed). hicCorrelate is "
              "better used on un-corrected matrices in order to exclude any changes introduced "
              "by the correction.");
    cli::ArgumentGroup& heatmap = parser.group("Heatmap arguments");
    heatmap.add({"--zMin", "-min"})
        .type("float")
        .help("Minimum value for the heatmap intensities. If not specified the value is set "
              "automatically.");
    heatmap.add({"--zMax", "-max"})
        .type("float")
        .help("Maximum value for the heatmap intensities.If not specified the value is set "
              "automatically.");
    heatmap.add({"--colorMap"})
        .default_value("jet")
        .metavar("")
        .help("Color map to use for the heatmap. Available values can be seen here: "
              "http://matplotlib.org/examples/color/colormaps_reference.html (Default: "
              "%(default)s).");
    heatmap.add({"--plotNumbers"})
        .action(cli::Action::StoreTrue)
        .required(false)
        .help("If set, then the correlation number is plotted on top of the heatmap.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--method"})
        .choices({"pearson", "spearman"})
        .default_value("pearson")
        .help("Correlation method to use (Default: %(default)s).");
    optional.add({"--log1p"})
        .action(cli::Action::StoreTrue)
        .help("If set, then the log1p of the matrix values is used.");
    optional.add({"--labels", "-l"})
        .metavar("sample1 sample2")
        .nargs("+")
        .help("User defined labels instead of default labels from file names.");
    optional.add({"--range"})
        .help("In bp with the format low_range:high_range, for example 1000000:2000000.");
    optional.add({"--outFileNameHeatmap", "-oh"})
        .default_value("heatmap.png")
        .required()
        .output({"png", "pdf", "svg"})
        .help("File name to save the resulting heatmap plot.");
    optional.add({"--outFileNameScatter", "-os"})
        .default_value("scatter.png")
        .required()
        .output({"png", "pdf", "svg"})
        .help("File name to save the resulting scatter plot.");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("List of chromosomes to be included in the correlation.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figures as JSON (and the vectors as FILE.npy), without "
                  "drawing them (cpp/PLAN.md tier 7).")
        .help("Write the data the figures are drawn from as JSON to this file and do not draw "
              "them.");

    const cli::Namespace ns = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("hicCorrelate", !ns.given("plotData")); refused != 0) {
        return refused;
    }
    const std::vector<std::string> matrices = ns.strs("matrices");
    std::vector<std::string> labels;
    if (ns.given("labels")) {
        labels = ns.strs("labels");
        if (labels.size() != matrices.size()) {
            std::fputs("The number of labels does not match the number of matrices. The "
                       "reference exits 0 here without writing anything; the C++ port exits 1 "
                       "(cpp/AGENTS_CONTRACT.md rule 7).\n",
                       stderr);
            return 1;
        }
    } else {
        for (const std::string& path : matrices) {
            labels.push_back(basename_of(path));
        }
    }
    const bool log1p = ns.flag("log1p");
    const std::string method = ns.str("method");
    const std::optional<std::vector<std::string>> chromosomes =
        ns.given("chromosomes") ? std::optional(ns.strs("chromosomes")) : std::nullopt;
    const std::optional<std::string> range = ns.opt_str("range");
    const std::optional<std::string> plot_data = ns.opt_str("plotData");

    try {
        const std::size_t num_files = matrices.size();
        std::vector<std::int64_t> keys;
        std::vector<double> big;
        std::vector<std::int64_t> all_nan;
        Kind kind = Kind::Int;
        std::int64_t size = 0;
        for (std::size_t s = 0; s < num_files; ++s) {
            SampleEntries sample = load_sample(matrices[s], chromosomes, range, log1p);
            if (s == 0) {
                size = sample.size;
            } else if (sample.size != size) {
                throw PythonError("ValueError: inconsistent shapes (the matrices differ in size)");
            }
            kind = s == 0 ? sample.kind : promote(kind, sample.kind);
            all_nan.insert(all_nan.end(), sample.nan_bins.begin(), sample.nan_bins.end());
            merge_sample(keys, big, s, num_files, sample);
        }

        // rows_keep: the positions outside the NaN bins of every matrix.
        std::sort(all_nan.begin(), all_nan.end());
        all_nan.erase(std::unique(all_nan.begin(), all_nan.end()), all_nan.end());
        std::vector<char> removed(static_cast<std::size_t>(size), 0);
        for (const std::int64_t bin : all_nan) {
            if (bin >= size) {
                throw PythonError("IndexError: index " + std::to_string(bin) +
                                  " is out of bounds for axis 0 with size " + std::to_string(size));
            }
            removed[static_cast<std::size_t>(bin)] = 1;
        }
        std::size_t entries = 0;
        for (std::size_t u = 0; u < keys.size(); ++u) {
            if (removed[static_cast<std::size_t>(keys[u] / size)] ||
                removed[static_cast<std::size_t>(keys[u] % size)]) {
                continue;
            }
            if (entries != u) {
                std::copy(big.begin() + static_cast<std::ptrdiff_t>(u * num_files),
                          big.begin() + static_cast<std::ptrdiff_t>((u + 1) * num_files),
                          big.begin() + static_cast<std::ptrdiff_t>(entries * num_files));
            }
            ++entries;
        }
        std::vector<std::int64_t>().swap(keys);
        big.resize(entries * num_files);
        big.shrink_to_fit();

        // all_mat is the sum of the matrices in command line order, and each
        // stored value becomes (v + a) - a. A stored value is never zero, and
        // adding the zero of a matrix that stores nothing is exact, so the
        // fold over the row equals scipy's sparse additions. Integers are
        // exact and keep their values.
        if (kind != Kind::Int) {
            for (std::size_t u = 0; u < entries; ++u) {
                double* row = big.data() + u * num_files;
                if (kind == Kind::F32) {
                    float a = static_cast<float>(row[0]);
                    for (std::size_t s = 1; s < num_files; ++s) {
                        a = a + static_cast<float>(row[s]);
                    }
                    for (std::size_t s = 0; s < num_files; ++s) {
                        if (row[s] != 0.0) {
                            row[s] = static_cast<double>((static_cast<float>(row[s]) + a) - a);
                        }
                    }
                } else {
                    double a = row[0];
                    for (std::size_t s = 1; s < num_files; ++s) {
                        a = a + row[s];
                    }
                    for (std::size_t s = 0; s < num_files; ++s) {
                        if (row[s] != 0.0) {
                            row[s] = (row[s] + a) - a;
                        }
                    }
                }
            }
        }
        for (const double value : big) {
            if (!std::isfinite(value)) {
                throw PythonError("a vector holds a non-finite value, which the reference masks "
                                  "with np.ma.masked_invalid; the port does not reproduce masked "
                                  "reductions and refuses such input");
            }
        }

        // The pairwise correlations.
        std::vector<double> results(num_files * num_files, 0.0);
        for (std::size_t row = 0; row < num_files; ++row) {
            results[row * num_files + row] = 1.0;
            for (std::size_t col = row + 1; col < num_files; ++col) {
                double r = 0.0;
                if (kind == Kind::F32) {
                    const auto [v1, v2] = pair_vectors<float>(big, entries, num_files, row, col);
                    r = method == "spearman" ? spearman(v1, v2) : pearson(v1, v2);
                } else {
                    const auto [v1, v2] = pair_vectors<double>(big, entries, num_files, row, col);
                    r = method == "spearman" ? spearman(v1, v2) : pearson(v1, v2);
                }
                results[row * num_files + col] = r;
            }
        }
        std::vector<double> symmetric = results;
        for (std::size_t row = 0; row < num_files; ++row) {
            for (std::size_t col = row + 1; col < num_files; ++col) {
                symmetric[col * num_files + row] = results[row * num_files + col];
            }
        }
        hicx::cluster::Samples observations;
        observations.samples = static_cast<std::int64_t>(num_files);
        observations.features = static_cast<std::int64_t>(num_files);
        observations.values = symmetric;
        const std::vector<double> linkage = hicx::cluster::detail::complete_linkage(observations);
        const std::vector<std::int64_t> leaves = hicx::cluster::detail::dendrogram_leaves(
            linkage, static_cast<std::int64_t>(num_files));

        // The vectors travel in a .npy next to the data.
        std::string npy_path;
        if (plot_data.has_value()) {
            npy_path = *plot_data + ".npy";
        } else {
            const char* tmpdir = std::getenv("TMPDIR");
            std::string pattern = std::string(tmpdir != nullptr && *tmpdir != '\0' ? tmpdir : "/tmp") +
                                  "/hicx-plot-XXXXXX";
            std::vector<char> name(pattern.begin(), pattern.end());
            name.push_back('\0');
            const int fd = ::mkstemp(name.data());
            if (fd < 0) {
                throw std::runtime_error("cannot create a temporary file in " + pattern);
            }
            ::close(fd);
            npy_path = name.data();
        }
        write_npy(npy_path, big, entries, num_files, kind);
        std::vector<double>().swap(big);

        std::vector<std::string> result_rows;
        for (std::size_t row = 0; row < num_files; ++row) {
            result_rows.push_back(plot::json_numbers(std::vector<double>(
                results.begin() + static_cast<std::ptrdiff_t>(row * num_files),
                results.begin() + static_cast<std::ptrdiff_t>((row + 1) * num_files))));
        }
        plot::JsonObject data;
        data.add("outFileNameHeatmap", plot::json_string(ns.str("outFileNameHeatmap")));
        data.add("outFileNameScatter", plot::json_string(ns.str("outFileNameScatter")));
        data.add("labels", plot::json_strings(labels));
        data.add("method", plot::json_string(method));
        data.add("log1p", plot::json_bool(log1p));
        const std::optional<double> z_min = ns.opt_real("zMin");
        const std::optional<double> z_max = ns.opt_real("zMax");
        data.add("zMin", z_min.has_value() ? plot::json_number(*z_min) : "null");
        data.add("zMax", z_max.has_value() ? plot::json_number(*z_max) : "null");
        data.add("colorMap", plot::json_string(ns.str("colorMap")));
        data.add("plotNumbers", plot::json_bool(ns.flag("plotNumbers")));
        data.add("results", plot::json_list(result_rows));
        data.add("linkage", plot::json_numbers(linkage));
        data.add("leaves", plot::json_ints(leaves));
        data.add("big_mat", plot::json_string(npy_path));
        data.add("temporary_files",
                 plot_data.has_value() ? "[]" : plot::json_strings({npy_path}));

        hicx::report_resource_usage("hicCorrelate");
        const int status = plot::draw("hicCorrelate", data.str(), plot_data);
        if (!plot_data.has_value()) {
            ::unlink(npy_path.c_str());  // reached only when the drawing could not start
        }
        return status;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicCorrelate: %s\n", error.what());
        return 1;
    }
}
