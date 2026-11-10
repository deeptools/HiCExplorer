// Writer for the cHi-C HDF5 files, the layout lib/viewpoint.py and the cHi-C
// tools create with h5py directly rather than through hicmatrix.
//
// This is a reusable core component. chicViewpoint writes the interaction file
// with it; chicSignificantInteractions, chicAggregateStatistic and
// chicDifferentialTest write files of the same kind (groups per matrix,
// chromosome and gene, a `genes` group of hard links, and per group a handful
// of scalar and one dimensional datasets).
//
// What h5py produces, and therefore what this writes:
//
//   str value (attribute or dataset)   variable length string, UTF-8, scalar
//   Python int / numpy int64           H5T_STD_I64LE
//   Python float / numpy float64       H5T_IEEE_F64LE
//   list of int, list of float         a one dimensional int64 or float64
//                                      array; an *empty* list becomes float64,
//                                      because np.array([]) is float64
//   list attribute of ints             int64 array attribute
//   create_dataset(..., compression='gzip', compression_opts=9)
//                                      chunked with h5py's guessed chunk
//                                      shape, deflate level 9, no shuffle
//   group[name] = other_group          a hard link
//   create_group('a/b')                intermediate groups created
//
// As in hicx::h5::FileWriter, object modification times are not recorded on
// groups and datasets, so that two runs produce identical bytes.

#ifndef HICX_CHIC_HDF5_HPP
#define HICX_CHIC_HDF5_HPP

#include <cstdint>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "hicx/chic_viewpoint.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/python_error.hpp"

namespace hicx::chic {

class Hdf5Writer {
  public:
    // h5py.File(path, 'w')
    explicit Hdf5Writer(const std::string& path);

    [[nodiscard]] bool exists(const std::string& object_path) const;

    // Throws when the group exists already, as h5py's create_group does
    // (ValueError: name already exists).
    void create_group(const std::string& path);

    void set_attribute(const std::string& object_path, const std::string& name,
                       const std::string& value);
    void set_attribute(const std::string& object_path, const std::string& name,
                       std::int64_t value);
    void set_attribute(const std::string& object_path, const std::string& name,
                       std::span<const std::int64_t> values);
    // A Python float: a float64 scalar attribute.
    void set_attribute(const std::string& object_path, const std::string& name, double value);
    // A Python bool: h5py's enumeration {FALSE: 0, TRUE: 1} over int8.
    void set_bool_attribute(const std::string& object_path, const std::string& name, bool value);

    void write_string(const std::string& path, const std::string& value);
    // create_dataset(path, data=[str, ...]): a one dimensional variable length
    // UTF-8 string array, contiguous, without a filter.
    void write_strings(const std::string& path, const std::vector<std::string>& values);
    void write_scalar(const std::string& path, std::int64_t value);
    void write_scalar(const std::string& path, double value);

    // One dimensional arrays, gzip compressed at `level` (0 = no filter).
    void write_array(const std::string& path, std::span<const std::int64_t> values, int level);
    void write_array(const std::string& path, std::span<const double> values, int level);

    // group[link_path] = file[target_path]. Returns false when the link cannot
    // be created, which h5py reports as an exception the cHi-C tools swallow.
    bool hard_link(const std::string& target_path, const std::string& link_path);

    void close() noexcept { file_.close(); }

  private:
    void write_numeric_array(const std::string& path, hid_t file_type, hid_t memory_type,
                             std::size_t length, std::size_t type_size, const void* data,
                             int level);

    std::string path_;
    h5::Handle file_;
};

// Viewpoint.writeInteractionFileHDF5's dataset writes into an existing group:
// chromosome, start_list, end_list, gene, sum_of_interactions,
// relative_position_list, interaction_data_list, pvalue, xfold, raw,
// reference_point_start, reference_point_end.
void write_interaction_datasets(Hdf5Writer& writer, const std::string& group_path,
                                const InteractionFileData& data,
                                std::int64_t reference_point_start,
                                std::int64_t reference_point_end);

// ---------------------------------------------------------------------------
// Reading

// h5py's `path in file`: every component of the path must exist.
[[nodiscard]] bool contains(const h5::File& file, const std::string& path);

// One entry of readInteractionFile's interaction_file_data:
// [chromosome, start, end, gene, sum_of_interactions, relative position,
//  relative interaction, p-value, x-fold, raw], in that order.
struct InteractionRecord {
    std::string chromosome;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::string gene;
    double sum_of_interactions = 0.0;
    std::int64_t relative_position = 0;
    double interaction = 0.0;
    double pvalue = 0.0;
    double xfold = 0.0;
    double raw = 0.0;
};

// Viewpoint.readInteractionFile(pFilePath, triplet), viewpoint.py:100-187.
//
// The two dicts it returns share their keys, the relative positions as numpy
// reads them; `keys` holds them in dict insertion order, so a repeated position
// keeps its first place and its last values. `reference_point` is the list the
// Python returns third: empty when the group is missing or the records cannot
// be built (the try block returns ({}, {}, [])), otherwise [start, end] with
// nullopt for a missing dataset (None).
struct InteractionTable {
    std::vector<double> keys;
    std::map<double, InteractionRecord> records;
    std::vector<std::optional<std::int64_t>> reference_point;
};
[[nodiscard]] InteractionTable read_interaction_table(const h5::File& file,
                                                      const std::vector<std::string>& triplet);

// ---------------------------------------------------------------------------
// h5py's reading side, for the public Python API (hicexplorer._hicx)
//
// The 3.7.6 functions of lib/viewpoint.py and the cHi-C modules open the files
// with h5py.File(path, 'r') and read `dataset[()]` or `dataset[:]`. H5Reader
// gives the same view: the stored dtype, the stored values, and h5py 3.12's
// exception types and messages (thrown as hicx::PythonError), so the
// bindings can build the objects h5py would have returned.

// One dataset as h5py reads it: numeric values in the file's type converted to
// the native one (numpy's dtype of the array), strings as their bytes (h5py
// returns bytes for variable length and numpy.bytes_ for fixed length ones).
struct H5Value {
    enum class Kind { Integer, Unsigned, Float, Bytes, Other };
    Kind kind = Kind::Other;
    // bytes of one numeric element
    std::size_t item_size = 0;
    bool scalar = false;
    bool variable_length = false;
    std::size_t length = 0;
    std::vector<unsigned char> raw;
    std::vector<std::string> strings;

    [[nodiscard]] double as_double(std::size_t i) const;
    [[nodiscard]] std::int64_t as_int64(std::size_t i) const;
    // numpy's dtype.str of the array, for example "<i8"; "|O" for variable
    // length strings and "|S<n>" for fixed length ones.
    [[nodiscard]] std::string dtype_str() const;
};

class H5Reader {
  public:
    // h5py.File(path, 'r'): FileNotFoundError (or the OSError subclass of the
    // errno) when the file cannot be opened, OSError "file signature not
    // found" when it is not an HDF5 file.
    explicit H5Reader(const std::string& path);

    enum class Kind { Missing, Group, Dataset };
    // `path in file`; a path through a dataset or a missing component is
    // Missing.
    [[nodiscard]] Kind kind(const std::string& path) const;
    [[nodiscard]] bool contains(const std::string& path) const { return kind(path) != Kind::Missing; }
    // list(file[path].keys()), in name order. KeyError for a missing group.
    [[nodiscard]] std::vector<std::string> keys(const std::string& path) const;
    // file[path][()] (or [:]). KeyError for a missing object; the dataset's
    // h5py message when `path` is a group is not needed by any caller, so a
    // group raises TypeError.
    [[nodiscard]] H5Value read(const std::string& path) const;
    // The KeyError message h5py gives for file[path] when it is missing.
    [[nodiscard]] std::string missing_message(const std::string& path) const;

    void close() noexcept { file_.reset(); }
    [[nodiscard]] const h5::File& file() const { return *file_; }

  private:
    std::optional<h5::File> file_;
};

// pybedtools.BedTool(lines).sort().merge(d=distance) of bedtools 2.31, the
// target regions chicSignificantInteractions writes: sorted by chromosome name
// and start, then every interval that starts at most `distance` after the end
// of the merged one before it is merged into it.
struct TargetRegion {
    std::string chromosome;
    std::int64_t start = 0;
    std::int64_t end = 0;
};
[[nodiscard]] std::vector<TargetRegion> sort_merge_regions(std::vector<TargetRegion> regions,
                                                           std::int64_t distance);

}  // namespace hicx::chic

#endif  // HICX_CHIC_HDF5_HPP
