// Thin RAII layer over the HDF5 C API plus the blosc filter.
//
// The HiCExplorer h5 format is written by PyTables with blosc compression
// (HDF5 filter id 32001). That filter is not built into libhdf5 and no filter
// plugin is installed next to it, so this file registers an implementation on
// top of libblosc itself, in both directions: decompression to read the
// existing corpus, compression to write files PyTables can read again.
//
// Verified on 2026-09-01 against the reference environment: a dataset written
// through this filter at complevel 5 with shuffle is opened by PyTables 3.10.1
// as a CArray with Filters(complevel=5, complib='blosc', shuffle=True), and
// the filter parameters HDF5 records are identical to the ones PyTables emits
// for the same data (cd_values 2, 2, typesize, chunk bytes, 5, 1). That closes
// open question 1 of cpp/STATUS.md.

#ifndef HICX_HDF5_UTIL_HPP
#define HICX_HDF5_UTIL_HPP

#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <hdf5.h>

namespace hicx::h5 {

class Error : public std::runtime_error {
  public:
    explicit Error(const std::string& what) : std::runtime_error(what) {}
};

// Registers the blosc filter with the HDF5 library. Idempotent and thread
// safe enough for our single threaded readers; called automatically by File.
void register_blosc_filter();

// Owning handle that closes with the right H5*close function.
class Handle {
  public:
    enum class Kind {
        File, Group, Dataset, DataType, DataSpace, Attribute, PropertyList
    };

    Handle() = default;
    Handle(hid_t id, Kind kind) : id_(id), kind_(kind) {}
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    Handle(Handle&& other) noexcept { swap(other); }
    Handle& operator=(Handle&& other) noexcept {
        if (this != &other) {
            close();
            swap(other);
        }
        return *this;
    }
    ~Handle() { close(); }

    [[nodiscard]] hid_t get() const noexcept { return id_; }
    [[nodiscard]] bool valid() const noexcept { return id_ >= 0; }
    void close() noexcept;

  private:
    void swap(Handle& other) noexcept {
        std::swap(id_, other.id_);
        std::swap(kind_, other.kind_);
    }
    hid_t id_ = -1;
    Kind kind_ = Kind::File;
};

// An HDF5 attribute value, in the three flavours the cool and h5 formats use.
using AttributeValue = std::variant<std::int64_t, double, std::string>;

class File {
  public:
    explicit File(const std::string& path);

    [[nodiscard]] const std::string& path() const noexcept { return path_; }
    [[nodiscard]] hid_t id() const noexcept { return file_.get(); }

    [[nodiscard]] bool exists(const std::string& object_path) const;

    // Attributes of a group or of the file root ("/").
    [[nodiscard]] std::map<std::string, AttributeValue> attributes(
        const std::string& object_path) const;

    // Names of the direct children of a group, in HDF5 name order.
    [[nodiscard]] std::vector<std::string> children(const std::string& group_path) const;

    [[nodiscard]] std::size_t dataset_length(const std::string& dataset_path) const;

    // True for a two dimensional dataset whose minor extent is one, the n-by-1
    // column layout that one h5 matrix in the corpus stores its correction
    // factors in.
    [[nodiscard]] bool dataset_is_column(const std::string& dataset_path) const;

    // Numeric datasets, converted by HDF5 to the requested C type.
    [[nodiscard]] std::vector<double> read_doubles(const std::string& dataset_path) const;
    [[nodiscard]] std::vector<std::int64_t> read_int64(const std::string& dataset_path) const;
    [[nodiscard]] std::vector<std::int32_t> read_int32(const std::string& dataset_path) const;

    // Fixed length or variable length string datasets.
    [[nodiscard]] std::vector<std::string> read_strings(const std::string& dataset_path) const;

    // The numpy style dtype name of a dataset, for example "int32" or
    // "float64". Needed because hicInfo prints minima and maxima with the
    // formatting of the stored type.
    [[nodiscard]] std::string dataset_dtype(const std::string& dataset_path) const;

  private:
    [[nodiscard]] Handle open_dataset(const std::string& dataset_path) const;

    std::string path_;
    Handle file_;
};

[[nodiscard]] bool is_hdf5(const std::string& path);

// --------------------------------------------------------------------------
// Writing

// The filter pipelines the writers use, named after the Python writer whose
// pipeline each one reproduces. cool files are written by coolercpp.
enum class Filter {
    // No filter at all.
    None,
    // PyTables Filters(complevel=5, complib='blosc') (hicmatrix/lib/h5.py:123).
    PyTablesBlosc,
};

// h5py's automatic chunk layout, h5py/_hl/filters.py guess_chunk, for one
// dimensional datasets. Output is compared including the chunk shape, so this
// has to be the same function and not merely a reasonable one. Returns the
// chunk length in elements.
[[nodiscard]] std::size_t guess_chunk(std::size_t length, std::size_t typesize);

// A fixed width byte string type, NUL padded and ASCII, which is what numpy
// writes for an S<width> array.
[[nodiscard]] Handle fixed_string_type(std::size_t width);

// h5py's file modes, of which the writers need two: the h5 writer truncates,
// and the cool adapter reopens a file coolercpp has written to set hicmatrix's
// root attributes.
enum class WriteMode {
    // h5py.File(path, 'w'): truncate whatever is there.
    Truncate,
    // h5py.File(path, 'a'): open for read and write, create when missing.
    Append,
};

class FileWriter {
  public:
    // Truncates an existing file, like h5py.File(path, 'w').
    explicit FileWriter(const std::string& path, WriteMode mode = WriteMode::Truncate);

    [[nodiscard]] hid_t id() const noexcept { return file_.get(); }
    void close() noexcept { file_.close(); }

    // Creates the group and every missing parent of it, like
    // h5py.Group.create_group with a multi component name.
    Handle create_group(const std::string& path);

    // A resizable dataset is created by passing a max_length larger than
    // length. chunk = 0 asks for h5py's guessed layout, which is derived from
    // the *initial* length exactly as h5py derives it. minor > 0 makes the
    // dataset two dimensional, which one matrix in the corpus needs for its
    // n-by-1 correction factor column.
    Handle create_dataset(const std::string& path, hid_t file_type,
                          std::size_t length, std::size_t max_length,
                          Filter filter, std::size_t chunk = 0,
                          std::size_t minor = 0);

    // Writes count rows at offset. For a two dimensional dataset a row is the
    // whole minor extent. The writers call this once per block so that no full
    // column ever exists in memory.
    static void write_block(hid_t dataset, hid_t mem_type, std::size_t offset,
                            std::size_t count, const void* data);

    void set_attribute(const std::string& object_path, const std::string& name,
                       const AttributeValue& value);

    // A fixed width byte string attribute. PyTables marks every node with
    // CLASS, VERSION and TITLE in this form; an empty leaf title is stored
    // with a null dataspace, which is what null_dataspace selects.
    void set_bytes_attribute(const std::string& object_path, const std::string& name,
                             const std::string& value, bool null_dataspace = false);


  private:
    std::string path_;
    Handle file_;
};

}  // namespace hicx::h5

#endif  // HICX_HDF5_UTIL_HPP
