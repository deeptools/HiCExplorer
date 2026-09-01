// Thin RAII layer over the HDF5 C API plus the blosc decompression filter.
//
// The HiCExplorer h5 format is written by PyTables with blosc compression
// (HDF5 filter id 32001). That filter is not built into libhdf5 and no filter
// plugin is installed next to it, so the reader registers a decompress only
// implementation on top of libblosc itself.

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
    enum class Kind { File, Group, Dataset, DataType, DataSpace, Attribute };

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

}  // namespace hicx::h5

#endif  // HICX_HDF5_UTIL_HPP
