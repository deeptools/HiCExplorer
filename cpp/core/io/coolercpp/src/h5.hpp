// Internal: a thin RAII layer over the HDF5 C API that behaves like the parts
// of h5py cooler relies on (file modes, link lookup, dataset slicing, dataset
// creation with h5py's filter pipeline order, attribute decoding).
//
// Derived from HiCExplorer v4 cpp/core/src/hdf5_util.cpp (the Handle class,
// guess_chunk, the string/enum type builders and the attribute reader), with
// the blosc filter removed and the h5py modes and dataset creation added.

#ifndef COOLERCPP_H5_HPP
#define COOLERCPP_H5_HPP

#include <hdf5.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "coolercpp/dtype.hpp"
#include "coolercpp/errors.hpp"
#include "coolercpp/json.hpp"
#include "coolercpp/table.hpp"

namespace coolercpp::h5 {

class Handle {
  public:
    enum class Kind { File, Group, Dataset, DataType, DataSpace, Attribute, PropertyList, Object };

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

// h5py.File modes.
enum class Mode { Read, ReadWrite, Append, Truncate, Exclusive };
[[nodiscard]] Mode parse_mode(std::string_view mode);

// h5py.is_hdf5.
[[nodiscard]] bool is_hdf5(const std::string& path);

// How h5py maps an HDF5 datatype to numpy.
struct TypeDesc {
    enum class Kind { Integer, Float, Bool, Enum, FixedString, VlenString, Other };
    Kind kind = Kind::Other;
    // Integer/Float: the numeric dtype. Enum: the base integer dtype.
    // Bool (h5py's int8 FALSE/TRUE enum): Bool. Strings: String.
    DType dtype = DType::Float64;
    std::size_t size = 0;
    std::vector<std::pair<std::string, std::int64_t>> enum_members;
};
[[nodiscard]] TypeDesc describe_type(hid_t type);

// A new copy of the little endian file type h5py creates for a numpy dtype;
// Bool is h5py's enum of int8 {FALSE: 0, TRUE: 1}.
[[nodiscard]] Handle file_type(DType dtype);
// The native memory type for a numeric dtype.
[[nodiscard]] Handle memory_type(DType dtype);
// numpy 'S<width>': NUL padded ASCII.
[[nodiscard]] Handle fixed_string_type(std::size_t width);
// h5py.special_dtype(enum=(base, mapping)) with members in id order.
[[nodiscard]] Handle enum_type(const std::vector<std::string>& names, DType base);

// h5py/_hl/filters.py guess_chunk for a one dimensional dataset.
[[nodiscard]] std::size_t guess_chunk(std::size_t length, std::size_t typesize);

// The dataset creation property list h5py builds (filters.fill_dcpl plus
// make_new_dset), in h5py's filter order: scale-offset, shuffle, compression,
// fletcher32.
struct DatasetCreate {
    std::size_t length = 0;
    // nullopt: fixed size (maxshape None); kUnlimited: H5S_UNLIMITED.
    std::optional<std::size_t> maxlength;
    static constexpr std::size_t kUnlimited = static_cast<std::size_t>(-1);
    // nullopt: contiguous layout.
    std::optional<std::size_t> chunk;
    std::optional<std::int64_t> scaleoffset;
    bool shuffle = false;
    enum class Compression { None, Gzip, Szip, FilterId };
    Compression compression = Compression::None;
    unsigned gzip_level = 0;
    unsigned szip_mask = 0;
    unsigned szip_pixels = 0;
    int filter_id = 0;
    std::vector<unsigned> filter_options;
    bool fletcher32 = false;
    // A numeric fill value, written in the dataset's own type.
    std::optional<json::Value> fillvalue;
    bool track_times = false;
};

class Dataset {
  public:
    Dataset() = default;
    Dataset(Handle handle, std::string path);

    [[nodiscard]] hid_t id() const noexcept { return handle_.get(); }
    [[nodiscard]] const std::string& path() const noexcept { return path_; }
    [[nodiscard]] std::size_t length() const;
    [[nodiscard]] TypeDesc type() const;

    // Reads rows [lo, lo + n) converted into `mem_type`.
    void read_raw(std::size_t lo, std::size_t n, hid_t mem_type, void* out) const;
    template <typename T>
    [[nodiscard]] std::vector<T> read(std::size_t lo, std::size_t hi) const {
        std::vector<T> out(hi > lo ? hi - lo : 0);
        if (!out.empty()) {
            const Handle mem = memory_type(dtype_v<T>);
            read_raw(lo, out.size(), mem.get(), out.data());
        }
        return out;
    }
    // The slice [lo, hi) in the dtype numpy would give it: numbers keep their
    // type, enums read as their base integer type, strings as String.
    [[nodiscard]] Column read_column(std::size_t lo, std::size_t hi) const;
    [[nodiscard]] std::vector<std::string> read_strings(std::size_t lo, std::size_t hi) const;

    // Writes n rows at lo from a buffer of `mem_type`; HDF5 converts to the
    // file type, as h5py leaves it to do.
    void write_raw(std::size_t lo, std::size_t n, hid_t mem_type, const void* data);
    // Writes a column at lo in its own dtype.
    void write_column(std::size_t lo, const Column& column);
    void resize(std::size_t length);

  private:
    Handle handle_;
    std::string path_;
};

class File {
  public:
    File(const std::string& path, Mode mode);

    [[nodiscard]] hid_t id() const noexcept { return file_.get(); }
    [[nodiscard]] const std::string& path() const noexcept { return path_; }

    // `path in file` in h5py: every component exists and resolves.
    [[nodiscard]] bool exists(const std::string& path) const;
    [[nodiscard]] bool is_group(const std::string& path) const;
    [[nodiscard]] bool is_dataset(const std::string& path) const;
    // Group member names in h5py's iteration order (name order, or creation
    // order when the group tracks it).
    [[nodiscard]] std::vector<std::string> keys(const std::string& group) const;

    // Throws KeyError when the dataset is missing.
    [[nodiscard]] Dataset open_dataset(const std::string& path) const;

    // Group.create_group: intermediate groups are created; an existing name
    // raises ValueError.
    void create_group(const std::string& path);
    // del group[name]; a missing name raises KeyError.
    void remove(const std::string& path);
    // Creates a one dimensional dataset. HDF5 failures raise ValueError, which
    // is what h5py raises for, for example, an enum header that is too large.
    Dataset create_dataset(const std::string& path, hid_t type, const DatasetCreate& options);

    // h5py attrs.items(): values as h5py returns them (str, numpy.bytes_,
    // numpy scalars, arrays), before any JSON decoding.
    [[nodiscard]] std::vector<std::pair<std::string, json::Value>> attributes(
        const std::string& path) const;
    // attrs[name] = value. Strings are stored as variable length UTF-8,
    // numbers in their dtype, booleans as h5py's bool enum.
    void set_attribute(const std::string& path, const std::string& name,
                       const json::Value& value);

    void flush();

  private:
    std::string path_;
    Handle file_;
};

}  // namespace coolercpp::h5

#endif  // COOLERCPP_H5_HPP
