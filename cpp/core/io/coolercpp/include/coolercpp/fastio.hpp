// Threaded HDF5 reading and writing of cooler files over the HDF5 C library:
// the chunked column writer, the chunk reader and in-place editing of files.
// No public header includes hdf5.h.
//
// Datasets are written the way h5py writes datasets
// (create_dataset with compression='gzip', compression_opts=6, shuffle=True
// and h5py's automatic chunk shape), but their chunks are filtered on several
// threads and stored with H5Dwrite_chunk. A chunk's bytes are what HDF5's own
// shuffle and deflate filters make of it, so the files are ordinary HDF5
// files, and they do not depend on the number of threads.
//
// Reading decodes chunks on several threads the same way when a dataset's
// filters are among shuffle, deflate and fletcher32, and falls back to
// H5Dread for any other layout or filter.

#ifndef COOLERCPP_FASTIO_HPP
#define COOLERCPP_FASTIO_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "coolercpp/parallel.hpp"

namespace coolercpp::fastio {

class Error : public std::runtime_error {
  public:
    explicit Error(const std::string& what) : std::runtime_error(what) {}
};

// The type of a dataset's elements in the file.
struct ColumnType {
    enum class Kind { Int32, Int64, Float64, FixedString, Enum };
    Kind kind = Kind::Int64;
    std::size_t width = 0;                // FixedString: bytes
    std::vector<std::string> enum_names;  // Enum: members in index order

    [[nodiscard]] static ColumnType int32() { return {Kind::Int32, 0, {}}; }
    [[nodiscard]] static ColumnType int64() { return {Kind::Int64, 0, {}}; }
    [[nodiscard]] static ColumnType float64() { return {Kind::Float64, 0, {}}; }
    // numpy S<width>: fixed width, NUL padded, ASCII.
    [[nodiscard]] static ColumnType fixed_string(std::size_t width) { return {Kind::FixedString, width, {}}; }
    // h5py.special_dtype(enum=(int32, {name: index})), members in index order.
    [[nodiscard]] static ColumnType enumeration(std::vector<std::string> names) {
        return {Kind::Enum, 0, std::move(names)};
    }
    // Bytes of one element.
    [[nodiscard]] std::size_t size() const;
};

// An HDF5 object identifier, closed when the handle goes out of scope.
class Handle {
  public:
    enum class Kind { File, Group, Dataset, DataType, DataSpace, Attribute, PropertyList, Object };

    Handle() = default;
    Handle(std::int64_t id, Kind kind) : id_(id), kind_(kind) {}
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    Handle(Handle&& other) noexcept : id_(other.id_), kind_(other.kind_) { other.id_ = -1; }
    Handle& operator=(Handle&& other) noexcept;
    ~Handle() { close(); }

    [[nodiscard]] std::int64_t get() const noexcept { return id_; }
    [[nodiscard]] bool valid() const noexcept { return id_ >= 0; }
    void close() noexcept;

  private:
    std::int64_t id_ = -1;
    Kind kind_ = Kind::File;
};

// A scalar attribute value: None (any type this layer does not read), an
// integer, a float or a string.
using Value = std::variant<std::monostate, std::int64_t, double, std::string>;

enum class Mode {
    Read,       // h5py.File(path, 'r')
    ReadWrite,  // h5py.File(path, 'r+')
    Create,     // h5py.File(path, 'a') on a missing file, 'w' otherwise
};

[[nodiscard]] bool is_hdf5(const std::string& path);

// h5py's automatic chunk length (h5py/_hl/filters.py guess_chunk) for a one
// dimensional dataset of `length` elements of `typesize` bytes.
[[nodiscard]] std::size_t guess_chunk(std::size_t length, std::size_t typesize);

class File {
  public:
    File(const std::string& path, Mode mode);

    [[nodiscard]] std::int64_t id() const noexcept { return file_.get(); }
    [[nodiscard]] const std::string& path() const noexcept { return path_; }
    void close() noexcept { file_.close(); }

    [[nodiscard]] bool exists(const std::string& object) const;
    [[nodiscard]] bool is_group(const std::string& object) const;
    [[nodiscard]] bool is_dataset(const std::string& object) const;
    // Direct children in name order, as h5py iterates a group.
    [[nodiscard]] std::vector<std::string> children(const std::string& group) const;

    // Attributes in name order, as h5py lists them.
    [[nodiscard]] std::vector<std::pair<std::string, Value>> attributes(const std::string& object) const;
    [[nodiscard]] std::optional<Value> attribute(const std::string& object, const std::string& name) const;
    // Python int -> int64, float -> float64, str -> variable length UTF-8;
    // replaces an attribute of the same name.
    void set_attribute(const std::string& object, const std::string& name, const Value& value);

    void create_group(const std::string& path);
    // del group[name]
    void unlink(const std::string& path);

    // A chunked dataset with shuffle and gzip level 6. chunk 0 takes h5py's
    // guess for `length`. resizable gives maxshape (None,).
    Handle create_dataset(const std::string& path, const ColumnType& type, std::size_t length, bool resizable,
                          std::size_t chunk = 0);
    [[nodiscard]] Handle open_dataset(const std::string& path) const;

    [[nodiscard]] std::size_t length(const std::string& dataset) const;

    // Whole datasets, converted by HDF5.
    [[nodiscard]] std::vector<double> read_doubles(const std::string& dataset) const;
    [[nodiscard]] std::vector<std::int64_t> read_int64(const std::string& dataset) const;
    [[nodiscard]] std::vector<std::string> read_strings(const std::string& dataset) const;
    // Overwrites a whole dataset of the same length (dset[:] = values).
    void write_doubles(const std::string& dataset, const std::vector<double>& values);

  private:
    std::string path_;
    Handle file_;
};

// Appends elements to a one dimensional dataset chunk by chunk. Elements are
// given as their bytes in the dataset's file type (little endian).
class ChunkedWriter {
  public:
    ChunkedWriter(Handle dataset, std::size_t element_size, bool resizable);

    void append(const void* elements, std::size_t count);
    [[nodiscard]] std::size_t buffered_chunks() const noexcept;
    [[nodiscard]] std::size_t rows() const noexcept { return written_ + buffer_.size() / element_size_; }

  private:
    friend void flush(ThreadPool& pool, const std::vector<ChunkedWriter*>& writers, bool final);

    Handle dataset_;
    std::size_t element_size_;
    std::size_t chunk_;
    bool resizable_;
    std::size_t written_ = 0;  // rows stored, always a multiple of chunk_ until the end
    std::vector<unsigned char> buffer_;
};

// Filters every complete chunk the writers hold, on the pool, and stores them
// in order; with final, the incomplete last chunk too (padded with zeros, as
// HDF5 fills it), after which the datasets have their final length.
void flush(ThreadPool& pool, const std::vector<ChunkedWriter*>& writers, bool final);

// Writes a whole fixed length dataset from its elements' bytes.
void write_dataset(ThreadPool& pool, File& file, const std::string& path, const ColumnType& type,
                   const void* elements, std::size_t count);

// Reads ranges of a one dimensional numeric dataset. The chunk index is read
// once (H5Dchunk_iter); chunks are then read from the file with pread and
// decoded on the pool, without going through the HDF5 library, which looks
// chunks up one at a time and serializes all access.
class ColumnReader {
  public:
    ColumnReader(const File& file, const std::string& path);
    ~ColumnReader();
    ColumnReader(const ColumnReader&) = delete;
    ColumnReader& operator=(const ColumnReader&) = delete;

    [[nodiscard]] std::size_t size() const noexcept { return length_; }
    void read(ThreadPool& pool, std::size_t lo, std::size_t hi, std::int64_t* out) const;
    void read(ThreadPool& pool, std::size_t lo, std::size_t hi, double* out) const;

  private:
    template <class T>
    void read_as(ThreadPool& pool, std::size_t lo, std::size_t hi, T* out) const;

    struct ChunkLocation {
        std::uint64_t address = 0;  // absolute file offset
        std::uint64_t size = 0;     // 0: not allocated, holds the fill value
        std::uint32_t mask = 0;
    };

    Handle dataset_;
    int fd_ = -1;
    std::vector<ChunkLocation> chunks_;  // by chunk number, from one H5Dchunk_iter
    std::size_t length_ = 0;
    std::size_t chunk_ = 0;
    std::size_t element_size_ = 0;
    int type_class_ = 0;  // H5T_class_t
    bool is_signed_ = true;
    bool direct_ = false;
    // Positions in the filter pipeline, -1 when absent.
    int shuffle_ = -1;
    int deflate_ = -1;
    int fletcher32_ = -1;
    std::vector<unsigned char> fill_;
};

}  // namespace coolercpp::fastio

#endif  // COOLERCPP_FASTIO_HPP
