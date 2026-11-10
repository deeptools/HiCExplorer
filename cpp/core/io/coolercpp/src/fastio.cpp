#include "coolercpp/fastio.hpp"

#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#include <zlib.h>

#include <fcntl.h>
#include <unistd.h>

namespace coolercpp::fastio {

namespace {

// Complete chunks filtered per round of flush; bounds the filtered bytes held.
constexpr std::size_t kMaxChunksPerFlush = 4096;

void silence_errors() {
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
}

// H5Z_filter_shuffle: byte j of every element goes to the j-th plane.
void shuffle(const unsigned char* in, unsigned char* out, std::size_t bytes, std::size_t element) {
    if (element <= 1) {
        std::memcpy(out, in, bytes);
        return;
    }
    const std::size_t n = bytes / element;
    for (std::size_t j = 0; j < element; ++j) {
        unsigned char* plane = out + j * n;
        for (std::size_t i = 0; i < n; ++i) {
            plane[i] = in[i * element + j];
        }
    }
    std::memcpy(out + n * element, in + n * element, bytes - n * element);
}

void unshuffle(const unsigned char* in, unsigned char* out, std::size_t bytes, std::size_t element) {
    if (element <= 1) {
        std::memcpy(out, in, bytes);
        return;
    }
    const std::size_t n = bytes / element;
    for (std::size_t j = 0; j < element; ++j) {
        const unsigned char* plane = in + j * n;
        for (std::size_t i = 0; i < n; ++i) {
            out[i * element + j] = plane[i];
        }
    }
    std::memcpy(out + n * element, in + n * element, bytes - n * element);
}

// The pipeline h5py builds for shuffle=True, compression='gzip',
// compression_opts=6: H5Z_filter_shuffle, then H5Z_filter_deflate, which
// calls compress2 at the given level.
std::vector<unsigned char> filter_chunk(const unsigned char* data, std::size_t bytes, std::size_t element) {
    std::vector<unsigned char> shuffled(bytes);
    shuffle(data, shuffled.data(), bytes, element);
    uLongf size = compressBound(static_cast<uLong>(bytes));
    std::vector<unsigned char> out(size);
    if (compress2(out.data(), &size, shuffled.data(), static_cast<uLong>(bytes), 6) != Z_OK) {
        throw Error("deflate failed on a chunk");
    }
    out.resize(size);
    return out;
}

}  // namespace

Handle& Handle::operator=(Handle&& other) noexcept {
    if (this != &other) {
        close();
        id_ = other.id_;
        kind_ = other.kind_;
        other.id_ = -1;
    }
    return *this;
}

void Handle::close() noexcept {
    if (id_ < 0) {
        return;
    }
    switch (kind_) {
        case Kind::File: H5Fclose(id_); break;
        case Kind::Group: H5Gclose(id_); break;
        case Kind::Dataset: H5Dclose(id_); break;
        case Kind::DataType: H5Tclose(id_); break;
        case Kind::DataSpace: H5Sclose(id_); break;
        case Kind::Attribute: H5Aclose(id_); break;
        case Kind::PropertyList: H5Pclose(id_); break;
        case Kind::Object: H5Oclose(id_); break;
    }
    id_ = -1;
}

bool is_hdf5(const std::string& path) {
    silence_errors();
    return H5Fis_accessible(path.c_str(), H5P_DEFAULT) > 0;
}

std::size_t guess_chunk(std::size_t length, std::size_t typesize) {
    // h5py/_hl/filters.py guess_chunk, one dimensional case.
    constexpr double kChunkBase = 16.0 * 1024.0;
    constexpr double kChunkMin = 8.0 * 1024.0;
    constexpr double kChunkMax = 1024.0 * 1024.0;
    double chunk = length != 0 ? static_cast<double>(length) : 1024.0;
    const double element = static_cast<double>(std::max<std::size_t>(typesize, 1));
    const double dataset_bytes = chunk * element;
    double target = kChunkBase * std::pow(2.0, std::log10(dataset_bytes / (1024.0 * 1024.0)));
    target = std::clamp(target, kChunkMin, kChunkMax);
    while (true) {
        const double chunk_bytes = chunk * element;
        if ((chunk_bytes < target || std::fabs(chunk_bytes - target) / target < 0.5) && chunk_bytes < kChunkMax) {
            break;
        }
        if (chunk == 1.0) {
            break;
        }
        chunk = std::ceil(chunk / 2.0);
    }
    return static_cast<std::size_t>(chunk);
}

Handle make_fixed_string_type(std::size_t width) {
    Handle type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
    if (!type.valid() || H5Tset_size(type.get(), std::max<std::size_t>(width, 1)) < 0 ||
        H5Tset_strpad(type.get(), H5T_STR_NULLPAD) < 0 || H5Tset_cset(type.get(), H5T_CSET_ASCII) < 0) {
        throw Error("cannot build a fixed width string type");
    }
    return type;
}

Handle make_enum_type(const std::vector<std::string>& names) {
    Handle type(H5Tenum_create(H5T_STD_I32LE), Handle::Kind::DataType);
    if (!type.valid()) {
        throw Error("cannot build an enumeration type");
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto value = static_cast<std::int32_t>(i);
        if (H5Tenum_insert(type.get(), names[i].c_str(), &value) < 0) {
            throw Error("cannot add " + names[i] + " to the chromosome enumeration");
        }
    }
    return type;
}

Handle make_type(const ColumnType& type) {
    switch (type.kind) {
        case ColumnType::Kind::Int32:
            return Handle(H5Tcopy(H5T_STD_I32LE), Handle::Kind::DataType);
        case ColumnType::Kind::Int64:
            return Handle(H5Tcopy(H5T_STD_I64LE), Handle::Kind::DataType);
        case ColumnType::Kind::Float64:
            return Handle(H5Tcopy(H5T_IEEE_F64LE), Handle::Kind::DataType);
        case ColumnType::Kind::FixedString:
            return make_fixed_string_type(type.width);
        case ColumnType::Kind::Enum:
            return make_enum_type(type.enum_names);
    }
    throw Error("unknown column type");
}

std::size_t ColumnType::size() const {
    return H5Tget_size(make_type(*this).get());
}

File::File(const std::string& path, Mode mode) : path_(path) {
    silence_errors();
    hid_t id = -1;
    switch (mode) {
        case Mode::Read: id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); break;
        case Mode::ReadWrite: id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT); break;
        case Mode::Create: id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); break;
    }
    if (id < 0) {
        throw Error("Unable to open file " + path);
    }
    file_ = Handle(id, Handle::Kind::File);
}

bool File::exists(const std::string& object) const {
    if (object.empty() || object == "/") {
        return true;
    }
    std::string partial;
    std::size_t start = object[0] == '/' ? 1 : 0;
    while (start <= object.size()) {
        const std::size_t slash = object.find('/', start);
        const std::size_t end = slash == std::string::npos ? object.size() : slash;
        if (end > start) {
            partial += "/" + object.substr(start, end - start);
            if (H5Lexists(file_.get(), partial.c_str(), H5P_DEFAULT) <= 0 ||
                H5Oexists_by_name(file_.get(), partial.c_str(), H5P_DEFAULT) <= 0) {
                return false;
            }
        }
        if (slash == std::string::npos) {
            break;
        }
        start = slash + 1;
    }
    return true;
}

namespace {

H5O_type_t object_type(hid_t file, const std::string& object) {
    H5O_info2_t info;
    if (H5Oget_info_by_name3(file, object.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) {
        return H5O_TYPE_UNKNOWN;
    }
    return info.type;
}

}  // namespace

bool File::is_group(const std::string& object) const {
    return exists(object) && object_type(file_.get(), object.empty() ? "/" : object) == H5O_TYPE_GROUP;
}

bool File::is_dataset(const std::string& object) const {
    return exists(object) && object_type(file_.get(), object) == H5O_TYPE_DATASET;
}

std::vector<std::string> File::children(const std::string& group) const {
    std::vector<std::string> names;
    H5G_info_t info;
    if (H5Gget_info_by_name(file_.get(), group.c_str(), &info, H5P_DEFAULT) < 0) {
        throw Error("cannot list group " + group + " in " + path_);
    }
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        const ssize_t size = H5Lget_name_by_idx(file_.get(), group.c_str(), H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0,
                                                H5P_DEFAULT);
        std::string name(static_cast<std::size_t>(std::max<ssize_t>(size, 0)), '\0');
        H5Lget_name_by_idx(file_.get(), group.c_str(), H5_INDEX_NAME, H5_ITER_INC, i, name.data(), name.size() + 1,
                           H5P_DEFAULT);
        names.push_back(name);
    }
    return names;
}

namespace {

Value read_attribute(hid_t attribute) {
    const Handle type(H5Aget_type(attribute), Handle::Kind::DataType);
    const Handle space(H5Aget_space(attribute), Handle::Kind::DataSpace);
    if (H5Sget_simple_extent_npoints(space.get()) != 1) {
        return std::monostate{};
    }
    switch (H5Tget_class(type.get())) {
        case H5T_INTEGER: {
            std::int64_t value = 0;
            return H5Aread(attribute, H5T_NATIVE_INT64, &value) < 0 ? Value{} : Value{value};
        }
        case H5T_FLOAT: {
            double value = 0;
            return H5Aread(attribute, H5T_NATIVE_DOUBLE, &value) < 0 ? Value{} : Value{value};
        }
        case H5T_STRING: {
            if (H5Tis_variable_str(type.get()) > 0) {
                const Handle memory(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
                H5Tset_size(memory.get(), H5T_VARIABLE);
                H5Tset_cset(memory.get(), H5Tget_cset(type.get()));
                char* text = nullptr;
                if (H5Aread(attribute, memory.get(), &text) < 0) {
                    return std::monostate{};
                }
                std::string value = text != nullptr ? text : "";
                H5free_memory(text);
                return value;
            }
            const std::size_t size = H5Tget_size(type.get());
            std::string value(size, '\0');
            if (H5Aread(attribute, type.get(), value.data()) < 0) {
                return std::monostate{};
            }
            value.resize(strnlen(value.c_str(), size));
            return value;
        }
        default: return std::monostate{};
    }
}

}  // namespace

std::vector<std::pair<std::string, Value>> File::attributes(const std::string& object) const {
    const Handle handle(H5Oopen(file_.get(), object.c_str(), H5P_DEFAULT), Handle::Kind::Object);
    if (!handle.valid()) {
        throw Error("cannot open " + object + " in " + path_);
    }
    H5O_info2_t info;
    H5Oget_info3(handle.get(), &info, H5O_INFO_NUM_ATTRS);
    std::vector<std::pair<std::string, Value>> out;
    for (hsize_t i = 0; i < info.num_attrs; ++i) {
        const Handle attribute(H5Aopen_by_idx(handle.get(), ".", H5_INDEX_NAME, H5_ITER_INC, i, H5P_DEFAULT,
                                              H5P_DEFAULT),
                               Handle::Kind::Attribute);
        const ssize_t size = H5Aget_name(attribute.get(), 0, nullptr);
        std::string name(static_cast<std::size_t>(std::max<ssize_t>(size, 0)), '\0');
        H5Aget_name(attribute.get(), name.size() + 1, name.data());
        out.emplace_back(name, read_attribute(attribute.get()));
    }
    return out;
}

std::optional<Value> File::attribute(const std::string& object, const std::string& name) const {
    if (H5Aexists_by_name(file_.get(), object.c_str(), name.c_str(), H5P_DEFAULT) <= 0) {
        return std::nullopt;
    }
    const Handle attribute(H5Aopen_by_name(file_.get(), object.c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT),
                           Handle::Kind::Attribute);
    return read_attribute(attribute.get());
}

void File::set_attribute(const std::string& object, const std::string& name, const Value& value) {
    const Handle handle(H5Oopen(file_.get(), object.c_str(), H5P_DEFAULT), Handle::Kind::Object);
    if (!handle.valid()) {
        throw Error("cannot open " + object + " in " + path_);
    }
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    Handle type;
    hid_t memory = -1;
    const void* buffer = nullptr;
    const char* text = nullptr;
    std::int64_t integer = 0;
    double number = 0;
    if (const auto* s = std::get_if<std::string>(&value)) {
        type = Handle(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
        H5Tset_size(type.get(), H5T_VARIABLE);
        H5Tset_cset(type.get(), H5T_CSET_UTF8);
        memory = type.get();
        text = s->c_str();
        buffer = &text;
    } else if (const auto* i = std::get_if<std::int64_t>(&value)) {
        type = Handle(H5Tcopy(H5T_STD_I64LE), Handle::Kind::DataType);
        memory = H5T_NATIVE_INT64;
        integer = *i;
        buffer = &integer;
    } else if (const auto* d = std::get_if<double>(&value)) {
        type = Handle(H5Tcopy(H5T_IEEE_F64LE), Handle::Kind::DataType);
        memory = H5T_NATIVE_DOUBLE;
        number = *d;
        buffer = &number;
    } else {
        throw Error("cannot write attribute " + name + " without a value");
    }
    if (H5Aexists(handle.get(), name.c_str()) > 0) {
        H5Adelete(handle.get(), name.c_str());
    }
    const Handle attribute(H5Acreate2(handle.get(), name.c_str(), type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
                           Handle::Kind::Attribute);
    if (!attribute.valid() || H5Awrite(attribute.get(), memory, buffer) < 0) {
        throw Error("cannot write attribute " + name + " on " + object + " in " + path_);
    }
}

void File::create_group(const std::string& path) {
    const Handle plist(H5Pcreate(H5P_LINK_CREATE), Handle::Kind::PropertyList);
    H5Pset_create_intermediate_group(plist.get(), 1);
    const Handle group(H5Gcreate2(file_.get(), path.c_str(), plist.get(), H5P_DEFAULT, H5P_DEFAULT),
                       Handle::Kind::Group);
    if (!group.valid()) {
        throw Error("cannot create group " + path + " in " + path_);
    }
}

void File::unlink(const std::string& path) {
    if (H5Ldelete(file_.get(), path.c_str(), H5P_DEFAULT) < 0) {
        throw Error("cannot delete " + path + " in " + path_);
    }
}

Handle File::create_dataset(const std::string& path, const ColumnType& type, std::size_t length, bool resizable,
                            std::size_t chunk) {
    const Handle file_type_handle = make_type(type);
    const hid_t file_type = file_type_handle.get();
    const hsize_t dims = length;
    const hsize_t maxdims = resizable ? H5S_UNLIMITED : static_cast<hsize_t>(length);
    const Handle space(H5Screate_simple(1, &dims, &maxdims), Handle::Kind::DataSpace);
    const Handle plist(H5Pcreate(H5P_DATASET_CREATE), Handle::Kind::PropertyList);
    if (chunk == 0) {
        chunk = guess_chunk(length, H5Tget_size(file_type));
    }
    const hsize_t chunk_dims = std::max<std::size_t>(chunk, 1);
    if (!space.valid() || !plist.valid() || H5Pset_chunk(plist.get(), 1, &chunk_dims) < 0 ||
        H5Pset_shuffle(plist.get()) < 0 || H5Pset_deflate(plist.get(), 6) < 0) {
        throw Error("cannot set up dataset " + path + " in " + path_);
    }
    const hid_t dataset =
        H5Dcreate2(file_.get(), path.c_str(), file_type, space.get(), H5P_DEFAULT, plist.get(), H5P_DEFAULT);
    if (dataset < 0) {
        throw Error("cannot create dataset " + path + " in " + path_);
    }
    return Handle(dataset, Handle::Kind::Dataset);
}

Handle File::open_dataset(const std::string& path) const {
    const hid_t dataset = H5Dopen2(file_.get(), path.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        throw Error("cannot open dataset " + path + " in " + path_);
    }
    return Handle(dataset, Handle::Kind::Dataset);
}

std::size_t File::length(const std::string& dataset) const {
    const Handle handle = open_dataset(dataset);
    const Handle space(H5Dget_space(handle.get()), Handle::Kind::DataSpace);
    const hssize_t n = H5Sget_simple_extent_npoints(space.get());
    return n < 0 ? 0 : static_cast<std::size_t>(n);
}

namespace {

template <class T>
std::vector<T> read_all(const File& file, const std::string& dataset, hid_t memory) {
    const Handle handle = file.open_dataset(dataset);
    std::vector<T> values(file.length(dataset));
    if (!values.empty() && H5Dread(handle.get(), memory, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
        throw Error("cannot read dataset " + dataset + " in " + file.path());
    }
    return values;
}

}  // namespace

std::vector<double> File::read_doubles(const std::string& dataset) const {
    return read_all<double>(*this, dataset, H5T_NATIVE_DOUBLE);
}

std::vector<std::int64_t> File::read_int64(const std::string& dataset) const {
    return read_all<std::int64_t>(*this, dataset, H5T_NATIVE_INT64);
}

std::vector<std::string> File::read_strings(const std::string& dataset) const {
    const Handle handle = open_dataset(dataset);
    const Handle type(H5Dget_type(handle.get()), Handle::Kind::DataType);
    const std::size_t n = length(dataset);
    std::vector<std::string> out;
    if (H5Tget_class(type.get()) != H5T_STRING) {
        throw Error("dataset " + dataset + " in " + path_ + " does not hold strings");
    }
    if (H5Tis_variable_str(type.get()) > 0) {
        std::vector<char*> texts(n, nullptr);
        const Handle memory(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
        H5Tset_size(memory.get(), H5T_VARIABLE);
        H5Tset_cset(memory.get(), H5Tget_cset(type.get()));
        if (n > 0 && H5Dread(handle.get(), memory.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, texts.data()) < 0) {
            throw Error("cannot read dataset " + dataset + " in " + path_);
        }
        for (char* text : texts) {
            out.emplace_back(text != nullptr ? text : "");
        }
        const Handle space(H5Dget_space(handle.get()), Handle::Kind::DataSpace);
        H5Treclaim(memory.get(), space.get(), H5P_DEFAULT, texts.data());
        return out;
    }
    const std::size_t width = H5Tget_size(type.get());
    std::vector<char> raw(n * width);
    if (n > 0 && H5Dread(handle.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data()) < 0) {
        throw Error("cannot read dataset " + dataset + " in " + path_);
    }
    for (std::size_t i = 0; i < n; ++i) {
        const char* begin = raw.data() + i * width;
        out.emplace_back(begin, strnlen(begin, width));
    }
    return out;
}

void File::write_doubles(const std::string& dataset, const std::vector<double>& values) {
    const Handle handle = open_dataset(dataset);
    if (!values.empty() && H5Dwrite(handle.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
        throw Error("cannot write dataset " + dataset + " in " + path_);
    }
}

// ------------------------------------------------------------------ writing

ChunkedWriter::ChunkedWriter(Handle dataset, std::size_t element_size, bool resizable)
    : dataset_(std::move(dataset)), element_size_(element_size), resizable_(resizable) {
    const Handle plist(H5Dget_create_plist(dataset_.get()), Handle::Kind::PropertyList);
    hsize_t chunk = 0;
    if (H5Pget_chunk(plist.get(), 1, &chunk) != 1) {
        throw Error("dataset is not chunked");
    }
    chunk_ = static_cast<std::size_t>(chunk);
}

void ChunkedWriter::append(const void* elements, std::size_t count) {
    const auto* bytes = static_cast<const unsigned char*>(elements);
    buffer_.insert(buffer_.end(), bytes, bytes + count * element_size_);
}

std::size_t ChunkedWriter::buffered_chunks() const noexcept {
    return buffer_.size() / (chunk_ * element_size_);
}

void flush(ThreadPool& pool, const std::vector<ChunkedWriter*>& writers, bool final) {
    struct Task {
        ChunkedWriter* writer;
        std::size_t offset;  // bytes into the writer's buffer
        std::size_t bytes;   // data bytes of the chunk (less than a chunk at the end)
        std::size_t index;   // chunk number in the dataset
        std::vector<unsigned char> filtered;
    };
    while (true) {
        std::vector<Task> tasks;
        bool more = false;
        for (ChunkedWriter* w : writers) {
            const std::size_t chunk_bytes = w->chunk_ * w->element_size_;
            const std::size_t full = w->buffer_.size() / chunk_bytes;
            const std::size_t take = std::min(full, kMaxChunksPerFlush);
            more = more || take < full;
            const std::size_t first = w->written_ / w->chunk_;
            for (std::size_t c = 0; c < take; ++c) {
                tasks.push_back(Task{w, c * chunk_bytes, chunk_bytes, first + c, {}});
            }
            const std::size_t rest = w->buffer_.size() - full * chunk_bytes;
            if (final && take == full && rest > 0) {
                tasks.push_back(Task{w, full * chunk_bytes, rest, first + full, {}});
            }
        }
        pool.run(tasks.size(), [&](std::size_t i) {
            Task& t = tasks[i];
            const std::size_t chunk_bytes = t.writer->chunk_ * t.writer->element_size_;
            if (t.bytes < chunk_bytes) {
                // HDF5 filters an edge chunk at its full size, the part past
                // the extent holding the fill value, 0.
                std::vector<unsigned char> padded(chunk_bytes, 0);
                std::memcpy(padded.data(), t.writer->buffer_.data() + t.offset, t.bytes);
                t.filtered = filter_chunk(padded.data(), chunk_bytes, t.writer->element_size_);
            } else {
                t.filtered = filter_chunk(t.writer->buffer_.data() + t.offset, chunk_bytes, t.writer->element_size_);
            }
        });
        for (ChunkedWriter* w : writers) {
            std::size_t consumed = 0;
            for (const Task& t : tasks) {
                consumed += t.writer == w ? t.bytes : 0;
            }
            if (consumed == 0) {
                continue;
            }
            const std::size_t rows_after = w->written_ + consumed / w->element_size_;
            if (w->resizable_) {
                const hsize_t extent = rows_after;
                if (H5Dset_extent(w->dataset_.get(), &extent) < 0) {
                    throw Error("cannot extend a dataset");
                }
            }
            for (const Task& t : tasks) {
                if (t.writer != w) {
                    continue;
                }
                const hsize_t offset = static_cast<hsize_t>(t.index) * w->chunk_;
                if (H5Dwrite_chunk(w->dataset_.get(), H5P_DEFAULT, 0, &offset, t.filtered.size(), t.filtered.data()) <
                    0) {
                    throw Error("cannot write a chunk");
                }
            }
            w->buffer_.erase(w->buffer_.begin(), w->buffer_.begin() + static_cast<std::ptrdiff_t>(consumed));
            w->written_ = rows_after;
        }
        if (!more) {
            break;
        }
    }
}

void write_dataset(ThreadPool& pool, File& file, const std::string& path, const ColumnType& type, const void* elements,
                   std::size_t count) {
    ChunkedWriter writer(file.create_dataset(path, type, count, false), type.size(), false);
    writer.append(elements, count);
    flush(pool, {&writer}, true);
}

// ------------------------------------------------------------------ reading

ColumnReader::ColumnReader(const File& file, const std::string& path) : dataset_(file.open_dataset(path)) {
    length_ = file.length(path);
    const Handle type(H5Dget_type(dataset_.get()), Handle::Kind::DataType);
    type_class_ = H5Tget_class(type.get());
    element_size_ = H5Tget_size(type.get());
    bool little = H5Tget_order(type.get()) == H5T_ORDER_LE || element_size_ == 1;
    if (type_class_ == H5T_ENUM) {
        const Handle base(H5Tget_super(type.get()), Handle::Kind::DataType);
        is_signed_ = H5Tget_sign(base.get()) == H5T_SGN_2;
        little = H5Tget_order(base.get()) == H5T_ORDER_LE || element_size_ == 1;
        type_class_ = H5T_INTEGER;
    } else if (type_class_ == H5T_INTEGER) {
        is_signed_ = H5Tget_sign(type.get()) == H5T_SGN_2;
    }
    const bool numeric = (type_class_ == H5T_INTEGER && element_size_ <= 8 && (element_size_ & (element_size_ - 1)) == 0) ||
                         (type_class_ == H5T_FLOAT && (element_size_ == 4 || element_size_ == 8));
    const Handle plist(H5Dget_create_plist(dataset_.get()), Handle::Kind::PropertyList);
    if (!little || !numeric || H5Pget_layout(plist.get()) != H5D_CHUNKED) {
        return;
    }
    hsize_t chunk = 0;
    H5Pget_chunk(plist.get(), 1, &chunk);
    const int nfilters = H5Pget_nfilters(plist.get());
    for (int i = 0; i < nfilters; ++i) {
        unsigned flags = 0;
        std::size_t nelmts = 8;
        unsigned values[8];
        char name[64];
        unsigned config = 0;
        const H5Z_filter_t id =
            H5Pget_filter2(plist.get(), static_cast<unsigned>(i), &flags, &nelmts, values, sizeof(name), name, &config);
        if (id == H5Z_FILTER_SHUFFLE && shuffle_ < 0) {
            shuffle_ = i;
        } else if (id == H5Z_FILTER_DEFLATE && deflate_ < 0) {
            deflate_ = i;
        } else if (id == H5Z_FILTER_FLETCHER32 && fletcher32_ < 0) {
            fletcher32_ = i;
        } else {
            return;  // a filter the decoder does not know: H5Dread
        }
    }
    fill_.assign(element_size_, 0);
    H5D_fill_value_t defined = H5D_FILL_VALUE_UNDEFINED;
    if (H5Pfill_value_defined(plist.get(), &defined) >= 0 && defined == H5D_FILL_VALUE_USER_DEFINED) {
        H5Pget_fill_value(plist.get(), type.get(), fill_.data());
    }
    chunk_ = static_cast<std::size_t>(chunk);
    if (chunk_ == 0) {
        return;
    }
    const Handle file_plist(H5Fget_create_plist(file.id()), Handle::Kind::PropertyList);
    hsize_t userblock = 0;
    H5Pget_userblock(file_plist.get(), &userblock);
    chunks_.assign((length_ + chunk_ - 1) / chunk_, ChunkLocation{});
    struct Context {
        std::vector<ChunkLocation>* chunks;
        std::size_t chunk;
        hsize_t base;
    } context{&chunks_, chunk_, userblock};
    const auto visit = [](const hsize_t* offset, unsigned mask, haddr_t address, hsize_t size, void* data) -> int {
        auto* c = static_cast<Context*>(data);
        const std::size_t index = static_cast<std::size_t>(offset[0]) / c->chunk;
        if (index < c->chunks->size() && address != HADDR_UNDEF) {
            (*c->chunks)[index] = ChunkLocation{static_cast<std::uint64_t>(address + c->base),
                                                static_cast<std::uint64_t>(size), mask};
        }
        return H5_ITER_CONT;
    };
    if (H5Dchunk_iter(dataset_.get(), H5P_DEFAULT, visit, &context) < 0) {
        chunks_.clear();
        return;
    }
    fd_ = ::open(file.path().c_str(), O_RDONLY | O_CLOEXEC);
    direct_ = fd_ >= 0;
}

ColumnReader::~ColumnReader() {
    if (fd_ >= 0) {
        ::close(fd_);
    }
}

namespace {

template <class T>
T element_as(const unsigned char* p, int type_class, std::size_t size, bool is_signed) {
    if (type_class == H5T_FLOAT) {
        if (size == 4) {
            float v;
            std::memcpy(&v, p, 4);
            return static_cast<T>(v);
        }
        double v;
        std::memcpy(&v, p, 8);
        return static_cast<T>(v);
    }
    std::uint64_t bits = 0;
    std::memcpy(&bits, p, size);
    if (is_signed && size < 8 && (bits >> (size * 8 - 1)) & 1u) {
        bits |= ~std::uint64_t{0} << (size * 8);
    }
    return is_signed ? static_cast<T>(static_cast<std::int64_t>(bits)) : static_cast<T>(bits);
}

}  // namespace

template <class T>
void ColumnReader::read_as(ThreadPool& pool, std::size_t lo, std::size_t hi, T* out) const {
    if (hi <= lo) {
        return;
    }
    if (hi > length_) {
        throw Error("read past the end of a dataset");
    }
    if (!direct_) {
        const Handle space(H5Dget_space(dataset_.get()), Handle::Kind::DataSpace);
        const hsize_t start = lo;
        const hsize_t count = hi - lo;
        H5Sselect_hyperslab(space.get(), H5S_SELECT_SET, &start, nullptr, &count, nullptr);
        const Handle memory_space(H5Screate_simple(1, &count, nullptr), Handle::Kind::DataSpace);
        const hid_t memory = std::is_same_v<T, double> ? H5T_NATIVE_DOUBLE : H5T_NATIVE_INT64;
        if (H5Dread(dataset_.get(), memory, memory_space.get(), space.get(), H5P_DEFAULT, out) < 0) {
            throw Error("cannot read a dataset range");
        }
        return;
    }
    const std::size_t chunk_bytes = chunk_ * element_size_;
    const std::size_t first = lo / chunk_;
    const std::size_t last = (hi - 1) / chunk_;
    pool.run(last - first + 1, [&](std::size_t i) {
        const std::size_t c = first + i;
        const ChunkLocation& location = chunks_[c];
        std::vector<unsigned char> data;
        if (location.size == 0) {
            data.resize(chunk_bytes);
            for (std::size_t k = 0; k < chunk_; ++k) {
                std::memcpy(data.data() + k * element_size_, fill_.data(), element_size_);
            }
        } else {
            data.resize(location.size);
            std::size_t done = 0;
            while (done < data.size()) {
                const ssize_t got = ::pread(fd_, data.data() + done, data.size() - done,
                                            static_cast<off_t>(location.address + done));
                if (got <= 0) {
                    throw Error("cannot read a chunk");
                }
                done += static_cast<std::size_t>(got);
            }
            // Undo the filters in reverse pipeline order; a set bit of the
            // mask marks a filter that was skipped for this chunk.
            std::vector<std::pair<int, int>> steps;  // position, kind
            if (shuffle_ >= 0) steps.emplace_back(shuffle_, 0);
            if (deflate_ >= 0) steps.emplace_back(deflate_, 1);
            if (fletcher32_ >= 0) steps.emplace_back(fletcher32_, 2);
            std::sort(steps.rbegin(), steps.rend());
            for (const auto& [position, kind] : steps) {
                if ((location.mask & (1u << static_cast<unsigned>(position))) != 0) {
                    continue;
                }
                if (kind == 2) {
                    if (data.size() < 4) {
                        throw Error("corrupt fletcher32 chunk");
                    }
                    data.resize(data.size() - 4);
                } else if (kind == 1) {
                    std::vector<unsigned char> inflated(chunk_bytes);
                    uLongf size = static_cast<uLongf>(chunk_bytes);
                    if (uncompress(inflated.data(), &size, data.data(), static_cast<uLong>(data.size())) != Z_OK ||
                        size != chunk_bytes) {
                        throw Error("cannot inflate a chunk");
                    }
                    data.swap(inflated);
                } else {
                    std::vector<unsigned char> plain(data.size());
                    unshuffle(data.data(), plain.data(), data.size(), element_size_);
                    data.swap(plain);
                }
            }
            if (data.size() < chunk_bytes) {
                throw Error("a chunk is shorter than its dataset's chunk size");
            }
        }
        const std::size_t chunk_lo = c * chunk_;
        const std::size_t from = std::max(lo, chunk_lo);
        const std::size_t to = std::min(hi, chunk_lo + chunk_);
        for (std::size_t k = from; k < to; ++k) {
            out[k - lo] =
                element_as<T>(data.data() + (k - chunk_lo) * element_size_, type_class_, element_size_, is_signed_);
        }
    });
}

void ColumnReader::read(ThreadPool& pool, std::size_t lo, std::size_t hi, std::int64_t* out) const {
    read_as<std::int64_t>(pool, lo, hi, out);
}

void ColumnReader::read(ThreadPool& pool, std::size_t lo, std::size_t hi, double* out) const {
    read_as<double>(pool, lo, hi, out);
}

}  // namespace coolercpp::fastio
