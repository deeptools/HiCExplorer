#include "h5.hpp"

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace coolercpp::h5 {

namespace {

// The HDF5 error stack is per thread in a threadsafe build; failures are
// reported through exceptions instead of being printed.
void quiet() { H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); }

herr_t collect_link(hid_t /*group*/, const char* name, const H5L_info2_t* /*info*/,
                    void* data) {
    static_cast<std::vector<std::string>*>(data)->emplace_back(name);
    return 0;
}

bool file_exists(const std::string& path, off_t* size) {
    struct stat info {};
    if (::stat(path.c_str(), &info) != 0) {
        return false;
    }
    if (size != nullptr) {
        *size = info.st_size;
    }
    return true;
}

Handle link_create_plist() {
    Handle plist(H5Pcreate(H5P_LINK_CREATE), Handle::Kind::PropertyList);
    if (!plist.valid() || H5Pset_create_intermediate_group(plist.get(), 1) < 0) {
        throw OSError("cannot configure HDF5 link creation");
    }
    return plist;
}

std::vector<std::string> split_path(const std::string& path) {
    std::vector<std::string> parts;
    std::size_t pos = 0;
    while (pos <= path.size()) {
        const std::size_t slash = path.find('/', pos);
        const std::size_t end = slash == std::string::npos ? path.size() : slash;
        if (end > pos) {
            parts.push_back(path.substr(pos, end - pos));
        }
        if (slash == std::string::npos) {
            break;
        }
        pos = slash + 1;
    }
    return parts;
}

DType integer_dtype(std::size_t size, bool is_signed) {
    switch (size) {
        case 1: return is_signed ? DType::Int8 : DType::UInt8;
        case 2: return is_signed ? DType::Int16 : DType::UInt16;
        case 4: return is_signed ? DType::Int32 : DType::UInt32;
        default: return is_signed ? DType::Int64 : DType::UInt64;
    }
}

template <typename T>
json::Value numeric_value(T value, DType dtype) {
    if (std::is_floating_point_v<T>) {
        return json::Value::number(static_cast<double>(value), dtype);
    }
    if (std::is_signed_v<T>) {
        return json::Value::integer(static_cast<std::int64_t>(value), dtype);
    }
    return json::Value::unsigned_integer(static_cast<std::uint64_t>(value), dtype);
}

// Decodes `count` elements of an attribute of type `type` into values.
std::vector<json::Value> read_attribute_values(hid_t attr, hid_t type, std::size_t count) {
    std::vector<json::Value> values;
    const TypeDesc desc = describe_type(type);
    switch (desc.kind) {
        case TypeDesc::Kind::VlenString: {
            std::vector<char*> raw(count, nullptr);
            const Handle mem(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
            H5Tset_size(mem.get(), H5T_VARIABLE);
            H5Tset_cset(mem.get(), H5Tget_cset(type));
            if (H5Aread(attr, mem.get(), raw.data()) < 0) {
                throw OSError("cannot read a string attribute");
            }
            for (char* item : raw) {
                values.emplace_back(std::string(item != nullptr ? item : ""));
            }
            const Handle space(H5Aget_space(attr), Handle::Kind::DataSpace);
            H5Treclaim(mem.get(), space.get(), H5P_DEFAULT, raw.data());
            return values;
        }
        case TypeDesc::Kind::FixedString: {
            std::vector<char> buffer(desc.size * count);
            if (H5Aread(attr, type, buffer.data()) < 0) {
                throw OSError("cannot read a string attribute");
            }
            for (std::size_t i = 0; i < count; ++i) {
                std::string item(buffer.data() + i * desc.size, desc.size);
                // numpy.bytes_ drops trailing NULs.
                while (!item.empty() && item.back() == '\0') {
                    item.pop_back();
                }
                values.push_back(json::Value::bytes(std::move(item)));
            }
            return values;
        }
        case TypeDesc::Kind::Bool:
        case TypeDesc::Kind::Enum:
        case TypeDesc::Kind::Integer: {
            const bool is_signed = is_signed_integer(desc.dtype) || desc.kind == TypeDesc::Kind::Bool;
            std::vector<std::int64_t> ints(count);
            std::vector<std::uint64_t> uints(count);
            if (desc.kind == TypeDesc::Kind::Integer) {
                if (is_signed) {
                    if (H5Aread(attr, H5T_NATIVE_INT64, ints.data()) < 0) {
                        throw OSError("cannot read an integer attribute");
                    }
                } else if (H5Aread(attr, H5T_NATIVE_UINT64, uints.data()) < 0) {
                    throw OSError("cannot read an integer attribute");
                }
            } else {
                // Enumerations cannot be converted to plain integers by HDF5;
                // read them in their native layout and widen here.
                const Handle native(H5Tget_native_type(type, H5T_DIR_ASCEND),
                                    Handle::Kind::DataType);
                std::vector<unsigned char> raw(desc.size * count);
                if (H5Aread(attr, native.get(), raw.data()) < 0) {
                    throw OSError("cannot read an enum attribute");
                }
                for (std::size_t i = 0; i < count; ++i) {
                    const unsigned char* p = raw.data() + i * desc.size;
                    std::int64_t v = 0;
                    switch (desc.size) {
                        case 1: v = static_cast<std::int8_t>(p[0]); break;
                        case 2: { std::int16_t x; std::memcpy(&x, p, 2); v = x; break; }
                        case 4: { std::int32_t x; std::memcpy(&x, p, 4); v = x; break; }
                        default: std::memcpy(&v, p, 8); break;
                    }
                    ints[i] = v;
                }
            }
            for (std::size_t i = 0; i < count; ++i) {
                if (desc.kind == TypeDesc::Kind::Bool) {
                    values.emplace_back(ints[i] != 0);
                } else if (is_signed) {
                    values.push_back(json::Value::integer(ints[i], desc.dtype));
                } else {
                    values.push_back(json::Value::unsigned_integer(uints[i], desc.dtype));
                }
            }
            return values;
        }
        case TypeDesc::Kind::Float: {
            std::vector<double> numbers(count);
            if (H5Aread(attr, H5T_NATIVE_DOUBLE, numbers.data()) < 0) {
                throw OSError("cannot read a float attribute");
            }
            for (double d : numbers) {
                values.push_back(json::Value::number(d, desc.dtype));
            }
            return values;
        }
        case TypeDesc::Kind::Other: break;
    }
    throw TypeError("unsupported attribute datatype");
}

}  // namespace

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

Mode parse_mode(std::string_view mode) {
    if (mode == "r") return Mode::Read;
    if (mode == "r+") return Mode::ReadWrite;
    if (mode == "a") return Mode::Append;
    if (mode == "w") return Mode::Truncate;
    if (mode == "w-" || mode == "x") return Mode::Exclusive;
    throw ValueError("Invalid mode; must be one of r, r+, w, w-, x, a");
}

bool is_hdf5(const std::string& path) {
    quiet();
    if (!file_exists(path, nullptr)) {
        return false;
    }
    return H5Fis_accessible(path.c_str(), H5P_DEFAULT) > 0;
}

TypeDesc describe_type(hid_t type) {
    TypeDesc desc;
    desc.size = H5Tget_size(type);
    switch (H5Tget_class(type)) {
        case H5T_INTEGER:
            desc.kind = TypeDesc::Kind::Integer;
            desc.dtype = integer_dtype(desc.size, H5Tget_sign(type) != H5T_SGN_NONE);
            break;
        case H5T_FLOAT:
            desc.kind = TypeDesc::Kind::Float;
            desc.dtype = desc.size == 4 ? DType::Float32 : DType::Float64;
            break;
        case H5T_STRING:
            desc.kind = H5Tis_variable_str(type) > 0 ? TypeDesc::Kind::VlenString
                                                     : TypeDesc::Kind::FixedString;
            desc.dtype = DType::String;
            break;
        case H5T_ENUM: {
            const Handle base(H5Tget_super(type), Handle::Kind::DataType);
            desc.dtype = integer_dtype(H5Tget_size(base.get()),
                                       H5Tget_sign(base.get()) != H5T_SGN_NONE);
            const int n = H5Tget_nmembers(type);
            for (int i = 0; i < n; ++i) {
                char* name = H5Tget_member_name(type, static_cast<unsigned>(i));
                std::int64_t value = 0;
                const Handle native(H5Tget_native_type(base.get(), H5T_DIR_ASCEND),
                                    Handle::Kind::DataType);
                std::vector<unsigned char> raw(H5Tget_size(native.get()));
                H5Tget_member_value(type, static_cast<unsigned>(i), raw.data());
                switch (raw.size()) {
                    case 1: value = static_cast<std::int8_t>(raw[0]); break;
                    case 2: { std::int16_t x; std::memcpy(&x, raw.data(), 2); value = x; break; }
                    case 4: { std::int32_t x; std::memcpy(&x, raw.data(), 4); value = x; break; }
                    default: std::memcpy(&value, raw.data(), 8); break;
                }
                desc.enum_members.emplace_back(name != nullptr ? name : "", value);
                H5free_memory(name);
            }
            // h5py reads an int8 enum of exactly FALSE=0, TRUE=1 as numpy bool.
            desc.kind = TypeDesc::Kind::Enum;
            if (desc.size == 1 && desc.enum_members.size() == 2 &&
                desc.enum_members[0] == std::pair<std::string, std::int64_t>{"FALSE", 0} &&
                desc.enum_members[1] == std::pair<std::string, std::int64_t>{"TRUE", 1}) {
                desc.kind = TypeDesc::Kind::Bool;
                desc.dtype = DType::Bool;
            }
            break;
        }
        default: desc.kind = TypeDesc::Kind::Other; break;
    }
    return desc;
}

Handle file_type(DType dtype) {
    hid_t base = -1;
    switch (dtype) {
        case DType::Bool: {
            Handle type(H5Tenum_create(H5T_NATIVE_INT8), Handle::Kind::DataType);
            std::int8_t v = 0;
            H5Tenum_insert(type.get(), "FALSE", &v);
            v = 1;
            H5Tenum_insert(type.get(), "TRUE", &v);
            return type;
        }
        case DType::Int8: base = H5T_STD_I8LE; break;
        case DType::Int16: base = H5T_STD_I16LE; break;
        case DType::Int32: base = H5T_STD_I32LE; break;
        case DType::Int64: base = H5T_STD_I64LE; break;
        case DType::UInt8: base = H5T_STD_U8LE; break;
        case DType::UInt16: base = H5T_STD_U16LE; break;
        case DType::UInt32: base = H5T_STD_U32LE; break;
        case DType::UInt64: base = H5T_STD_U64LE; break;
        case DType::Float32: base = H5T_IEEE_F32LE; break;
        case DType::Float64: base = H5T_IEEE_F64LE; break;
        default: throw TypeError("no HDF5 file type for dtype " + std::string(dtype_name(dtype)));
    }
    return Handle(H5Tcopy(base), Handle::Kind::DataType);
}

Handle memory_type(DType dtype) {
    hid_t base = -1;
    switch (dtype) {
        case DType::Bool: return file_type(DType::Bool);
        case DType::Int8: base = H5T_NATIVE_INT8; break;
        case DType::Int16: base = H5T_NATIVE_INT16; break;
        case DType::Int32: base = H5T_NATIVE_INT32; break;
        case DType::Int64: base = H5T_NATIVE_INT64; break;
        case DType::UInt8: base = H5T_NATIVE_UINT8; break;
        case DType::UInt16: base = H5T_NATIVE_UINT16; break;
        case DType::UInt32: base = H5T_NATIVE_UINT32; break;
        case DType::UInt64: base = H5T_NATIVE_UINT64; break;
        case DType::Float32: base = H5T_NATIVE_FLOAT; break;
        case DType::Float64: base = H5T_NATIVE_DOUBLE; break;
        default:
            throw TypeError("no HDF5 memory type for dtype " + std::string(dtype_name(dtype)));
    }
    return Handle(H5Tcopy(base), Handle::Kind::DataType);
}

Handle fixed_string_type(std::size_t width) {
    Handle type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
    if (!type.valid() || H5Tset_size(type.get(), std::max<std::size_t>(width, 1)) < 0 ||
        H5Tset_strpad(type.get(), H5T_STR_NULLPAD) < 0 ||
        H5Tset_cset(type.get(), H5T_CSET_ASCII) < 0) {
        throw OSError("cannot build a fixed width string type");
    }
    return type;
}

Handle enum_type(const std::vector<std::string>& names, DType base) {
    const Handle base_type = file_type(base);
    Handle type(H5Tenum_create(base_type.get()), Handle::Kind::DataType);
    if (!type.valid()) {
        throw OSError("cannot build an enumeration type");
    }
    // h5py inserts the members of an enum dtype in the sorted order of their
    // names, not in the order of their values.
    std::vector<std::size_t> order(names.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&names](std::size_t a, std::size_t b) { return names[a] < names[b]; });
    const std::size_t width = H5Tget_size(base_type.get());
    for (const std::size_t i : order) {
        const auto value = static_cast<std::int64_t>(i);
        std::vector<unsigned char> raw(width);
        const auto v32 = static_cast<std::int32_t>(value);
        const auto v16 = static_cast<std::int16_t>(value);
        const auto v8 = static_cast<std::int8_t>(value);
        switch (width) {
            case 1: std::memcpy(raw.data(), &v8, 1); break;
            case 2: std::memcpy(raw.data(), &v16, 2); break;
            case 4: std::memcpy(raw.data(), &v32, 4); break;
            default: std::memcpy(raw.data(), &value, 8); break;
        }
        if (H5Tenum_insert(type.get(), names[i].c_str(), raw.data()) < 0) {
            throw ValueError("cannot add '" + names[i] + "' to an enumeration");
        }
    }
    return type;
}

std::size_t guess_chunk(std::size_t length, std::size_t typesize) {
    // h5py/_hl/filters.py guess_chunk, one dimensional case.
    constexpr double kChunkBase = 16.0 * 1024.0;
    constexpr double kChunkMin = 8.0 * 1024.0;
    constexpr double kChunkMax = 1024.0 * 1024.0;
    double chunk = length != 0 ? static_cast<double>(length) : 1024.0;
    const double element = static_cast<double>(typesize);
    const double dataset_bytes = chunk * element;
    double target = kChunkBase * std::pow(2.0, std::log10(dataset_bytes / (1024.0 * 1024.0)));
    if (target > kChunkMax) {
        target = kChunkMax;
    } else if (target < kChunkMin) {
        target = kChunkMin;
    }
    while (true) {
        const double chunk_bytes = chunk * element;
        if ((chunk_bytes < target || std::fabs(chunk_bytes - target) / target < 0.5) &&
            chunk_bytes < kChunkMax) {
            break;
        }
        if (chunk == 1.0) {
            break;
        }
        chunk = std::ceil(chunk / 2.0);
    }
    return static_cast<std::size_t>(chunk);
}

// ---------------------------------------------------------------------------
// Dataset

Dataset::Dataset(Handle handle, std::string path)
    : handle_(std::move(handle)), path_(std::move(path)) {}

std::size_t Dataset::length() const {
    const Handle space(H5Dget_space(id()), Handle::Kind::DataSpace);
    hsize_t dims[H5S_MAX_RANK];
    const int rank = H5Sget_simple_extent_dims(space.get(), dims, nullptr);
    if (rank < 0) {
        throw OSError("cannot read the extent of " + path_);
    }
    if (rank == 0) {
        return 1;
    }
    return static_cast<std::size_t>(dims[0]);
}

TypeDesc Dataset::type() const {
    const Handle type(H5Dget_type(id()), Handle::Kind::DataType);
    return describe_type(type.get());
}

void Dataset::read_raw(std::size_t lo, std::size_t n, hid_t mem_type, void* out) const {
    if (n == 0) {
        return;
    }
    const Handle space(H5Dget_space(id()), Handle::Kind::DataSpace);
    const hsize_t start = lo;
    const hsize_t count = n;
    if (H5Sselect_hyperslab(space.get(), H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0) {
        throw OSError("cannot select rows of " + path_);
    }
    const Handle mem_space(H5Screate_simple(1, &count, nullptr), Handle::Kind::DataSpace);
    if (H5Dread(id(), mem_type, mem_space.get(), space.get(), H5P_DEFAULT, out) < 0) {
        throw OSError("Can't read data (" + path_ + ")");
    }
}

std::vector<std::string> Dataset::read_strings(std::size_t lo, std::size_t hi) const {
    const std::size_t n = hi > lo ? hi - lo : 0;
    std::vector<std::string> out;
    out.reserve(n);
    if (n == 0) {
        return out;
    }
    const Handle type(H5Dget_type(id()), Handle::Kind::DataType);
    if (H5Tget_class(type.get()) != H5T_STRING) {
        throw TypeError(path_ + " is not a string dataset");
    }
    if (H5Tis_variable_str(type.get()) > 0) {
        std::vector<char*> raw(n, nullptr);
        const Handle mem(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
        H5Tset_size(mem.get(), H5T_VARIABLE);
        H5Tset_cset(mem.get(), H5Tget_cset(type.get()));
        read_raw(lo, n, mem.get(), raw.data());
        for (char* item : raw) {
            out.emplace_back(item != nullptr ? item : "");
        }
        const hsize_t count = n;
        const Handle mem_space(H5Screate_simple(1, &count, nullptr), Handle::Kind::DataSpace);
        H5Treclaim(mem.get(), mem_space.get(), H5P_DEFAULT, raw.data());
        return out;
    }
    const std::size_t item = H5Tget_size(type.get());
    std::vector<char> buffer(item * n);
    read_raw(lo, n, type.get(), buffer.data());
    for (std::size_t i = 0; i < n; ++i) {
        const char* p = buffer.data() + i * item;
        std::size_t used = item;
        // numpy.bytes_ strips trailing NULs only.
        while (used > 0 && p[used - 1] == '\0') {
            --used;
        }
        out.emplace_back(p, used);
    }
    return out;
}

Column Dataset::read_column(std::size_t lo, std::size_t hi) const {
    const std::size_t n = hi > lo ? hi - lo : 0;
    const Handle type(H5Dget_type(id()), Handle::Kind::DataType);
    const TypeDesc desc = describe_type(type.get());
    switch (desc.kind) {
        case TypeDesc::Kind::Integer:
        case TypeDesc::Kind::Float: {
            Column column = Column::empty(desc.dtype, n);
            const Handle mem = memory_type(desc.dtype);
            std::visit(
                [&](auto& v) {
                    using V = std::decay_t<decltype(v)>;
                    if constexpr (!std::is_same_v<V, CategoricalData> &&
                                  !std::is_same_v<V, std::vector<std::string>>) {
                        read_raw(lo, n, mem.get(), v.data());
                    }
                },
                column.data());
            return column;
        }
        case TypeDesc::Kind::Bool:
        case TypeDesc::Kind::Enum: {
            Column column = Column::empty(desc.dtype, n);
            const Handle native(H5Tget_native_type(type.get(), H5T_DIR_ASCEND),
                                Handle::Kind::DataType);
            std::visit(
                [&](auto& v) {
                    using V = std::decay_t<decltype(v)>;
                    if constexpr (!std::is_same_v<V, CategoricalData> &&
                                  !std::is_same_v<V, std::vector<std::string>>) {
                        read_raw(lo, n, native.get(), v.data());
                    }
                },
                column.data());
            return column;
        }
        case TypeDesc::Kind::FixedString:
        case TypeDesc::Kind::VlenString: return Column(read_strings(lo, hi));
        case TypeDesc::Kind::Other: break;
    }
    throw TypeError("unsupported datatype in " + path_);
}

void Dataset::write_raw(std::size_t lo, std::size_t n, hid_t mem_type, const void* data) {
    if (n == 0) {
        return;
    }
    const Handle space(H5Dget_space(id()), Handle::Kind::DataSpace);
    const hsize_t start = lo;
    const hsize_t count = n;
    if (H5Sselect_hyperslab(space.get(), H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0) {
        throw OSError("cannot select rows of " + path_);
    }
    const Handle mem_space(H5Screate_simple(1, &count, nullptr), Handle::Kind::DataSpace);
    if (H5Dwrite(id(), mem_type, mem_space.get(), space.get(), H5P_DEFAULT, data) < 0) {
        throw OSError("Can't write data (" + path_ + ")");
    }
}

void Dataset::write_column(std::size_t lo, const Column& column) {
    const std::size_t n = column.size();
    if (column.dtype() == DType::String) {
        const auto& strings = column.values<std::string>();
        const Handle type(H5Dget_type(id()), Handle::Kind::DataType);
        const std::size_t item = H5Tget_size(type.get());
        std::vector<char> buffer(item * n, '\0');
        for (std::size_t i = 0; i < n; ++i) {
            std::memcpy(buffer.data() + i * item, strings[i].data(),
                        std::min(item, strings[i].size()));
        }
        write_raw(lo, n, type.get(), buffer.data());
        return;
    }
    if (column.dtype() == DType::Categorical) {
        const auto& codes = column.categorical().codes;
        const Handle type(H5Dget_type(id()), Handle::Kind::DataType);
        if (H5Tget_class(type.get()) == H5T_ENUM) {
            const Handle native(H5Tget_native_type(type.get(), H5T_DIR_ASCEND),
                                Handle::Kind::DataType);
            if (H5Tget_size(native.get()) != 4) {
                throw TypeError("categorical codes need an int32 enumeration");
            }
            write_raw(lo, n, native.get(), codes.data());
        } else {
            write_raw(lo, n, H5T_NATIVE_INT32, codes.data());
        }
        return;
    }
    const Handle mem = memory_type(column.dtype());
    std::visit(
        [&](const auto& v) {
            using V = std::decay_t<decltype(v)>;
            if constexpr (!std::is_same_v<V, CategoricalData> &&
                          !std::is_same_v<V, std::vector<std::string>>) {
                write_raw(lo, n, mem.get(), v.data());
            }
        },
        column.data());
}

void Dataset::resize(std::size_t length) {
    const hsize_t dims = length;
    if (H5Dset_extent(id(), &dims) < 0) {
        throw ValueError("Unable to set extent of dataset (" + path_ + ")");
    }
}

// ---------------------------------------------------------------------------
// File

File::File(const std::string& path, Mode mode) : path_(path) {
    quiet();
    hid_t id = -1;
    const auto missing = [&path]() {
        return OSError("Unable to synchronously open file (unable to open file: name = '" +
                       path + "', errno = 2, error message = 'No such file or directory')");
    };
    switch (mode) {
        case Mode::Read:
            id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (id < 0) {
                if (!file_exists(path, nullptr)) {
                    throw missing();
                }
                throw OSError("Unable to synchronously open file (file signature not found)");
            }
            break;
        case Mode::ReadWrite:
            id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (id < 0) {
                if (!file_exists(path, nullptr)) {
                    throw missing();
                }
                throw OSError("Unable to synchronously open file (file signature not found)");
            }
            break;
        case Mode::Append: {
            id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (id < 0) {
                off_t size = 0;
                const bool exists = file_exists(path, &size);
                if (!exists) {
                    id = H5Fcreate(path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
                } else if (size == 0) {
                    // h5py opens an empty file in mode 'a' as a new HDF5 file;
                    // cooler's unordered writer relies on it for the temporary
                    // file tempfile.NamedTemporaryFile has already created.
                    id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                } else {
                    throw OSError("Unable to synchronously open file (file signature not found)");
                }
            }
            break;
        }
        case Mode::Truncate:
            id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            break;
        case Mode::Exclusive:
            id = H5Fcreate(path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
            if (id < 0 && file_exists(path, nullptr)) {
                throw OSError("Unable to synchronously create file (file exists)");
            }
            break;
    }
    if (id < 0) {
        throw OSError("Unable to synchronously create file (unable to open file: name = '" +
                      path + "')");
    }
    file_ = Handle(id, Handle::Kind::File);
}

bool File::exists(const std::string& path) const {
    if (path.empty() || path == "/") {
        return true;
    }
    std::string prefix;
    for (const std::string& part : split_path(path)) {
        prefix += "/" + part;
        if (H5Lexists(id(), prefix.c_str(), H5P_DEFAULT) <= 0) {
            return false;
        }
    }
    return H5Oexists_by_name(id(), prefix.c_str(), H5P_DEFAULT) > 0;
}

bool File::is_group(const std::string& path) const {
    if (!exists(path)) {
        return false;
    }
    H5O_info2_t info{};
    if (H5Oget_info_by_name3(id(), path.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) {
        return false;
    }
    return info.type == H5O_TYPE_GROUP;
}

bool File::is_dataset(const std::string& path) const {
    if (!exists(path)) {
        return false;
    }
    H5O_info2_t info{};
    if (H5Oget_info_by_name3(id(), path.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) {
        return false;
    }
    return info.type == H5O_TYPE_DATASET;
}

std::vector<std::string> File::keys(const std::string& group) const {
    const Handle handle(H5Gopen2(id(), group.c_str(), H5P_DEFAULT), Handle::Kind::Group);
    if (!handle.valid()) {
        throw KeyError("Unable to synchronously open object (component not found)");
    }
    H5_index_t index = H5_INDEX_NAME;
    const Handle gcpl(H5Gget_create_plist(handle.get()), Handle::Kind::PropertyList);
    unsigned flags = 0;
    if (gcpl.valid() && H5Pget_link_creation_order(gcpl.get(), &flags) >= 0 &&
        (flags & H5P_CRT_ORDER_TRACKED) != 0) {
        index = H5_INDEX_CRT_ORDER;
    }
    std::vector<std::string> names;
    hsize_t position = 0;
    if (H5Literate2(handle.get(), index, H5_ITER_INC, &position, collect_link, &names) < 0 &&
        index == H5_INDEX_CRT_ORDER) {
        names.clear();
        position = 0;
        H5Literate2(handle.get(), H5_INDEX_NAME, H5_ITER_INC, &position, collect_link, &names);
    }
    return names;
}

Dataset File::open_dataset(const std::string& path) const {
    if (!exists(path)) {
        throw KeyError("Unable to synchronously open object (component not found): " + path);
    }
    Handle handle(H5Dopen2(id(), path.c_str(), H5P_DEFAULT), Handle::Kind::Dataset);
    if (!handle.valid()) {
        throw KeyError("Unable to open dataset " + path);
    }
    return Dataset(std::move(handle), path);
}

void File::create_group(const std::string& path) {
    if (exists(path)) {
        throw ValueError("Unable to synchronously create group (name already exists)");
    }
    const Handle lcpl = link_create_plist();
    const Handle group(H5Gcreate2(id(), path.c_str(), lcpl.get(), H5P_DEFAULT, H5P_DEFAULT),
                       Handle::Kind::Group);
    if (!group.valid()) {
        throw ValueError("Unable to synchronously create group (" + path + ")");
    }
}

void File::remove(const std::string& path) {
    if (!exists(path)) {
        throw KeyError("Couldn't delete link (" + path + ")");
    }
    if (H5Ldelete(id(), path.c_str(), H5P_DEFAULT) < 0) {
        throw KeyError("Couldn't delete link (" + path + ")");
    }
}

Dataset File::create_dataset(const std::string& path, hid_t type,
                             const DatasetCreate& options) {
    const hsize_t dims = options.length;
    hsize_t maxdims = dims;
    if (options.maxlength.has_value()) {
        maxdims = *options.maxlength == DatasetCreate::kUnlimited
                      ? H5S_UNLIMITED
                      : static_cast<hsize_t>(*options.maxlength);
    }
    const Handle space(H5Screate_simple(1, &dims, &maxdims), Handle::Kind::DataSpace);
    if (!space.valid()) {
        throw ValueError("Unable to create dataspace for " + path);
    }
    const Handle dcpl(H5Pcreate(H5P_DATASET_CREATE), Handle::Kind::PropertyList);
    const auto fail = [&path](const char* what) {
        return ValueError(std::string("Unable to synchronously create dataset (") + what +
                          ") " + path);
    };
    if (options.chunk.has_value()) {
        const hsize_t chunk = *options.chunk;
        if (H5Pset_chunk(dcpl.get(), 1, &chunk) < 0) {
            throw fail("chunk");
        }
    }
    if (options.scaleoffset.has_value()) {
        const H5T_class_t cls = H5Tget_class(type);
        const H5Z_SO_scale_type_t scale =
            (cls == H5T_FLOAT) ? H5Z_SO_FLOAT_DSCALE : H5Z_SO_INT;
        if (H5Pset_scaleoffset(dcpl.get(), scale, static_cast<int>(*options.scaleoffset)) < 0) {
            throw fail("scaleoffset");
        }
    }
    if (options.shuffle && H5Pset_shuffle(dcpl.get()) < 0) {
        throw fail("shuffle");
    }
    switch (options.compression) {
        case DatasetCreate::Compression::None: break;
        case DatasetCreate::Compression::Gzip:
            if (H5Pset_deflate(dcpl.get(), options.gzip_level) < 0) {
                throw fail("gzip");
            }
            break;
        case DatasetCreate::Compression::Szip:
            if (H5Pset_szip(dcpl.get(), options.szip_mask, options.szip_pixels) < 0) {
                throw fail("szip");
            }
            break;
        case DatasetCreate::Compression::FilterId:
            if (H5Pset_filter(dcpl.get(), static_cast<H5Z_filter_t>(options.filter_id),
                              H5Z_FLAG_OPTIONAL, options.filter_options.size(),
                              options.filter_options.data()) < 0) {
                throw fail("filter");
            }
            break;
    }
    if (options.fletcher32 && H5Pset_fletcher32(dcpl.get()) < 0) {
        throw fail("fletcher32");
    }
    if (options.fillvalue.has_value()) {
        const TypeDesc desc = describe_type(type);
        const json::Value& fill = *options.fillvalue;
        if (desc.kind == TypeDesc::Kind::Float) {
            const double v = fill.as_double();
            if (H5Pset_fill_value(dcpl.get(), H5T_NATIVE_DOUBLE, &v) < 0) {
                throw fail("fillvalue");
            }
        } else if (desc.kind == TypeDesc::Kind::Integer || desc.kind == TypeDesc::Kind::Enum ||
                   desc.kind == TypeDesc::Kind::Bool) {
            if (desc.kind == TypeDesc::Kind::Integer) {
                const std::int64_t v = fill.as_int();
                if (H5Pset_fill_value(dcpl.get(), H5T_NATIVE_INT64, &v) < 0) {
                    throw fail("fillvalue");
                }
            } else {
                const Handle native(H5Tget_native_type(type, H5T_DIR_ASCEND),
                                    Handle::Kind::DataType);
                std::vector<unsigned char> raw(H5Tget_size(native.get()));
                const std::int64_t v = fill.as_int();
                std::memcpy(raw.data(), &v, std::min<std::size_t>(raw.size(), 8));
                if (H5Pset_fill_value(dcpl.get(), native.get(), raw.data()) < 0) {
                    throw fail("fillvalue");
                }
            }
        }
    }
    if (H5Pset_obj_track_times(dcpl.get(), options.track_times ? 1 : 0) < 0) {
        throw fail("track_times");
    }
    const Handle lcpl = link_create_plist();
    Handle dataset(H5Dcreate2(id(), path.c_str(), type, space.get(), lcpl.get(), dcpl.get(),
                              H5P_DEFAULT),
                   Handle::Kind::Dataset);
    if (!dataset.valid()) {
        throw ValueError("Unable to synchronously create dataset (" + path + ")");
    }
    return Dataset(std::move(dataset), path);
}

std::vector<std::pair<std::string, json::Value>> File::attributes(
    const std::string& path) const {
    const Handle object(H5Oopen(id(), path.c_str(), H5P_DEFAULT), Handle::Kind::Object);
    if (!object.valid()) {
        throw KeyError("Unable to synchronously open object (component not found)");
    }
    H5_index_t index = H5_INDEX_NAME;
    {
        H5O_info2_t info{};
        H5Oget_info3(object.get(), &info, H5O_INFO_BASIC);
        Handle plist;
        if (info.type == H5O_TYPE_GROUP) {
            plist = Handle(H5Gget_create_plist(object.get()), Handle::Kind::PropertyList);
        } else if (info.type == H5O_TYPE_DATASET) {
            plist = Handle(H5Dget_create_plist(object.get()), Handle::Kind::PropertyList);
        }
        unsigned flags = 0;
        if (plist.valid() && H5Pget_attr_creation_order(plist.get(), &flags) >= 0 &&
            (flags & H5P_CRT_ORDER_TRACKED) != 0) {
            index = H5_INDEX_CRT_ORDER;
        }
    }
    H5O_info2_t info{};
    if (H5Oget_info3(object.get(), &info, H5O_INFO_NUM_ATTRS) < 0) {
        throw OSError("cannot count the attributes of " + path);
    }
    std::vector<std::pair<std::string, json::Value>> out;
    for (hsize_t i = 0; i < info.num_attrs; ++i) {
        Handle attr(H5Aopen_by_idx(object.get(), ".", index, H5_ITER_INC, i, H5P_DEFAULT,
                                   H5P_DEFAULT),
                    Handle::Kind::Attribute);
        if (!attr.valid() && index == H5_INDEX_CRT_ORDER) {
            attr = Handle(H5Aopen_by_idx(object.get(), ".", H5_INDEX_NAME, H5_ITER_INC, i,
                                         H5P_DEFAULT, H5P_DEFAULT),
                          Handle::Kind::Attribute);
        }
        if (!attr.valid()) {
            continue;
        }
        const ssize_t name_length = H5Aget_name(attr.get(), 0, nullptr);
        std::string name(static_cast<std::size_t>(std::max<ssize_t>(name_length, 0)), '\0');
        H5Aget_name(attr.get(), name.size() + 1, name.data());
        const Handle type(H5Aget_type(attr.get()), Handle::Kind::DataType);
        const Handle space(H5Aget_space(attr.get()), Handle::Kind::DataSpace);
        const H5S_class_t space_class = H5Sget_simple_extent_type(space.get());
        if (space_class == H5S_NULL) {
            out.emplace_back(std::move(name), json::Value());
            continue;
        }
        const hssize_t points = H5Sget_simple_extent_npoints(space.get());
        std::vector<json::Value> values =
            read_attribute_values(attr.get(), type.get(), static_cast<std::size_t>(points));
        if (space_class == H5S_SCALAR) {
            out.emplace_back(std::move(name), std::move(values.at(0)));
        } else {
            out.emplace_back(std::move(name), json::Value::array(std::move(values)));
        }
    }
    return out;
}

void File::set_attribute(const std::string& path, const std::string& name,
                         const json::Value& value) {
    const Handle object(H5Oopen(id(), path.c_str(), H5P_DEFAULT), Handle::Kind::Object);
    if (!object.valid()) {
        throw KeyError("Unable to synchronously open object (component not found)");
    }
    if (H5Aexists(object.get(), name.c_str()) > 0) {
        H5Adelete(object.get(), name.c_str());
    }
    const Handle scalar(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    const auto create = [&](hid_t type, hid_t space) {
        Handle attr(H5Acreate2(object.get(), name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT),
                    Handle::Kind::Attribute);
        if (!attr.valid()) {
            throw OSError("Unable to synchronously create attribute (" + name + ")");
        }
        return attr;
    };
    switch (value.type()) {
        case json::Type::String: {
            Handle type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
            H5Tset_size(type.get(), H5T_VARIABLE);
            H5Tset_cset(type.get(), H5T_CSET_UTF8);
            const Handle attr = create(type.get(), scalar.get());
            const char* text = value.as_string().c_str();
            if (H5Awrite(attr.get(), type.get(), &text) < 0) {
                throw OSError("cannot write attribute " + name);
            }
            return;
        }
        case json::Type::Bytes: {
            const Handle type = fixed_string_type(std::max<std::size_t>(value.as_string().size(), 1));
            H5Tset_strpad(type.get(), H5T_STR_NULLPAD);
            const Handle attr = create(type.get(), scalar.get());
            std::string buffer = value.as_string();
            buffer.resize(std::max<std::size_t>(buffer.size(), 1), '\0');
            if (H5Awrite(attr.get(), type.get(), buffer.data()) < 0) {
                throw OSError("cannot write attribute " + name);
            }
            return;
        }
        case json::Type::Bool: {
            const Handle type = file_type(DType::Bool);
            const Handle attr = create(type.get(), scalar.get());
            const std::int8_t v = value.as_bool() ? 1 : 0;
            if (H5Awrite(attr.get(), type.get(), &v) < 0) {
                throw OSError("cannot write attribute " + name);
            }
            return;
        }
        case json::Type::Int:
        case json::Type::UInt:
        case json::Type::Double: {
            DType dtype = value.dtype();
            if (!is_numeric(dtype) || dtype == DType::Bool) {
                dtype = value.type() == json::Type::Double ? DType::Float64 : DType::Int64;
            }
            const Handle type = file_type(dtype);
            const Handle attr = create(type.get(), scalar.get());
            herr_t status = -1;
            if (value.type() == json::Type::Double) {
                const double v = value.as_double();
                status = H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &v);
            } else if (value.type() == json::Type::Int) {
                const std::int64_t v = value.as_int();
                status = H5Awrite(attr.get(), H5T_NATIVE_INT64, &v);
            } else {
                const std::uint64_t v = value.as_uint();
                status = H5Awrite(attr.get(), H5T_NATIVE_UINT64, &v);
            }
            if (status < 0) {
                throw OSError("cannot write attribute " + name);
            }
            return;
        }
        case json::Type::Array: {
            const auto& items = value.as_array();
            const hsize_t n = items.size();
            const Handle space(H5Screate_simple(1, &n, nullptr), Handle::Kind::DataSpace);
            if (!items.empty() && items[0].is_number()) {
                DType dtype = items[0].dtype();
                const bool floating = items[0].type() == json::Type::Double;
                const Handle type = file_type(dtype);
                const Handle attr = create(type.get(), space.get());
                if (floating) {
                    std::vector<double> v;
                    for (const auto& item : items) v.push_back(item.as_double());
                    H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, v.data());
                } else {
                    std::vector<std::int64_t> v;
                    for (const auto& item : items) v.push_back(item.as_int());
                    H5Awrite(attr.get(), H5T_NATIVE_INT64, v.data());
                }
                return;
            }
            throw TypeError("only numeric array attributes are supported");
        }
        case json::Type::Null:
        case json::Type::Object: break;
    }
    throw TypeError("Object dtype dtype('O') has no native HDF5 equivalent");
}

void File::flush() {
    if (H5Fflush(id(), H5F_SCOPE_LOCAL) < 0) {
        throw OSError("cannot flush " + path_);
    }
}

}  // namespace coolercpp::h5
