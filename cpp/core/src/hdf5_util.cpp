#include "hicx/hdf5_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>

#include <blosc.h>

namespace hicx::h5 {

namespace {

constexpr H5Z_filter_t kBloscFilterId = 32001;
// The filter revision that PyTables writes into cd_values[0]. Keeping it makes
// the parameter block of a dataset written here identical to PyTables'.
constexpr unsigned int kBloscFilterVersion = 2;
constexpr int kPyTablesComplevel = 5;
constexpr int kPyTablesShuffle = 1;

// Both halves of the reference hdf5-blosc filter.
//
//   cd_values[0] filter revision      cd_values[3] uncompressed chunk bytes
//   cd_values[1] blosc format version cd_values[4] compression level
//   cd_values[2] type size            cd_values[5] shuffle flag
//
// The first four are filled in by blosc_set_local below, so a caller only has
// to supply the level and the shuffle flag, exactly as PyTables does.
size_t blosc_filter_impl(unsigned int flags, size_t cd_nelmts,
                         const unsigned int cd_values[], size_t nbytes,
                         size_t* buf_size, void** buf) {
    size_t outbuf_size = 0;
    if ((flags & H5Z_FLAG_REVERSE) == 0) {
        const size_t typesize = cd_nelmts > 2 ? cd_values[2] : 8;
        const int clevel =
            cd_nelmts >= 5 ? static_cast<int>(cd_values[4]) : kPyTablesComplevel;
        const int doshuffle =
            cd_nelmts >= 6 ? static_cast<int>(cd_values[5]) : kPyTablesShuffle;
        outbuf_size = *buf_size;
        void* outbuf = std::malloc(outbuf_size + BLOSC_MAX_OVERHEAD);
        if (outbuf == nullptr) {
            return 0;
        }
        // blosclz is what PyTables selects for complib='blosc'.
        if (blosc_set_compressor("blosclz") < 0) {
            std::free(outbuf);
            return 0;
        }
        const int status = blosc_compress(clevel, doshuffle, typesize, nbytes, *buf,
                                          outbuf, nbytes + BLOSC_MAX_OVERHEAD);
        if (status <= 0) {
            // A negative status is an error; zero means the data did not
            // compress, and returning zero makes HDF5 fall back to storing the
            // chunk uncompressed because the filter is optional.
            std::free(outbuf);
            return 0;
        }
        std::free(*buf);
        *buf = outbuf;
        *buf_size = outbuf_size;
        return static_cast<size_t>(status);
    }

    (void)cd_nelmts;
    (void)cd_values;
    (void)nbytes;
    size_t cbytes = 0;
    size_t blocksize = 0;
    blosc_cbuffer_sizes(*buf, &outbuf_size, &cbytes, &blocksize);
    if (outbuf_size == 0) {
        return 0;
    }
    void* outbuf = std::malloc(outbuf_size);
    if (outbuf == nullptr) {
        return 0;
    }
    const int status = blosc_decompress(*buf, outbuf, outbuf_size);
    if (status <= 0) {
        std::free(outbuf);
        return 0;
    }
    std::free(*buf);
    *buf = outbuf;
    *buf_size = outbuf_size;
    return static_cast<size_t>(status);
}

// Fills in the parameters that depend on the dataset rather than on the
// caller: the type size and the uncompressed size of one chunk.
herr_t blosc_set_local(hid_t dcpl, hid_t type, hid_t /*space*/) {
    unsigned int flags = 0;
    size_t nelmts = 8;
    unsigned int values[8] = {0};
    if (H5Pget_filter_by_id2(dcpl, kBloscFilterId, &flags, &nelmts, values, 0, nullptr,
                             nullptr) < 0) {
        return -1;
    }
    hsize_t chunk_dims[H5S_MAX_RANK];
    const int rank = H5Pget_chunk(dcpl, H5S_MAX_RANK, chunk_dims);
    if (rank < 0) {
        return -1;
    }
    const size_t typesize = H5Tget_size(type);
    size_t chunk_bytes = typesize;
    for (int i = 0; i < rank; ++i) {
        chunk_bytes *= static_cast<size_t>(chunk_dims[i]);
    }
    if (nelmts < 5) {
        values[4] = kPyTablesComplevel;
    }
    if (nelmts < 6) {
        values[5] = kPyTablesShuffle;
    }
    values[0] = kBloscFilterVersion;
    values[1] = BLOSC_VERSION_FORMAT;
    values[2] = static_cast<unsigned int>(typesize);
    values[3] = static_cast<unsigned int>(chunk_bytes);
    return H5Pmodify_filter(dcpl, kBloscFilterId, flags, 6, values);
}

const H5Z_class2_t kBloscClass = {
    H5Z_CLASS_T_VERS,
    kBloscFilterId,
    1,  // encoder_present
    1,  // decoder_present
    "blosc",
    nullptr,
    blosc_set_local,
    blosc_filter_impl,
};

std::string dtype_name(hid_t type_id) {
    const H5T_class_t cls = H5Tget_class(type_id);
    const size_t size = H5Tget_size(type_id);
    if (cls == H5T_INTEGER) {
        const H5T_sign_t sign = H5Tget_sign(type_id);
        const std::string prefix = (sign == H5T_SGN_NONE) ? "uint" : "int";
        return prefix + std::to_string(size * 8);
    }
    if (cls == H5T_FLOAT) {
        return "float" + std::to_string(size * 8);
    }
    if (cls == H5T_STRING) {
        return "string";
    }
    return "unknown";
}

herr_t collect_name(hid_t /*group*/, const char* name, const H5L_info2_t* /*info*/,
                    void* op_data) {
    auto* names = static_cast<std::vector<std::string>*>(op_data);
    names->emplace_back(name);
    return 0;
}

}  // namespace

void register_blosc_filter() {
    static std::once_flag once;
    std::call_once(once, [] {
        blosc_init();
        if (H5Zfilter_avail(kBloscFilterId) <= 0) {
            if (H5Zregister(&kBloscClass) < 0) {
                throw Error("failed to register the blosc HDF5 filter");
            }
        }
    });
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
    }
    id_ = -1;
}

bool is_hdf5(const std::string& path) {
    H5E_auto2_t old_func = nullptr;
    void* old_data = nullptr;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    const htri_t result = H5Fis_accessible(path.c_str(), H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
    return result > 0;
}

File::File(const std::string& path) : path_(path) {
    register_blosc_filter();
    // Keep the HDF5 error stack quiet; we translate failures into exceptions.
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    const hid_t id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (id < 0) {
        throw Error("cannot open HDF5 file: " + path);
    }
    file_ = Handle(id, Handle::Kind::File);
}

bool File::exists(const std::string& object_path) const {
    if (object_path == "/" || object_path.empty()) {
        return true;
    }
    return H5Lexists(file_.get(), object_path.c_str(), H5P_DEFAULT) > 0;
}

std::map<std::string, AttributeValue> File::attributes(const std::string& object_path) const {
    std::map<std::string, AttributeValue> result;
    const hid_t obj = H5Oopen(file_.get(), object_path.c_str(), H5P_DEFAULT);
    if (obj < 0) {
        throw Error("cannot open object " + object_path + " in " + path_);
    }
    const Handle obj_handle(obj, Handle::Kind::Group);

    H5O_info2_t info{};
    if (H5Oget_info3(obj, &info, H5O_INFO_NUM_ATTRS) < 0) {
        throw Error("cannot read attribute count of " + object_path);
    }
    for (hsize_t i = 0; i < info.num_attrs; ++i) {
        const hid_t attr = H5Aopen_by_idx(obj, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                                          H5P_DEFAULT, H5P_DEFAULT);
        if (attr < 0) {
            continue;
        }
        const Handle attr_handle(attr, Handle::Kind::Attribute);
        char name[512];
        const ssize_t name_length = H5Aget_name(attr, sizeof(name), name);
        if (name_length <= 0) {
            continue;
        }
        const Handle type(H5Aget_type(attr), Handle::Kind::DataType);
        const H5T_class_t cls = H5Tget_class(type.get());
        if (cls == H5T_STRING) {
            if (H5Tis_variable_str(type.get()) > 0) {
                char* value = nullptr;
                const Handle mem_type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
                H5Tset_size(mem_type.get(), H5T_VARIABLE);
                H5Tset_cset(mem_type.get(), H5Tget_cset(type.get()));
                if (H5Aread(attr, mem_type.get(), &value) >= 0 && value != nullptr) {
                    result.emplace(name, std::string(value));
                    H5free_memory(value);
                }
            } else {
                const size_t size = H5Tget_size(type.get());
                std::string buffer(size, '\0');
                if (H5Aread(attr, type.get(), buffer.data()) >= 0) {
                    const std::size_t nul = buffer.find('\0');
                    if (nul != std::string::npos) {
                        buffer.resize(nul);
                    }
                    result.emplace(name, buffer);
                }
            }
        } else if (cls == H5T_INTEGER) {
            std::int64_t value = 0;
            if (H5Aread(attr, H5T_NATIVE_INT64, &value) >= 0) {
                result.emplace(name, value);
            }
        } else if (cls == H5T_FLOAT) {
            double value = 0.0;
            if (H5Aread(attr, H5T_NATIVE_DOUBLE, &value) >= 0) {
                result.emplace(name, value);
            }
        }
    }
    return result;
}

std::vector<std::string> File::children(const std::string& group_path) const {
    std::vector<std::string> names;
    const hid_t group = H5Gopen2(file_.get(), group_path.c_str(), H5P_DEFAULT);
    if (group < 0) {
        throw Error("cannot open group " + group_path + " in " + path_);
    }
    const Handle group_handle(group, Handle::Kind::Group);
    hsize_t index = 0;
    H5Literate2(group, H5_INDEX_NAME, H5_ITER_INC, &index, collect_name, &names);
    return names;
}

Handle File::open_dataset(const std::string& dataset_path) const {
    const hid_t dataset = H5Dopen2(file_.get(), dataset_path.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        throw Error("cannot open dataset " + dataset_path + " in " + path_);
    }
    return Handle(dataset, Handle::Kind::Dataset);
}

std::size_t File::dataset_length(const std::string& dataset_path) const {
    const Handle dataset = open_dataset(dataset_path);
    const Handle space(H5Dget_space(dataset.get()), Handle::Kind::DataSpace);
    const hssize_t points = H5Sget_simple_extent_npoints(space.get());
    if (points < 0) {
        throw Error("cannot determine the size of " + dataset_path);
    }
    return static_cast<std::size_t>(points);
}

bool File::dataset_is_column(const std::string& dataset_path) const {
    const Handle dataset = open_dataset(dataset_path);
    const Handle space(H5Dget_space(dataset.get()), Handle::Kind::DataSpace);
    hsize_t extent[2] = {0, 0};
    const int rank = H5Sget_simple_extent_dims(space.get(), extent, nullptr);
    return rank == 2 && extent[1] == 1;
}

std::string File::dataset_dtype(const std::string& dataset_path) const {
    const Handle dataset = open_dataset(dataset_path);
    const Handle type(H5Dget_type(dataset.get()), Handle::Kind::DataType);
    return dtype_name(type.get());
}

namespace {

template <typename T>
std::vector<T> read_typed(const File& file, const std::string& dataset_path,
                          hid_t mem_type, std::size_t length) {
    std::vector<T> values(length);
    if (length == 0) {
        return values;
    }
    const hid_t dataset = H5Dopen2(file.id(), dataset_path.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        throw Error("cannot open dataset " + dataset_path);
    }
    const Handle dataset_handle(dataset, Handle::Kind::Dataset);
    if (H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
        throw Error("cannot read dataset " + dataset_path +
                    " (missing HDF5 filter, or type conversion failed)");
    }
    return values;
}

}  // namespace

std::vector<double> File::read_doubles(const std::string& dataset_path) const {
    return read_typed<double>(*this, dataset_path, H5T_NATIVE_DOUBLE,
                              dataset_length(dataset_path));
}

std::vector<std::int64_t> File::read_int64(const std::string& dataset_path) const {
    return read_typed<std::int64_t>(*this, dataset_path, H5T_NATIVE_INT64,
                                    dataset_length(dataset_path));
}

std::vector<std::int32_t> File::read_int32(const std::string& dataset_path) const {
    return read_typed<std::int32_t>(*this, dataset_path, H5T_NATIVE_INT32,
                                    dataset_length(dataset_path));
}

// --------------------------------------------------------------------------
// Writing

std::size_t guess_chunk(std::size_t length, std::size_t typesize) {
    // h5py/_hl/filters.py guess_chunk, one dimensional case. The constants and
    // the loop are copied rather than approximated because the resulting chunk
    // shape is compared byte for byte against the Python written cool files.
    constexpr double kChunkBase = 16.0 * 1024.0;
    constexpr double kChunkMin = 8.0 * 1024.0;
    constexpr double kChunkMax = 1024.0 * 1024.0;

    double chunk = length != 0 ? static_cast<double>(length) : 1024.0;
    const double element = static_cast<double>(typesize);
    const double dataset_bytes = chunk * element;
    double target = kChunkBase * std::pow(2.0, std::log10(dataset_bytes / (1024.0 * 1024.0)));
    target = std::min(target, kChunkMax);
    target = std::max(target, kChunkMin);

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

Handle fixed_string_type(std::size_t width) {
    Handle type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
    if (!type.valid() || H5Tset_size(type.get(), std::max<std::size_t>(width, 1)) < 0 ||
        H5Tset_strpad(type.get(), H5T_STR_NULLPAD) < 0 ||
        H5Tset_cset(type.get(), H5T_CSET_ASCII) < 0) {
        throw Error("cannot build a fixed width string type");
    }
    return type;
}

Handle enum_type(const std::vector<std::string>& names, hid_t base) {
    Handle type(H5Tenum_create(base), Handle::Kind::DataType);
    if (!type.valid()) {
        throw Error("cannot build an enumeration type");
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        const std::int32_t value = static_cast<std::int32_t>(i);
        if (H5Tenum_insert(type.get(), names[i].c_str(), &value) < 0) {
            throw Error("cannot add " + names[i] + " to the chromosome enumeration");
        }
    }
    return type;
}

FileWriter::FileWriter(const std::string& path, WriteMode mode) : path_(path) {
    register_blosc_filter();
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    hid_t id = -1;
    if (mode == WriteMode::Append) {
        // h5py.File(path, 'a'): open an existing file for writing, create it
        // when it is not there.
        id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (id < 0) {
            id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
    } else {
        id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (id < 0) {
        throw Error("cannot create HDF5 file: " + path);
    }
    file_ = Handle(id, Handle::Kind::File);
}

Handle FileWriter::create_group(const std::string& path) {
    // Missing parents are created too, which is what an mcool needs: the
    // second resolution finds /resolutions already there, the first does not.
    const Handle link_plist(H5Pcreate(H5P_LINK_CREATE), Handle::Kind::PropertyList);
    if (!link_plist.valid() ||
        H5Pset_create_intermediate_group(link_plist.get(), 1) < 0) {
        throw Error("cannot configure the link creation of group " + path);
    }
    // Object headers carry an optional modification time (message 0x12), which
    // HDF5 writes by default. Two runs a second apart then produce different
    // bytes for the same content, which makes the reproducibility check of
    // cpp/OPTIMIZATION.md section 3 a coin toss. Nothing in the corpus reads
    // the field and no comparator compares it, so it is cleared.
    const Handle group_plist(H5Pcreate(H5P_GROUP_CREATE), Handle::Kind::PropertyList);
    if (!group_plist.valid() || H5Pset_obj_track_times(group_plist.get(), 0) < 0) {
        throw Error("cannot disable object time tracking on group " + path);
    }
    const hid_t group = H5Gcreate2(file_.get(), path.c_str(), link_plist.get(),
                                   group_plist.get(), H5P_DEFAULT);
    if (group < 0) {
        throw Error("cannot create group " + path + " in " + path_);
    }
    return Handle(group, Handle::Kind::Group);
}

bool FileWriter::exists(const std::string& object_path) const {
    if (object_path.empty() || object_path == "/") {
        return true;
    }
    // H5Lexists only answers for one link at a time, so every component of the
    // path has to be probed in turn.
    std::string prefix;
    std::size_t position = 0;
    while (position < object_path.size()) {
        const std::size_t slash = object_path.find('/', position);
        const std::size_t end = slash == std::string::npos ? object_path.size() : slash;
        const std::string component = object_path.substr(position, end - position);
        position = end + 1;
        if (component.empty()) {
            continue;
        }
        prefix += "/" + component;
        if (H5Lexists(file_.get(), prefix.c_str(), H5P_DEFAULT) <= 0) {
            return false;
        }
    }
    return true;
}

void FileWriter::unlink(const std::string& object_path) {
    if (!exists(object_path) || object_path == "/") {
        return;
    }
    if (H5Ldelete(file_.get(), object_path.c_str(), H5P_DEFAULT) < 0) {
        throw Error("cannot remove " + object_path + " from " + path_);
    }
}

Handle FileWriter::create_dataset(const std::string& path, hid_t file_type,
                                  std::size_t length, std::size_t max_length,
                                  Filter filter, std::size_t chunk,
                                  std::size_t minor) {
    const int rank = minor > 0 ? 2 : 1;
    const hsize_t dims[2] = {static_cast<hsize_t>(length),
                             static_cast<hsize_t>(minor)};
    const hsize_t maxdims[2] = {max_length == kUnlimited
                                    ? H5S_UNLIMITED
                                    : static_cast<hsize_t>(std::max(max_length, length)),
                                static_cast<hsize_t>(minor)};
    const Handle space(H5Screate_simple(rank, dims, maxdims), Handle::Kind::DataSpace);
    if (!space.valid()) {
        throw Error("cannot create the dataspace of " + path);
    }

    const Handle plist(H5Pcreate(H5P_DATASET_CREATE), Handle::Kind::PropertyList);
    if (!plist.valid()) {
        throw Error("cannot create the property list of " + path);
    }
    // See create_group: the object modification time is the only part of a
    // dataset that changes between two identical runs.
    if (H5Pset_obj_track_times(plist.get(), 0) < 0) {
        throw Error("cannot disable object time tracking on " + path);
    }
    if (chunk == 0) {
        chunk = guess_chunk(length, H5Tget_size(file_type));
    }
    const hsize_t chunk_dims[2] = {static_cast<hsize_t>(std::max<std::size_t>(chunk, 1)),
                                   static_cast<hsize_t>(std::max<std::size_t>(minor, 1))};
    if (H5Pset_chunk(plist.get(), rank, chunk_dims) < 0) {
        throw Error("cannot set the chunk shape of " + path);
    }
    switch (filter) {
        case Filter::None:
            break;
        case Filter::CoolerDefault:
            if (H5Pset_shuffle(plist.get()) < 0 || H5Pset_deflate(plist.get(), 6) < 0) {
                throw Error("cannot set the shuffle and gzip filters of " + path);
            }
            break;
        case Filter::CoolerColumn:
            if (H5Pset_deflate(plist.get(), 6) < 0) {
                throw Error("cannot set the gzip filter of " + path);
            }
            break;
        case Filter::PyTablesBlosc: {
            // Only the level and the shuffle flag are given here; the type size
            // and the chunk size are filled in by blosc_set_local.
            const unsigned int cd_values[6] = {0, 0, 0, 0, kPyTablesComplevel,
                                               kPyTablesShuffle};
            if (H5Pset_filter(plist.get(), kBloscFilterId, H5Z_FLAG_OPTIONAL, 6,
                              cd_values) < 0) {
                throw Error("cannot set the blosc filter of " + path);
            }
            break;
        }
    }

    const hid_t dataset = H5Dcreate2(file_.get(), path.c_str(), file_type, space.get(),
                                     H5P_DEFAULT, plist.get(), H5P_DEFAULT);
    if (dataset < 0) {
        throw Error("cannot create dataset " + path + " in " + path_);
    }
    return Handle(dataset, Handle::Kind::Dataset);
}

void FileWriter::resize(hid_t dataset, std::size_t length) {
    const hsize_t dims = static_cast<hsize_t>(length);
    if (H5Dset_extent(dataset, &dims) < 0) {
        throw Error("cannot resize a dataset");
    }
}

void FileWriter::write_block(hid_t dataset, hid_t mem_type, std::size_t offset,
                             std::size_t count, const void* data) {
    if (count == 0) {
        return;
    }
    const Handle file_space(H5Dget_space(dataset), Handle::Kind::DataSpace);
    hsize_t extent[2] = {0, 0};
    const int rank = H5Sget_simple_extent_dims(file_space.get(), extent, nullptr);
    if (rank < 1 || rank > 2) {
        throw Error("only one and two dimensional datasets are written");
    }
    const hsize_t start[2] = {static_cast<hsize_t>(offset), 0};
    const hsize_t block[2] = {static_cast<hsize_t>(count), extent[1]};
    if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, block,
                            nullptr) < 0) {
        throw Error("cannot select the destination of a dataset write");
    }
    const Handle mem_space(H5Screate_simple(rank, block, block),
                           Handle::Kind::DataSpace);
    if (H5Dwrite(dataset, mem_type, mem_space.get(), file_space.get(), H5P_DEFAULT,
                 data) < 0) {
        char name[512] = {0};
        H5Iget_name(dataset, name, sizeof(name));
        throw Error(std::string("cannot write to dataset ") + name);
    }
}

void FileWriter::set_attribute(const std::string& object_path, const std::string& name,
                               const AttributeValue& value) {
    const hid_t object = H5Oopen(file_.get(), object_path.c_str(), H5P_DEFAULT);
    if (object < 0) {
        throw Error("cannot open object " + object_path + " in " + path_);
    }
    const Handle object_handle(object, Handle::Kind::Group);
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);

    Handle type;
    const void* buffer = nullptr;
    const char* text = nullptr;
    std::int64_t integer = 0;
    double number = 0.0;
    if (std::holds_alternative<std::string>(value)) {
        // h5py writes a Python str as a variable length UTF-8 string.
        type = Handle(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
        H5Tset_size(type.get(), H5T_VARIABLE);
        H5Tset_cset(type.get(), H5T_CSET_UTF8);
        text = std::get<std::string>(value).c_str();
        buffer = &text;
    } else if (std::holds_alternative<std::int64_t>(value)) {
        type = Handle(H5Tcopy(H5T_STD_I64LE), Handle::Kind::DataType);
        integer = std::get<std::int64_t>(value);
        buffer = &integer;
    } else {
        type = Handle(H5Tcopy(H5T_IEEE_F64LE), Handle::Kind::DataType);
        number = std::get<double>(value);
        buffer = &number;
    }

    // h5py assignment replaces an existing attribute, and hicmatrix relies on
    // that: cooler writes 'format', 'format-url' and 'generated-by' and
    // Cool.save overwrites them (hicmatrix/lib/cool.py:422-426).
    if (H5Aexists(object, name.c_str()) > 0) {
        H5Adelete(object, name.c_str());
    }
    const hid_t attribute = H5Acreate2(object, name.c_str(), type.get(), space.get(),
                                       H5P_DEFAULT, H5P_DEFAULT);
    if (attribute < 0) {
        throw Error("cannot create attribute " + name + " on " + object_path);
    }
    const Handle attribute_handle(attribute, Handle::Kind::Attribute);
    Handle mem_type(H5Tcopy(type.get()), Handle::Kind::DataType);
    if (std::holds_alternative<std::int64_t>(value)) {
        mem_type = Handle(H5Tcopy(H5T_NATIVE_INT64), Handle::Kind::DataType);
    } else if (std::holds_alternative<double>(value)) {
        mem_type = Handle(H5Tcopy(H5T_NATIVE_DOUBLE), Handle::Kind::DataType);
    }
    if (H5Awrite(attribute, mem_type.get(), buffer) < 0) {
        throw Error("cannot write attribute " + name + " on " + object_path);
    }
}

void FileWriter::set_bytes_attribute(const std::string& object_path,
                                     const std::string& name, const std::string& value,
                                     bool null_dataspace) {
    const hid_t object = H5Oopen(file_.get(), object_path.c_str(), H5P_DEFAULT);
    if (object < 0) {
        throw Error("cannot open object " + object_path + " in " + path_);
    }
    const Handle object_handle(object, Handle::Kind::Group);
    const Handle type = fixed_string_type(std::max<std::size_t>(value.size(), 1));
    // numpy byte strings are NUL terminated when they are shorter than the
    // type, which is the padding PyTables writes for these markers.
    H5Tset_strpad(type.get(), H5T_STR_NULLTERM);
    const Handle space(H5Screate(null_dataspace ? H5S_NULL : H5S_SCALAR),
                       Handle::Kind::DataSpace);
    // h5py assignment replaces an existing attribute, and hicmatrix relies on
    // that: cooler writes 'format', 'format-url' and 'generated-by' and
    // Cool.save overwrites them (hicmatrix/lib/cool.py:422-426).
    if (H5Aexists(object, name.c_str()) > 0) {
        H5Adelete(object, name.c_str());
    }
    const hid_t attribute = H5Acreate2(object, name.c_str(), type.get(), space.get(),
                                       H5P_DEFAULT, H5P_DEFAULT);
    if (attribute < 0) {
        throw Error("cannot create attribute " + name + " on " + object_path);
    }
    const Handle attribute_handle(attribute, Handle::Kind::Attribute);
    if (!null_dataspace && H5Awrite(attribute, type.get(), value.data()) < 0) {
        throw Error("cannot write attribute " + name + " on " + object_path);
    }
}

std::vector<std::string> File::read_strings(const std::string& dataset_path) const {
    const std::size_t length = dataset_length(dataset_path);
    std::vector<std::string> out;
    out.reserve(length);
    if (length == 0) {
        return out;
    }
    const Handle dataset = open_dataset(dataset_path);
    const Handle type(H5Dget_type(dataset.get()), Handle::Kind::DataType);
    if (H5Tget_class(type.get()) != H5T_STRING) {
        throw Error(dataset_path + " is not a string dataset");
    }
    if (H5Tis_variable_str(type.get()) > 0) {
        std::vector<char*> raw(length, nullptr);
        const Handle mem_type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
        H5Tset_size(mem_type.get(), H5T_VARIABLE);
        if (H5Dread(dataset.get(), mem_type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    raw.data()) < 0) {
            throw Error("cannot read string dataset " + dataset_path);
        }
        for (char* item : raw) {
            out.emplace_back(item != nullptr ? item : "");
        }
        const Handle space(H5Dget_space(dataset.get()), Handle::Kind::DataSpace);
        H5Dvlen_reclaim(mem_type.get(), space.get(), H5P_DEFAULT, raw.data());
        return out;
    }
    const size_t item_size = H5Tget_size(type.get());
    std::vector<char> buffer(item_size * length, '\0');
    if (H5Dread(dataset.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                buffer.data()) < 0) {
        throw Error("cannot read string dataset " + dataset_path);
    }
    for (std::size_t i = 0; i < length; ++i) {
        const char* start = buffer.data() + i * item_size;
        // numpy fixed width byte strings are NUL padded.
        const std::size_t used = ::strnlen(start, item_size);
        out.emplace_back(start, used);
    }
    return out;
}

}  // namespace hicx::h5
