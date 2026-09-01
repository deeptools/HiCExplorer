#include "hicx/hdf5_util.hpp"

#include <cstdlib>
#include <cstring>
#include <mutex>

#include <blosc.h>

namespace hicx::h5 {

namespace {

constexpr H5Z_filter_t kBloscFilterId = 32001;

// Decompression half of the reference hdf5-blosc filter. Compression is not
// implemented; version 4 does not write PyTables files yet.
size_t blosc_filter_impl(unsigned int flags, size_t cd_nelmts,
                         const unsigned int cd_values[], size_t nbytes,
                         size_t* buf_size, void** buf) {
    (void)cd_nelmts;
    (void)cd_values;
    (void)nbytes;
    if ((flags & H5Z_FLAG_REVERSE) == 0) {
        return 0;  // write path unsupported
    }
    size_t outbuf_size = 0;
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

const H5Z_class2_t kBloscClass = {
    H5Z_CLASS_T_VERS,
    kBloscFilterId,
    1,  // encoder_present
    1,  // decoder_present
    "blosc",
    nullptr,
    nullptr,
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
