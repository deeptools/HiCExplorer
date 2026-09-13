// See hicx/chic_hdf5.hpp.

#include "hicx/chic_hdf5.hpp"

#include <hdf5.h>

#include <vector>

namespace hicx::chic {

namespace {

using h5::Error;
using h5::Handle;

Handle variable_string_type() {
    Handle type(H5Tcopy(H5T_C_S1), Handle::Kind::DataType);
    if (!type.valid() || H5Tset_size(type.get(), H5T_VARIABLE) < 0 ||
        H5Tset_cset(type.get(), H5T_CSET_UTF8) < 0 ||
        H5Tset_strpad(type.get(), H5T_STR_NULLTERM) < 0) {
        throw Error("cannot create a variable length string type");
    }
    return type;
}

Handle link_create_plist() {
    Handle plist(H5Pcreate(H5P_LINK_CREATE), Handle::Kind::PropertyList);
    if (!plist.valid() || H5Pset_create_intermediate_group(plist.get(), 1) < 0) {
        throw Error("cannot create a link creation property list");
    }
    return plist;
}

Handle dataset_create_plist() {
    Handle plist(H5Pcreate(H5P_DATASET_CREATE), Handle::Kind::PropertyList);
    if (!plist.valid() || H5Pset_obj_track_times(plist.get(), 0) < 0) {
        throw Error("cannot create a dataset creation property list");
    }
    return plist;
}

Handle open_object(hid_t file, const std::string& path) {
    const hid_t id = H5Oopen(file, path.c_str(), H5P_DEFAULT);
    if (id < 0) {
        throw Error("cannot open " + path);
    }
    // H5Oopen hands back a group or a dataset; Handle closes either through
    // H5Oclose when told it is a group, which H5Oclose accepts for both.
    return Handle(id, H5Iget_type(id) == H5I_DATASET ? Handle::Kind::Dataset
                                                     : Handle::Kind::Group);
}

}  // namespace

Hdf5Writer::Hdf5Writer(const std::string& path) : path_(path) {
    const hid_t id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (id < 0) {
        throw Error("cannot create HDF5 file: " + path);
    }
    file_ = Handle(id, Handle::Kind::File);
}

bool Hdf5Writer::exists(const std::string& object_path) const {
    if (object_path.empty() || object_path == "/") {
        return true;
    }
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

void Hdf5Writer::create_group(const std::string& path) {
    if (exists(path)) {
        throw Error("Unable to create group (name already exists): " + path);
    }
    const Handle link_plist = link_create_plist();
    const Handle group_plist(H5Pcreate(H5P_GROUP_CREATE), Handle::Kind::PropertyList);
    if (!group_plist.valid() || H5Pset_obj_track_times(group_plist.get(), 0) < 0) {
        throw Error("cannot create a group creation property list");
    }
    const Handle group(H5Gcreate2(file_.get(), path.c_str(), link_plist.get(),
                                  group_plist.get(), H5P_DEFAULT),
                       Handle::Kind::Group);
    if (!group.valid()) {
        throw Error("cannot create group " + path + " in " + path_);
    }
}

void Hdf5Writer::set_attribute(const std::string& object_path, const std::string& name,
                               const std::string& value) {
    const Handle object = open_object(file_.get(), object_path);
    const Handle type = variable_string_type();
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    const Handle attribute(
        H5Acreate2(object.get(), name.c_str(), type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
        Handle::Kind::Attribute);
    const char* pointer = value.c_str();
    if (!attribute.valid() || H5Awrite(attribute.get(), type.get(), &pointer) < 0) {
        throw Error("cannot write attribute " + name + " on " + object_path);
    }
}

void Hdf5Writer::set_attribute(const std::string& object_path, const std::string& name,
                               std::int64_t value) {
    const Handle object = open_object(file_.get(), object_path);
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    const Handle attribute(H5Acreate2(object.get(), name.c_str(), H5T_STD_I64LE, space.get(),
                                      H5P_DEFAULT, H5P_DEFAULT),
                           Handle::Kind::Attribute);
    if (!attribute.valid() || H5Awrite(attribute.get(), H5T_NATIVE_INT64, &value) < 0) {
        throw Error("cannot write attribute " + name + " on " + object_path);
    }
}

void Hdf5Writer::set_attribute(const std::string& object_path, const std::string& name,
                               std::span<const std::int64_t> values) {
    const Handle object = open_object(file_.get(), object_path);
    const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    const Handle space(H5Screate_simple(1, dims, nullptr), Handle::Kind::DataSpace);
    const Handle attribute(H5Acreate2(object.get(), name.c_str(), H5T_STD_I64LE, space.get(),
                                      H5P_DEFAULT, H5P_DEFAULT),
                           Handle::Kind::Attribute);
    if (!attribute.valid() ||
        H5Awrite(attribute.get(), H5T_NATIVE_INT64, values.data()) < 0) {
        throw Error("cannot write attribute " + name + " on " + object_path);
    }
}

void Hdf5Writer::write_string(const std::string& path, const std::string& value) {
    const Handle type = variable_string_type();
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    const Handle link_plist = link_create_plist();
    const Handle plist = dataset_create_plist();
    const Handle dataset(H5Dcreate2(file_.get(), path.c_str(), type.get(), space.get(),
                                    link_plist.get(), plist.get(), H5P_DEFAULT),
                         Handle::Kind::Dataset);
    const char* pointer = value.c_str();
    if (!dataset.valid() ||
        H5Dwrite(dataset.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointer) < 0) {
        throw Error("cannot write dataset " + path + " in " + path_);
    }
}

void Hdf5Writer::write_scalar(const std::string& path, std::int64_t value) {
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    const Handle link_plist = link_create_plist();
    const Handle plist = dataset_create_plist();
    const Handle dataset(H5Dcreate2(file_.get(), path.c_str(), H5T_STD_I64LE, space.get(),
                                    link_plist.get(), plist.get(), H5P_DEFAULT),
                         Handle::Kind::Dataset);
    if (!dataset.valid() || H5Dwrite(dataset.get(), H5T_NATIVE_INT64, H5S_ALL, H5S_ALL,
                                     H5P_DEFAULT, &value) < 0) {
        throw Error("cannot write dataset " + path + " in " + path_);
    }
}

void Hdf5Writer::write_scalar(const std::string& path, double value) {
    const Handle space(H5Screate(H5S_SCALAR), Handle::Kind::DataSpace);
    const Handle link_plist = link_create_plist();
    const Handle plist = dataset_create_plist();
    const Handle dataset(H5Dcreate2(file_.get(), path.c_str(), H5T_IEEE_F64LE, space.get(),
                                    link_plist.get(), plist.get(), H5P_DEFAULT),
                         Handle::Kind::Dataset);
    if (!dataset.valid() || H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                     H5P_DEFAULT, &value) < 0) {
        throw Error("cannot write dataset " + path + " in " + path_);
    }
}

void Hdf5Writer::write_numeric_array(const std::string& path, hid_t file_type,
                                     hid_t memory_type, std::size_t length,
                                     std::size_t type_size, const void* data, int level) {
    const hsize_t dims[1] = {static_cast<hsize_t>(length)};
    const Handle space(H5Screate_simple(1, dims, nullptr), Handle::Kind::DataSpace);
    const Handle link_plist = link_create_plist();
    const Handle plist = dataset_create_plist();
    if (level > 0) {
        const hsize_t chunk[1] = {static_cast<hsize_t>(h5::guess_chunk(length, type_size))};
        if (H5Pset_chunk(plist.get(), 1, chunk) < 0 || H5Pset_deflate(plist.get(), level) < 0) {
            throw Error("cannot set the gzip filter of " + path);
        }
    }
    const Handle dataset(H5Dcreate2(file_.get(), path.c_str(), file_type, space.get(),
                                    link_plist.get(), plist.get(), H5P_DEFAULT),
                         Handle::Kind::Dataset);
    if (!dataset.valid()) {
        throw Error("cannot create dataset " + path + " in " + path_);
    }
    if (length > 0 &&
        H5Dwrite(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
        throw Error("cannot write dataset " + path + " in " + path_);
    }
}

void Hdf5Writer::write_array(const std::string& path, std::span<const std::int64_t> values,
                             int level) {
    write_numeric_array(path, H5T_STD_I64LE, H5T_NATIVE_INT64, values.size(), 8,
                        values.data(), level);
}

void Hdf5Writer::write_array(const std::string& path, std::span<const double> values,
                             int level) {
    write_numeric_array(path, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, values.size(), 8,
                        values.data(), level);
}

bool Hdf5Writer::hard_link(const std::string& target_path, const std::string& link_path) {
    if (exists(link_path)) {
        return false;
    }
    const Handle link_plist = link_create_plist();
    return H5Lcreate_hard(file_.get(), target_path.c_str(), file_.get(), link_path.c_str(),
                          link_plist.get(), H5P_DEFAULT) >= 0;
}

void write_interaction_datasets(Hdf5Writer& writer, const std::string& group_path,
                                const InteractionFileData& data,
                                std::int64_t reference_point_start,
                                std::int64_t reference_point_end) {
    constexpr int kGzip = 9;
    const std::string prefix = group_path + "/";
    // A list of ints that happens to be empty is np.array([]), float64.
    const auto write_int_list = [&](const std::string& name,
                                    const std::vector<std::int64_t>& values) {
        if (values.empty()) {
            writer.write_array(prefix + name, std::span<const double>(), kGzip);
        } else {
            writer.write_array(prefix + name, std::span<const std::int64_t>(values), kGzip);
        }
    };
    writer.write_string(prefix + "chromosome", data.chromosome);
    write_int_list("start_list", data.starts);
    write_int_list("end_list", data.ends);
    writer.write_string(prefix + "gene", data.gene);
    writer.write_scalar(prefix + "sum_of_interactions", data.sum_of_interactions);
    write_int_list("relative_position_list", data.relative_positions);
    writer.write_array(prefix + "interaction_data_list",
                       std::span<const double>(data.interaction_data), kGzip);
    writer.write_array(prefix + "pvalue", std::span<const double>(data.pvalues), kGzip);
    writer.write_array(prefix + "xfold", std::span<const double>(data.xfold), kGzip);
    writer.write_array(prefix + "raw", std::span<const double>(data.raw), kGzip);
    writer.write_scalar(prefix + "reference_point_start", reference_point_start);
    writer.write_scalar(prefix + "reference_point_end", reference_point_end);
}

}  // namespace hicx::chic
