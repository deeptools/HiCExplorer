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
#include <span>
#include <string>

#include "hicx/chic_viewpoint.hpp"
#include "hicx/hdf5_util.hpp"

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

    void write_string(const std::string& path, const std::string& value);
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

}  // namespace hicx::chic

#endif  // HICX_CHIC_HDF5_HPP
