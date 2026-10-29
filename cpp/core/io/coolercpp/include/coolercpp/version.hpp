#ifndef COOLERCPP_VERSION_HPP
#define COOLERCPP_VERSION_HPP

// Vendored from the coolercpp 0.4.0 CMake template (version.hpp.in): the
// configure-time project version is written out.
namespace coolercpp {

inline constexpr const char* kVersion = "0.4.0";
// The cooler release whose API and file format this library reproduces.
inline constexpr const char* kCoolerVersion = "0.10.2";
// cooler.create's format constants.
inline constexpr int kFormatVersion = 3;
inline constexpr const char* kMagic = "HDF5::Cooler";
inline constexpr const char* kMagicScool = "HDF5::SCOOL";
inline constexpr const char* kMagicMcool = "HDF5::MCOOL";
inline constexpr const char* kFormatUrl = "https://github.com/open2c/cooler";

}  // namespace coolercpp

#endif  // COOLERCPP_VERSION_HPP
