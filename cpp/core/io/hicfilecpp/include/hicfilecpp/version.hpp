#ifndef HICFILECPP_VERSION_HPP
#define HICFILECPP_VERSION_HPP

// Vendored from version.hpp.in: the configure-time project version is written out.
namespace hicfilecpp {

inline constexpr const char* kVersion = "0.4.0";
// The hicstraw release whose reading API this library reproduces.
inline constexpr const char* kHicstrawVersion = "1.3.1";
// The .hic versions this library reads and writes.
inline constexpr int kMinVersion = 6;
inline constexpr int kMaxVersion = 9;

}  // namespace hicfilecpp

#endif  // HICFILECPP_VERSION_HPP
