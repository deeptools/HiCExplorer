#ifndef HICFILECPP_ERRORS_HPP
#define HICFILECPP_ERRORS_HPP

#include <stdexcept>
#include <string>

namespace hicfilecpp {

// Every failure hicfilecpp reports: an unreadable or malformed file, an
// unsupported version, an unknown chromosome, a missing zoom level or vector,
// or invalid writer input.
class HicError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

}  // namespace hicfilecpp

#endif  // HICFILECPP_ERRORS_HPP
