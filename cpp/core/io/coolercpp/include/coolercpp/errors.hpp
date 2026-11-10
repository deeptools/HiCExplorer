// Exception types that mirror the Python exceptions cooler raises.
//
// A C++ caller catches coolercpp::Error (or std::exception); code that ports a
// Python `except ValueError:` catches coolercpp::ValueError. python_type()
// names the Python class so that the equivalence harness can compare error
// behaviour case by case.

#ifndef COOLERCPP_ERRORS_HPP
#define COOLERCPP_ERRORS_HPP

#include <functional>
#include <stdexcept>
#include <string>

namespace coolercpp {

class Error : public std::runtime_error {
  public:
    explicit Error(const std::string& what) : std::runtime_error(what) {}
    [[nodiscard]] virtual const char* python_type() const noexcept { return "Exception"; }
};

#define COOLERCPP_DEFINE_ERROR(Name, Base)                                         \
    class Name : public Base {                                                     \
      public:                                                                      \
        explicit Name(const std::string& what) : Base(what) {}                     \
        [[nodiscard]] const char* python_type() const noexcept override {          \
            return #Name;                                                          \
        }                                                                          \
    };

COOLERCPP_DEFINE_ERROR(KeyError, Error)
COOLERCPP_DEFINE_ERROR(ValueError, Error)
COOLERCPP_DEFINE_ERROR(IndexError, Error)
COOLERCPP_DEFINE_ERROR(TypeError, Error)
COOLERCPP_DEFINE_ERROR(AttributeError, Error)
COOLERCPP_DEFINE_ERROR(NotImplementedError, Error)
// h5py raises OSError (or its subclass FileNotFoundError) for files that
// cannot be opened; HDF5 library failures during reads and writes map here.
COOLERCPP_DEFINE_ERROR(OSError, Error)
// cooler.create.BadInputError, a ValueError subclass raised by the pixel
// validation pipeline.
COOLERCPP_DEFINE_ERROR(BadInputError, ValueError)

#undef COOLERCPP_DEFINE_ERROR

// Python's warnings.warn. The default handler prints
// "UserWarning: <message>" to stderr; a caller may redirect or silence it.
using WarningHandler = std::function<void(const std::string& message)>;
void set_warning_handler(WarningHandler handler);
void warn(const std::string& message);

}  // namespace coolercpp

#endif  // COOLERCPP_ERRORS_HPP
