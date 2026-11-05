// An error that the Python API must raise as a specific Python exception.
//
// The public Python functions of HiCExplorer 4 are thin wrappers over C++
// (cpp/python/API_BINDINGS.md). Where the 3.7.6 function raised, say, a
// NameError or an argparse.ArgumentTypeError, the C++ implementation throws a
// PythonError naming that type, and the extension's exception translator
// (cpp/python/src/api/support/errors.cpp) raises it in Python with the same
// message. The core stays free of Python headers: the type is only a name.
//
// A C++ executable that calls the same code sees an ordinary
// std::runtime_error whose what() is the message.

#ifndef HICX_PYTHON_ERROR_HPP
#define HICX_PYTHON_ERROR_HPP

#include <stdexcept>
#include <string>
#include <utility>

namespace hicx {

class PythonError : public std::runtime_error {
  public:
    // `type` is a builtin exception name ("ValueError") or a dotted path to
    // an exception class ("argparse.ArgumentTypeError").
    PythonError(std::string type, const std::string& message)
        : std::runtime_error(message), type_(std::move(type)) {}

    [[nodiscard]] const std::string& type() const noexcept { return type_; }

  private:
    std::string type_;
};

}  // namespace hicx

#endif  // HICX_PYTHON_ERROR_HPP
