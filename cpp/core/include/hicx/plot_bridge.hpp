// The bridge from a C++ plotting tool to its drawing module (cpp/PLAN.md tier
// 7, option (a), settled by the GUI request of tier 10).
//
// A plotting tool computes the data of its figure in C++ and writes the
// tool's data outputs itself. The figure is drawn by
// plot/hicexplorer_plot/<tool>.py with the matplotlib calls of the Python
// tool, from a JSON document the C++ side builds with the helpers below.
//
// draw() either writes that document to the file of the C++-only option
// --plotData and returns without drawing (the explicit data mode), or writes
// it to a temporary file and replaces the process with
//
//     $HICX_PLOT_PYTHON -m hicexplorer_plot <tool> <file> --remove-data
//
// (python3 when the variable is unset). Replacing the process means the C++
// working set is released before matplotlib starts, so the peak RSS of the
// pair is the larger of the two, not their sum. When a directory `plot` with
// the package exists next to the directory of the binary (the build tree
// copies it there), it is put in front of PYTHONPATH.
//
// Numbers in the document are written as CPython's float repr, which
// json.loads reads back to the identical double; NaN and the infinities as
// NaN, Infinity and -Infinity, which json.loads accepts.

#ifndef HICX_PLOT_BRIDGE_HPP
#define HICX_PLOT_BRIDGE_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hicx::plot {

[[nodiscard]] std::string json_string(const std::string& text);
[[nodiscard]] std::string json_number(double value);
[[nodiscard]] std::string json_int(std::int64_t value);
[[nodiscard]] std::string json_bool(bool value);
[[nodiscard]] std::string json_numbers(const std::vector<double>& values);
[[nodiscard]] std::string json_ints(const std::vector<std::int64_t>& values);
[[nodiscard]] std::string json_strings(const std::vector<std::string>& values);
// "[a, b, ...]" from JSON texts.
[[nodiscard]] std::string json_list(const std::vector<std::string>& items);

// An object whose values are JSON texts, in insertion order.
class JsonObject {
  public:
    JsonObject& add(const std::string& key, std::string json_text);
    [[nodiscard]] std::string str() const;

  private:
    std::vector<std::pair<std::string, std::string>> fields_;
};

// A new, empty file for data that is too large for the JSON document, in
// $TMPDIR (else /tmp). The document lists it under "temporary_files" so the
// drawing process removes it.
[[nodiscard]] std::string temporary_file();

// A float64, C order .npy file of rows by cols values, which numpy.load
// reads as a two dimensional array.
void write_npy_float64(const std::string& path, const std::vector<double>& values,
                       std::int64_t rows, std::int64_t cols);

// Hands the figure to the drawing layer, see the file comment. Returns the
// exit status the tool should return: 0 after writing the --plotData file, 1
// when the data could not be written or the interpreter could not be started.
// On success without data_file it does not return.
[[nodiscard]] int draw(const std::string& tool, const std::string& data_json,
                       const std::optional<std::string>& data_file);

}  // namespace hicx::plot

#endif  // HICX_PLOT_BRIDGE_HPP
