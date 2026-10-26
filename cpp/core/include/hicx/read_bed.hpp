// The C++ implementation behind hicexplorer.readBed.ReadBed (HiCExplorer
// 3.7.6, upstream c2ad8630), which reads BED3, BED6, BED9 and BED12 lines
// into namedtuples. The Python class keeps its attributes and methods; each
// method body is one call into the extension (cpp/python/API_BINDINGS.md),
// which uses the rules here.
//
// Reproduced as the Python has them:
//  * the file type is guessed from the tab separated fields of the first line
//    that is not a comment, a track or browser line or blank, *without*
//    stripping the line, so a trailing newline stays in the last field;
//  * a bed12 file is never checked for 12 fields, because the Python compares
//    the file handle, not the file type, with 'bed12';
//  * every field of a line is converted by its column, also beyond the 12th:
//    columns 0 and 3 stay strings, 5 is the strand, 1, 2, 6, 7 and 9 are
//    integers (a failure makes the line an empty dict), 8 is the colour,
//    10 and 11 the block lists, everything else a float when it parses.
// Integer and float parsing is Python's int() and float(), which the
// extension takes from the interpreter.

#ifndef HICX_READ_BED_HPP
#define HICX_READ_BED_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace hicx::readbed {

// guess_file_type by the number of fields. `assumed` is set to "bed6" or
// "bed3" when the count is not a standard one and the type is assumed.
[[nodiscard]] std::string guess_file_type(std::size_t field_count, std::string* assumed = nullptr);

// The namedtuple fields of a file type: 12 for bed12, 9 for bed9, 6 otherwise.
[[nodiscard]] std::vector<std::string> interval_fields(std::string_view file_type);

// The field count get_bed_interval asserts for a file type, if any.
[[nodiscard]] std::optional<std::size_t> required_field_count(std::string_view file_type);

enum class Column { Text, Strand, Integer, Rgb, Blocks, Float };
[[nodiscard]] Column column_kind(std::size_t index) noexcept;

// The strand column: '+', '-' and '.' stay, '1' is '+', '-1' is '-', anything
// else is '.' with `invalid` set.
[[nodiscard]] std::string normalize_strand(std::string_view value, bool& invalid);

// Comment lines: '#', 'track', 'browser' or only whitespace.
[[nodiscard]] bool is_comment_line(std::string_view line) noexcept;

// line.strip().split('\t') for get_bed_interval; line.split('\t') for the
// type guess (no strip).
[[nodiscard]] std::vector<std::string> split_line(std::string_view line, bool strip);

// The parts of the rgb column (r.split(',')), and of a block column without
// the empty parts.
[[nodiscard]] std::vector<std::string> split_commas(std::string_view value, bool drop_empty);

// The messages of the Python's assertions and warnings; `line_number` is
// str(self.line_number).
[[nodiscard]] std::string field_count_message(std::string_view file_type, std::size_t count,
                                              std::string_view line_number, std::string_view line);
[[nodiscard]] std::string start_end_message(std::string_view line_number, std::string_view line);
[[nodiscard]] std::string unsorted_message(std::string_view file_name, std::string_view previous_line,
                                           std::string_view line);
[[nodiscard]] std::string invalid_strand_message(std::string_view value, std::string_view line,
                                                 std::string_view line_number);
[[nodiscard]] std::string not_integer_message(std::string_view value, std::size_t column,
                                              std::string_view line_number);

}  // namespace hicx::readbed

#endif  // HICX_READ_BED_HPP
