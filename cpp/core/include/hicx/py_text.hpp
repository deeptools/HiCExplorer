// Python's str semantics on UTF-8 text, for the C++ bodies of the public
// Python API (cpp/python/API_BINDINGS.md).
//
// A Python str reaches C++ as UTF-8 encoded with the "surrogatepass" error
// handler, so a lone surrogate arrives as its three byte form and goes back
// unchanged. Every function here works on code points in that encoding:
//
//  * whitespace is Python's str.isspace() set, not the C locale's: the ASCII
//    controls 0x09-0x0D and 0x1C-0x1F, the space, U+0085, U+00A0, U+1680,
//    U+2000-U+200A, U+2028, U+2029, U+202F, U+205F and U+3000;
//  * strip(), split() and split(sep) follow str.strip(), str.split() and
//    str.split(sep) exactly, including the empty fields split(sep) keeps;
//  * text files are read as open(path, 'r') reads them in a UTF-8 locale:
//    strict UTF-8 (no surrogates, no overlong forms) and universal newlines,
//    so "\r\n" and a lone "\r" both end a line and become "\n".

#ifndef HICX_PY_TEXT_HPP
#define HICX_PY_TEXT_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace hicx::pytext {

// The byte length of the whitespace code point starting at text[pos], or 0
// when the code point there is not whitespace. `pos` must start a code point.
[[nodiscard]] std::size_t whitespace_length(std::string_view text, std::size_t pos) noexcept;

// The byte length of the code point that ends at text[end - 1], for walking
// backwards. `end` must be the end of a code point and greater than zero.
[[nodiscard]] std::size_t previous_code_point_length(std::string_view text,
                                                     std::size_t end) noexcept;

// The number of code points, len() of the str.
[[nodiscard]] std::size_t code_point_count(std::string_view text) noexcept;

// str.strip() with no argument.
[[nodiscard]] std::string_view strip(std::string_view text) noexcept;

// str.split() with no argument: runs of whitespace separate, leading and
// trailing whitespace yields no empty field.
[[nodiscard]] std::vector<std::string_view> split_whitespace(std::string_view text);

// str.split(sep) for a non-empty separator: every occurrence separates, so
// n occurrences give n + 1 fields, empty ones included.
[[nodiscard]] std::vector<std::string_view> split(std::string_view text, std::string_view sep);

// ''.join(text.split()): the text with every whitespace code point removed.
[[nodiscard]] std::string remove_whitespace(std::string_view text);

// The offset of the first byte that strict UTF-8 decoding (Python's 'utf-8'
// codec with errors='strict') rejects, or nullopt when the whole text decodes.
// A sequence cut off by the end of the text is rejected at its first byte.
[[nodiscard]] std::optional<std::size_t> first_invalid_utf8(std::string_view bytes) noexcept;

// The lines file.readlines() returns for a text file opened with newline=None:
// "\r\n" and "\r" are translated to "\n", every line but possibly the last
// ends with "\n", and an empty file has no line.
[[nodiscard]] std::vector<std::string> universal_newline_lines(std::string_view text);

}  // namespace hicx::pytext

#endif  // HICX_PY_TEXT_HPP
