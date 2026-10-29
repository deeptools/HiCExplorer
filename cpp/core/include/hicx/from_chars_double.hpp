// std::from_chars for double, also where the C++ library does not offer it.
//
// conda-forge's libc++ marks the floating-point overloads of std::from_chars
// as unavailable on macOS, so code calling them does not compile there. On
// macOS with libc++ this parses with strtod_l in the C locale instead and
// applies the rules of std::from_chars with chars_format::general: no leading
// whitespace, no '+' sign, no hexadecimal form; "inf", "infinity" and "nan"
// in any case are accepted. strtod and std::from_chars both round correctly,
// so the value is the same. Everywhere else std::from_chars is called.

#ifndef HICX_FROM_CHARS_DOUBLE_HPP
#define HICX_FROM_CHARS_DOUBLE_HPP

#include <charconv>
#include <system_error>

#if defined(__APPLE__) && defined(_LIBCPP_VERSION)
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <string>
#include <xlocale.h>
#define HICX_FROM_CHARS_DOUBLE_STRTOD 1
#endif

namespace hicx {

inline std::from_chars_result from_chars_double(const char* first, const char* last, double& value) {
#ifdef HICX_FROM_CHARS_DOUBLE_STRTOD
    const char* body = first;
    if (body != last && *body == '-') {
        ++body;
    }
    const bool rejected = body == last || *body == ' ' || *body == '\t' || *body == '\n' ||
                          *body == '\r' || *body == '\f' || *body == '\v' || *body == '+' || *body == '-' ||
                          (last - body >= 2 && body[0] == '0' && (body[1] == 'x' || body[1] == 'X'));
    if (rejected) {
        return {first, std::errc::invalid_argument};
    }
    static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", nullptr);
    char small[64];
    std::string large;
    const auto length = static_cast<std::size_t>(last - first);
    const char* text;
    if (length < sizeof(small)) {
        std::char_traits<char>::copy(small, first, length);
        small[length] = '\0';
        text = small;
    } else {
        large.assign(first, last);
        text = large.c_str();
    }
    char* end = nullptr;
    errno = 0;
    const double parsed = strtod_l(text, &end, c_locale);
    if (end == text) {
        return {first, std::errc::invalid_argument};
    }
    const char* ptr = first + (end - text);
    if (errno == ERANGE && (parsed == 0.0 || std::isinf(parsed))) {
        return {ptr, std::errc::result_out_of_range};
    }
    value = parsed;
    return {ptr, std::errc()};
#else
    return std::from_chars(first, last, value);
#endif
}

}  // namespace hicx

#endif
