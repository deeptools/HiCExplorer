// A minimal reader and writer for NumPy .npy files (format 1.0, C order,
// little endian), used by the stage dumps of hicDetectStripes' oracle
// comparison (cpp/scripts/stripenn_oracle) and by the unit tests that read
// fixtures produced by the Stripenn oracle. Not a general npy implementation.

#ifndef HICX_TOOLS_STRIPES_NPY_HPP
#define HICX_TOOLS_STRIPES_NPY_HPP

#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hicx::stripes::npy {

struct Array {
    std::string descr;  // e.g. "<f8"
    std::vector<std::size_t> shape;
    std::vector<char> bytes;

    [[nodiscard]] std::size_t count() const {
        std::size_t n = 1;
        for (const std::size_t s : shape) {
            n *= s;
        }
        return n;
    }
};

template <class T>
inline const char* descr_of();
template <>
inline const char* descr_of<double>() { return "<f8"; }
template <>
inline const char* descr_of<float>() { return "<f4"; }
template <>
inline const char* descr_of<std::int64_t>() { return "<i8"; }
template <>
inline const char* descr_of<std::int32_t>() { return "<i4"; }
template <>
inline const char* descr_of<std::uint8_t>() { return "|u1"; }

template <class T>
inline void save(const std::string& path, const T* data, const std::vector<std::size_t>& shape) {
    std::string header = std::string("{'descr': '") + descr_of<T>() + "', 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        header += std::to_string(shape[i]);
        header += shape.size() == 1 ? "," : (i + 1 < shape.size() ? ", " : "");
    }
    header += "), }";
    while ((10 + header.size() + 1) % 64 != 0) {
        header += ' ';
    }
    header += '\n';
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    const char magic[] = "\x93NUMPY\x01\x00";
    out.write(magic, 8);
    const auto length = static_cast<std::uint16_t>(header.size());
    out.write(reinterpret_cast<const char*>(&length), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    std::size_t n = 1;
    for (const std::size_t s : shape) {
        n *= s;
    }
    out.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(n * sizeof(T)));
}

template <class T>
inline void save(const std::string& path, const std::vector<T>& data, const std::vector<std::size_t>& shape) {
    save<T>(path, data.data(), shape);
}

inline Array load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot read " + path);
    }
    char magic[8];
    in.read(magic, 8);
    if (!in || std::memcmp(magic, "\x93NUMPY", 6) != 0) {
        throw std::runtime_error("not an npy file: " + path);
    }
    std::uint32_t header_length = 0;
    if (magic[6] == 1) {
        std::uint16_t h16 = 0;
        in.read(reinterpret_cast<char*>(&h16), 2);
        header_length = h16;
    } else {
        in.read(reinterpret_cast<char*>(&header_length), 4);
    }
    std::string header(header_length, '\0');
    in.read(header.data(), header_length);
    Array array;
    const std::size_t d = header.find("'descr'");
    const std::size_t q1 = header.find('\'', header.find(':', d) + 1);
    const std::size_t q2 = header.find('\'', q1 + 1);
    array.descr = header.substr(q1 + 1, q2 - q1 - 1);
    if (header.find("'fortran_order': True") != std::string::npos) {
        throw std::runtime_error("fortran order not supported: " + path);
    }
    const std::size_t s = header.find("'shape'");
    const std::size_t open = header.find('(', s);
    const std::size_t close = header.find(')', open);
    std::string dims = header.substr(open + 1, close - open - 1);
    std::stringstream ss(dims);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.find_first_of("0123456789") != std::string::npos) {
            array.shape.push_back(static_cast<std::size_t>(std::stoull(token)));
        }
    }
    std::size_t item = 8;
    if (array.descr == "<f4" || array.descr == "<i4") {
        item = 4;
    } else if (array.descr == "|u1" || array.descr == "|b1") {
        item = 1;
    }
    array.bytes.resize(array.count() * item);
    in.read(array.bytes.data(), static_cast<std::streamsize>(array.bytes.size()));
    return array;
}

template <class T>
inline std::vector<T> as(const Array& array) {
    std::vector<T> out(array.count());
    if (sizeof(T) * out.size() != array.bytes.size()) {
        throw std::runtime_error("npy element size mismatch for " + array.descr);
    }
    std::memcpy(out.data(), array.bytes.data(), array.bytes.size());
    return out;
}

}  // namespace hicx::stripes::npy

#endif  // HICX_TOOLS_STRIPES_NPY_HPP
