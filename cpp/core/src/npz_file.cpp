#include "hicx/npz_file.hpp"

#include <zlib.h>

#include <cstring>
#include <fstream>
#include <stdexcept>

namespace hicx::npz {

namespace {

void put_u16(std::string& out, std::uint16_t value) {
    out.push_back(static_cast<char>(value & 0xFF));
    out.push_back(static_cast<char>((value >> 8) & 0xFF));
}

void put_u32(std::string& out, std::uint32_t value) {
    for (int i = 0; i < 4; ++i) {
        out.push_back(static_cast<char>((value >> (8 * i)) & 0xFF));
    }
}

std::uint16_t get_u16(const std::string& text, std::size_t offset) {
    return static_cast<std::uint16_t>(
        static_cast<unsigned char>(text[offset]) |
        (static_cast<unsigned char>(text[offset + 1]) << 8));
}

std::uint32_t get_u32(const std::string& text, std::size_t offset) {
    std::uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        value |= static_cast<std::uint32_t>(
                     static_cast<unsigned char>(text[offset + static_cast<std::size_t>(i)]))
                 << (8 * i);
    }
    return value;
}

std::string deflate_raw(const std::string& input) {
    // Python's zipfile compresses with zlib.compressobj(-1, DEFLATED, -15),
    // that is a raw deflate stream at the default level.
    z_stream stream{};
    if (deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK) {
        throw std::runtime_error("deflateInit2 failed");
    }
    std::string out;
    out.resize(deflateBound(&stream, static_cast<uLong>(input.size())) + 64);
    stream.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    stream.avail_in = static_cast<uInt>(input.size());
    stream.next_out = reinterpret_cast<Bytef*>(out.data());
    stream.avail_out = static_cast<uInt>(out.size());
    const int status = deflate(&stream, Z_FINISH);
    if (status != Z_STREAM_END) {
        deflateEnd(&stream);
        throw std::runtime_error("deflate failed");
    }
    out.resize(stream.total_out);
    deflateEnd(&stream);
    return out;
}

std::string inflate_raw(const std::string& input, std::size_t expected_size) {
    z_stream stream{};
    if (inflateInit2(&stream, -15) != Z_OK) {
        throw std::runtime_error("inflateInit2 failed");
    }
    std::string out;
    out.resize(expected_size);
    stream.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    stream.avail_in = static_cast<uInt>(input.size());
    stream.next_out = reinterpret_cast<Bytef*>(out.data());
    stream.avail_out = static_cast<uInt>(out.size());
    const int status = inflate(&stream, Z_FINISH);
    inflateEnd(&stream);
    if (status != Z_STREAM_END && status != Z_OK) {
        throw std::runtime_error("inflate failed");
    }
    out.resize(stream.total_out);
    return out;
}

std::uint32_t crc32_of(const std::string& text) {
    return static_cast<std::uint32_t>(
        crc32(0L, reinterpret_cast<const Bytef*>(text.data()),
              static_cast<uInt>(text.size())));
}

}  // namespace

std::string npy_header(const std::string& dtype,
                       const std::vector<std::int64_t>& shape) {
    // The keys are written in sorted order, which for these three is descr,
    // fortran_order, shape, and every entry ends in ", " including the last.
    std::string body = "{'descr': '" + dtype + "', 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        body += std::to_string(shape[i]);
        // A one element tuple keeps its trailing comma, a longer one separates
        // with ", ", which is what repr(tuple) produces.
        if (shape.size() == 1) {
            body += ",";
        } else if (i + 1 != shape.size()) {
            body += ", ";
        }
    }
    body += "), }";

    // MAGIC_LEN is 8 (six magic bytes plus two version bytes) and the length
    // field is two more, so the padding makes 10 + len(body) + padding + 1 a
    // multiple of 64.
    std::size_t total = 10 + body.size() + 1;
    const std::size_t pad = 64 - (total % 64);
    body.append(pad == 64 ? 0 : pad, ' ');
    body.push_back('\n');

    std::string header;
    header += '\x93';
    header += "NUMPY";
    header.push_back('\x01');
    header.push_back('\x00');
    put_u16(header, static_cast<std::uint16_t>(body.size()));
    header += body;
    return header;
}

void write_npz(const std::string& path, const std::vector<Array>& arrays,
               bool compressed) {
    std::string file;
    std::string central;
    std::uint32_t entries = 0;

    for (const Array& array : arrays) {
        const std::string name = array.name + ".npy";
        const std::string payload = npy_header(array.dtype, array.shape) + array.data;
        const std::string stored = compressed ? deflate_raw(payload) : payload;
        const std::uint32_t crc = crc32_of(payload);
        const std::uint32_t offset = static_cast<std::uint32_t>(file.size());

        // Local file header.
        put_u32(file, 0x04034b50);
        put_u16(file, 20);                              // version needed
        put_u16(file, 0);                               // flags
        put_u16(file, compressed ? 8 : 0);              // method
        put_u16(file, 0);                               // time, the ZIP epoch
        put_u16(file, 0x0021);                          // date, 1980-01-01
        put_u32(file, crc);
        put_u32(file, static_cast<std::uint32_t>(stored.size()));
        put_u32(file, static_cast<std::uint32_t>(payload.size()));
        put_u16(file, static_cast<std::uint16_t>(name.size()));
        put_u16(file, 0);                               // extra length
        file += name;
        file += stored;

        // Central directory entry.
        put_u32(central, 0x02014b50);
        put_u16(central, 20 | (3 << 8));                // version made by, unix
        put_u16(central, 20);
        put_u16(central, 0);
        put_u16(central, compressed ? 8 : 0);
        put_u16(central, 0);
        put_u16(central, 0x0021);
        put_u32(central, crc);
        put_u32(central, static_cast<std::uint32_t>(stored.size()));
        put_u32(central, static_cast<std::uint32_t>(payload.size()));
        put_u16(central, static_cast<std::uint16_t>(name.size()));
        put_u16(central, 0);                            // extra
        put_u16(central, 0);                            // comment
        put_u16(central, 0);                            // disk number
        put_u16(central, 0);                            // internal attributes
        put_u32(central, 0600U << 16);                  // external attributes
        put_u32(central, offset);
        central += name;
        ++entries;
    }

    const std::uint32_t central_offset = static_cast<std::uint32_t>(file.size());
    file += central;
    put_u32(file, 0x06054b50);
    put_u16(file, 0);
    put_u16(file, 0);
    put_u16(file, static_cast<std::uint16_t>(entries));
    put_u16(file, static_cast<std::uint16_t>(entries));
    put_u32(file, static_cast<std::uint32_t>(central.size()));
    put_u32(file, central_offset);
    put_u16(file, 0);

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("cannot write '" + path + "'");
    }
    out.write(file.data(), static_cast<std::streamsize>(file.size()));
    if (!out) {
        throw std::runtime_error("failed to write '" + path + "'");
    }
}

std::vector<Array> read_npz(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open '" + path + "'");
    }
    const std::string file((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());

    // Find the end of central directory record, scanning back over the
    // optional comment.
    if (file.size() < 22) {
        throw std::runtime_error("'" + path + "' is not a zip archive");
    }
    std::size_t eocd = std::string::npos;
    for (std::size_t back = 0; back + 22 <= file.size() && back < 65558; ++back) {
        const std::size_t position = file.size() - 22 - back;
        if (get_u32(file, position) == 0x06054b50) {
            eocd = position;
            break;
        }
    }
    if (eocd == std::string::npos) {
        throw std::runtime_error("'" + path + "' has no zip end of central directory");
    }
    const std::uint16_t count = get_u16(file, eocd + 10);
    std::size_t cursor = get_u32(file, eocd + 16);

    std::vector<Array> arrays;
    for (std::uint16_t entry = 0; entry < count; ++entry) {
        if (get_u32(file, cursor) != 0x02014b50) {
            throw std::runtime_error("corrupt zip central directory in '" + path + "'");
        }
        const std::uint16_t method = get_u16(file, cursor + 10);
        const std::uint32_t compressed_size = get_u32(file, cursor + 20);
        const std::uint32_t raw_size = get_u32(file, cursor + 24);
        const std::uint16_t name_length = get_u16(file, cursor + 28);
        const std::uint16_t extra_length = get_u16(file, cursor + 30);
        const std::uint16_t comment_length = get_u16(file, cursor + 32);
        const std::uint32_t local_offset = get_u32(file, cursor + 42);
        std::string name = file.substr(cursor + 46, name_length);
        cursor += 46 + name_length + extra_length + comment_length;

        const std::uint16_t local_name = get_u16(file, local_offset + 26);
        const std::uint16_t local_extra = get_u16(file, local_offset + 28);
        const std::size_t data_offset = local_offset + 30 + local_name + local_extra;
        const std::string stored = file.substr(data_offset, compressed_size);
        const std::string payload =
            method == 0 ? stored : inflate_raw(stored, raw_size);

        if (payload.size() < 10 || payload.compare(0, 6, "\x93NUMPY") != 0) {
            throw std::runtime_error("entry '" + name + "' is not a .npy array");
        }
        const unsigned char major = static_cast<unsigned char>(payload[6]);
        std::size_t header_length = 0;
        std::size_t header_start = 0;
        if (major == 1) {
            header_length = get_u16(payload, 8);
            header_start = 10;
        } else if (major == 2) {
            header_length = get_u32(payload, 8);
            header_start = 12;
        } else {
            throw std::runtime_error("unsupported .npy version in '" + name + "'");
        }
        const std::string header = payload.substr(header_start, header_length);

        Array array;
        if (name.size() > 4 && name.compare(name.size() - 4, 4, ".npy") == 0) {
            name.erase(name.size() - 4);
        }
        array.name = name;

        const std::size_t descr = header.find("'descr': '");
        if (descr == std::string::npos) {
            throw std::runtime_error("entry '" + name +
                                     "' has no plain dtype; a pickled object array is "
                                     "not supported");
        }
        const std::size_t descr_end = header.find('\'', descr + 10);
        array.dtype = header.substr(descr + 10, descr_end - descr - 10);
        if (header.find("'fortran_order': True") != std::string::npos) {
            throw std::runtime_error("entry '" + name + "' is Fortran ordered");
        }

        const std::size_t shape_open = header.find("'shape': (");
        const std::size_t shape_close = header.find(')', shape_open);
        const std::string shape_text =
            header.substr(shape_open + 10, shape_close - shape_open - 10);
        std::string number;
        for (char c : shape_text + ",") {
            if (c >= '0' && c <= '9') {
                number.push_back(c);
                continue;
            }
            if (!number.empty()) {
                array.shape.push_back(std::stoll(number));
                number.clear();
            }
        }
        array.data = payload.substr(header_start + header_length);
        arrays.push_back(std::move(array));
    }
    return arrays;
}

void save_csr_npz(const std::string& path, std::int64_t rows, std::int64_t cols,
                  const std::vector<std::int32_t>& indptr,
                  const std::vector<std::int32_t>& indices,
                  const std::vector<double>& data) {
    const auto bytes_of = [](const void* source, std::size_t length) {
        std::string out;
        out.resize(length);
        std::memcpy(out.data(), source, length);
        return out;
    };

    std::vector<Array> arrays;
    // scipy.sparse.save_npz fills its dictionary in this order and numpy's
    // savez writes the entries in dictionary order.
    arrays.push_back(Array{"indices", "<i4",
                           {static_cast<std::int64_t>(indices.size())},
                           bytes_of(indices.data(), indices.size() * 4)});
    arrays.push_back(Array{"indptr", "<i4",
                           {static_cast<std::int64_t>(indptr.size())},
                           bytes_of(indptr.data(), indptr.size() * 4)});
    arrays.push_back(Array{"format", "|S3", {}, "csr"});
    const std::int64_t shape[2] = {rows, cols};
    arrays.push_back(Array{"shape", "<i8", {2}, bytes_of(shape, 16)});
    arrays.push_back(Array{"data", "<f8",
                           {static_cast<std::int64_t>(data.size())},
                           bytes_of(data.data(), data.size() * 8)});
    write_npz(path, arrays, true);
}

}  // namespace hicx::npz
