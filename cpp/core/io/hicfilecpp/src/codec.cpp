#include "codec.hpp"

#include <zlib.h>

#include <algorithm>
#include <string>

#include "hicfilecpp/errors.hpp"

namespace hicfilecpp::detail {

namespace {

// One inflate and one deflate state per thread, reset between streams.
struct InflateState {
    z_stream stream{};
    bool ready = false;
    InflateState() {
        if (inflateInit(&stream) != Z_OK) {
            throw HicError("zlib inflateInit failed");
        }
        ready = true;
    }
    ~InflateState() {
        if (ready) {
            inflateEnd(&stream);
        }
    }
};

}  // namespace

void inflateZlib(const char* data, size_t size, std::vector<char>& out) {
    thread_local InflateState state;
    z_stream& zs = state.stream;
    if (inflateReset(&zs) != Z_OK) {
        throw HicError("zlib inflateReset failed");
    }
    out.resize(std::max<size_t>(size * 4, 1024));
    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(data));
    zs.avail_in = static_cast<uInt>(size);
    size_t produced = 0;
    while (true) {
        zs.next_out = reinterpret_cast<Bytef*>(out.data() + produced);
        zs.avail_out = static_cast<uInt>(out.size() - produced);
        const int status = inflate(&zs, Z_NO_FLUSH);
        produced = out.size() - zs.avail_out;
        if (status == Z_STREAM_END) {
            break;
        }
        if (status == Z_OK) {
            if (zs.avail_out == 0) {
                out.resize(out.size() * 2);
            }
            continue;
        }
        if (status == Z_BUF_ERROR) {
            if (zs.avail_out == 0) {
                out.resize(out.size() * 2);
                continue;
            }
            break;  // input ended before the stream did
        }
        throw HicError("zlib inflate failed with status " + std::to_string(status));
    }
    out.resize(produced);
}

void deflateZlib(const char* data, size_t size, int level, std::vector<char>& out) {
    uLongf bound = compressBound(static_cast<uLong>(size));
    out.resize(bound);
    const int status = compress2(reinterpret_cast<Bytef*>(out.data()), &bound,
                                 reinterpret_cast<const Bytef*>(data), static_cast<uLong>(size), level);
    if (status != Z_OK) {
        throw HicError("zlib compress2 failed with status " + std::to_string(status));
    }
    out.resize(bound);
}

}  // namespace hicfilecpp::detail
