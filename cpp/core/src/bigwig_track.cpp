#include "hicx/bigwig_track.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <stdexcept>

extern "C" {
#include "bigWig.h"
}

namespace hicx::bigwig {

namespace {

void initialise() {
    static std::once_flag once;
    std::call_once(once, [] {
        if (bwInit(1 << 17) != 0) {
            throw std::runtime_error("libBigWig initialisation failed");
        }
    });
}

bigWigFile_t* handle(void* file) {
    return static_cast<bigWigFile_t*>(file);
}

}  // namespace

Track::Track(const std::string& path) {
    initialise();
    bigWigFile_t* file = bwOpen(const_cast<char*>(path.c_str()), nullptr, "r");
    if (file == nullptr) {
        throw std::runtime_error("Received an error during file opening!");
    }
    file_ = file;
    if (file->cl != nullptr) {
        for (std::int64_t k = 0; k < file->cl->nKeys; ++k) {
            chroms_.emplace_back(file->cl->chrom[k], static_cast<std::int64_t>(file->cl->len[k]));
        }
    }
}

Track::~Track() {
    if (file_ != nullptr) {
        bwClose(handle(file_));
    }
}

std::optional<std::int64_t> Track::length(const std::string& chrom) const {
    for (const auto& [name, size] : chroms_) {
        if (name == chrom) {
            return size;
        }
    }
    return std::nullopt;
}

std::vector<double> Track::mean(const std::string& chrom, std::int64_t start, std::int64_t end,
                                std::int64_t bins) const {
    const std::optional<std::int64_t> size = length(chrom);
    if (!size || start < 0 || end > *size || start >= end) {
        throw std::invalid_argument("Invalid interval bounds!");
    }
    if (bins < 1) {
        bins = 1;
    }
    double* values = bwStats(handle(file_), const_cast<char*>(chrom.c_str()), static_cast<std::uint32_t>(start),
                             static_cast<std::uint32_t>(end), static_cast<std::uint32_t>(bins), ::mean);
    if (values == nullptr) {
        throw std::runtime_error("An error was encountered while fetching statistics.");
    }
    std::vector<double> out(values, values + bins);
    std::free(values);
    return out;
}

}  // namespace hicx::bigwig
