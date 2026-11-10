#include "hic_header.hpp"

#include <algorithm>
#include <fstream>

namespace hic2cool::detail {

std::string lowercase(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

void check_hic_magic(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw ExitError("[Errno 2] No such file or directory: " + path);
    }
    char magic[4] = {0, 0, 0, 0};
    in.read(magic, 4);
    if (in.gcount() < 4 || std::string(magic, 3) != "HIC") {
        throw ExitError("... This does not appear to be a HiC file; magic string is incorrect");
    }
}

HicHeader read_hic_header(const hicfilecpp::HiCFile& hic, const Console& console) {
    HicHeader header;
    for (const auto& chromosome : hic.getChromosomes()) {
        if (!chromosome.name.empty() && chromosome.length != 0) {
            header.used.push_back(HicChrom{chromosome.index, chromosome.name, chromosome.length});
        }
    }
    header.resolutions = hic.getResolutions();
    header.genome = hic.getGenomeID();
    for (const auto& [key, value] : hic.attributes()) {
        const auto it = std::find_if(header.metadata.begin(), header.metadata.end(),
                                     [&](const auto& entry) { return entry.first == key; });
        if (it == header.metadata.end()) {
            header.metadata.emplace_back(key, value);
        } else {
            it->second = value;
        }
    }
    if (!hic.hasNormalizedExpectedSection()) {
        console.err("!!! WARNING. No normalization vectors found in the hic file.");
        return header;
    }
    for (const auto& key : hic.expectedValuesKeys()) {
        if (key.normalization != "NONE" &&
            std::find(header.norms.begin(), header.norms.end(), key.normalization) == header.norms.end()) {
            header.norms.push_back(key.normalization);
        }
    }
    return header;
}

}  // namespace hic2cool::detail
