// hic2cool_extractnorms (hic2cool_utils.py of hic2cool 1.0.1), with
// cooler.fileops.list_coolers, cooler.binnify and cooler.create.append
// reduced to what it uses of them.

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include <hicfilecpp/hicfilecpp.hpp>

#include <coolercpp/coolercpp.hpp>

#include "io.hpp"
#include "hic2cool/hic2cool.hpp"
#include "hic_header.hpp"
#include "pyformat.hpp"

namespace hic2cool {

void hic2cool_extractnorms(const std::string& infile, const std::string& outfile, bool exclude_mt,
                           bool show_warnings, bool silent, const Console& console) {
    using namespace detail;
    bool warn = false;
    check_hic_magic(infile);
    const hicfilecpp::HiCFile hic(infile);
    HicHeader header = read_hic_header(hic, console);

    std::vector<std::string> chr_names;
    for (const auto& chrom : header.used) {
        chr_names.push_back(chrom.name);
    }
    if (!silent) {
        console.out("################################");
        console.out("### hic2cool / extract-norms ###");
        console.out("################################");
        console.out("Header info from hic:");
        console.out("... Chromosomes:  " + py_repr(chr_names));
        console.out("... Resolutions:  " + py_repr(header.resolutions));
        console.out("... Normalizations:  " + py_repr(header.norms));
        console.out("... Genome:  " + header.genome);
    }

    if (exclude_mt) {
        std::vector<std::size_t> found;
        for (std::size_t i = 0; i < header.used.size(); ++i) {
            const std::string name = lowercase(header.used[i].name);
            if (name == "m" || name == "mt" || name == "chrm" || name == "chrmt") {
                found.push_back(i);
            }
        }
        if (found.size() == 1) {
            const HicChrom excluded = header.used[found[0]];
            header.used.erase(header.used.begin() + static_cast<std::ptrdiff_t>(found[0]));
            if (!silent) {
                console.out("... Excluding chromosome " + excluded.name + " with index " +
                            std::to_string(excluded.index));
            }
        }
        // As in hic2cool: the else belongs to the second test, so a single
        // match prints both messages.
        if (found.size() > 1) {
            throw ExitError("ERROR. More than one chromosome was found when attempting to exclude MT. Found "
                            "chromosomes: " +
                            py_repr(chr_names));
        } else if (!silent) {
            console.out("... No chromosome found when attempting to exclude MT.");
        }
    }

    std::vector<const HicChrom*> chromosomes;
    for (const auto& chrom : header.used) {
        if (lowercase(chrom.name) != "all") {
            chromosomes.push_back(&chrom);
        }
    }

    ThreadPool pool(1);
    const std::vector<std::string> cooler_paths = coolercpp::list_coolers(outfile);
    h5::File file(outfile, h5::Mode::ReadWrite);
    std::vector<std::pair<std::int64_t, std::string>> cooler_groups;
    for (const auto& path : cooler_paths) {
        const auto size = file.attribute(path, "bin-size");
        if (!size || !std::holds_alternative<std::int64_t>(*size)) {
            throw ExitError("!!! ERROR. Cooler " + path + " in " + outfile + " has no integer bin-size");
        }
        const std::int64_t binsize = std::get<std::int64_t>(*size);
        const auto it = std::find_if(cooler_groups.begin(), cooler_groups.end(),
                                     [&](const auto& entry) { return entry.first == binsize; });
        if (it == cooler_groups.end()) {
            cooler_groups.emplace_back(binsize, path);
        } else {
            it->second = path;
        }
    }
    if (!silent) {
        console.out("### Found cooler contents:");
        console.out("... " + py_repr(cooler_groups));
    }

    for (const auto& norm : header.norms) {
        for (const std::int32_t binsize : header.resolutions) {
            const auto group = std::find_if(cooler_groups.begin(), cooler_groups.end(),
                                            [&](const auto& entry) { return entry.first == binsize; });
            if (group == cooler_groups.end()) {
                if (!silent) {
                    console.out("... Skip resolution " + std::to_string(binsize) + "; it is not in cooler file");
                }
                continue;
            }
            if (!silent) {
                console.out("... Extracting " + norm + " normalization vector at " + std::to_string(binsize) + " BP");
            }
            std::vector<double> column;
            std::vector<BinsPreviewRow> preview;
            std::size_t n_bins = 0;
            for (const HicChrom* chrom : chromosomes) {
                const auto chr_bins = static_cast<std::size_t>((chrom->length + binsize - 1) / binsize);
                for (std::size_t b = 0; b < chr_bins && preview.size() < 5; ++b) {
                    preview.push_back(BinsPreviewRow{chrom->name, static_cast<std::int64_t>(b) * binsize,
                                                     std::min<std::int64_t>((static_cast<std::int64_t>(b) + 1) * binsize,
                                                                            chrom->length),
                                                     0.0});
                }
                n_bins += chr_bins;
                const auto vector = hic.readNormVector(norm, chrom->index, "BP", binsize);
                if (!vector) {
                    warn = true;
                    if (show_warnings && !silent) {
                        console.err("!!! WARNING. Normalization vector " + norm + " does not exist for " +
                                    chrom->name + ".");
                    }
                    column.insert(column.end(), chr_bins, std::numeric_limits<double>::quiet_NaN());
                } else {
                    column.insert(column.end(), vector->begin(),
                                  vector->begin() + static_cast<std::ptrdiff_t>(std::min(vector->size(), chr_bins)));
                }
            }
            for (std::size_t i = 0; i < preview.size() && i < column.size(); ++i) {
                preview[i].value = column[i];
            }
            if (column.size() != n_bins) {
                throw ExitError("ValueError: Length of values (" + std::to_string(column.size()) +
                                ") does not match length of index (" + std::to_string(n_bins) + ")");
            }
            if (!silent) {
                console.out("... Writing to cool file ...");
                console.out(pandas_bins_head(norm, preview) + "\n... Truncated ...");
            }
            // cooler.create.append(uri, 'bins', {norm: values}, force=True)
            const std::string path = (group->second == "/" ? std::string("/") : group->second + "/") + "bins/" + norm;
            if (file.exists(path)) {
                file.unlink(path);
            }
            h5::ChunkedWriter writer(file.create_dataset(path, h5::ColumnType::float64(), column.size(), true), 8, true);
            writer.append(column.data(), column.size());
            h5::flush(pool, {&writer}, true);
        }
    }
    file.close();
    if (!silent) {
        if (warn && !show_warnings) {
            console.out("... Warnings were found in this run. Run with -v to display them.");
        }
        console.out("### Finished! Output written to: " + outfile);
    }
}

}  // namespace hic2cool
