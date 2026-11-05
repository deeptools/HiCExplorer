// hic2cool_update (hic2cool_utils.py and hic2cool_updates.py of hic2cool
// 1.0.1): upgrades cool and mcool files written by hic2cool before 0.7.1.

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <functional>
#include <limits>

#include "io.hpp"
#include "hic2cool/hic2cool.hpp"
#include "pyformat.hpp"

namespace hic2cool {

namespace {

using namespace detail;

struct Update {
    std::string title;
    std::string effect;
    std::string detail;
    std::function<void(const std::string&, const Console&)> run;
};

std::vector<std::string> resolution_groups(const h5::File& file) {
    return file.is_group("/resolutions") ? file.children("/resolutions") : std::vector<std::string>{};
}

void invert_weights(const std::string& writefile, const Console& console) {
    h5::File file(writefile, h5::Mode::ReadWrite);
    const auto invert = [&](const std::string& bins, const std::string& res) {
        std::vector<std::string> found;
        for (const auto& name : file.children(bins)) {
            if (name != "chrom" && name != "start" && name != "end") {
                found.push_back(name);
            }
        }
        for (const auto& name : found) {
            std::vector<double> values = file.read_doubles(bins + "/" + name);
            for (double& v : values) {
                v = v != 0.0 ? 1.0 / v : std::numeric_limits<double>::quiet_NaN();
            }
            file.write_doubles(bins + "/" + name, values);
        }
        if (!res.empty()) {
            console.out("... For resolution " + res + ", inverted following weights: " + py_repr(found));
        } else {
            console.out("... Inverted following weights: " + py_repr(found));
        }
    };
    if (file.is_group("/resolutions")) {
        for (const auto& res : resolution_groups(file)) {
            invert("/resolutions/" + res + "/bins", res);
        }
    } else {
        invert("/bins", "");
    }
}

void cooler_schema_v3(const std::string& writefile, const Console& console) {
    h5::File file(writefile, h5::Mode::ReadWrite);
    const auto add = [&](const std::string& group, const std::string& res) {
        file.set_attribute(group, "format-version", std::int64_t{3});
        file.set_attribute(group, "storage-mode", std::string("symmetric-upper"));
        if (!res.empty()) {
            console.out("... For resolution " + res + ", added format-version and storage-mode attributes");
        } else {
            console.out("... Added format-version and storage-mode attributes");
        }
    };
    if (file.is_group("/resolutions")) {
        for (const auto& res : resolution_groups(file)) {
            add("/resolutions/" + res, res);
        }
    } else {
        add("/", "");
    }
}

void mcool_schema_v2(const std::string& writefile, const Console& console) {
    h5::File file(writefile, h5::Mode::ReadWrite);
    if (file.is_group("/resolutions")) {
        file.set_attribute("/", "format", std::string("HDF5::MCOOL"));
        file.set_attribute("/", "format-version", std::int64_t{2});
        console.out("... Added format and format-version attributes for the mcool");
    } else {
        console.out("... Not a multi-res file, so will not add mcool schema attributes");
    }
}

std::vector<Update> prepare_updates(const std::vector<int>& v) {
    std::vector<Update> updates;
    if (v[0] == 0 && v[1] < 5) {
        updates.push_back(Update{
            "Invert weights", "Invert cooler weights so that they match original hic normalization values",
            "cooler uses multiplicative weights and hic uses divisive weights. Before version 0.5.0, hic2cool "
            "inverted normalization vectors for consistency with cooler behavior, but now that is no longer done "
            "for consistency with 4DN analysis pipelines.",
            invert_weights});
    }
    if (v[0] == 0 && v[1] < 6) {
        updates.push_back(Update{"Add cooler schema version", "Add a couple important cooler schema attributes",
                                 "Adds format-version and storage-mode attributes to hdf5 for compatibility with "
                                 "cooler schema v3.",
                                 cooler_schema_v3});
    }
    if (v[0] == 0 && ((v[1] == 7 && v[2] < 1) || v[1] < 7)) {
        updates.push_back(Update{"Add mcool schema attributes",
                                 "Adds missing schema attributes if this is a multi-resolution cooler",
                                 "Adds format and format-version attributes to the \"/\" hdf5 collection for mcool "
                                 "schema v2.",
                                 mcool_schema_v2});
    }
    return updates;
}

// datetime.strptime(text, '%Y-%m-%dT%H:%M:%S.%f') formatted as
// '%Y-%m-%d at %H:%M:%S'. isoformat() leaves out the fraction when it is
// zero, which strptime rejects; that form is accepted here.
std::string timestamp(const std::string& text) {
    std::tm tm{};
    const char* rest = strptime(text.c_str(), "%Y-%m-%dT%H:%M:%S", &tm);
    if (rest == nullptr || (*rest != '\0' && *rest != '.')) {
        throw ExitError("ValueError: time data '" + text + "' does not match format '%Y-%m-%dT%H:%M:%S.%f'");
    }
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d at %H:%M:%S", &tm);
    return buffer;
}

std::string text_attribute(const std::vector<std::pair<std::string, h5::Value>>& attrs, const std::string& key,
                           bool& present) {
    for (const auto& [name, value] : attrs) {
        if (name == key) {
            present = true;
            if (const auto* s = std::get_if<std::string>(&value)) {
                return *s;
            }
            if (const auto* i = std::get_if<std::int64_t>(&value)) {
                return std::to_string(*i);
            }
            return {};
        }
    }
    present = false;
    return {};
}

}  // namespace

void hic2cool_update(const std::string& infile, const std::string& outfile, bool /*show_warnings*/, bool silent,
                     const Console& console) {
    if (!silent) {
        console.out("#########################");
        console.out("### hic2cool / update ###");
        console.out("#########################");
    }
    std::vector<std::pair<std::string, h5::Value>> attrs;
    std::string resolutions_text;
    {
        h5::File file(infile, h5::Mode::ReadWrite);
        if (file.is_group("/resolutions")) {
            const auto groups = resolution_groups(file);
            // hic2cool compares 'generated_by', which no file has, so
            // resolutions never disagree; the first one's attributes are used.
            if (!groups.empty()) {
                attrs = file.attributes("/resolutions/" + groups.front());
            }
            resolutions_text = py_repr(groups);
        } else {
            attrs = file.attributes("/");
            bool present = false;
            const std::string size = text_attribute(attrs, "bin-size", present);
            if (!present) {
                throw ExitError("KeyError: 'bin-size'");
            }
            resolutions_text = "[" + size + "]";
        }
    }
    bool present = false;
    const std::string generated_by = text_attribute(attrs, "generated-by", present);
    if (!present) {
        throw ExitError("!!! ERROR. Input file doesn't seem to made by hic2cool. Exiting.");
    }
    if (generated_by.rfind("hic2cool-", 0) != 0) {
        throw ExitError("!!! ERROR. Malformed 'generated-by' attribute: " + generated_by);
    }
    std::vector<int> version;
    {
        std::string part;
        const std::string numbers = generated_by.substr(9) + ".";
        bool malformed = false;
        for (char c : numbers) {
            if (c == '.') {
                if (part.empty() || !std::all_of(part.begin(), part.end(), [](unsigned char d) { return std::isdigit(d); })) {
                    malformed = true;
                } else {
                    version.push_back(std::stoi(part));
                }
                part.clear();
            } else {
                part += c;
            }
        }
        if (malformed || version.size() != 3) {
            throw ExitError("!!! ERROR. Malformed 'generated-by' attribute: " + generated_by);
        }
    }
    const std::string created = timestamp(text_attribute(attrs, "creation-date", present));
    const std::string update_date = text_attribute(attrs, "update-date", present);
    const std::string updated = present && !update_date.empty() ? timestamp(update_date) : "None";
    if (!silent) {
        console.out("### Info from input file");
        console.out("... File path: " + infile);
        console.out("... Found creation timestamp: " + created);
        console.out("... Found update timestamp: " + updated);
        console.out("... Found resolutions: " + resolutions_text);
        console.out("... Found version: " + generated_by.substr(9));
        console.out(std::string("... Target version: ") + kVersion);
    }
    const std::string writefile = outfile.empty() ? infile : outfile;
    const std::vector<Update> updates = prepare_updates(version);
    if (updates.empty()) {
        console.out("### No updates found!\n... Exiting");
        return;
    }
    if (!silent) {
        console.out(std::string("### Updates found. Will upgrade hic2cool file to version ") + kVersion);
        console.out("### This is what will change:");
        for (const auto& u : updates) {
            console.out("- " + u.title + "\n    - Effect: " + u.effect + "\n    - Detail: " + u.detail);
        }
        console.out("### Will write to " + writefile + "\n### Continue? [y/n]");
        std::string response;
        const auto accepted = [](std::string r) {
            std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c) { return std::tolower(c); });
            return r;
        };
        while (response.empty() || (accepted(response) != "y" && accepted(response) != "n" &&
                                    accepted(response) != "yes" && accepted(response) != "no")) {
            response = console.input("[y/n] ");
        }
        if (accepted(response) == "n" || accepted(response) == "no") {
            return;
        }
    }
    // run_hic2cool_updates
    if (infile != writefile) {
        std::filesystem::copy_file(infile, writefile, std::filesystem::copy_options::overwrite_existing);
    }
    console.out("### Updating...");
    for (const auto& u : updates) {
        console.out("... Running: " + u.title);
        u.run(writefile, console);
        console.out("... Finished: " + u.title);
    }
    const std::string generated = std::string("hic2cool-") + kVersion;
    const std::string now = utcnow_isoformat();
    {
        h5::File file(writefile, h5::Mode::ReadWrite);
        if (file.is_group("/resolutions")) {
            for (const auto& res : resolution_groups(file)) {
                file.set_attribute("/resolutions/" + res, "generated-by", generated);
                file.set_attribute("/resolutions/" + res, "update-date", now);
                console.out("... Updated metadata for resolution " + res);
            }
        } else {
            file.set_attribute("/", "generated-by", generated);
            file.set_attribute("/", "update-date", now);
            console.out("... Updated metadata");
        }
    }
    console.out("### Finished! Output written to: " + writefile);
}

}  // namespace hic2cool
