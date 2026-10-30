// The three behaviours of cooler 0.12.0/0.12.1 that coolercpp adopts, checked
// against files cooler 0.12.1 wrote (fixtures/cooler012/cases.txt, made by
// cpp/scripts/gen_cooler012_fixture.py, which also documents the case format):
//
//   1. chroms/length, bins/start and bins/end are int64 when a chromosome is
//      longer than 2147483647 bp, int32 otherwise (cooler #480),
//
// Every case writes the same calls with coolercpp and compares a structural
// dump of the file (groups, links, attributes except creation-date and
// generated-by, dataset dtype, shape, chunking, maxshape, filters, values)
// with the dump of the fixture, line by line. The failure tests inject an
// invalid pixel, which makes create_cooler throw after the datasets exist.

#include <doctest/doctest.h>
#include <hdf5.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <optional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "coolercpp/create.hpp"
#include "coolercpp/errors.hpp"

namespace {

namespace fs = std::filesystem;
using namespace coolercpp;

const std::string kFixtures = std::string(HICX_CPYTHON_FIXTURE_DIR) + "/../cooler012/";
const std::string kOutput = std::string(HICX_TEST_OUTPUT_DIR) + "/cooler012";

// ---------------------------------------------------------------------------
// The structural dump (same text as dump() of the generator script)

std::string hex64(std::uint64_t v) {
    char buffer[24];
    std::snprintf(buffer, sizeof buffer, "%016llx", static_cast<unsigned long long>(v));
    return buffer;
}

std::string attribute_text(hid_t object, const char* name) {
    const hid_t attr = H5Aopen(object, name, H5P_DEFAULT);
    const hid_t type = H5Aget_type(attr);
    std::string out;
    switch (H5Tget_class(type)) {
        case H5T_INTEGER: {
            std::int64_t v = 0;
            H5Aread(attr, H5T_NATIVE_INT64, &v);
            out = "i" + std::to_string(v);
            break;
        }
        case H5T_FLOAT: {
            double v = 0;
            H5Aread(attr, H5T_NATIVE_DOUBLE, &v);
            std::uint64_t bits;
            std::memcpy(&bits, &v, 8);
            out = "f" + hex64(bits);
            break;
        }
        case H5T_STRING: {
            if (H5Tis_variable_str(type) > 0) {
                char* text = nullptr;
                const hid_t mem = H5Tcopy(H5T_C_S1);
                H5Tset_size(mem, H5T_VARIABLE);
                H5Tset_cset(mem, H5Tget_cset(type));
                H5Aread(attr, mem, &text);
                out = "s" + std::string(text != nullptr ? text : "");
                H5free_memory(text);
                H5Tclose(mem);
            } else {
                std::vector<char> buffer(H5Tget_size(type) + 1, '\0');
                H5Aread(attr, type, buffer.data());
                out = "s" + std::string(buffer.data());
            }
            break;
        }
        default:
            out = "?";
    }
    H5Tclose(type);
    H5Aclose(attr);
    return out;
}

std::string attributes_text(hid_t object) {
    std::map<std::string, std::string> items;
    H5O_info2_t info{};
    H5Oget_info3(object, &info, H5O_INFO_NUM_ATTRS);
    for (hsize_t i = 0; i < info.num_attrs; ++i) {
        const ssize_t length = H5Aget_name_by_idx(object, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0,
                                                  H5P_DEFAULT);
        std::string name(static_cast<std::size_t>(length), '\0');
        H5Aget_name_by_idx(object, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(),
                           name.size() + 1, H5P_DEFAULT);
        if (name == "creation-date" || name == "generated-by") {
            continue;
        }
        items[name] = attribute_text(object, name.c_str());
    }
    std::string out;
    for (const auto& [name, value] : items) {
        out += (out.empty() ? "" : "|") + name + "=" + value;
    }
    return out;
}

std::string links_text(hid_t group) {
    std::string out;
    H5G_info_t info{};
    H5Gget_info(group, &info);
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        const ssize_t length =
            H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
        std::string name(static_cast<std::size_t>(length), '\0');
        H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(), name.size() + 1,
                           H5P_DEFAULT);
        out += (out.empty() ? "" : ",") + name;
    }
    return out;
}

std::string dtype_text(hid_t type) {
    const H5T_class_t cls = H5Tget_class(type);
    const auto base = [](hid_t t) {
        const char kind = H5Tget_class(t) == H5T_FLOAT ? 'f' : (H5Tget_sign(t) == H5T_SGN_NONE ? 'u' : 'i');
        return std::string(1, kind) + std::to_string(H5Tget_size(t));
    };
    if (cls == H5T_ENUM) {
        const hid_t super = H5Tget_super(type);
        std::map<long long, std::string> members;
        const int n = H5Tget_nmembers(type);
        for (int i = 0; i < n; ++i) {
            char* name = H5Tget_member_name(type, static_cast<unsigned>(i));
            long long value = 0;
            H5Tget_member_value(type, static_cast<unsigned>(i), &value);
            members[value] = name;
            H5free_memory(name);
        }
        std::string out = "enum:" + base(super) + "{";
        bool first = true;
        for (const auto& [value, name] : members) {
            out += (first ? "" : ",") + name + "=" + std::to_string(value);
            first = false;
        }
        H5Tclose(super);
        return out + "}";
    }
    if (cls == H5T_STRING) {
        return "S" + std::to_string(H5Tget_size(type));
    }
    return base(type);
}

std::string dataset_text(hid_t dataset, const std::string& path) {
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims = 0;
    hsize_t maxdims = 0;
    H5Sget_simple_extent_dims(space, &dims, &maxdims);
    const hid_t type = H5Dget_type(dataset);
    const hid_t plist = H5Dget_create_plist(dataset);
    std::string chunks = "-";
    if (H5Pget_layout(plist) == H5D_CHUNKED) {
        hsize_t chunk = 0;
        H5Pget_chunk(plist, 1, &chunk);
        chunks = std::to_string(chunk);
    }
    std::string filters;
    for (int i = 0; i < H5Pget_nfilters(plist); ++i) {
        unsigned flags = 0;
        size_t cd_n = 8;
        unsigned cd[8] = {0};
        char name[64];
        const H5Z_filter_t id = H5Pget_filter2(plist, static_cast<unsigned>(i), &flags, &cd_n, cd, 64,
                                               name, nullptr);
        filters += (i ? ";" : "") + std::to_string(id) + ":";
        for (size_t k = 0; k < cd_n; ++k) {
            filters += (k ? "," : "") + std::to_string(cd[k]);
        }
    }
    // values: as the generator prints them
    std::vector<std::string> items;
    bool integer = false;
    std::int64_t sum = 0;
    const std::size_t n = static_cast<std::size_t>(dims);
    if (n > 0) {
        const H5T_class_t cls = H5Tget_class(type);
        if (cls == H5T_STRING) {
            const std::size_t width = H5Tget_size(type);
            std::vector<char> buffer(n * width);
            H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
            for (std::size_t i = 0; i < n; ++i) {
                items.emplace_back(buffer.data() + i * width, strnlen(buffer.data() + i * width, width));
            }
        } else if (cls == H5T_FLOAT) {
            std::vector<double> values(n);
            H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
            for (const double v : values) {
                std::uint64_t bits;
                std::memcpy(&bits, &v, 8);
                items.push_back(hex64(bits));
            }
        } else {
            integer = true;
            std::vector<std::int64_t> values(n);
            H5Dread(dataset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
            for (const std::int64_t v : values) {
                items.push_back(std::to_string(v));
                sum += v;
            }
        }
    }
    std::string values;
    if (items.size() <= 40) {
        for (std::size_t i = 0; i < items.size(); ++i) {
            values += (i ? "," : "") + items[i];
        }
    } else {
        for (std::size_t i = 0; i < 32; ++i) {
            values += (i ? "," : "") + items[i];
        }
        values += ",...";
        for (std::size_t i = items.size() - 8; i < items.size(); ++i) {
            values += "," + items[i];
        }
    }
    if (integer && n > 0) {
        values += " sum=" + std::to_string(sum);
    }
    std::string out = "D " + path + " dtype=" + dtype_text(type) + " shape=" + std::to_string(dims) +
                      " chunks=" + chunks + " maxshape=" +
                      (maxdims == H5S_UNLIMITED ? std::string("inf") : std::to_string(maxdims)) +
                      " filters=" + filters + " values=" + values;
    H5Pclose(plist);
    H5Tclose(type);
    H5Sclose(space);
    return out;
}

void dump_group(hid_t file, const std::string& path, std::vector<std::string>& lines) {
    const hid_t group = H5Gopen2(file, path.c_str(), H5P_DEFAULT);
    lines.push_back("G " + path + " links=" + links_text(group) + " attrs=" + attributes_text(group));
    H5G_info_t info{};
    H5Gget_info(group, &info);
    std::vector<std::string> names;
    {
        const std::string joined = links_text(group);
        std::stringstream stream(joined);
        std::string name;
        while (std::getline(stream, name, ',')) {
            names.push_back(name);
        }
    }
    for (const std::string& name : names) {
        const std::string child = (path == "/" ? "" : path) + "/" + name;
        H5O_info2_t oinfo{};
        H5Oget_info_by_name3(file, child.c_str(), &oinfo, H5O_INFO_BASIC, H5P_DEFAULT);
        if (oinfo.type == H5O_TYPE_GROUP) {
            dump_group(file, child, lines);
        } else {
            const hid_t dataset = H5Dopen2(file, child.c_str(), H5P_DEFAULT);
            lines.push_back(dataset_text(dataset, child));
            H5Dclose(dataset);
        }
    }
    H5Gclose(group);
}

std::vector<std::string> dump_file(const std::string& path) {
    std::vector<std::string> lines;
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file >= 0);
    dump_group(file, "/", lines);
    H5Fclose(file);
    return lines;
}

// ---------------------------------------------------------------------------
// Cases

struct Create {
    std::string uri, mode, input, chroms, pixels;
    std::int64_t binsize = 0;
    int chunks = 1;
};

struct Case {
    std::string name;
    std::vector<Create> creates;
    std::vector<std::string> expected;
};

std::vector<std::string> split(const std::string& text, char sep) {
    std::vector<std::string> out;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, sep)) {
        out.push_back(item);
    }
    return out;
}

std::vector<Case> read_cases(const std::string& file) {
    std::ifstream in(kFixtures + file);
    REQUIRE_MESSAGE(in.good(), "missing fixture " << kFixtures << file);
    std::vector<Case> cases;
    std::string line;
    bool in_dump = false;
    while (std::getline(in, line)) {
        if (line.rfind("case ", 0) == 0) {
            cases.push_back(Case{line.substr(5), {}, {}});
            in_dump = false;
        } else if (line == "dump") {
            in_dump = true;
        } else if (line == "end") {
            in_dump = false;
        } else if (in_dump) {
            cases.back().expected.push_back(line);
        } else if (line.rfind("create ", 0) == 0) {
            const std::vector<std::string> fields = split(line, ' ');
            Create c;
            c.uri = fields[1];
            for (std::size_t i = 2; i < fields.size(); ++i) {
                const std::size_t eq = fields[i].find('=');
                const std::string key = fields[i].substr(0, eq);
                const std::string value = fields[i].substr(eq + 1);
                if (key == "mode") c.mode = value;
                else if (key == "input") c.input = value;
                else if (key == "chroms") c.chroms = value;
                else if (key == "binsize") c.binsize = std::stoll(value);
                else if (key == "chunks") c.chunks = std::stoi(value);
                else if (key == "pixels") c.pixels = value;
            }
            cases.back().creates.push_back(c);
        }
    }
    return cases;
}

// cooler.util.binnify
Table make_bins(const std::string& chroms, std::int64_t binsize) {
    std::vector<std::string> chrom;
    std::vector<std::int64_t> start;
    std::vector<std::int64_t> end;
    for (const std::string& item : split(chroms, ',')) {
        const std::size_t colon = item.find(':');
        const std::string name = item.substr(0, colon);
        const std::int64_t length = std::stoll(item.substr(colon + 1));
        for (std::int64_t s = 0; s < length; s += binsize) {
            chrom.push_back(name);
            start.push_back(s);
            end.push_back(std::min(s + binsize, length));
        }
    }
    Table bins;
    bins.set("chrom", Column(chrom));
    bins.set("start", Column(start));
    bins.set("end", Column(end));
    return bins;
}

Table make_pixels(const std::string& text, std::size_t lo, std::size_t hi) {
    std::vector<std::int64_t> b1, b2, count;
    const std::vector<std::string> items = text.empty() ? std::vector<std::string>{} : split(text, ',');
    for (std::size_t i = lo; i < hi; ++i) {
        const std::vector<std::string> f = split(items[i], ':');
        b1.push_back(std::stoll(f[0]));
        b2.push_back(std::stoll(f[1]));
        count.push_back(std::stoll(f[2]));
    }
    Table t;
    t.set("bin1_id", Column(b1));
    t.set("bin2_id", Column(b2));
    t.set("count", Column(count));
    return t;
}

std::size_t pixel_count(const std::string& text) { return text.empty() ? 0 : split(text, ',').size(); }

void run_create(const std::string& directory, const Create& c) {
    const Table bins = make_bins(c.chroms, c.binsize);
    CreateOptions options;
    options.mode = c.mode;
    options.creation_date = "2000-01-01T00:00:00";
    const std::string target = directory + "/" + c.uri;
    const std::size_t total = pixel_count(c.pixels);
    if (c.input == "table") {
        create_cooler(target, bins, make_pixels(c.pixels, 0, total), options);
        return;
    }
    options.ordered = c.input == "ordered";
    std::vector<Table> parts;
    for (int i = 0; i < c.chunks; ++i) {
        parts.push_back(make_pixels(c.pixels, total * i / c.chunks, total * (i + 1) / c.chunks));
    }
    create_cooler(target, bins, iterate_chunks(std::move(parts)), options);
}

std::string fresh_directory(const std::string& name) {
    const std::string directory = kOutput + "/" + name;
    fs::remove_all(directory);
    fs::create_directories(directory);
    return directory;
}

std::vector<std::string> directory_listing(const std::string& directory) {
    std::vector<std::string> names;
    for (const auto& entry : fs::directory_iterator(directory)) {
        names.push_back(entry.path().filename().string());
    }
    std::sort(names.begin(), names.end());
    return names;
}

std::string file_bytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

}  // namespace

// Writes every case of a fixture file with coolercpp and compares the dumps.
void check_cases(const std::string& file, std::size_t minimum) {
    const std::vector<Case> cases = read_cases(file);
    REQUIRE(cases.size() >= minimum);
    for (const Case& c : cases) {
        const std::string directory = fresh_directory(c.name);
        for (const Create& create : c.creates) {
            run_create(directory, create);
        }
        // No staging file or group is left behind.
        CHECK_MESSAGE(directory_listing(directory) == std::vector<std::string>{"main.cool"}, c.name);
        const std::vector<std::string> got = dump_file(directory + "/main.cool");
        REQUIRE_MESSAGE(got.size() == c.expected.size(), c.name);
        for (std::size_t i = 0; i < got.size(); ++i) {
            CHECK_MESSAGE(got[i] == c.expected[i], c.name << " line " << i);
        }
    }
}

TEST_CASE("ordinary coolers, mcool resolutions and root appends are written as before") {
    check_cases("ordinary.txt", 6);
}

TEST_CASE("chromosomes beyond the int32 maximum make chroms/length, bins/start and bins/end int64") {
    check_cases("int64.txt", 8);
}

TEST_CASE("a cooler without pixels has empty pixel datasets") {
    check_cases("empty.txt", 7);
}
