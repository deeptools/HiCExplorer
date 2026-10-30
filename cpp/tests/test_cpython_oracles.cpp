// Own implementations of behaviour that the Python interpreter and two Python
// packages define, compared with fixtures that the real Python wrote (the
// generator scripts are named at each test). No test here depends on Python
// being installed.

#include <doctest/doctest.h>

#include <cstdint>
#include <cmath>
#include <cstring>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/py_intervaltree.hpp"
#include "hicx/py_set_order.hpp"
#include "hicx/py_text.hpp"
#include "../core/io/hic2cool/src/pyformat.hpp"

namespace {

const std::string kFixtures = HICX_CPYTHON_FIXTURE_DIR;

namespace po = hicx::pyorder;

constexpr std::uint64_t kMask64 = ~std::uint64_t{0};

std::uint64_t splitmix(std::uint64_t& state) {
    state += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// A key of a fixture case: up to three integers, mapped to an identifier.
struct Key {
    std::int64_t part[3] = {0, 0, 0};
    bool operator<(const Key& o) const {
        for (int i = 0; i < 3; ++i) {
            if (part[i] != o.part[i]) return part[i] < o.part[i];
        }
        return false;
    }
};

struct Universe {
    int arity = 1;
    std::map<Key, std::int64_t> ids;
    std::vector<Key> keys;

    std::int64_t id_of(const Key& key) {
        auto [it, inserted] = ids.emplace(key, static_cast<std::int64_t>(keys.size()));
        if (inserted) keys.push_back(key);
        return it->second;
    }
    std::int64_t hash_of(const Key& key) const {
        if (arity == 1) return po::hash_of_int(key.part[0]);
        if (arity == 2) return po::hash_of_int_pair(key.part[0], key.part[1]);
        return po::hash_of_tuple({po::hash_of_int(key.part[0]), po::hash_of_int(key.part[1]),
                                  po::hash_of_int(key.part[2])});
    }
};

Key read_key(std::istream& in, int arity) {
    Key key;
    for (int i = 0; i < arity; ++i) in >> key.part[i];
    return key;
}

std::uint64_t fnv_digest(const std::vector<Key>& order, int arity) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (const Key& key : order) {
        for (int i = 0; i < arity; ++i) {
            const auto v = static_cast<std::uint64_t>(key.part[i]);
            for (int b = 0; b < 8; ++b) {
                h = (h ^ ((v >> (8 * b)) & 0xFF)) * 0x100000001B3ULL;
            }
        }
    }
    return h;
}

std::vector<Key> bulk_keys(int arity, std::uint64_t seed, std::size_t count) {
    std::vector<Key> keys;
    std::uint64_t state = seed;
    for (std::size_t j = 0; j < count; ++j) {
        const std::uint64_t x = splitmix(state);
        Key key;
        if (arity == 1) {
            const std::uint64_t sel = x & 3;
            const std::uint64_t four_n = 4 * count;
            if (sel == 0) key.part[0] = static_cast<std::int64_t>(x);
            else if (sel == 1) key.part[0] = static_cast<std::int64_t>((x >> 2) % four_n);
            else if (sel == 2) key.part[0] = -static_cast<std::int64_t>((x >> 2) % four_n) - 1;
            else key.part[0] = static_cast<std::int64_t>(((x >> 2) % (1ULL << 40)) << 21);
        } else if (arity == 2) {
            key.part[0] = static_cast<std::int64_t>((x >> 8) % 100003) - 50000;
            key.part[1] = static_cast<std::int64_t>((x >> 40) % 777);
        } else {
            key.part[0] = static_cast<std::int64_t>((x >> 8) % 1009) - 500;
            key.part[1] = static_cast<std::int64_t>((x >> 24) % 211);
            key.part[2] = static_cast<std::int64_t>((x >> 40) % 61) - 30;
        }
        keys.push_back(key);
    }
    return keys;
}

int arity_of(const std::string& kind) { return kind == "i" ? 1 : (kind == "t2" ? 2 : 3); }

}  // namespace

TEST_CASE("OrderedSet iterates like a Python set (replay of set_order.txt)") {
    std::ifstream in(kFixtures + "/set_order.txt");
    REQUIRE(in.good());
    std::string word;
    in >> word;
    REQUIRE(word == "version");
    in >> word;
    CHECK(word == "3.12.7");

    Universe uni;
    po::OrderedSet set;
    std::size_t cases = 0;
    std::size_t bulk = 0;
    const auto add = [&](po::OrderedSet& s, const Key& k) { s.add(uni.id_of(k), uni.hash_of(k)); };
    const auto read_set = [&](std::istream& is) {
        std::size_t n = 0;
        is >> n;
        std::vector<Key> keys;
        for (std::size_t i = 0; i < n; ++i) keys.push_back(read_key(is, uni.arity));
        return keys;
    };
    std::string name;
    while (in >> word) {
        if (word == "case") {
            std::string kind;
            in >> name >> kind;
            uni = Universe{};
            uni.arity = arity_of(kind);
            set = po::OrderedSet{};
        } else if (word == "a") {
            add(set, read_key(in, uni.arity));
        } else if (word == "d") {
            const Key k = read_key(in, uni.arity);
            set.discard(uni.id_of(k), uni.hash_of(k));
        } else if (word == "u") {
            for (const Key& k : read_set(in)) add(set, k);
        } else if (word == "U") {
            po::OrderedSet other;
            for (const Key& k : read_set(in)) add(other, k);
            set.update(other);
        } else if (word == "x") {
            po::OrderedSet other;
            for (const Key& k : read_set(in)) add(other, k);
            set.difference_update(other);
        } else if (word == "c") {
            po::OrderedSet copy;
            copy.update(set);
            set = copy;
        } else if (word == "order") {
            const std::vector<Key> expected = read_set(in);
            std::vector<std::int64_t> expected_ids;
            for (const Key& k : expected) expected_ids.push_back(uni.id_of(k));
            INFO("case " << name);
            CHECK(set.items() == expected_ids);
            ++cases;
        } else if (word == "bulk") {
            std::string kind;
            std::uint64_t seed = 0;
            std::size_t count = 0;
            std::uint64_t digest = 0;
            in >> name >> kind >> seed >> count >> digest;
            Universe u;
            u.arity = arity_of(kind);
            const std::vector<Key> keys = bulk_keys(u.arity, seed, count);
            po::OrderedSet s;
            const auto put = [&](po::OrderedSet& t, const Key& k) { t.add(u.id_of(k), u.hash_of(k)); };
            for (const Key& k : keys) put(s, k);
            for (std::size_t j = 0; j < keys.size(); j += 3) s.discard(u.id_of(keys[j]), u.hash_of(keys[j]));
            for (std::size_t j = 0; j < keys.size(); j += 5) put(s, keys[j]);
            po::OrderedSet seventh;
            for (std::size_t j = 0; j < keys.size(); j += 7) put(seventh, keys[j]);
            s.update(seventh);
            po::OrderedSet eleventh;
            for (std::size_t j = 0; j < keys.size(); j += 11) put(eleventh, keys[j]);
            s.difference_update(eleventh);
            std::vector<Key> order;
            for (const std::int64_t id : s.items()) order.push_back(u.keys[static_cast<std::size_t>(id)]);
            INFO("bulk " << name);
            CHECK(fnv_digest(order, u.arity) == digest);
            ++bulk;
        } else {
            FAIL("unknown record " << word);
        }
    }
    CHECK(cases > 80);
    CHECK(bulk == 6);
}

TEST_CASE("hash_of_int and hash_of_tuple follow the interpreter's rules") {
    CHECK(po::hash_of_int(0) == 0);
    CHECK(po::hash_of_int(-1) == -2);
    CHECK(po::hash_of_int(-2) == -2);
    CHECK(po::hash_of_int((std::int64_t{1} << 61) - 1) == 0);
    CHECK(po::hash_of_int(std::int64_t{1} << 61) == 1);
    CHECK(po::hash_of_int(-(std::int64_t{1} << 61)) == -1 * 1 - 1);  // -1 becomes -2
}

TEST_CASE("Tree returns overlap sets in intervaltree's order (replay of intervaltree.txt)") {
    std::ifstream in(kFixtures + "/intervaltree.txt");
    REQUIRE(in.good());
    std::string word;
    std::string python_version;
    std::string package_version;
    in >> word >> python_version >> package_version;
    REQUIRE(word == "version");
    CHECK(python_version == "3.12.7");
    CHECK(package_version == "3.1.0");

    hicx::pytree::Tree tree;
    std::string name;
    std::size_t next = 0;
    std::size_t queries = 0;
    while (in >> word) {
        if (word == "case") {
            in >> name;
            tree = hicx::pytree::Tree();
            next = 0;
        } else if (word == "a") {
            std::int64_t b = 0;
            std::int64_t e = 0;
            in >> b >> e;
            tree.add(b, e, next++);
        } else if (word == "r") {
            std::int32_t key = 0;
            in >> key;
            CHECK(tree.remove(key));
        } else if (word == "q") {
            std::int64_t b = 0;
            std::int64_t e = 0;
            std::size_t n = 0;
            in >> b >> e >> n;
            std::vector<std::int32_t> expected(n);
            for (auto& key : expected) in >> key;
            INFO("case " << name << " query " << b << " " << e);
            CHECK(tree.overlap(b, e) == expected);
            ++queries;
        } else {
            FAIL("unknown record " << word);
        }
    }
    CHECK(queries > 4000);
}

namespace {

std::string from_hex(const std::string& hex) {
    std::string out;
    if (hex == "-") return out;
    for (std::size_t i = 0; i + 1 < hex.size(); i += 2) {
        out.push_back(static_cast<char>(std::stoi(hex.substr(i, 2), nullptr, 16)));
    }
    return out;
}

std::string encode_code_point(std::uint32_t c) {  // UTF-8 with surrogatepass
    std::string out;
    if (c < 0x80) {
        out.push_back(static_cast<char>(c));
    } else if (c < 0x800) {
        out.push_back(static_cast<char>(0xC0 | (c >> 6)));
        out.push_back(static_cast<char>(0x80 | (c & 0x3F)));
    } else if (c < 0x10000) {
        out.push_back(static_cast<char>(0xE0 | (c >> 12)));
        out.push_back(static_cast<char>(0x80 | ((c >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (c & 0x3F)));
    } else {
        out.push_back(static_cast<char>(0xF0 | (c >> 18)));
        out.push_back(static_cast<char>(0x80 | ((c >> 12) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | ((c >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (c & 0x3F)));
    }
    return out;
}

std::vector<std::string> read_hex_list(std::istream& in) {
    std::size_t n = 0;
    in >> n;
    std::vector<std::string> out(n);
    for (auto& item : out) {
        std::string hex;
        in >> hex;
        item = from_hex(hex);
    }
    return out;
}

std::vector<std::string> to_strings(const std::vector<std::string_view>& views) {
    return {views.begin(), views.end()};
}

}  // namespace

TEST_CASE("py_text follows str, strict UTF-8 and universal newlines (replay of py_text.txt)") {
    std::ifstream in(kFixtures + "/py_text.txt");
    REQUIRE(in.good());
    std::string word;
    in >> word;
    REQUIRE(word == "version");
    in >> word;
    CHECK(word == "3.12.7");
    in >> word;
    REQUIRE(word == "spaces");
    std::size_t n = 0;
    in >> n;
    std::vector<bool> is_space(0x110000, false);
    for (std::size_t i = 0; i < n; ++i) {
        std::uint32_t c = 0;
        in >> c;
        is_space[c] = true;
    }
    // Every code point, surrogates included.
    for (std::uint32_t c = 0; c < 0x110000; ++c) {
        const std::string text = encode_code_point(c);
        const bool whitespace = hicx::pytext::whitespace_length(text, 0) == text.size();
        if (whitespace != is_space[c]) {
            FAIL("whitespace classification differs at U+" << std::hex << c);
        }
    }

    std::string text;
    std::size_t strings = 0;
    std::size_t decodings = 0;
    std::size_t files = 0;
    while (in >> word) {
        if (word == "t") {
            std::string hex;
            in >> hex;
            text = from_hex(hex);
            ++strings;
        } else if (word == "s") {
            std::string hex;
            in >> hex;
            CHECK(std::string(hicx::pytext::strip(text)) == from_hex(hex));
        } else if (word == "w") {
            CHECK(to_strings(hicx::pytext::split_whitespace(text)) == read_hex_list(in));
        } else if (word == "p") {
            std::string sep;
            in >> sep;
            CHECK(to_strings(hicx::pytext::split(text, from_hex(sep))) == read_hex_list(in));
        } else if (word == "r") {
            std::string hex;
            in >> hex;
            CHECK(hicx::pytext::remove_whitespace(text) == from_hex(hex));
        } else if (word == "v") {
            std::string hex;
            long expected = 0;
            in >> hex >> expected;
            const auto got = hicx::pytext::first_invalid_utf8(from_hex(hex));
            INFO("bytes " << hex);
            CHECK((got ? static_cast<long>(*got) : -1L) == expected);
            ++decodings;
        } else if (word == "e") {
            std::string hex;
            std::string repr_hex;
            in >> hex >> repr_hex;
            INFO("bytes " << hex);
            CHECK(hicx::pytext::repr(from_hex(hex)) == from_hex(repr_hex));
        } else if (word == "reprall") {
            std::uint64_t expected = 0;
            in >> expected;
            std::uint64_t h = 0xCBF29CE484222325ULL;
            for (std::uint32_t c = 0; c < 0x110000; ++c) {
                if (c >= 0xD800 && c < 0xE000) continue;
                const std::string r = hicx::pytext::repr(encode_code_point(c)) + '\0';
                for (const char ch : r) h = (h ^ static_cast<unsigned char>(ch)) * 0x100000001B3ULL;
            }
            CHECK(h == expected);
        } else if (word == "n") {
            std::string hex;
            in >> hex;
            const std::string file = from_hex(hex);
            CHECK(hicx::pytext::universal_newline_lines(file) == read_hex_list(in));
            ++files;
        } else {
            FAIL("unknown record " << word);
        }
    }
    CHECK(strings == 4000);
    CHECK(decodings > 6000);
    CHECK(files == 2000);
}

TEST_CASE("hic2cool's Python text matches repr and pandas (replay of pyformat.txt)") {
    namespace hd = hic2cool::detail;
    std::ifstream in(kFixtures + "/pyformat.txt");
    REQUIRE(in.good());
    std::string word;
    std::string python_version;
    std::string pandas_version;
    in >> word >> python_version >> pandas_version;
    REQUIRE(word == "version");
    CHECK(python_version == "3.12.7");
    CHECK(pandas_version == "2.2.3");
    std::size_t floats = 0;
    std::size_t frames = 0;
    while (in >> word) {
        if (word == "f") {
            std::string bits;
            std::string expected;
            in >> bits >> expected;
            const std::uint64_t raw = std::stoull(bits, nullptr, 16);
            double value = 0;
            std::memcpy(&value, &raw, sizeof(value));
            INFO("bits " << bits);
            CHECK(hd::py_repr(value) == from_hex(expected));
            ++floats;
        } else if (word == "l") {
            std::size_t n = 0;
            in >> n;
            std::vector<std::int64_t> values(n);
            for (auto& v : values) in >> v;
            std::string expected;
            in >> expected;
            CHECK(hd::py_repr(values) == from_hex(expected));
        } else if (word == "ls") {
            const std::vector<std::string> values = read_hex_list(in);
            std::string expected;
            in >> expected;
            CHECK(hd::py_repr(values) == from_hex(expected));
        } else if (word == "d") {
            std::size_t n = 0;
            in >> n;
            std::vector<std::pair<std::int64_t, std::string>> values(n);
            for (auto& [k, v] : values) {
                std::string hex;
                in >> k >> hex;
                v = from_hex(hex);
            }
            std::string expected;
            in >> expected;
            CHECK(hd::py_repr(values) == from_hex(expected));
        } else if (word == "df") {
            std::string norm;
            std::size_t n = 0;
            in >> norm >> n;
            std::vector<hd::BinsPreviewRow> rows(n);
            for (auto& row : rows) {
                std::string chrom;
                std::string bits;
                in >> chrom >> row.start >> row.end >> bits;
                row.chrom = from_hex(chrom);
                const std::uint64_t raw = std::stoull(bits, nullptr, 16);
                std::memcpy(&row.value, &raw, sizeof(row.value));
            }
            std::string expected;
            in >> expected;
            INFO("frame " << frames);
            CHECK(hd::pandas_bins_head(from_hex(norm), rows) == from_hex(expected));
            ++frames;
        } else {
            FAIL("unknown record " << word);
        }
    }
    CHECK(floats == 3000);
    CHECK(frames == 1500);
}

TEST_CASE("utcnow_isoformat has the shape of datetime.isoformat") {
    // YYYY-MM-DDTHH:MM:SS, plus .ffffff unless the microsecond is zero.
    const std::regex shape(R"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d{6})?)");
    for (int i = 0; i < 50; ++i) {
        CHECK(std::regex_match(hic2cool::detail::utcnow_isoformat(), shape));
    }
}

TEST_CASE("python_int and python_float accept what int() and float() accept (replay of py_number.txt)") {
    std::ifstream in(kFixtures + "/py_number.txt");
    REQUIRE(in.good());
    std::string word;
    std::string version;
    in >> word >> version;
    REQUIRE(word == "version");
    std::size_t ints = 0;
    std::size_t floats = 0;
    while (in >> word) {
        std::string hex;
        std::string expected;
        in >> hex >> expected;
        const std::string text = from_hex(hex);
        INFO("text " << hex);
        if (word == "i") {
            std::int64_t value = 0;
            const bool ok = hicx::cli::python_int(text, &value);
            if (expected == "X" || expected == "B") {
                CHECK_FALSE(ok);
            } else {
                CHECK(ok);
                CHECK(value == std::stoll(expected));
            }
            ++ints;
        } else if (word == "g") {
            double value = 0;
            const bool ok = hicx::cli::python_float(text, &value);
            if (expected == "X") {
                CHECK_FALSE(ok);
            } else {
                REQUIRE(ok);
                const std::uint64_t raw = std::stoull(expected, nullptr, 16);
                double want = 0;
                std::memcpy(&want, &raw, sizeof(want));
                if (std::isnan(want)) {
                    CHECK(std::isnan(value));
                } else {
                    CHECK(value == want);
                    CHECK(std::signbit(value) == std::signbit(want));
                }
            }
            ++floats;
        } else {
            FAIL("unknown record " << word);
        }
    }
    CHECK(ints == 40000);
    CHECK(floats == 40000);
}
