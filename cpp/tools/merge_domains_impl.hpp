// The computational core of hicMergeDomains, a port of
// hicexplorer/hicMergeDomains.py.
//
// It lives in its own translation unit, as detect_loops_impl does, so that
// cpp/tests/test_merge_domains.cpp links against exactly the code the tool
// runs. Every function here is a literal port of the Python function of the
// same name. The Python works on lists of lists of strings with index loops,
// and several of its results depend on details of that representation, so the
// port keeps the representation instead of reinterpreting the algorithm:
//
//  * A row is the list of tab separated fields of one line, and a domain list
//    is a list of *references* to rows. merge_list can put the same row object
//    into the merged list twice, and add_id then numbers that object twice, so
//    both slots print the later ID. RowPool and index lists reproduce that
//    aliasing (pinned by the reversed file order case, where one TAD is
//    written twice under the same ID).
//  * `x in list` and `list.remove(x)` compare rows by content, not identity.
//  * Every list access that can raise IndexError in the Python raises
//    PythonError here, in the same evaluation order, so an input the
//    reference rejects is rejected too, with exit status 1.
//
// Nothing in the tool clusters anything: the Python imports
// scipy.cluster.hierarchy.linkage and dendrogram at hicMergeDomains.py:2 and
// never calls either.

#ifndef HICX_TOOLS_MERGE_DOMAINS_IMPL_HPP
#define HICX_TOOLS_MERGE_DOMAINS_IMPL_HPP

#include <cstddef>
#include <cstdint>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hicx::merge_domains {

// An exception the Python reference raises on this input. The tool reports
// it and exits 1, the status of an uncaught Python exception.
class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// An input on which the Python reference never terminates. The port refuses
// it with exit status 1 instead of reproducing an endless loop.
class ReferenceNeverTerminates : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

using Row = std::vector<std::string>;

// Every row read from a domain file is stored once; lists refer to rows by
// index, the way Python list slots refer to list objects.
struct RowPool {
    std::vector<Row> rows;
};

// ---------------------------------------------------------------------------
// Python string and number semantics

// str.rstrip() over the ASCII whitespace str.isspace() accepts.
[[nodiscard]] std::string py_rstrip(const std::string& text);
// str.split('\t').
[[nodiscard]] std::vector<std::string> split_tab(const std::string& text);
// [line.rstrip() for line in open(path)]: universal newlines, so '\n', '\r'
// and '\r\n' all end a line.
[[nodiscard]] std::vector<std::string> read_lines(const std::string& path);
// int(text) for base 10, including surrounding whitespace, a sign and single
// underscores between digits. Throws PythonError(ValueError) otherwise.
[[nodiscard]] std::int64_t py_int(const std::string& text);
// float(text). Throws PythonError(ValueError) when CPython would.
[[nodiscard]] double py_float(const std::string& text);

// ---------------------------------------------------------------------------
// The functions of hicMergeDomains.py

struct DomainList {
    std::vector<std::size_t> rows;
    std::int64_t bin_size = 10000000;
};

// create_list_of_file (:90-115). The "bin size" is taken from the start
// coordinate of each of the first 20 lines by keeping only its last non zero
// digit followed by its trailing zeros (710000 gives 10000, 1250000 gives
// 50000), and the smallest of those wins.
[[nodiscard]] DomainList create_list_of_file(RowPool& pool, const std::string& path);

// read_protein (:322-335): the first three fields of every line, grouped into
// runs of consecutive lines with the same chromosome.
using ProteinList = std::vector<std::vector<Row>>;
[[nodiscard]] ProteinList read_protein(const std::string& path);

struct ProteinBin {
    std::string chrom;
    std::int64_t left = 0;
    std::int64_t right = 0;
    std::int64_t count = 0;
};
using MergedProtein = std::vector<std::vector<ProteinBin>>;

// merge_protein (:338-363). The last bin of every chromosome is never
// appended, because the append only happens when a later peak opens a new bin.
[[nodiscard]] MergedProtein merge_protein(const ProteinList& proteins,
                                          std::int64_t bin_size,
                                          std::int64_t min_peak);

// compare_boundaries_protein (:366-391). Modifies b_list in place, as the
// Python does.
void compare_boundaries_protein(const RowPool& pool, std::vector<std::size_t>& b_list,
                                const MergedProtein& c_list, double para_score = 0.2);

// merge_list (:118-182).
[[nodiscard]] std::vector<std::size_t> merge_list(const RowPool& pool,
                                                  const std::vector<std::size_t>& d1,
                                                  const std::vector<std::size_t>& d2,
                                                  std::int64_t p_value);

// add_id (:185-191): overwrites field 4 of every referenced row, in list order.
void add_id(RowPool& pool, const std::vector<std::size_t>& list);

struct Relation {
    std::string chrom;
    std::string parent;
    std::vector<std::string> children;
};

// create_relationsship_list (:194-222) together with add_relation_to_list.
[[nodiscard]] std::vector<Relation> create_relationship_list(
    const RowPool& pool, const std::vector<std::size_t>& list, double percent);

// write_in_file (:242-264), the two layouts.
void write_domain_list(std::ostream& out, const RowPool& pool,
                       const std::vector<std::size_t>& list);
void write_relation_list(std::ostream& out, const std::vector<Relation>& relations);

// One graphviz.Digraph at the moment create_tree renders it.
struct TreeGraph {
    std::string chrom;     // the chromosome printed in "Saved relation tree of"
    std::string filename;  // prefix + '_' + chromosome
    std::string source;    // Digraph.source, byte for byte
};

// create_tree (:267-303) and create_small_list (:306-319). `render` is called
// at each point the Python calls g.render, after its print, and may throw.
void create_tree(const std::vector<Relation>& relations, const RowPool& pool,
                 const std::vector<std::size_t>& list, const std::string& prefix,
                 const std::function<void(const TreeGraph&)>& render);

// main (:403-423) up to add_id: every domain file read, filtered against the
// protein peaks when a protein file is given, merged in order and numbered.
// Merging only happens with more than one domain file, as in the Python.
struct MergedDomains {
    RowPool pool;
    std::vector<std::size_t> merged;
};
[[nodiscard]] MergedDomains merge_domain_files(const std::vector<std::string>& domain_files,
                                               const ProteinList* proteins,
                                               std::int64_t minimum_number_of_peaks,
                                               std::int64_t value);

// graphviz.parameters.verify_format: the lower cased format must be one of
// graphviz.parameters.FORMATS (graphviz 0.20.3).
[[nodiscard]] bool is_graphviz_format(const std::string& format);

}  // namespace hicx::merge_domains

#endif  // HICX_TOOLS_MERGE_DOMAINS_IMPL_HPP
