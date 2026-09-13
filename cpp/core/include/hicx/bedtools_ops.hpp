// The three bedtools operations pybedtools reaches for in the tier 2 tools:
// sort, merge and intersect -c. cpp/PLAN.md 3.5 decides to reimplement them
// rather than depend on the bedtools binary at runtime.
//
// The only hard part is sort, and it is hard for an unexpected reason.
// bedtools' sortBed loads the records into a std::map<string, vector<BED>> in
// file order and then, per chromosome, calls
//
//     sort(v.begin(), v.end(), sortByStart)
//
// where sortByStart compares the start coordinate and nothing else. That
// comparator is not a total order over records that share a start, and
// std::sort is not stable, so the relative order of records with equal starts
// is whatever libstdc++'s introsort leaves behind. It is deterministic for a
// given input sequence but it is neither the input order nor an order derived
// from the other fields.
//
// That order is observable: hicValidateLocations writes out a subset of the
// sorted loop table in its sorted order, and loops_1.bedgraph has 3,730 rows
// involved in a start tie. Reproducing it means running the same algorithm:
// the same std::sort, with the same comparator, over the same input sequence.
// The permutation introsort produces depends only on the length of the range
// and on the outcomes of the comparisons, not on the element type, so sorting
// (start, row) pairs here gives the same permutation bedtools gets sorting its
// BED structs. Verified on the real loop file: the 11,723 sorted rows match
// `bedtools sort` byte for byte, ties included.
//
// This does tie the port to libstdc++'s std::sort. It is the same standard
// library on the same machine, and the algorithm has been unchanged for many
// years, but a change there would show up as a reordered output file rather
// than as wrong data, and the harness would catch it.

#ifndef HICX_BEDTOOLS_OPS_HPP
#define HICX_BEDTOOLS_OPS_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace hicx::bedtools {

// A set of intervals in file order, the minimum bedtools needs.
struct Intervals {
    std::vector<std::string> chrom;
    std::vector<std::int64_t> start;
    std::vector<std::int64_t> end;

    [[nodiscard]] std::size_t size() const { return chrom.size(); }
    void push_back(std::string c, std::int64_t s, std::int64_t e) {
        chrom.push_back(std::move(c));
        start.push_back(s);
        end.push_back(e);
    }
};

// `bedtools sort`: the row order of the output, as indices into the input.
// Chromosomes come out in std::map order, that is lexicographic by byte, and
// the records of a chromosome come out in the order the start-only comparator
// and std::sort leave them in.
[[nodiscard]] std::vector<std::size_t> sort_order(const std::vector<std::string>& chrom,
                                                  const std::vector<std::int64_t>& start);

// `bedtools merge` with the default -d 0, which merges overlapping and
// book-ended intervals. The input must already be sorted, as bedtools requires;
// an unsorted input throws rather than producing a silently wrong answer,
// because bedtools aborts with "you have a problem with your input file".
[[nodiscard]] Intervals merge(const Intervals& input);

// `bedtools intersect -a A -b B -c`: for every interval of A, the number of
// intervals of B that overlap it by at least one base. The result is in A's
// order, which is what pybedtools' to_dataframe then indexes by.
//
// B is expected to be the output of merge(), that is disjoint and sorted per
// chromosome, which is how both call sites in hicValidateLocations use it;
// that case is answered with two binary searches per query. A B that is not
// disjoint falls back to a linear scan of the chromosome, which is correct but
// quadratic, and is never taken by the ported tools.
[[nodiscard]] std::vector<std::int64_t> intersect_count(const Intervals& a,
                                                        const Intervals& b);

}  // namespace hicx::bedtools

#endif  // HICX_BEDTOOLS_OPS_HPP
