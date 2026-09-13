// Streaming FASTA reader and IUPAC reverse complement, replacing Bio.SeqIO
// and Bio.Seq for hicFindRestSite (cpp/PLAN.md 3.5).
//
// Two details of the Biopython behaviour are reproduced because they are
// observable in hicFindRestSite's output:
//
//   * `record.name` is the header line up to the first whitespace, not the
//     whole description. That string becomes the chromosome column.
//   * `Seq.reverse_complement()` complements through the IUPAC table
//     Bio.Seq._dna_complement_table and leaves every other character alone,
//     then reverses. Since hicFindRestSite hands it a regular expression
//     rather than a sequence (hicFindRestSite.py:102), characters such as '.',
//     '[' and '^' pass through and end up mirrored: the reverse complement of
//     "A[CG]T" is "A]CG[T". That is a defect in the reference, but a pattern
//     with a character class is rare and the behaviour is what it is, so the
//     port reproduces it rather than guessing at an intent.

#ifndef HICX_FASTA_READER_HPP
#define HICX_FASTA_READER_HPP

#include <functional>
#include <string>

// The functions live in hicx::fasta rather than in hicx because
// core/src/build_matrix.cpp already defines a hicx::reverse_complement for the
// restriction-fragment path of hicBuildMatrix. That one is a switch over the
// uppercase IUPAC codes plus the four lowercase bases and leaves the lowercase
// ambiguity codes alone, while Bio.Seq's table complements those too, so the
// two are not interchangeable. Keeping them apart is the additive change; a
// later commit can decide which one hicBuildMatrix should be using.
namespace hicx::fasta {

// Bio.SeqIO.parse(handle, 'fasta'), one record at a time. `visit` receives the
// record name and its sequence with spaces and carriage returns removed, as
// SimpleFastaParser produces it. `gzipped` selects the gzip reader; the caller
// decides, because the Python decides from the file name alone through
// mimetypes.guess_type and not from the file's magic bytes.
void read_fasta(const std::string& path, bool gzipped,
                const std::function<void(const std::string& name,
                                         const std::string& sequence)>& visit);

// mimetypes.guess_type(name)[1] == 'gzip', which is true only for a name
// ending in '.gz'. A '.bz2' name reports 'bzip2', which the Python compares
// against 'gzip' and then opens as plain text.
[[nodiscard]] bool looks_gzipped(const std::string& path);

// str(Seq(text).reverse_complement()).
[[nodiscard]] std::string reverse_complement(const std::string& text);

}  // namespace hicx::fasta

#endif  // HICX_FASTA_READER_HPP
