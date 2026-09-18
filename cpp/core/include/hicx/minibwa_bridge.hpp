// The external-process bridge to minibwa (cpp/PLAN.md tier 13), used by
// hicBuildIndex and hicAlignReads.
//
// Neither tool reimplements alignment: minibwa is a real, already installed
// conda-packaged aligner (see the tools' own file comments for the exact
// path), and both tools are thin wrappers that build its command line, run
// it and, for hicAlignReads, hand the SAM it wrote to htslib to become a
// BAM. This mirrors the fork/exec/waitpid shape hicx::plot::preflight and
// hicx::plot::draw already use to run the Python drawing interpreter
// (core/src/plot_bridge.cpp): the same reasons apply here (the child's exit
// status must be inspected before the C++ tool can do its own next step,
// so this cannot be a plain execvp replacement the way the Python plotting
// bridge's draw() is for its last step).

#ifndef HICX_MINIBWA_BRIDGE_HPP
#define HICX_MINIBWA_BRIDGE_HPP

#include <optional>
#include <string>
#include <vector>

namespace hicx::minibwa {

// HICX_MINIBWA_BIN, or "minibwa" on PATH when unset (the same pattern as
// HICX_PLOT_PYTHON for the drawing interpreter).
[[nodiscard]] std::string binary_path();

// Runs "<binary> <args...>", waits for it and returns its exit status
// (Python's os.WEXITSTATUS(status), i.e. 127 when the binary itself could
// not be started, matching the shell's own convention and what
// hicx::plot::preflight does on the same failure). Nothing is captured:
// minibwa's own stdout/stderr go to the caller's, so a user sees exactly
// what running it by hand would print. tool_name is used only in the
// wrapper's own error messages ("hicAlignReads: ...").
[[nodiscard]] int run(const std::string& tool_name, const std::vector<std::string>& args);

// Runs "<binary> map <args...> -o <sam_path> <index_prefix> <fastq_path>",
// then converts the SAM file it wrote into a BAM at bam_path with htslib
// (hicx::BamReader / hicx::BamWriter, cpp/core/src/bam_file.cpp), which is
// the same htslib `samtools view -b --no-PG` uses, so the alignment records
// the two routes produce are byte-identical (cpp/PLAN.md 5.1 class E0; see
// hicAlignReads.cpp for the measured result and the one header field, the
// minibwa @PG CL line's own temporary-file path, that differs by
// construction and is not part of that comparison). The temporary SAM file
// is removed afterwards, on success or failure. Returns
// minibwa's own exit status, or a wrapper-specific nonzero status when the
// SAM-to-BAM conversion itself fails (minibwa exited 0 but its output could
// not be read or the BAM could not be written).
[[nodiscard]] int run_map_to_bam(const std::string& tool_name, const std::vector<std::string>& map_args,
                                 const std::string& index_prefix, const std::string& fastq_path,
                                 const std::string& bam_path);

}  // namespace hicx::minibwa

#endif  // HICX_MINIBWA_BRIDGE_HPP
