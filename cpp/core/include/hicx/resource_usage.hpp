// Self reported resource usage.
//
// Reduced memory use is a design goal of version 4, so every tool must be able
// to say what it actually used. The equivalence harness measures both
// implementations from the outside with /usr/bin/time, which is the number to
// compare across implementations; this is the cheap in-process check that does
// not need a wrapper.
//
// The report goes to stderr only when HICX_REPORT_RSS is set to a non empty
// value other than "0", so that the tools stay byte identical to the Python
// implementation by default.

#ifndef HICX_RESOURCE_USAGE_HPP
#define HICX_RESOURCE_USAGE_HPP

#include <cstdint>
#include <string>

namespace hicx {

// Peak resident set size in kilobytes, from /proc/self/status VmHWM with
// getrusage as the fallback. Returns 0 when neither is available.
[[nodiscard]] std::int64_t peak_rss_kb();

// Elapsed wall clock seconds since the process started running hicx code.
[[nodiscard]] double elapsed_seconds();

// Writes "tool: peak RSS 123.4 MB, 0.567 s" to stderr when HICX_REPORT_RSS is
// enabled, otherwise does nothing.
void report_resource_usage(const std::string& tool_name);

}  // namespace hicx

#endif  // HICX_RESOURCE_USAGE_HPP
