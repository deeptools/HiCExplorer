#include "hicx/resource_usage.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include <sys/resource.h>

namespace hicx {

namespace {

const std::chrono::steady_clock::time_point kStart = std::chrono::steady_clock::now();

bool reporting_enabled() {
    const char* value = std::getenv("HICX_REPORT_RSS");
    return value != nullptr && value[0] != '\0' && std::strcmp(value, "0") != 0;
}

}  // namespace

std::int64_t peak_rss_kb() {
    std::ifstream status("/proc/self/status");
    if (status) {
        std::string line;
        while (std::getline(status, line)) {
            if (line.rfind("VmHWM:", 0) == 0) {
                const std::size_t digits = line.find_first_of("0123456789");
                if (digits != std::string::npos) {
                    return std::strtoll(line.c_str() + digits, nullptr, 10);
                }
            }
        }
    }
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        return static_cast<std::int64_t>(usage.ru_maxrss);  // kilobytes on Linux
    }
    return 0;
}

double elapsed_seconds() {
    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - kStart;
    return elapsed.count();
}

void report_resource_usage(const std::string& tool_name) {
    if (!reporting_enabled()) {
        return;
    }
    std::fprintf(stderr, "%s: peak RSS %.1f MB, %.3f s\n", tool_name.c_str(),
                 static_cast<double>(peak_rss_kb()) / 1024.0, elapsed_seconds());
}

}  // namespace hicx
