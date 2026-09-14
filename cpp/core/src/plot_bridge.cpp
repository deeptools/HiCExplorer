#include "hicx/plot_bridge.hpp"

#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"

namespace hicx::plot {

std::string json_string(const std::string& text) {
    std::string out = "\"";
    for (const char ch : text) {
        const auto c = static_cast<unsigned char>(ch);
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            default:
                if (c < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof buffer, "\\u%04x", static_cast<unsigned>(c));
                    out += buffer;
                } else {
                    out += ch;
                }
        }
    }
    out += "\"";
    return out;
}

std::string json_number(double value) {
    if (std::isnan(value)) {
        return "NaN";
    }
    if (std::isinf(value)) {
        return value > 0 ? "Infinity" : "-Infinity";
    }
    return npy::float_repr(value);
}

std::string json_int(std::int64_t value) { return std::to_string(value); }

std::string json_bool(bool value) { return value ? "true" : "false"; }

std::string json_list(const std::vector<std::string>& items) {
    std::string out = "[";
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += items[i];
    }
    out += "]";
    return out;
}

std::string json_numbers(const std::vector<double>& values) {
    std::vector<std::string> items;
    items.reserve(values.size());
    for (const double value : values) {
        items.push_back(json_number(value));
    }
    return json_list(items);
}

std::string json_ints(const std::vector<std::int64_t>& values) {
    std::vector<std::string> items;
    items.reserve(values.size());
    for (const std::int64_t value : values) {
        items.push_back(json_int(value));
    }
    return json_list(items);
}

std::string json_strings(const std::vector<std::string>& values) {
    std::vector<std::string> items;
    items.reserve(values.size());
    for (const std::string& value : values) {
        items.push_back(json_string(value));
    }
    return json_list(items);
}

JsonObject& JsonObject::add(const std::string& key, std::string json_text) {
    fields_.emplace_back(key, std::move(json_text));
    return *this;
}

std::string JsonObject::str() const {
    std::string out = "{";
    for (std::size_t i = 0; i < fields_.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += json_string(fields_[i].first) + ": " + fields_[i].second;
    }
    out += "}";
    return out;
}

namespace {

bool write_all(int fd, const std::string& text) {
    std::size_t done = 0;
    while (done < text.size()) {
        const ssize_t n = ::write(fd, text.data() + done, text.size() - done);
        if (n < 0) {
            if (errno == EINTR) {
                continue;
            }
            return false;
        }
        done += static_cast<std::size_t>(n);
    }
    return true;
}

// <directory of the binary>/../plot, when it holds the package.
void prepend_package_path() {
    char exe[PATH_MAX];
    const ssize_t n = ::readlink("/proc/self/exe", exe, sizeof exe - 1);
    if (n <= 0) {
        return;
    }
    exe[n] = '\0';
    const std::string path(exe);
    const std::string::size_type slash = path.rfind('/');
    if (slash == std::string::npos) {
        return;
    }
    const std::string binary_dir = path.substr(0, slash);
    const std::string::size_type parent_slash = binary_dir.rfind('/');
    const std::string root =
        parent_slash == std::string::npos ? std::string(".") : binary_dir.substr(0, parent_slash);
    const std::string package = root + "/plot";
    struct stat info {};
    if (::stat((package + "/hicexplorer_plot/__main__.py").c_str(), &info) != 0) {
        return;
    }
    const char* old = std::getenv("PYTHONPATH");
    std::string value = package;
    if (old != nullptr && *old != '\0') {
        value += ":";
        value += old;
    }
    ::setenv("PYTHONPATH", value.c_str(), 1);
}

}  // namespace

int draw(const std::string& tool, const std::string& data_json,
         const std::optional<std::string>& data_file) {
    // The peak of the C++ step alone, for the harness's memory gate: VmHWM
    // belongs to the address space, which the exec below replaces, while
    // /usr/bin/time reports the larger of the two processes.
    if (const char* rss_file = std::getenv("HICX_COMPUTE_RSS_FILE");
        rss_file != nullptr && *rss_file != '\0') {
        std::ofstream out(rss_file, std::ios::trunc);
        out << peak_rss_kb() << "\n";
    }
    if (data_file.has_value()) {
        std::ofstream out(*data_file, std::ios::binary | std::ios::trunc);
        out << data_json;
        out.close();
        if (!out) {
            std::fprintf(stderr, "%s: cannot write the plot data to %s\n", tool.c_str(),
                         data_file->c_str());
            return 1;
        }
        return 0;
    }

    const char* tmpdir = std::getenv("TMPDIR");
    std::string pattern = std::string(tmpdir != nullptr && *tmpdir != '\0' ? tmpdir : "/tmp") +
                          "/hicx-plot-XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    const int fd = ::mkstemp(name.data());
    if (fd < 0) {
        std::fprintf(stderr, "%s: cannot create a temporary file for the plot data in %s: %s\n",
                     tool.c_str(), pattern.c_str(), std::strerror(errno));
        return 1;
    }
    const std::string path(name.data());
    const bool written = write_all(fd, data_json);
    ::close(fd);
    if (!written) {
        std::fprintf(stderr, "%s: cannot write the plot data to %s\n", tool.c_str(), path.c_str());
        ::unlink(path.c_str());
        return 1;
    }

    const char* configured = std::getenv("HICX_PLOT_PYTHON");
    const std::string python =
        configured != nullptr && *configured != '\0' ? configured : "python3";
    prepend_package_path();

    std::vector<std::string> args = {python, "-m", "hicexplorer_plot", tool, path, "--remove-data"};
    std::vector<char*> argv;
    argv.reserve(args.size() + 1);
    for (std::string& arg : args) {
        argv.push_back(arg.data());
    }
    argv.push_back(nullptr);
    std::fflush(nullptr);
    ::execvp(python.c_str(), argv.data());

    const int error = errno;
    ::unlink(path.c_str());
    std::fprintf(stderr,
                 "%s: cannot start the drawing interpreter '%s': %s. Set HICX_PLOT_PYTHON to a "
                 "Python with matplotlib 3.8.4 and the hicexplorer_plot package (plot/ in the "
                 "repository).\n",
                 tool.c_str(), python.c_str(), std::strerror(error));
    return 1;
}

}  // namespace hicx::plot
