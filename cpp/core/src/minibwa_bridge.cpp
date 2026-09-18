#include "hicx/minibwa_bridge.hpp"

#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "hicx/bam_file.hpp"

namespace hicx::minibwa {

std::string binary_path() {
    const char* configured = std::getenv("HICX_MINIBWA_BIN");
    return configured != nullptr && *configured != '\0' ? std::string(configured) : std::string("minibwa");
}

int run(const std::string& tool_name, const std::vector<std::string>& args) {
    const std::string binary = binary_path();
    std::vector<std::string> full_args = args;
    full_args.insert(full_args.begin(), binary);
    std::vector<char*> argv;
    argv.reserve(full_args.size() + 1);
    for (std::string& arg : full_args) {
        argv.push_back(arg.data());
    }
    argv.push_back(nullptr);

    std::fflush(nullptr);
    const pid_t child = ::fork();
    if (child < 0) {
        std::fprintf(stderr, "%s: cannot start %s: %s\n", tool_name.c_str(), binary.c_str(),
                     std::strerror(errno));
        return 1;
    }
    if (child == 0) {
        ::execvp(binary.c_str(), argv.data());
        std::fprintf(stderr,
                     "%s: cannot start '%s': %s. Set HICX_MINIBWA_BIN to the minibwa executable "
                     "(for example ~/miniconda3/envs/__minibwa@0.5/bin/minibwa).\n",
                     tool_name.c_str(), binary.c_str(), std::strerror(errno));
        ::_exit(127);
    }
    int status = 0;
    while (::waitpid(child, &status, 0) < 0) {
        if (errno != EINTR) {
            std::fprintf(stderr, "%s: %s was lost: %s\n", tool_name.c_str(), binary.c_str(),
                         std::strerror(errno));
            return 1;
        }
    }
    if (WIFEXITED(status)) {
        return WEXITSTATUS(status);
    }
    if (WIFSIGNALED(status)) {
        std::fprintf(stderr, "%s: %s was killed by signal %d\n", tool_name.c_str(), binary.c_str(),
                     WTERMSIG(status));
    }
    return 1;
}

namespace {

// A private temporary path next to the requested BAM output, so the SAM
// intermediate lands on the same filesystem (a same-filesystem unlink is
// cheap and avoids a cross-device rename anywhere downstream that might
// later want to move it).
std::string make_temp_sam_path(const std::string& bam_path) {
    // mkstemp requires "XXXXXX" as the literal last six characters, so the
    // temporary file carries no ".sam" extension; BamReader detects the SAM
    // format from its content, not from the name, so this does not matter.
    std::string pattern = bam_path + ".hicx-minibwa-XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    const int fd = ::mkstemp(name.data());
    if (fd < 0) {
        throw std::runtime_error("cannot create a temporary SAM file next to " + bam_path + ": " +
                                 std::strerror(errno));
    }
    ::close(fd);
    return std::string(name.data());
}

}  // namespace

int run_map_to_bam(const std::string& tool_name, const std::vector<std::string>& map_args,
                   const std::string& index_prefix, const std::string& fastq_path,
                   const std::string& bam_path) {
    std::string sam_path;
    try {
        sam_path = make_temp_sam_path(bam_path);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: %s\n", tool_name.c_str(), error.what());
        return 1;
    }

    std::vector<std::string> args = {"map"};
    args.insert(args.end(), map_args.begin(), map_args.end());
    args.push_back("-o");
    args.push_back(sam_path);
    args.push_back(index_prefix);
    args.push_back(fastq_path);

    const int status = run(tool_name, args);
    if (status != 0) {
        ::unlink(sam_path.c_str());
        return status;
    }

    try {
        hicx::BamReader reader(sam_path);
        hicx::BamWriter writer(bam_path, reader);
        hicx::BamRecord record;
        while (reader.read(record)) {
            writer.write(record);
        }
        writer.close();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: cannot convert minibwa's SAM output to BAM: %s\n", tool_name.c_str(),
                     error.what());
        ::unlink(sam_path.c_str());
        return 1;
    }
    ::unlink(sam_path.c_str());
    return 0;
}

}  // namespace hicx::minibwa
