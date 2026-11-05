// The last line of error reporting, linked into every tool executable.
//
// Each tool reports the errors it expects (bad arguments, unreadable input)
// with a message and exit status 1. What escapes that, an exception thrown
// outside main's own try block or in a worker thread, or a crash such as a
// segmentation fault, would otherwise end the process with no explanation:
// std::terminate aborts silently under libc++ and with a bare
// "terminate called after ..." under libstdc++, and a signal leaves only the
// shell's exit status. Here every such end prints the tool name and the reason
// on stderr. An escaped exception exits with status 1; a crash keeps its
// signal, so the exit status still says which one it was (128 + signal).
//
// The handlers are installed by a static object before main runs, so the
// tools need no change to have them.

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <string>

#include <unistd.h>

#ifdef __APPLE__
#include <stdlib.h>  // getprogname
#else
#include <errno.h>  // program_invocation_short_name
#endif

namespace {

const char* program_name() {
#ifdef __APPLE__
    return getprogname();
#else
    return program_invocation_short_name;
#endif
}

const char* const kReport =
    "Please report this at https://github.com/deeptools/HiCExplorer/issues "
    "with the command line and, if possible, the input files.\n";

[[noreturn]] void on_terminate() {
    std::string reason = "the tool stopped on an unexpected error";
    const char* report = kReport;
    if (const std::exception_ptr current = std::current_exception()) {
        try {
            std::rethrow_exception(current);
        } catch (const std::bad_alloc&) {
            reason = "out of memory; the input needs more memory than is available";
            report = "";
        } catch (const std::exception& error) {
            reason = error.what();
        } catch (...) {
            reason = "the tool stopped on an unknown exception";
        }
    }
    std::fprintf(stderr, "%s: error: %s\n%s", program_name(), reason.c_str(), report);
    std::fflush(stderr);
    std::_Exit(1);
}

// Only async-signal-safe calls: write, signal, raise.
void write_stderr(const char* text) {
    const std::size_t length = std::strlen(text);
    std::size_t written = 0;
    while (written < length) {
        const ssize_t n = ::write(STDERR_FILENO, text + written, length - written);
        if (n <= 0) {
            return;
        }
        written += static_cast<std::size_t>(n);
    }
}

extern "C" void on_crash(int signal_number) {
    const char* what = "a fatal signal";
    switch (signal_number) {
        case SIGSEGV: what = "a segmentation fault"; break;
        case SIGBUS: what = "a bus error"; break;
        case SIGFPE: what = "a floating-point exception"; break;
        case SIGILL: what = "an illegal instruction (is the binary built for this CPU?)"; break;
        case SIGABRT: what = "an abort"; break;
        default: break;
    }
    write_stderr(program_name());
    write_stderr(": error: the tool crashed with ");
    write_stderr(what);
    write_stderr(".\n");
    write_stderr(kReport);
    std::signal(signal_number, SIG_DFL);
    std::raise(signal_number);
}

struct InstallHandlers {
    InstallHandlers() {
        std::set_terminate(on_terminate);
        for (const int signal_number : {SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT}) {
            std::signal(signal_number, on_crash);
        }
    }
};

const InstallHandlers install_handlers;

}  // namespace
