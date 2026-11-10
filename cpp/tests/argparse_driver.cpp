// A test driver for hicx::cli (core/include/hicx/argparse.hpp). It declares
// the parsers of cpp/scripts/test_argparse_layer.py, which builds the same
// parsers with Python's argparse, parses its argv and prints the resulting
// namespace as JSON (json.dumps(vars(ns), sort_keys=True)). The test runs
// both on the same command lines and compares exit status, stdout and stderr.
//
//   HICX_ARGDRIVER_MODE=main|sub hicx_argparse_driver ARGS...

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/numpy_compat.hpp"

namespace {

namespace cli = hicx::cli;
using hicx::json::Value;

std::string quote(const std::string& text) {
    std::string out = "\"";
    for (const char c : text) {
        if (c == '"' || c == '\\') {
            out += '\\';
        }
        out += c;
    }
    return out + "\"";
}

std::string number(double value) {
    if (std::isnan(value)) {
        return "NaN";
    }
    if (std::isinf(value)) {
        return value > 0 ? "Infinity" : "-Infinity";
    }
    std::string text = hicx::npy::float_repr(value);
    if (text.find_first_of(".eE") == std::string::npos) {
        text += ".0";
    }
    return text;
}

std::string print(const std::map<std::string, std::string>& fields) {
    std::string out = "{";
    bool first = true;
    for (const auto& [key, value] : fields) {
        out += (first ? "" : ", ") + quote(key) + ": " + value;
        first = false;
    }
    return out + "}";
}

std::string opt_int(const cli::Namespace& ns, const std::string& dest) {
    const auto v = ns.opt_integer(dest);
    return v ? std::to_string(*v) : "null";
}
std::string opt_float(const cli::Namespace& ns, const std::string& dest) {
    const auto v = ns.opt_real(dest);
    return v ? number(*v) : "null";
}
std::string opt_string(const cli::Namespace& ns, const std::string& dest) {
    const auto v = ns.opt_str(dest);
    return v ? quote(*v) : "null";
}
std::string boolean(const cli::Namespace& ns, const std::string& dest) {
    return ns.flag(dest) ? "true" : "false";
}
std::string list(const cli::Namespace& ns, const std::string& dest, const std::string& type) {
    if (!ns.given(dest) && ns.default_of(dest).is_null()) {
        return "null";
    }
    std::string out = "[";
    const std::vector<std::string> values = ns.strs(dest);
    for (std::size_t i = 0; i < values.size(); ++i) {
        out += i ? ", " : "";
        if (type == "int") {
            std::int64_t v = 0;
            (void)cli::python_int(values[i], &v);
            out += std::to_string(v);
        } else if (type == "float") {
            double v = 0;
            (void)cli::python_float(values[i], &v);
            out += number(v);
        } else {
            out += quote(values[i]);
        }
    }
    return out + "]";
}

int main_mode(int argc, char** argv) {
    cli::Parser parser("argdriver", "A parser with every feature the tools use.");
    parser.set_usage("usage: argdriver [options]\n").set_help("\nhelp text\n").set_version_string("1.0");
    cli::ArgumentGroup& g = parser.group("Options");
    g.add({"--int", "-i"}).type("int").default_value(3);
    g.add({"--float"}).type("float").default_value(0.5);
    g.add({"--name"}).choices({"a", "b"}).default_value("a");
    g.add({"--list"}).type("int").nargs("+");
    g.add({"--star"}).nargs("*");
    g.add({"--pair"}).type("float").nargs(2);
    g.add({"--single"}).nargs(1);
    g.add({"--opt"}).nargs("?").const_value(Value::string("C")).default_value("D");
    g.add({"--flag"}).action(cli::Action::StoreTrue);
    g.add({"--no-thing"}).action(cli::Action::StoreFalse);
    g.add({"--append"}).action(cli::Action::Append);
    g.add({"--choice-int"}).type("int").choices({Value::integer(1), Value::integer(2)});
    g.add({"--required-when"}).default_value("x");
    cli::MutuallyExclusiveGroup& m = parser.mutually_exclusive(g);
    m.add({"--left"}).action(cli::Action::StoreTrue);
    m.add({"--right"}).action(cli::Action::StoreTrue);
    g.add({"pos"}).nargs("?");
    g.add({"--help", "-h"}).action(cli::Action::Help);
    g.add({"--version", "-v"}).version("%(prog)s 1.0");
    const cli::Namespace ns = parser.parse(argc, argv);
    std::map<std::string, std::string> f;
    f["int"] = opt_int(ns, "int");
    f["float"] = opt_float(ns, "float");
    f["name"] = opt_string(ns, "name");
    f["list"] = list(ns, "list", "int");
    f["star"] = list(ns, "star", "str");
    f["pair"] = list(ns, "pair", "float");
    f["single"] = list(ns, "single", "str");
    f["opt"] = opt_string(ns, "opt");
    f["flag"] = boolean(ns, "flag");
    f["no_thing"] = boolean(ns, "no_thing");
    f["append"] = list(ns, "append", "str");
    f["choice_int"] = opt_int(ns, "choice_int");
    f["required_when"] = opt_string(ns, "required_when");
    f["left"] = boolean(ns, "left");
    f["right"] = boolean(ns, "right");
    f["pos"] = opt_string(ns, "pos");
    std::printf("%s\n", print(f).c_str());
    return 0;
}

int sub_mode(int argc, char** argv) {
    cli::Parser parser("argdriver", "Subcommands.");
    parser.set_usage("usage: argdriver {run,show} ...\n").set_help("\nhelp text\n");
    cli::ArgumentGroup& top = parser.group("Options");
    top.add({"--verbose"}).action(cli::Action::StoreTrue);
    parser.subcommands("command", true);
    cli::Parser& run = parser.add_subcommand("run", "Run.");
    run.set_usage("usage: argdriver run [options]\n").set_help("\nrun help\n");
    run.set_prog("argdriver {run,show} ... run");
    cli::ArgumentGroup& r = run.group("Run options");
    r.add({"--x"}).type("int").required();
    r.add({"--mode"}).choices({"fast", "slow"}).default_value("fast");
    r.add({"--help", "-h"}).action(cli::Action::Help);
    r.add({"--version"}).version("%(prog)s 2.0");
    cli::Parser& show = parser.add_subcommand("show", "Show.");
    show.set_usage("usage: argdriver show [options]\n").set_help("\nshow help\n");
    show.set_prog("argdriver {run,show} ... show");
    cli::ArgumentGroup& s = show.group("Show options");
    s.add({"--y"}).action(cli::Action::StoreTrue);
    const cli::Namespace ns = parser.parse(argc, argv);
    std::map<std::string, std::string> f;
    f["verbose"] = boolean(ns, "verbose");
    f["command"] = quote(ns.command());
    if (ns.command() == "run") {
        f["x"] = opt_int(ns, "x");
        f["mode"] = opt_string(ns, "mode");
    } else {
        f["y"] = boolean(ns, "y");
    }
    std::printf("%s\n", print(f).c_str());
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    const char* mode = std::getenv("HICX_ARGDRIVER_MODE");
    if (mode != nullptr && std::string(mode) == "sub") {
        return sub_mode(argc, argv);
    }
    return main_mode(argc, argv);
}
