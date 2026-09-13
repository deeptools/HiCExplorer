// A small argparse emulation shared by the cHi-C tools.
//
// The cHi-C parsers are plain argparse.ArgumentParser objects, and two of
// argparse's habits are part of the command line contract:
//
//   * a long option may be abbreviated to any unique prefix. The Python test
//     suite itself calls chicViewpoint with --backgroundModel for
//     --backgroundModelFile, so this is not hypothetical;
//   * an option with nargs='+' or nargs=2 consumes the following arguments up
//     to the next one that looks like an option.
//
// Errors print the usage and "<prog>: error: <message>" and exit with 2, as
// argparse does. The help text is a summary, not a byte copy of argparse's.
//
// Beyond argparse, the parser records whether an option was given explicitly
// or left at its default, which is what the project rule for figure outputs
// needs: a tool refuses when a figure file was explicitly requested, and only
// skips a figure the user never named.

#ifndef HICX_TOOLS_CHIC_ARGUMENTS_HPP
#define HICX_TOOLS_CHIC_ARGUMENTS_HPP

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "hicx/chic_viewpoint.hpp"

namespace hicx::chic_cli {

enum class Arity { Flag, One, Two, OneOrMore };
enum class Type { String, Int, Float };

struct Option {
    std::vector<std::string> flags;  // long form first, as argparse lists them
    std::string dest;
    Arity arity = Arity::One;
    Type type = Type::String;
    bool required = false;
    std::vector<std::string> defaults;
    std::string help;
};

class Parsed {
  public:
    [[nodiscard]] const std::vector<std::string>& list(const std::string& dest) const {
        return values_.at(dest);
    }
    [[nodiscard]] const std::string& str(const std::string& dest) const {
        return values_.at(dest).at(0);
    }
    [[nodiscard]] std::int64_t integer(const std::string& dest, std::size_t index = 0) const {
        return chic::python_int(values_.at(dest).at(index));
    }
    [[nodiscard]] double real(const std::string& dest) const {
        return chic::python_float(values_.at(dest).at(0));
    }
    [[nodiscard]] bool flag(const std::string& dest) const {
        return !values_.at(dest).empty() && values_.at(dest)[0] == "1";
    }
    [[nodiscard]] bool explicitly_given(const std::string& dest) const {
        return given_.count(dest) > 0;
    }

  private:
    friend class Parser;
    std::map<std::string, std::vector<std::string>> values_;
    std::set<std::string> given_;
};

class Parser {
  public:
    Parser(std::string prog, std::string usage, std::string description, std::string version)
        : prog_(std::move(prog)),
          usage_(std::move(usage)),
          description_(std::move(description)),
          version_(std::move(version)) {}

    void add(Option option) { options_.push_back(std::move(option)); }

    [[nodiscard]] Parsed parse(int argc, char** argv) const {
        Parsed parsed;
        for (const Option& option : options_) {
            if (option.arity == Arity::Flag) {
                parsed.values_[option.dest] = {"0"};
            } else {
                parsed.values_[option.dest] = option.defaults;
            }
        }
        std::vector<std::string> tokens(argv + 1, argv + argc);
        std::vector<std::string> unrecognized;
        for (std::size_t i = 0; i < tokens.size(); ++i) {
            std::string token = tokens[i];
            if (token == "-h" || token == "--help" || (token.size() > 2 && token.rfind("--h", 0) == 0 &&
                                                       std::string("--help").rfind(token, 0) == 0)) {
                print_help();
                std::exit(0);
            }
            if (token == "--version" || (token.size() > 3 && token.rfind("--v", 0) == 0 &&
                                         std::string("--version").rfind(token, 0) == 0)) {
                std::printf("%s %s\n", prog_.c_str(), version_.c_str());
                std::exit(0);
            }
            if (!looks_like_option(token)) {
                unrecognized.push_back(token);
                continue;
            }
            std::string explicit_value;
            bool has_explicit = false;
            const Option* option = match(token, explicit_value, has_explicit);
            if (option == nullptr) {
                unrecognized.push_back(token);
                continue;
            }
            const std::string name = display_name(*option);
            std::vector<std::string> values;
            const auto take_following = [&](std::size_t maximum) {
                while (values.size() < maximum && i + 1 < tokens.size() &&
                       !looks_like_option(tokens[i + 1])) {
                    values.push_back(tokens[++i]);
                }
            };
            switch (option->arity) {
                case Arity::Flag:
                    if (has_explicit) {
                        error("argument " + name + ": ignored explicit argument '" +
                              explicit_value + "'");
                    }
                    values = {"1"};
                    break;
                case Arity::One:
                    if (has_explicit) {
                        values = {explicit_value};
                    } else {
                        take_following(1);
                        if (values.empty()) {
                            error("argument " + name + ": expected one argument");
                        }
                    }
                    break;
                case Arity::Two:
                    if (has_explicit) {
                        error("argument " + name + ": expected 2 arguments");
                    }
                    take_following(2);
                    if (values.size() != 2) {
                        error("argument " + name + ": expected 2 arguments");
                    }
                    break;
                case Arity::OneOrMore:
                    if (has_explicit) {
                        values = {explicit_value};
                    }
                    take_following(static_cast<std::size_t>(-1));
                    if (values.empty()) {
                        error("argument " + name + ": expected at least one argument");
                    }
                    break;
            }
            for (const std::string& value : values) {
                check_type(*option, name, value);
            }
            parsed.values_[option->dest] = values;
            parsed.given_.insert(option->dest);
        }

        std::string missing;
        for (const Option& option : options_) {
            if (option.required && parsed.given_.count(option.dest) == 0) {
                missing += (missing.empty() ? "" : ", ") + display_name(option);
            }
        }
        if (!missing.empty()) {
            error("the following arguments are required: " + missing);
        }
        if (!unrecognized.empty()) {
            std::string joined;
            for (const std::string& token : unrecognized) {
                joined += (joined.empty() ? "" : " ") + token;
            }
            error("unrecognized arguments: " + joined);
        }
        return parsed;
    }

    [[noreturn]] void error(const std::string& message) const {
        std::fprintf(stderr, "%s%s: error: %s\n", usage_.c_str(), prog_.c_str(),
                     message.c_str());
        std::exit(2);
    }

  private:
    // argparse treats "-5" or "-.5" as a value when no option looks like a
    // negative number, which holds for every cHi-C parser.
    static bool looks_like_option(const std::string& token) {
        if (token.size() < 2 || token[0] != '-') {
            return false;
        }
        std::size_t i = 1;
        bool digits = false;
        bool dot = false;
        for (; i < token.size(); ++i) {
            if (token[i] >= '0' && token[i] <= '9') {
                digits = true;
            } else if (token[i] == '.' && !dot) {
                dot = true;
            } else {
                break;
            }
        }
        const bool negative_number = i == token.size() && digits;
        return !negative_number;
    }

    const Option* match(const std::string& token, std::string& explicit_value,
                        bool& has_explicit) const {
        std::string name = token;
        const std::size_t equals = token.find('=');
        if (token.rfind("--", 0) == 0 && equals != std::string::npos) {
            name = token.substr(0, equals);
            explicit_value = token.substr(equals + 1);
            has_explicit = true;
        }
        for (const Option& option : options_) {
            for (const std::string& flag : option.flags) {
                if (flag == name) {
                    return &option;
                }
            }
        }
        if (name.rfind("--", 0) == 0) {
            const Option* found = nullptr;
            std::vector<std::string> candidates;
            for (const Option& option : options_) {
                for (const std::string& flag : option.flags) {
                    if (flag.rfind("--", 0) == 0 && flag.rfind(name, 0) == 0) {
                        if (found != &option) {
                            candidates.push_back(flag);
                        }
                        found = &option;
                    }
                }
            }
            if (candidates.size() > 1) {
                std::string joined;
                for (const std::string& candidate : candidates) {
                    joined += (joined.empty() ? "" : ", ") + candidate;
                }
                error("ambiguous option: " + name + " could match " + joined);
            }
            return found;
        }
        // A single dash: "-mvalue" is -m with an attached value.
        if (name.size() > 2) {
            const std::string short_flag = name.substr(0, 2);
            for (const Option& option : options_) {
                for (const std::string& flag : option.flags) {
                    if (flag == short_flag && option.arity == Arity::One) {
                        explicit_value = name.substr(2);
                        has_explicit = true;
                        return &option;
                    }
                }
            }
        }
        return nullptr;
    }

    static std::string display_name(const Option& option) {
        std::string joined;
        for (const std::string& flag : option.flags) {
            joined += (joined.empty() ? "" : "/") + flag;
        }
        return joined;
    }

    void check_type(const Option& option, const std::string& name,
                    const std::string& value) const {
        try {
            if (option.type == Type::Int) {
                (void)chic::python_int(value);
            } else if (option.type == Type::Float) {
                (void)chic::python_float(value);
            }
        } catch (const chic::ValueError&) {
            error("argument " + name + ": invalid " +
                  std::string(option.type == Type::Int ? "int" : "float") + " value: '" +
                  value + "'");
        }
    }

    void print_help() const {
        std::printf("%s\n%s\n\noptions:\n", usage_.c_str(), description_.c_str());
        for (const Option& option : options_) {
            std::printf("  %s\n", display_name(option).c_str());
            if (!option.help.empty()) {
                std::printf("                        %s\n", option.help.c_str());
            }
        }
        std::printf("  --help, -h            show this help message and exit\n"
                    "  --version             show program's version number and exit\n");
    }

    std::string prog_;
    std::string usage_;
    std::string description_;
    std::string version_;
    std::vector<Option> options_;
};

}  // namespace hicx::chic_cli

#endif  // HICX_TOOLS_CHIC_ARGUMENTS_HPP
