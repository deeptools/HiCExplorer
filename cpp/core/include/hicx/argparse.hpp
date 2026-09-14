// The argument layer every C++ tool parses its command line with, and the
// single source of the tool's machine-readable specification (--help-json,
// cpp/PLAN.md tier 10, item 10.1).
//
// A tool declares its arguments the way the Python tool declares them for
// argparse: a parser with a description, argument groups, arguments with
// flags, dest, type, nargs, choices, default, required and action, mutually
// exclusive groups, and subcommands. On top of argparse an argument carries
// its file role (input or output, accepted formats) and whether it exists
// only in the C++ port. parse() reproduces the argparse behaviour the tools
// rely on:
//
//  * long options by their full name, `--name=value`, and unique prefixes
//    (allow_abbrev); short options, also with the value attached (`-m3`);
//  * a token that looks like a negative number (`-1`, `-.5`) is a value;
//  * nargs None, '?', '*', '+' and N; store, store_true, store_false,
//    store_const, append, help and version actions; repeated options: the
//    last one wins (append collects);
//  * Python's int() and float() for the int and float types, choices checked
//    on the converted value;
//  * positionals, and subcommands (argparse sub-parsers: the first positional
//    token selects the command and every later token belongs to it);
//  * the errors argparse prints, as "usage text" then "prog: error: message"
//    on stderr with exit status 2; --help prints the tool's help text and
//    exits 0.
//
// Help texts stay hand-written per tool (set_help), because several of them
// are compared byte for byte against the Python's argparse rendering.
//
// Every parser also accepts --help-json, which prints the specification (the
// format is documented in spec_json) and exits 0.

#ifndef HICX_ARGPARSE_HPP
#define HICX_ARGPARSE_HPP

#include <cstdint>
#include <deque>
#include <functional>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "hicx/json_lite.hpp"

namespace hicx::cli {

enum class Action { Store, StoreTrue, StoreFalse, StoreConst, Append, Help, Version };

// argparse.SUPPRESS as a default.
inline constexpr const char* kSuppress = "==SUPPRESS==";

struct Argument {
    std::vector<std::string> flags;  // empty for a positional
    std::string dest_name;
    Action action_kind = Action::Store;
    std::string type_name = "str";  // "str", "int", "float", "FileType" or a custom type name
    std::optional<std::string> nargs_value;  // nullopt: argparse's None; "?", "*", "+" or a count
    std::optional<std::vector<json::Value>> choice_values;
    json::Value default_json;  // null: None
    std::optional<json::Value> const_json;
    bool is_required = false;
    std::optional<std::string> metavar_text;
    std::string help_text;
    std::string version_template;  // for Action::Version: "%(prog)s <version>"
    std::string file_role;          // "", "input" or "output"
    std::vector<std::string> file_formats;
    std::string file_kind = "file";  // "file", "directory" or "prefix"
    bool only_cpp = false;
    std::optional<std::string> note_text;
    // A custom type: returns the error text after "argument X: " when the
    // value is invalid (for example "invalid genomicRegion value: 'x'").
    std::function<std::optional<std::string>(const std::string&)> validator;
    std::optional<std::string> file_type_mode;  // FileType('r') opens the file while parsing

    Argument& dest(std::string value);
    Argument& type(std::string value);
    // FileType(mode): with mode "r" a file that cannot be opened is an
    // argparse error, as in Python; other modes are not opened.
    Argument& file_type(std::string mode);
    Argument& nargs(std::string value);
    Argument& nargs(int count);
    Argument& choices(std::vector<json::Value> values);
    Argument& choices(std::initializer_list<const char*> values);
    Argument& default_value(json::Value value);
    Argument& default_value(const char* value);
    Argument& default_value(std::int64_t value);
    Argument& default_value(int value);
    Argument& default_value(double value);
    Argument& default_value(bool value);
    Argument& const_value(json::Value value);
    Argument& required(bool value = true);
    Argument& metavar(std::string value);
    Argument& help(std::string value);
    Argument& action(Action value);
    Argument& version(std::string template_text);
    Argument& input(std::vector<std::string> formats = {}, std::string kind = "file");
    Argument& output(std::vector<std::string> formats = {}, std::string kind = "file");
    Argument& cpp_only(std::string note = {});
    Argument& note(std::string value);
    Argument& check(std::function<std::optional<std::string>(const std::string&)> validator);

    [[nodiscard]] bool positional() const noexcept { return flags.empty(); }
    // argparse's _get_action_name.
    [[nodiscard]] std::string display_name() const;
};

class Parser;

class ArgumentGroup {
  public:
    // flags {"--matrix", "-m"}, or {"name"} for a positional. The dest follows
    // argparse: the first long flag without dashes, else the first short one.
    Argument& add(std::vector<std::string> flags);
    [[nodiscard]] const std::string& title() const noexcept { return title_; }

  private:
    friend class Parser;
    ArgumentGroup(Parser* parser, std::string title) : parser_(parser), title_(std::move(title)) {}
    Parser* parser_;
    std::string title_;
    std::vector<Argument*> arguments_;
};

class MutuallyExclusiveGroup {
  public:
    Argument& add(std::vector<std::string> flags);

  private:
    friend class Parser;
    MutuallyExclusiveGroup(ArgumentGroup* group, bool required) : group_(group), required_(required) {}
    ArgumentGroup* group_;
    bool required_;
    std::vector<Argument*> arguments_;
};

class Namespace {
  public:
    // Whether the option or positional was given on the command line.
    [[nodiscard]] bool given(const std::string& dest) const;
    // The value of a store action (the default when not given). Throws
    // std::logic_error when the dest is unknown or the value is None.
    [[nodiscard]] std::string str(const std::string& dest) const;
    [[nodiscard]] std::optional<std::string> opt_str(const std::string& dest) const;
    [[nodiscard]] std::int64_t integer(const std::string& dest) const;
    [[nodiscard]] std::optional<std::int64_t> opt_integer(const std::string& dest) const;
    [[nodiscard]] double real(const std::string& dest) const;
    [[nodiscard]] std::optional<double> opt_real(const std::string& dest) const;
    // The values of an nargs '+', '*' or N argument, or of an append action
    // (the default list when not given; empty for a None default).
    [[nodiscard]] std::vector<std::string> strs(const std::string& dest) const;
    [[nodiscard]] std::vector<std::int64_t> integers(const std::string& dest) const;
    [[nodiscard]] std::vector<double> reals(const std::string& dest) const;
    // store_true / store_false / store_const presence semantics: the stored
    // boolean (the default when not given).
    [[nodiscard]] bool flag(const std::string& dest) const;
    // The raw default of a dest, as declared.
    [[nodiscard]] const json::Value& default_of(const std::string& dest) const;
    // The subcommand that was selected, "" when the parser has none.
    [[nodiscard]] const std::string& command() const noexcept { return command_; }

  private:
    friend class Parser;
    struct Entry {
        const Argument* argument = nullptr;
        bool given = false;
        bool flag_value = false;
        std::vector<std::string> values;  // raw tokens (store/append)
    };
    const Entry& entry(const std::string& dest) const;
    std::map<std::string, Entry> entries_;
    std::string command_;
};

class Parser {
  public:
    // prog: the tool name ("hicInfo"); for a subcommand parser the caller
    // passes nothing, the name becomes "<tool> <command>".
    Parser(std::string prog, std::string description);
    Parser(const Parser&) = delete;
    Parser& operator=(const Parser&) = delete;

    ArgumentGroup& group(std::string title);
    MutuallyExclusiveGroup& mutually_exclusive(ArgumentGroup& group, bool required = false);

    // Verbatim texts: usage goes before every error message and before the
    // help text on --help.
    Parser& set_usage(std::string usage);
    Parser& set_help(std::string help);

    // argparse add_subparsers(dest=..., metavar=...): dest may be "" for
    // argparse's SUPPRESS; metavar nullopt for None.
    Parser& subcommands(std::string dest, bool required, std::optional<std::string> metavar = std::nullopt);
    Parser& add_subcommand(const std::string& name, std::string description);

    // Parses argv[1..]; on --help, --version, --help-json and on errors it
    // prints and exits as argparse does.
    [[nodiscard]] Namespace parse(int argc, char** argv) const;
    [[nodiscard]] Namespace parse(const std::vector<std::string>& tokens) const;

    // The specification (schema "hicexplorer-tool-spec", version 1):
    // {"schema", "schema_version", "tool", "version", "description",
    //  "groups": [{"title", "arguments": [ARG]}],
    //  "mutually_exclusive": [{"required", "dests"}],
    //  "subcommands": null | {"dest", "required", "commands": [{"name",
    //    "description", "groups", "mutually_exclusive"}]}}
    // ARG: {"dest", "flags", "positional", "action", "type", "nargs",
    //  "choices", "default", "required", "metavar", "help", "file",
    //  "cpp_only", "note"}; --help-json itself is an argument with action
    //  "help_json".
    [[nodiscard]] std::string spec_json(const std::string& version) const;

    // The version string --help-json reports (the tools pass hicx::kVersion).
    Parser& set_version_string(std::string version);

    // The prog that error messages name. A sub-parser's prog in argparse is
    // built from the parent's usage ("hicCorrectMatrix correct" without a
    // custom usage; "<the usage text> correct" with one); set it when it is
    // not "<tool> <command>". The spec keeps the tool name.
    Parser& set_prog(std::string prog);

  private:
    struct Subcommands {
        std::string dest;
        bool required = false;
        std::optional<std::string> metavar;
        std::vector<std::pair<std::string, std::unique_ptr<Parser>>> commands;
    };

    [[noreturn]] void error(const std::string& message) const;
    void parse_into(const std::vector<std::string>& tokens, std::size_t begin, Namespace& ns,
                    const Parser& root) const;
    [[nodiscard]] std::string spec_body(const std::string& extra_group) const;

    std::string prog_;
    std::string description_;
    std::string usage_;
    std::string help_;
    std::string version_string_;
    std::string error_prog_;
    std::deque<Argument> arguments_;
    std::deque<ArgumentGroup> groups_;
    std::deque<MutuallyExclusiveGroup> mutex_groups_;
    std::optional<Subcommands> subcommands_;
    friend class ArgumentGroup;
    friend class MutuallyExclusiveGroup;
};

// Python's int(text) and float(text) for a str argument.
[[nodiscard]] bool python_int(const std::string& text, std::int64_t* value);
[[nodiscard]] bool python_float(const std::string& text, double* value);

}  // namespace hicx::cli

#endif  // HICX_ARGPARSE_HPP
