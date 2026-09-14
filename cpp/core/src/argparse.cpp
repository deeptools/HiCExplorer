#include "hicx/argparse.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <set>
#include <stdexcept>

#include "hicx/numpy_compat.hpp"

namespace hicx::cli {

namespace {

std::string strip(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin])) != 0) {
        ++begin;
    }
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }
    return text.substr(begin, end - begin);
}

// Digits with single underscores between them, as Python's numeric literals
// in int() and float() allow. Returns the digits without underscores.
std::optional<std::string> digit_run(const std::string& text) {
    if (text.empty() || text.front() == '_' || text.back() == '_') {
        return std::nullopt;
    }
    std::string out;
    char previous = 0;
    for (const char c : text) {
        if (c == '_') {
            if (previous == '_') {
                return std::nullopt;
            }
        } else if (std::isdigit(static_cast<unsigned char>(c)) != 0) {
            out += c;
        } else {
            return std::nullopt;
        }
        previous = c;
    }
    return out;
}

std::string lower(std::string text) {
    for (char& c : text) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return text;
}

// Python's repr of a str.
std::string py_repr(const std::string& text) {
    const bool has_single = text.find('\'') != std::string::npos;
    const bool has_double = text.find('"') != std::string::npos;
    const char quote = has_single && !has_double ? '"' : '\'';
    std::string out(1, quote);
    for (const char c : text) {
        if (c == quote || c == '\\') {
            out += '\\';
            out += c;
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\t') {
            out += "\\t";
        } else {
            out += c;
        }
    }
    out += quote;
    return out;
}

std::string py_repr(const json::Value& value) {
    switch (value.type()) {
        case json::Type::Null: return "None";
        case json::Type::Bool: return value.as_bool() ? "True" : "False";
        case json::Type::Int: return std::to_string(value.as_int());
        case json::Type::Double: return npy::float_repr(value.as_double());
        case json::Type::String: return py_repr(value.as_string());
        default: return "?";
    }
}

std::string join(const std::vector<std::string>& parts, const std::string& separator) {
    std::string out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        out += (i ? separator : "") + parts[i];
    }
    return out;
}

// argparse's _negative_number_matcher: '^-\d+$|^-\d*\.\d+$'.
bool looks_like_negative_number(const std::string& token) {
    if (token.size() < 2 || token[0] != '-') {
        return false;
    }
    std::size_t i = 1;
    std::size_t digits_before = 0;
    while (i < token.size() && std::isdigit(static_cast<unsigned char>(token[i])) != 0) {
        ++i;
        ++digits_before;
    }
    if (i == token.size()) {
        return digits_before > 0;
    }
    if (token[i] != '.') {
        return false;
    }
    ++i;
    std::size_t digits_after = 0;
    while (i < token.size() && std::isdigit(static_cast<unsigned char>(token[i])) != 0) {
        ++i;
        ++digits_after;
    }
    return i == token.size() && digits_after > 0;
}

// ------------------------------------------------------------------ JSON

void json_string(std::string& out, const std::string& text) {
    out += '"';
    for (const unsigned char c : text) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x", c);
                    out += buffer;
                } else {
                    out += static_cast<char>(c);
                }
        }
    }
    out += '"';
}

void json_value(std::string& out, const json::Value& value) {
    switch (value.type()) {
        case json::Type::Null: out += "null"; break;
        case json::Type::Bool: out += value.as_bool() ? "true" : "false"; break;
        case json::Type::Int: out += std::to_string(value.as_int()); break;
        case json::Type::Double: {
            const double d = value.as_double();
            if (std::isfinite(d)) {
                std::string text = npy::float_repr(d);
                if (text.find_first_of(".eE") == std::string::npos) {
                    text += ".0";
                }
                out += text;
            } else {
                // JSON has no inf/nan; Python's json module writes these.
                out += std::isnan(d) ? "NaN" : (d > 0 ? "Infinity" : "-Infinity");
            }
            break;
        }
        case json::Type::String: json_string(out, value.as_string()); break;
        case json::Type::Array: {
            out += '[';
            bool first = true;
            for (const json::Value& item : value.as_array()) {
                if (!first) {
                    out += ", ";
                }
                first = false;
                json_value(out, item);
            }
            out += ']';
            break;
        }
        case json::Type::Object: {
            out += '{';
            bool first = true;
            for (const auto& [key, item] : value.as_object()) {
                if (!first) {
                    out += ", ";
                }
                first = false;
                json_string(out, key);
                out += ": ";
                json_value(out, item);
            }
            out += '}';
            break;
        }
    }
}

const char* action_name(Action action) {
    switch (action) {
        case Action::Store: return "store";
        case Action::StoreTrue: return "store_true";
        case Action::StoreFalse: return "store_false";
        case Action::StoreConst: return "store_const";
        case Action::Append: return "append";
        case Action::Help: return "help";
        case Action::Version: return "version";
    }
    return "store";
}

bool takes_values(Action action) {
    return action == Action::Store || action == Action::Append;
}

std::string argument_json(const Argument& a) {
    std::string out = "{\"dest\": ";
    json_string(out, a.dest_name);
    out += ", \"flags\": [";
    for (std::size_t i = 0; i < a.flags.size(); ++i) {
        if (i) {
            out += ", ";
        }
        json_string(out, a.flags[i]);
    }
    out += "], \"positional\": ";
    out += a.positional() ? "true" : "false";
    out += ", \"action\": ";
    json_string(out, action_name(a.action_kind));
    out += ", \"type\": ";
    if (takes_values(a.action_kind)) {
        json_string(out, a.type_name);
    } else {
        out += "null";
    }
    out += ", \"nargs\": ";
    if (!takes_values(a.action_kind) || !a.nargs_value.has_value()) {
        out += "null";
    } else if (*a.nargs_value == "?" || *a.nargs_value == "*" || *a.nargs_value == "+") {
        json_string(out, *a.nargs_value);
    } else {
        out += *a.nargs_value;
    }
    out += ", \"choices\": ";
    if (a.choice_values.has_value()) {
        json_value(out, json::Value::array(*a.choice_values));
    } else {
        out += "null";
    }
    out += ", \"default\": ";
    json_value(out, a.default_json);
    out += ", \"required\": ";
    out += a.is_required ? "true" : "false";
    out += ", \"metavar\": ";
    if (a.metavar_text.has_value()) {
        json_string(out, *a.metavar_text);
    } else {
        out += "null";
    }
    out += ", \"help\": ";
    json_string(out, a.help_text);
    out += ", \"file\": ";
    if (a.file_role.empty()) {
        out += "null";
    } else {
        out += "{\"role\": ";
        json_string(out, a.file_role);
        out += ", \"formats\": [";
        for (std::size_t i = 0; i < a.file_formats.size(); ++i) {
            if (i) {
                out += ", ";
            }
            json_string(out, a.file_formats[i]);
        }
        out += "], \"kind\": ";
        json_string(out, a.file_kind);
        out += "}";
    }
    out += ", \"cpp_only\": ";
    out += a.only_cpp ? "true" : "false";
    out += ", \"note\": ";
    if (a.note_text.has_value()) {
        json_string(out, *a.note_text);
    } else {
        out += "null";
    }
    out += "}";
    return out;
}

std::string dest_from_flags(const std::vector<std::string>& flags) {
    if (flags.size() == 1 && flags[0][0] != '-') {
        return flags[0];
    }
    std::string chosen;
    for (const std::string& flag : flags) {
        if (flag.rfind("--", 0) == 0) {
            chosen = flag.substr(2);
            break;
        }
    }
    if (chosen.empty() && !flags.empty()) {
        chosen = flags[0].substr(flags[0].find_first_not_of('-'));
    }
    std::replace(chosen.begin(), chosen.end(), '-', '_');
    return chosen;
}

}  // namespace

// --------------------------------------------------------------- numbers

bool python_int(const std::string& text, std::int64_t* value) {
    std::string t = strip(text);
    if (t.empty()) {
        return false;
    }
    bool negative = false;
    if (t[0] == '+' || t[0] == '-') {
        negative = t[0] == '-';
        t = t.substr(1);
    }
    const std::optional<std::string> digits = digit_run(t);
    if (!digits.has_value() || digits->empty()) {
        return false;
    }
    errno = 0;
    const long long parsed = std::strtoll(((negative ? "-" : "") + *digits).c_str(), nullptr, 10);
    if (errno == ERANGE) {
        return false;
    }
    *value = parsed;
    return true;
}

bool python_float(const std::string& text, double* value) {
    std::string t = strip(text);
    if (t.empty()) {
        return false;
    }
    std::string sign;
    if (t[0] == '+' || t[0] == '-') {
        sign = t.substr(0, 1);
        t = t.substr(1);
    }
    const std::string l = lower(t);
    if (l == "inf" || l == "infinity") {
        *value = sign == "-" ? -HUGE_VAL : HUGE_VAL;
        return true;
    }
    if (l == "nan") {
        *value = std::nan("");
        return true;
    }
    // digits [. digits] [e [sign] digits], or . digits
    std::string mantissa = t;
    std::string exponent;
    const std::size_t e = l.find('e');
    if (e != std::string::npos) {
        mantissa = t.substr(0, e);
        exponent = t.substr(e + 1);
        if (exponent.empty()) {
            return false;
        }
        std::string exponent_sign;
        if (exponent[0] == '+' || exponent[0] == '-') {
            exponent_sign = exponent.substr(0, 1);
            exponent = exponent.substr(1);
        }
        const std::optional<std::string> d = digit_run(exponent);
        if (!d.has_value() || d->empty()) {
            return false;
        }
        exponent = exponent_sign + *d;
    }
    std::string integer_part = mantissa;
    std::string fraction_part;
    const std::size_t dot = mantissa.find('.');
    bool has_dot = false;
    if (dot != std::string::npos) {
        has_dot = true;
        integer_part = mantissa.substr(0, dot);
        fraction_part = mantissa.substr(dot + 1);
    }
    std::string clean;
    if (!integer_part.empty()) {
        const std::optional<std::string> d = digit_run(integer_part);
        if (!d.has_value()) {
            return false;
        }
        clean += *d;
    }
    if (has_dot) {
        clean += '.';
        if (!fraction_part.empty()) {
            const std::optional<std::string> d = digit_run(fraction_part);
            if (!d.has_value()) {
                return false;
            }
            clean += *d;
        }
    }
    if (integer_part.empty() && fraction_part.empty()) {
        return false;
    }
    if (!exponent.empty()) {
        clean += "e" + exponent;
    }
    *value = std::strtod((sign + clean).c_str(), nullptr);
    return true;
}

// -------------------------------------------------------------- Argument

Argument& Argument::dest(std::string value) {
    dest_name = std::move(value);
    return *this;
}
Argument& Argument::type(std::string value) {
    type_name = std::move(value);
    return *this;
}
Argument& Argument::file_type(std::string mode) {
    type_name = "FileType";
    file_type_mode = std::move(mode);
    return *this;
}
Argument& Argument::nargs(std::string value) {
    nargs_value = std::move(value);
    return *this;
}
Argument& Argument::nargs(int count) {
    nargs_value = std::to_string(count);
    return *this;
}
Argument& Argument::choices(std::vector<json::Value> values) {
    choice_values = std::move(values);
    return *this;
}
Argument& Argument::choices(std::initializer_list<const char*> values) {
    std::vector<json::Value> out;
    for (const char* v : values) {
        out.push_back(json::Value::string(v));
    }
    choice_values = std::move(out);
    return *this;
}
Argument& Argument::default_value(json::Value value) {
    default_json = std::move(value);
    return *this;
}
Argument& Argument::default_value(const char* value) {
    default_json = json::Value::string(value);
    return *this;
}
Argument& Argument::default_value(std::int64_t value) {
    default_json = json::Value::integer(value);
    return *this;
}
Argument& Argument::default_value(int value) {
    default_json = json::Value::integer(value);
    return *this;
}
Argument& Argument::default_value(double value) {
    default_json = json::Value::number(value);
    return *this;
}
Argument& Argument::default_value(bool value) {
    default_json = json::Value::boolean(value);
    return *this;
}
Argument& Argument::const_value(json::Value value) {
    const_json = std::move(value);
    return *this;
}
Argument& Argument::required(bool value) {
    is_required = value;
    return *this;
}
Argument& Argument::metavar(std::string value) {
    metavar_text = std::move(value);
    return *this;
}
Argument& Argument::help(std::string value) {
    help_text = std::move(value);
    return *this;
}
Argument& Argument::action(Action value) {
    action_kind = value;
    if (value == Action::StoreTrue) {
        default_json = json::Value::boolean(false);
    } else if (value == Action::StoreFalse) {
        default_json = json::Value::boolean(true);
    } else if (value == Action::Help || value == Action::Version) {
        default_json = json::Value::string(kSuppress);
    }
    return *this;
}
Argument& Argument::version(std::string template_text) {
    action(Action::Version);
    version_template = std::move(template_text);
    return *this;
}
Argument& Argument::input(std::vector<std::string> formats, std::string kind) {
    file_role = "input";
    file_formats = std::move(formats);
    file_kind = std::move(kind);
    return *this;
}
Argument& Argument::output(std::vector<std::string> formats, std::string kind) {
    file_role = "output";
    file_formats = std::move(formats);
    file_kind = std::move(kind);
    return *this;
}
Argument& Argument::cpp_only(std::string note) {
    only_cpp = true;
    if (!note.empty()) {
        note_text = std::move(note);
    }
    return *this;
}
Argument& Argument::note(std::string value) {
    note_text = std::move(value);
    return *this;
}
Argument& Argument::check(std::function<std::optional<std::string>(const std::string&)> v) {
    validator = std::move(v);
    return *this;
}

std::string Argument::display_name() const {
    if (!flags.empty()) {
        return join(flags, "/");
    }
    if (metavar_text.has_value() && *metavar_text != kSuppress) {
        return *metavar_text;
    }
    return dest_name;
}

Argument& ArgumentGroup::add(std::vector<std::string> flags) {
    parser_->arguments_.emplace_back();
    Argument& argument = parser_->arguments_.back();
    argument.dest_name = dest_from_flags(flags);
    if (flags.size() == 1 && flags[0][0] != '-') {
        flags.clear();
        argument.is_required = true;  // a positional with nargs None is required
    }
    argument.flags = std::move(flags);
    arguments_.push_back(&argument);
    return argument;
}

Argument& MutuallyExclusiveGroup::add(std::vector<std::string> flags) {
    Argument& argument = group_->add(std::move(flags));
    arguments_.push_back(&argument);
    return argument;
}

// ------------------------------------------------------------- Namespace

const Namespace::Entry& Namespace::entry(const std::string& dest) const {
    const auto it = entries_.find(dest);
    if (it == entries_.end()) {
        throw std::logic_error("no argument with dest '" + dest + "'");
    }
    return it->second;
}

bool Namespace::given(const std::string& dest) const {
    return entry(dest).given;
}

const json::Value& Namespace::default_of(const std::string& dest) const {
    return entry(dest).argument->default_json;
}

std::optional<std::string> Namespace::opt_str(const std::string& dest) const {
    const Entry& e = entry(dest);
    if (e.given) {
        if (e.values.empty()) {
            if (e.argument->const_json.has_value() && e.argument->const_json->is_string()) {
                return e.argument->const_json->as_string();
            }
            return std::nullopt;
        }
        return e.values.back();
    }
    const json::Value& d = e.argument->default_json;
    if (d.is_null()) {
        return std::nullopt;
    }
    if (d.is_string()) {
        return d.as_string();
    }
    if (d.type() == json::Type::Int) {
        return std::to_string(d.as_int());
    }
    if (d.type() == json::Type::Double) {
        return npy::float_repr(d.as_double());
    }
    throw std::logic_error("default of '" + dest + "' is not a scalar");
}

std::string Namespace::str(const std::string& dest) const {
    std::optional<std::string> value = opt_str(dest);
    if (!value.has_value()) {
        throw std::logic_error("argument '" + dest + "' is None");
    }
    return *value;
}

std::optional<std::int64_t> Namespace::opt_integer(const std::string& dest) const {
    const Entry& e = entry(dest);
    if (!e.given) {
        const json::Value& d = e.argument->default_json;
        if (d.is_null()) {
            return std::nullopt;
        }
        if (d.type() == json::Type::Int) {
            return d.as_int();
        }
    }
    const std::optional<std::string> text = opt_str(dest);
    if (!text.has_value()) {
        return std::nullopt;
    }
    std::int64_t value = 0;
    if (!python_int(*text, &value)) {
        throw std::logic_error("argument '" + dest + "' is not an int: " + *text);
    }
    return value;
}

std::int64_t Namespace::integer(const std::string& dest) const {
    const std::optional<std::int64_t> value = opt_integer(dest);
    if (!value.has_value()) {
        throw std::logic_error("argument '" + dest + "' is None");
    }
    return *value;
}

std::optional<double> Namespace::opt_real(const std::string& dest) const {
    const Entry& e = entry(dest);
    if (!e.given) {
        const json::Value& d = e.argument->default_json;
        if (d.is_null()) {
            return std::nullopt;
        }
        if (d.is_number()) {
            return d.as_double();
        }
    }
    const std::optional<std::string> text = opt_str(dest);
    if (!text.has_value()) {
        return std::nullopt;
    }
    double value = 0.0;
    if (!python_float(*text, &value)) {
        throw std::logic_error("argument '" + dest + "' is not a float: " + *text);
    }
    return value;
}

double Namespace::real(const std::string& dest) const {
    const std::optional<double> value = opt_real(dest);
    if (!value.has_value()) {
        throw std::logic_error("argument '" + dest + "' is None");
    }
    return *value;
}

std::vector<std::string> Namespace::strs(const std::string& dest) const {
    const Entry& e = entry(dest);
    if (e.given) {
        return e.values;
    }
    std::vector<std::string> out;
    const json::Value& d = e.argument->default_json;
    if (d.type() == json::Type::Array) {
        for (const json::Value& item : d.as_array()) {
            if (item.is_string()) {
                out.push_back(item.as_string());
            } else if (item.type() == json::Type::Int) {
                out.push_back(std::to_string(item.as_int()));
            } else if (item.type() == json::Type::Double) {
                out.push_back(npy::float_repr(item.as_double()));
            }
        }
    } else if (!d.is_null()) {
        out.push_back(opt_str(dest).value_or(""));
    }
    return out;
}

std::vector<std::int64_t> Namespace::integers(const std::string& dest) const {
    std::vector<std::int64_t> out;
    for (const std::string& text : strs(dest)) {
        std::int64_t value = 0;
        if (!python_int(text, &value)) {
            throw std::logic_error("argument '" + dest + "' holds a non-int: " + text);
        }
        out.push_back(value);
    }
    return out;
}

std::vector<double> Namespace::reals(const std::string& dest) const {
    std::vector<double> out;
    for (const std::string& text : strs(dest)) {
        double value = 0.0;
        if (!python_float(text, &value)) {
            throw std::logic_error("argument '" + dest + "' holds a non-float: " + text);
        }
        out.push_back(value);
    }
    return out;
}

bool Namespace::flag(const std::string& dest) const {
    const Entry& e = entry(dest);
    if (e.given) {
        return e.flag_value;
    }
    const json::Value& d = e.argument->default_json;
    return d.type() == json::Type::Bool ? d.as_bool() : !d.is_null();
}

// ---------------------------------------------------------------- Parser

Parser::Parser(std::string prog, std::string description)
    : prog_(std::move(prog)), description_(std::move(description)) {}

ArgumentGroup& Parser::group(std::string title) {
    groups_.push_back(ArgumentGroup(this, std::move(title)));
    return groups_.back();
}

MutuallyExclusiveGroup& Parser::mutually_exclusive(ArgumentGroup& group, bool required) {
    mutex_groups_.push_back(MutuallyExclusiveGroup(&group, required));
    return mutex_groups_.back();
}

Parser& Parser::set_usage(std::string usage) {
    usage_ = std::move(usage);
    return *this;
}

Parser& Parser::set_help(std::string help) {
    help_ = std::move(help);
    return *this;
}

Parser& Parser::set_version_string(std::string version) {
    version_string_ = std::move(version);
    return *this;
}

Parser& Parser::set_prog(std::string prog) {
    error_prog_ = std::move(prog);
    return *this;
}

Parser& Parser::subcommands(std::string dest, bool required, std::optional<std::string> metavar) {
    subcommands_.emplace();
    subcommands_->dest = std::move(dest);
    subcommands_->required = required;
    subcommands_->metavar = std::move(metavar);
    return *this;
}

Parser& Parser::add_subcommand(const std::string& name, std::string description) {
    if (!subcommands_.has_value()) {
        throw std::logic_error("add_subcommand before subcommands()");
    }
    subcommands_->commands.emplace_back(
        name, std::make_unique<Parser>(prog_ + " " + name, std::move(description)));
    return *subcommands_->commands.back().second;
}

void Parser::error(const std::string& message) const {
    std::fputs(usage_.c_str(), stderr);
    std::fprintf(stderr, "%s: error: %s\n", (error_prog_.empty() ? prog_ : error_prog_).c_str(),
                 message.c_str());
    std::exit(2);
}

Namespace Parser::parse(int argc, char** argv) const {
    std::vector<std::string> tokens;
    for (int i = 1; i < argc; ++i) {
        tokens.emplace_back(argv[i]);
    }
    return parse(tokens);
}

Namespace Parser::parse(const std::vector<std::string>& tokens) const {
    Namespace ns;
    parse_into(tokens, 0, ns, *this);
    return ns;
}

void Parser::parse_into(const std::vector<std::string>& tokens, std::size_t begin, Namespace& ns,
                        const Parser& root) const {
    for (const Argument& a : arguments_) {
        Namespace::Entry e;
        e.argument = &a;
        ns.entries_[a.dest_name] = e;
    }

    // Option lookup (exact, "--x=value", attached short value, unique prefix).
    struct Match {
        const Argument* argument = nullptr;
        std::optional<std::string> explicit_value;
    };
    const auto lookup = [&](const std::string& token) -> std::optional<Match> {
        for (const Argument& a : arguments_) {
            if (std::find(a.flags.begin(), a.flags.end(), token) != a.flags.end()) {
                return Match{&a, std::nullopt};
            }
        }
        if (token == "--help-json") {
            return Match{nullptr, std::nullopt};
        }
        const std::size_t equals = token.find('=');
        if (equals != std::string::npos) {
            const std::string name = token.substr(0, equals);
            for (const Argument& a : arguments_) {
                if (std::find(a.flags.begin(), a.flags.end(), name) != a.flags.end()) {
                    return Match{&a, token.substr(equals + 1)};
                }
            }
        }
        std::vector<Match> matches;
        std::vector<std::string> matched_names;
        if (token.rfind("--", 0) == 0) {
            const std::string prefix = equals != std::string::npos ? token.substr(0, equals) : token;
            for (const Argument& a : arguments_) {
                for (const std::string& flag : a.flags) {
                    if (flag.rfind(prefix, 0) == 0) {
                        matches.push_back(Match{&a, equals != std::string::npos
                                                        ? std::optional<std::string>(token.substr(equals + 1))
                                                        : std::nullopt});
                        matched_names.push_back(flag);
                    }
                }
            }
        } else if (token.size() > 2) {
            const std::string short_prefix = token.substr(0, 2);
            for (const Argument& a : arguments_) {
                for (const std::string& flag : a.flags) {
                    if (flag == short_prefix) {
                        matches.push_back(Match{&a, token.substr(2)});
                        matched_names.push_back(flag);
                    } else if (flag.rfind(token, 0) == 0) {
                        matches.push_back(Match{&a, std::nullopt});
                        matched_names.push_back(flag);
                    }
                }
            }
        }
        if (matches.size() > 1) {
            std::vector<std::string> names = matched_names;
            // argparse 3.12 names the last string of the argument list here,
            // not the ambiguous one (its loop variable leaks into the message).
            error("ambiguous option: " + tokens.back() + " could match " + join(names, ", "));
        }
        if (matches.size() == 1) {
            return matches[0];
        }
        return std::nullopt;
    };
    const auto is_option_like = [&](const std::string& token) {
        if (token.empty() || token[0] != '-' || token == "-") {
            return false;
        }
        for (const Argument& a : arguments_) {
            if (std::find(a.flags.begin(), a.flags.end(), token) != a.flags.end()) {
                return true;
            }
        }
        if (token.size() == 1) {
            return false;
        }
        if (token.find('=') != std::string::npos) {
            const std::string name = token.substr(0, token.find('='));
            for (const Argument& a : arguments_) {
                if (std::find(a.flags.begin(), a.flags.end(), name) != a.flags.end()) {
                    return true;
                }
            }
        }
        if (looks_like_negative_number(token)) {
            return false;
        }
        if (token.find(' ') != std::string::npos) {
            return false;
        }
        return true;
    };

    std::map<const Argument*, const MutuallyExclusiveGroup*> mutex_of;
    for (const MutuallyExclusiveGroup& g : mutex_groups_) {
        for (const Argument* a : g.arguments_) {
            mutex_of[a] = &g;
        }
    }
    std::map<const MutuallyExclusiveGroup*, const Argument*> mutex_seen;

    const auto convert = [&](const Argument& a, const std::string& value) {
        std::optional<json::Value> converted;
        if (a.type_name == "int") {
            std::int64_t v = 0;
            if (!python_int(value, &v)) {
                error("argument " + a.display_name() + ": invalid int value: " + py_repr(value));
            }
            converted = json::Value::integer(v);
        } else if (a.type_name == "float") {
            double v = 0.0;
            if (!python_float(value, &v)) {
                error("argument " + a.display_name() + ": invalid float value: " + py_repr(value));
            }
            converted = json::Value::number(v);
        } else {
            if (a.validator) {
                if (std::optional<std::string> problem = a.validator(value)) {
                    error("argument " + a.display_name() + ": " + *problem);
                }
            }
            if (a.file_type_mode.has_value() && a.file_type_mode->find('r') != std::string::npos &&
                value != "-") {
                std::ifstream probe(value);
                if (!probe) {
                    const int code = errno;
                    error("argument " + a.display_name() + ": can't open " + py_repr(value) + ": [Errno " +
                          std::to_string(code) + "] " + std::strerror(code) + ": " + py_repr(value));
                }
            }
            converted = json::Value::string(value);
        }
        if (a.choice_values.has_value()) {
            bool found = false;
            for (const json::Value& choice : *a.choice_values) {
                if (choice.type() == converted->type() ||
                    (choice.is_number() && converted->is_number())) {
                    if (choice.is_string() ? choice.as_string() == converted->as_string()
                                           : choice.as_double() == converted->as_double()) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                std::vector<std::string> names;
                for (const json::Value& choice : *a.choice_values) {
                    names.push_back(py_repr(choice));
                }
                error("argument " + a.display_name() + ": invalid choice: " + py_repr(*converted) +
                      " (choose from " + join(names, ", ") + ")");
            }
        }
    };

    const auto mark_seen = [&](const Argument& a) {
        const auto it = mutex_of.find(&a);
        if (it != mutex_of.end()) {
            const auto seen = mutex_seen.find(it->second);
            if (seen != mutex_seen.end() && seen->second != &a) {
                error("argument " + a.display_name() + ": not allowed with argument " +
                      seen->second->display_name());
            }
            mutex_seen[it->second] = &a;
        }
    };

    std::vector<std::string> loose;  // tokens for positionals, in order
    std::vector<std::string> unrecognized;
    bool only_positionals = false;
    std::size_t i = begin;
    const std::vector<const Argument*> positionals = [&] {
        std::vector<const Argument*> out;
        for (const Argument& a : arguments_) {
            if (a.positional()) {
                out.push_back(&a);
            }
        }
        return out;
    }();

    // Positional tokens following an option that takes a variable number of
    // values belong to it; "consumable" tells whether the next token is a value.
    const auto consumable = [&](std::size_t k) {
        return k < tokens.size() && (only_positionals || !is_option_like(tokens[k]) ) && tokens[k] != "--";
    };

    // argparse classifies every string before consuming any, so an ambiguous
    // abbreviation is reported before value errors that precede it.
    for (std::size_t k = begin; k < tokens.size() && tokens[k] != "--"; ++k) {
        if (is_option_like(tokens[k])) {
            (void)lookup(tokens[k]);
        } else if (subcommands_.has_value()) {
            break;
        }
    }

    while (i < tokens.size()) {
        const std::string& token = tokens[i];
        if (!only_positionals && token == "--") {
            only_positionals = true;
            ++i;
            continue;
        }
        if (only_positionals || !is_option_like(token)) {
            if (subcommands_.has_value()) {
                const auto& commands = subcommands_->commands;
                const auto it = std::find_if(commands.begin(), commands.end(),
                                             [&](const auto& c) { return c.first == token; });
                if (it == commands.end()) {
                    std::vector<std::string> names;
                    for (const auto& c : commands) {
                        names.push_back(py_repr(c.first));
                    }
                    const std::string name = subcommands_->metavar.has_value()
                                                 ? *subcommands_->metavar
                                                 : subcommands_->dest;
                    error("argument " + name + ": invalid choice: " + py_repr(token) + " (choose from " +
                          join(names, ", ") + ")");
                }
                ns.command_ = token;
                it->second->parse_into(tokens, i + 1, ns, root);
                i = tokens.size();
                break;
            }
            loose.push_back(token);
            ++i;
            continue;
        }
        const std::optional<Match> match = lookup(token);
        if (!match.has_value()) {
            unrecognized.push_back(token);
            ++i;
            continue;
        }
        if (match->argument == nullptr) {
            std::string spec = root.spec_json(root.version_string_);
            std::fputs(spec.c_str(), stdout);
            std::fputs("\n", stdout);
            std::exit(0);
        }
        const Argument& a = *match->argument;
        Namespace::Entry& e = ns.entries_[a.dest_name];
        ++i;
        switch (a.action_kind) {
            case Action::Help:
                std::fputs(usage_.c_str(), stdout);
                std::fputs(help_.c_str(), stdout);
                std::exit(0);
            case Action::Version: {
                std::string text = a.version_template;
                // %(prog)s is the prog of the parser that holds the action:
                // "tool command" inside a subcommand, as argparse prints it.
                const std::size_t at = text.find("%(prog)s");
                if (at != std::string::npos) {
                    text.replace(at, 8, error_prog_.empty() ? prog_ : error_prog_);
                }
                std::fputs((text + "\n").c_str(), stdout);
                std::exit(0);
            }
            case Action::StoreTrue:
            case Action::StoreFalse:
            case Action::StoreConst:
                if (match->explicit_value.has_value()) {
                    error("argument " + a.display_name() + ": ignored explicit argument " +
                          py_repr(*match->explicit_value));
                }
                mark_seen(a);
                e.given = true;
                e.flag_value = a.action_kind != Action::StoreFalse;
                e.values.clear();
                if (a.action_kind == Action::StoreConst && a.const_json.has_value() &&
                    a.const_json->is_string()) {
                    e.values.push_back(a.const_json->as_string());
                }
                break;
            case Action::Store:
            case Action::Append: {
                std::vector<std::string> values;
                const std::string n = a.nargs_value.value_or("");
                if (match->explicit_value.has_value()) {
                    if (n == "?" || n.empty() || n == "*" || n == "+" || n == "1") {
                        values.push_back(*match->explicit_value);
                    } else {
                        error("argument " + a.display_name() + ": expected " + n + " arguments");
                    }
                } else if (n.empty()) {
                    if (!consumable(i)) {
                        error("argument " + a.display_name() + ": expected one argument");
                    }
                    values.push_back(tokens[i++]);
                } else if (n == "?") {
                    if (consumable(i)) {
                        values.push_back(tokens[i++]);
                    }
                } else if (n == "*" || n == "+") {
                    while (consumable(i)) {
                        values.push_back(tokens[i++]);
                    }
                    if (n == "+" && values.empty()) {
                        error("argument " + a.display_name() + ": expected at least one argument");
                    }
                } else {
                    // Python 3.12.7 matches a fixed count of values with
                    // [AO]{N}: the next N strings are taken even when they
                    // look like options; only "--" stops them.
                    const int count = std::stoi(n);
                    for (int k = 0; k < count; ++k) {
                        if (i >= tokens.size() || (!only_positionals && tokens[i] == "--")) {
                            error("argument " + a.display_name() + ": expected " + n + " argument" +
                                  (count == 1 ? "" : "s"));
                        }
                        values.push_back(tokens[i++]);
                    }
                }
                for (const std::string& v : values) {
                    convert(a, v);
                }
                mark_seen(a);
                e.given = true;
                if (a.action_kind == Action::Append) {
                    e.values.insert(e.values.end(), values.begin(), values.end());
                } else {
                    e.values = std::move(values);
                }
                break;
            }
        }
    }

    // Positionals, in declaration order.
    std::size_t next = 0;
    for (const Argument* p : positionals) {
        Namespace::Entry& e = ns.entries_[p->dest_name];
        const std::string n = p->nargs_value.value_or("");
        std::vector<std::string> values;
        if (n.empty()) {
            if (next < loose.size()) {
                values.push_back(loose[next++]);
            }
        } else if (n == "?") {
            if (next < loose.size()) {
                values.push_back(loose[next++]);
            }
        } else if (n == "*" || n == "+") {
            while (next < loose.size()) {
                values.push_back(loose[next++]);
            }
        } else {
            const std::size_t count = static_cast<std::size_t>(std::stoi(n));
            for (std::size_t k = 0; k < count && next < loose.size(); ++k) {
                values.push_back(loose[next++]);
            }
        }
        if (!values.empty()) {
            for (const std::string& v : values) {
                convert(*p, v);
            }
            e.given = true;
            e.values = std::move(values);
        }
    }
    while (next < loose.size()) {
        unrecognized.push_back(loose[next++]);
    }

    std::vector<std::string> missing;
    for (const Argument& a : arguments_) {
        const Namespace::Entry& e = ns.entries_[a.dest_name];
        bool satisfied = e.given;
        if (a.positional() && !satisfied) {
            const std::string n = a.nargs_value.value_or("");
            satisfied = n == "?" || n == "*";
        }
        if ((a.is_required || a.positional()) && !satisfied) {
            missing.push_back(a.positional() && a.metavar_text.has_value() ? *a.metavar_text : a.display_name());
        }
    }
    if (subcommands_.has_value() && subcommands_->required && ns.command_.empty()) {
        missing.push_back(subcommands_->metavar.has_value() ? *subcommands_->metavar : subcommands_->dest);
    }
    if (!missing.empty()) {
        error("the following arguments are required: " + join(missing, ", "));
    }
    for (const MutuallyExclusiveGroup& g : mutex_groups_) {
        if (g.required_ && mutex_seen.find(&g) == mutex_seen.end()) {
            std::vector<std::string> names;
            for (const Argument* a : g.arguments_) {
                names.push_back(a->display_name());
            }
            error("one of the arguments " + join(names, " ") + " is required");
        }
    }
    if (!unrecognized.empty()) {
        if (&root == this || ns.command_.empty()) {
            error("unrecognized arguments: " + join(unrecognized, " "));
        }
        root.error("unrecognized arguments: " + join(unrecognized, " "));
    }
}

std::string Parser::spec_body(const std::string& extra_group) const {
    std::string out = "\"description\": ";
    json_string(out, description_);
    out += ", \"groups\": [";
    bool first_group = true;
    for (const ArgumentGroup& g : groups_) {
        if (!first_group) {
            out += ", ";
        }
        first_group = false;
        out += "{\"title\": ";
        json_string(out, g.title_);
        out += ", \"arguments\": [";
        for (std::size_t k = 0; k < g.arguments_.size(); ++k) {
            if (k) {
                out += ", ";
            }
            out += argument_json(*g.arguments_[k]);
        }
        out += "]}";
    }
    if (!extra_group.empty()) {
        out += (first_group ? "" : ", ") + extra_group;
    }
    out += "], \"mutually_exclusive\": [";
    for (std::size_t k = 0; k < mutex_groups_.size(); ++k) {
        if (k) {
            out += ", ";
        }
        out += "{\"required\": ";
        out += mutex_groups_[k].required_ ? "true" : "false";
        out += ", \"dests\": [";
        for (std::size_t m = 0; m < mutex_groups_[k].arguments_.size(); ++m) {
            if (m) {
                out += ", ";
            }
            json_string(out, mutex_groups_[k].arguments_[m]->dest_name);
        }
        out += "]}";
    }
    out += "]";
    return out;
}

std::string Parser::spec_json(const std::string& version) const {
    std::string out = "{\"schema\": \"hicexplorer-tool-spec\", \"schema_version\": 1, \"tool\": ";
    json_string(out, prog_);
    out += ", \"version\": ";
    json_string(out, version);
    // --help-json is an argument of every top-level parser, listed in a group
    // of its own with action "help_json".
    Argument help_json;
    help_json.flags = {"--help-json"};
    help_json.dest_name = "help_json";
    help_json.action_kind = Action::Help;
    help_json.default_json = json::Value::string(kSuppress);
    help_json.help_text = "Print this tool's argument specification as JSON and exit.";
    help_json.only_cpp = true;
    std::string argument = argument_json(help_json);
    const std::string help_action = "\"action\": \"help\"";
    argument.replace(argument.find(help_action), help_action.size(), "\"action\": \"help_json\"");
    out += ", " + spec_body("{\"title\": \"C++ specification\", \"arguments\": [" + argument + "]}");
    out += ", \"subcommands\": ";
    if (!subcommands_.has_value()) {
        out += "null";
    } else {
        out += "{\"dest\": ";
        json_string(out, subcommands_->dest);
        out += ", \"required\": ";
        out += subcommands_->required ? "true" : "false";
        out += ", \"commands\": [";
        for (std::size_t k = 0; k < subcommands_->commands.size(); ++k) {
            if (k) {
                out += ", ";
            }
            out += "{\"name\": ";
            json_string(out, subcommands_->commands[k].first);
            out += ", " + subcommands_->commands[k].second->spec_body({}) + "}";
        }
        out += "]}";
    }
    out += "}";
    return out;
}

}  // namespace hicx::cli
