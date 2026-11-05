#include "hicx/read_bed.hpp"

#include "hicx/py_text.hpp"

namespace hicx::readbed {

std::string guess_file_type(std::size_t field_count, std::string* assumed) {
    switch (field_count) {
        case 3: return "bed3";
        case 4: return "bedgraph";
        case 6: return "bed6";
        case 12: return "bed12";
        case 9: return "bed9";
        default: break;
    }
    const std::string type = field_count > 6 ? "bed6" : "bed3";
    if (assumed != nullptr) {
        *assumed = type;
    }
    return type;
}

std::vector<std::string> interval_fields(std::string_view file_type) {
    std::vector<std::string> fields = {"chromosome", "start", "end", "name", "score", "strand",
                                       "thick_start", "thick_end", "rgb", "block_count",
                                       "block_sizes", "block_starts"};
    if (file_type == "bed12") {
        return fields;
    }
    fields.resize(file_type == "bed9" ? 9 : 6);
    return fields;
}

std::optional<std::size_t> required_field_count(std::string_view file_type) {
    // `if self.file_handle == 'bed12'` is never true, so bed12 is not checked.
    if (file_type == "bed3") {
        return 3;
    }
    if (file_type == "bed6") {
        return 6;
    }
    return std::nullopt;
}

Column column_kind(std::size_t index) noexcept {
    switch (index) {
        case 0: case 3: return Column::Text;
        case 5: return Column::Strand;
        case 1: case 2: case 6: case 7: case 9: return Column::Integer;
        case 8: return Column::Rgb;
        case 10: case 11: return Column::Blocks;
        default: return Column::Float;
    }
}

std::string normalize_strand(std::string_view value, bool& invalid) {
    invalid = false;
    if (value == "+" || value == "-" || value == ".") {
        return std::string(value);
    }
    if (value == "1") {
        return "+";
    }
    if (value == "-1") {
        return "-";
    }
    invalid = true;
    return ".";
}

bool is_comment_line(std::string_view line) noexcept {
    return line.substr(0, 1) == "#" || line.substr(0, 5) == "track" ||
           line.substr(0, 7) == "browser" || pytext::strip(line).empty();
}

std::vector<std::string> split_line(std::string_view line, bool strip) {
    const std::string_view text = strip ? pytext::strip(line) : line;
    std::vector<std::string> out;
    for (const std::string_view field : pytext::split(text, "\t")) {
        out.emplace_back(field);
    }
    return out;
}

std::vector<std::string> split_commas(std::string_view value, bool drop_empty) {
    std::vector<std::string> out;
    for (const std::string_view part : pytext::split(value, ",")) {
        if (drop_empty && part.empty()) {
            continue;
        }
        out.emplace_back(part);
    }
    return out;
}

std::string field_count_message(std::string_view file_type, std::size_t count,
                                std::string_view line_number, std::string_view line) {
    return "File type detected is " + std::string(file_type) + " but line " +
           std::string(line_number) + ": " + std::string(line) + " does not have " +
           std::to_string(count) + " fields.";
}

std::string start_end_message(std::string_view line_number, std::string_view line) {
    return "Start position larger or equal than end for line #" + std::string(line_number) +
           ":\n" + std::string(line) + "\n";
}

std::string unsorted_message(std::string_view file_name, std::string_view previous_line,
                             std::string_view line) {
    return "Bed file not sorted. Please use a sorted bed file.\nFile: " + std::string(file_name) +
           "\nPrevious line: " + std::string(previous_line) + "\n Current line" +
           std::string(line) + " ";
}

std::string invalid_strand_message(std::string_view value, std::string_view line,
                                   std::string_view line_number) {
    // The Python passes the line where it means the line number and the
    // other way round.
    return "*Warning, invalid strand value found " + std::string(value) + " for line #" +
           std::string(line) + ":\n" + std::string(line_number) + "\n Setting strand to '.'\n";
}

std::string not_integer_message(std::string_view value, std::size_t column,
                                std::string_view line_number) {
    return "Value: " + std::string(value) + " in field " + std::to_string(column) + " at line " +
           std::string(line_number) + " is not an integer\n";
}

}  // namespace hicx::readbed
