#include "hicx/qc_report_api.hpp"

#include "hicx/python_error.hpp"

namespace hicx::qcapi {

namespace {

const char* type_of(const QcStylerCell& cell) {
    return cell.is_int ? "int" : (cell.is_float ? "float" : "str");
}

double number(const QcStylerCell& cell) {
    return cell.is_int ? static_cast<double>(cell.integer) : cell.real;
}

[[noreturn]] void refuse(const char* left, const char* right) {
    throw PythonError("TypeError", std::string("unsupported operand type(s) for /: '") + left +
                                       "' and '" + right + "'");
}

}  // namespace

std::vector<double> divide_by_scalar(const Column& column, double scalar) {
    std::vector<double> out;
    out.reserve(column.size());
    for (const QcStylerCell& cell : column) {
        if (!cell.is_int && !cell.is_float) {
            refuse("str", "float");
        }
        out.push_back(number(cell) / scalar);
    }
    return out;
}

std::vector<double> divide_columns(const Column& numerators, const Column& denominators) {
    std::vector<double> out;
    out.reserve(numerators.size());
    for (std::size_t i = 0; i < numerators.size(); ++i) {
        const QcStylerCell& a = numerators[i];
        const QcStylerCell& b = denominators[i];
        if ((!a.is_int && !a.is_float) || (!b.is_int && !b.is_float)) {
            refuse(type_of(a), type_of(b));
        }
        out.push_back(number(a) / number(b));
    }
    return out;
}

Column sum_columns(const std::vector<Column>& columns) {
    Column out;
    if (columns.empty()) {
        return out;
    }
    bool all_int = true;
    for (const Column& column : columns) {
        for (const QcStylerCell& cell : column) {
            all_int = all_int && cell.is_int;
        }
    }
    for (std::size_t row = 0; row < columns[0].size(); ++row) {
        QcStylerCell total;
        if (all_int) {
            total.is_int = true;
            std::uint64_t sum = 0;
            for (const Column& column : columns) {
                sum += static_cast<std::uint64_t>(column[row].integer);
            }
            total.integer = static_cast<std::int64_t>(sum);
        } else {
            // numpy adds the objects of an object column in order
            total = columns[0][row];
            for (std::size_t c = 1; c < columns.size(); ++c) {
                const QcStylerCell& cell = columns[c][row];
                const bool left_text = !total.is_int && !total.is_float;
                const bool right_text = !cell.is_int && !cell.is_float;
                if (left_text && right_text) {
                    total.text += cell.text;
                } else if (left_text) {
                    throw PythonError("TypeError", std::string("can only concatenate str (not \"") +
                                                       type_of(cell) + "\") to str");
                } else if (right_text) {
                    throw PythonError("TypeError",
                                      std::string("unsupported operand type(s) for +: '") +
                                          type_of(total) + "' and 'str'");
                } else if (total.is_int && cell.is_int) {
                    total.integer = static_cast<std::int64_t>(
                        static_cast<std::uint64_t>(total.integer) +
                        static_cast<std::uint64_t>(cell.integer));
                } else {
                    total.real = number(total) + number(cell);
                    total.is_int = false;
                    total.is_float = true;
                }
            }
        }
        out.push_back(total);
    }
    return out;
}

}  // namespace hicx::qcapi
