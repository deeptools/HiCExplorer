#include <iostream>
#include <string>

#include <coolercpp/parallel.hpp>

#include "hic2cool/hic2cool.hpp"

namespace hic2cool {

int available_threads() { return coolercpp::available_threads(); }

Console Console::standard() {
    Console console;
    console.out = [](const std::string& line) { std::cout << line << '\n'; };
    console.err = [](const std::string& line) {
        std::cout.flush();
        std::cerr << line << std::endl;
    };
    console.input = [](const std::string& prompt) {
        std::cout << prompt << std::flush;
        std::string line;
        if (!std::getline(std::cin, line)) {
            throw ExitError("EOFError: EOF when reading a line");
        }
        return line;
    };
    return console;
}

}  // namespace hic2cool
