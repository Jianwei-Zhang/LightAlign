#pragma once
#include <climits>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
namespace lightalign_threads {
inline unsigned int requested = 12;
inline unsigned int positive_integer(const std::string& text) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("expected a positive integer");
    const auto value = std::stoull(text);
    if (!value || value > INT_MAX) throw std::out_of_range("invalid thread count");
    return static_cast<unsigned int>(value);
}
inline unsigned int worker_count() {
    static const unsigned int workers = [] {
        const unsigned int selected = requested ? requested : 1;
        std::cout << "Thread configuration: hardware=" << std::thread::hardware_concurrency()
                  << ", requested=" << requested << ", workers=" << selected << std::endl;
        return selected;
    }();
    return workers;
}
}
