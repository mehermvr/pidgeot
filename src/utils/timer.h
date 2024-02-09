#pragma once
#include <chrono>
#include <iostream>
#include <string>

namespace utils {
class Timer {
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

    std::string _label;
    TimePoint _start_time;

public:
    Timer() : _label("Execution"), _start_time(Clock::now()) {}
    explicit Timer(const std::string &label) : _label(label), _start_time(Clock::now()) {}

    /* default other stuff */
    Timer(const Timer &) = default;
    Timer(Timer &&) = default;
    Timer &operator=(const Timer &) = default;
    Timer &operator=(Timer &&) = default;

    ~Timer() {
        auto end_time = Clock::now();
        Duration duration = end_time - _start_time;

        std::cout << _label << " took " << duration.count() << " seconds.\n";
    }
};
}  // namespace utils
