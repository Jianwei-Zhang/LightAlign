#pragma once
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#ifdef __linux__
#include <sched.h>
#endif

namespace lightalign_threads {
inline unsigned int requested = 0;

inline unsigned int positive_integer(const std::string& text) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("expected a positive integer");
    const auto value = std::stoull(text);
    if (value == 0 || value > INT_MAX)
        throw std::out_of_range("thread count must be in 1..INT_MAX");
    return static_cast<unsigned int>(value);
}

// Conservative rounding for fractional CPU quotas: 1.5 CPUs -> 1 worker.
inline unsigned int quota_threads(long long quota, long long period) {
    if (quota <= 0 || period <= 0) return 0;
    return static_cast<unsigned int>(std::min<long long>(INT_MAX,
        std::max<long long>(1, quota / period)));
}

inline unsigned int affinity_threads() {
#ifdef __linux__
    // Grow the buffer rather than assuming the kernel has at most 1024 CPUs.
    for (int count = 128; count <= 1048576; count *= 2) {
        cpu_set_t* mask = CPU_ALLOC(count);
        if (!mask) break;
        const size_t bytes = CPU_ALLOC_SIZE(count);
        CPU_ZERO_S(bytes, mask);
        const int status = sched_getaffinity(0, bytes, mask);
        const int saved_errno = errno;
        const int result = status == 0 ? CPU_COUNT_S(bytes, mask) : 0;
        CPU_FREE(mask);
        if (result > 0) return static_cast<unsigned int>(result);
        if (saved_errno != EINVAL) break;
    }
#endif
    return 0;
}

inline std::string unescape_mount(const std::string& input) {
    std::string result;
    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i] == '\\' && i + 3 < input.size() &&
            input[i+1] >= '0' && input[i+1] <= '7' &&
            input[i+2] >= '0' && input[i+2] <= '7' &&
            input[i+3] >= '0' && input[i+3] <= '7') {
            result += static_cast<char>((input[i+1]-'0')*64 +
                (input[i+2]-'0')*8 + input[i+3]-'0');
            i += 3;
        } else result += input[i];
    }
    return result;
}

inline bool has_controller(const std::string& list, const std::string& name) {
    return ("," + list + ",").find("," + name + ",") != std::string::npos;
}

// Resolve the process cgroup through mountinfo, including nonstandard mount
// points, then inspect all visible ancestors (a parent may have a tighter cap).
inline unsigned int cgroup_threads() {
#ifdef __linux__
    std::string v2_path, v1_path, line;
    std::ifstream groups("/proc/self/cgroup");
    while (std::getline(groups, line)) {
        const auto first = line.find(':');
        const auto second = line.find(':', first == std::string::npos ? 0 : first+1);
        if (first == std::string::npos || second == std::string::npos) continue;
        const auto controllers = line.substr(first+1, second-first-1);
        if (controllers.empty()) v2_path = line.substr(second+1);
        if (has_controller(controllers, "cpu")) v1_path = line.substr(second+1);
    }
    unsigned int limit = 0;
    std::ifstream mounts("/proc/self/mountinfo");
    while (std::getline(mounts, line)) {
        const auto separator = line.find(" - ");
        if (separator == std::string::npos) continue;
        std::istringstream left(line.substr(0, separator));
        std::istringstream right(line.substr(separator+3));
        std::string id, parent, device, root, mount, type, source, options;
        if (!(left >> id >> parent >> device >> root >> mount) ||
            !(right >> type >> source >> options)) continue;
        const bool v2 = type == "cgroup2";
        if (!v2 && !(type == "cgroup" && has_controller(options, "cpu"))) continue;
        const std::string group = v2 ? v2_path : v1_path;
        if (group.empty()) continue;
        root = unescape_mount(root);
        mount = unescape_mount(mount);
        std::string relative;
        if (root == "/") relative = group.substr(1);
        else if (group == root) relative = "";
        else if (group.compare(0, root.size()+1, root+"/") == 0)
            relative = group.substr(root.size()+1);
        else if (group == "/") relative = ""; // cgroup namespace root
        else continue;
        const auto boundary = std::filesystem::path(mount).lexically_normal();
        auto current = (boundary / relative).lexically_normal();
        const auto rel = current.lexically_relative(boundary);
        if (!rel.empty() && *rel.begin() == "..") continue;
        while (true) {
            long long quota = -1, period = 0;
            if (v2) {
                std::string value;
                std::ifstream file(current / "cpu.max");
                if (file >> value >> period) {
                    if (value != "max") {
                        try { quota = std::stoll(value); } catch (...) {}
                    }
                }
            } else {
                std::ifstream(current / "cpu.cfs_quota_us") >> quota;
                std::ifstream(current / "cpu.cfs_period_us") >> period;
            }
            const auto cap = quota_threads(quota, period);
            if (cap) limit = limit ? std::min(limit, cap) : cap;
            if (current == boundary || current == current.parent_path()) break;
            current = current.parent_path();
        }
    }
    return limit;
#else
    return 0;
#endif
}

inline unsigned int worker_count() {
    // One shared, cached decision for every pipeline stage.
    static const unsigned int selected = [] {
        const auto hardware = std::thread::hardware_concurrency();
        const auto affinity = affinity_threads();
        const auto quota = cgroup_threads();
        // An affinity mask is more specific than the hardware hint.
        unsigned int limit = affinity ? affinity : (hardware ? hardware : 1);
        unsigned int scheduler = 0;
        if (const char* value = std::getenv("SLURM_CPUS_PER_TASK")) {
            try { scheduler = positive_integer(value); }
            catch (...) { std::cerr << "Warning: ignoring invalid SLURM_CPUS_PER_TASK.\n"; }
        }
        if (scheduler) limit = std::min(limit, scheduler);
        if (quota) limit = std::min(limit, quota);
        // Explicit -t is the user's worker count, independent of detected caps.
        const unsigned int result = requested ? requested : limit;
        std::cout << "Thread configuration: hardware=" << hardware
                  << ", affinity=" << affinity << ", slurm=" << scheduler
                  << ", quota=" << quota << ", requested=" << requested
                  << ", workers=" << result << std::endl;
        return result;
    }();
    return selected;
}
} // namespace lightalign_threads
