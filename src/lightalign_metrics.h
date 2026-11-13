#ifndef LIGHTALIGN_WINDOWS_METRICS_H
#define LIGHTALIGN_WINDOWS_METRICS_H

// Optional Windows-only resource reporting.  It is enabled only when the
// LIGHTALIGN_METRICS_FILE environment variable names an output JSON file.
// The peak is queried inside the process at normal termination, not estimated
// from periodic samples.  An abnormal termination may not produce this file.

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef _MSC_VER
#pragma comment(lib, "psapi.lib")
#endif

namespace lightalign_metrics {

struct State {
    std::wstring output_path;
    LARGE_INTEGER started{};
    LARGE_INTEGER frequency{};
    bool initialized = false;
};

inline State& state() {
    static State value;
    return value;
}

inline double filetime_seconds(const FILETIME& value) {
    ULARGE_INTEGER ticks{};
    ticks.LowPart = value.dwLowDateTime;
    ticks.HighPart = value.dwHighDateTime;
    return static_cast<double>(ticks.QuadPart) / 10000000.0;
}

inline void write_final_metrics() noexcept {
    try {
        State& value = state();
        PROCESS_MEMORY_COUNTERS_EX memory{};
        memory.cb = sizeof(memory);
        const BOOL memory_ok = GetProcessMemoryInfo(
            GetCurrentProcess(),
            reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&memory),
            sizeof(memory));
        const DWORD memory_error = memory_ok ? 0 : GetLastError();

        FILETIME created{}, exited{}, kernel{}, user{};
        const BOOL cpu_ok = GetProcessTimes(
            GetCurrentProcess(), &created, &exited, &kernel, &user);
        LARGE_INTEGER ended{};
        QueryPerformanceCounter(&ended);
        const double elapsed = value.frequency.QuadPart > 0
            ? static_cast<double>(ended.QuadPart - value.started.QuadPart)
                / static_cast<double>(value.frequency.QuadPart)
            : 0.0;

        char peak_working_set[64] = "null";
        char peak_commit[64] = "null";
        char user_seconds[64] = "null";
        char kernel_seconds[64] = "null";
        if (memory_ok) {
            std::snprintf(peak_working_set, sizeof(peak_working_set), "%llu",
                static_cast<unsigned long long>(memory.PeakWorkingSetSize));
            std::snprintf(peak_commit, sizeof(peak_commit), "%llu",
                static_cast<unsigned long long>(memory.PeakPagefileUsage));
        }
        if (cpu_ok) {
            std::snprintf(user_seconds, sizeof(user_seconds), "%.9f", filetime_seconds(user));
            std::snprintf(kernel_seconds, sizeof(kernel_seconds), "%.9f", filetime_seconds(kernel));
        }

        char json[2048]{};
        const int length = std::snprintf(json, sizeof(json),
            "{\n"
            "  \"schema_version\": 1,\n"
            "  \"process_id\": %lu,\n"
            "  \"elapsed_seconds\": %.9f,\n"
            "  \"cpu_user_seconds\": %s,\n"
            "  \"cpu_kernel_seconds\": %s,\n"
            "  \"memory_query_succeeded\": %s,\n"
            "  \"memory_query_error\": %lu,\n"
            "  \"peak_working_set_bytes\": %s,\n"
            "  \"peak_commit_bytes\": %s,\n"
            "  \"peak_source\": \"Windows GetProcessMemoryInfo at normal termination\"\n"
            "}\n",
            static_cast<unsigned long>(GetCurrentProcessId()), elapsed,
            user_seconds, kernel_seconds, memory_ok ? "true" : "false",
            static_cast<unsigned long>(memory_error), peak_working_set, peak_commit);
        if (length <= 0 || static_cast<size_t>(length) >= sizeof(json)) {
            std::fprintf(stderr, "Resource reporting: JSON formatting failed.\n");
            return;
        }

        HANDLE output = CreateFileW(value.output_path.c_str(), GENERIC_WRITE,
            FILE_SHARE_READ, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (output == INVALID_HANDLE_VALUE) {
            std::fprintf(stderr, "Resource reporting: cannot create metrics file (error %lu).\n",
                static_cast<unsigned long>(GetLastError()));
            return;
        }
        DWORD written = 0;
        const BOOL wrote = WriteFile(output, json, static_cast<DWORD>(length), &written, nullptr);
        CloseHandle(output);
        if (!wrote || written != static_cast<DWORD>(length)) {
            std::fprintf(stderr, "Resource reporting: incomplete metrics file.\n");
        }
    } catch (...) {
        std::fprintf(stderr, "Resource reporting: failed to write final metrics.\n");
    }
}

inline void initialize() {
    const wchar_t* configured = _wgetenv(L"LIGHTALIGN_METRICS_FILE");
    if (configured == nullptr || configured[0] == L'\0') {
        return;
    }
    State& value = state();
    if (value.initialized) {
        return;
    }
    value.output_path = configured;
    QueryPerformanceFrequency(&value.frequency);
    QueryPerformanceCounter(&value.started);
    if (std::atexit(write_final_metrics) == 0) {
        value.initialized = true;
    } else {
        std::fprintf(stderr, "Resource reporting: cannot register termination callback.\n");
    }
}

} // namespace lightalign_metrics
#else
namespace lightalign_metrics {
inline void initialize() {}
}
#endif

#endif
