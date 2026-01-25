#ifndef LAHUTA_CLI_TIME_UTILS_HPP
#define LAHUTA_CLI_TIME_UTILS_HPP

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

namespace lahuta::cli {

inline std::string current_timestamp_string() {
  using namespace std::chrono;
  const auto now = system_clock::now();
  const auto tt  = system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &tt);
#else
  localtime_r(&tt, &tm);
#endif
  const auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % milliseconds(1000);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S") << '_' << std::setw(3) << std::setfill('0') << ms.count();
  return oss.str();
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_TIME_UTILS_HPP
