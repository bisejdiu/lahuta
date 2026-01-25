#ifndef LAHUTA_CLI_EXTENSION_UTILS_HPP
#define LAHUTA_CLI_EXTENSION_UTILS_HPP

#include <cctype>
#include <string>
#include <string_view>
#include <vector>

namespace lahuta::cli {

namespace detail {

inline std::string_view trim_view(std::string_view sv) {
  while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front()))) {
    sv.remove_prefix(1);
  }
  while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back()))) {
    sv.remove_suffix(1);
  }
  return sv;
}

inline void split_argument_list(std::string_view raw, bool keep_empty_tokens, std::vector<std::string> &out) {
  if (raw.empty()) {
    if (keep_empty_tokens) out.emplace_back();
    return;
  }

  std::size_t pos = 0;
  while (pos < raw.size()) {
    std::size_t end = raw.find(',', pos);
    if (end == std::string_view::npos) end = raw.size();

    std::string_view token = trim_view(raw.substr(pos, end - pos));
    if (!token.empty()) {
      out.emplace_back(token);
    } else if (keep_empty_tokens) {
      out.emplace_back();
    }

    if (end == raw.size()) break;
    pos = end + 1;
  }
}

} // namespace detail

inline void parse_extension_argument(std::string_view raw, std::vector<std::string> &out) {
  detail::split_argument_list(raw, true, out);
}

inline void parse_file_argument(std::string_view raw, std::vector<std::string> &out) {
  detail::split_argument_list(raw, false, out);
}

inline std::string describe_extensions(const std::vector<std::string> &extensions) {
  if (extensions.empty()) return "(none)";
  std::string buf;
  for (std::size_t i = 0; i < extensions.size(); ++i) {
    if (i > 0) buf.append(", ");
    buf.append(extensions[i]);
  }
  return buf;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_EXTENSION_UTILS_HPP
