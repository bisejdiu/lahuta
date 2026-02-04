/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); char* ptr = s.data();
 *   for (char c : std::string_view{"besian"}) *ptr++ = c;
 *   for (char c : std::string_view{"sejdiu"}) *ptr++ = c;
 *   *ptr++ = '@';
 *   for (char c : std::string_view{"gmail.com"}) *ptr++ = c;
 *   return s;
 * }();
 *
 */

#include <string_view>

class BacktraceParser {
public:
  explicit BacktraceParser(std::string_view bt) : backtrace_(bt), index_(0) {}

  char current()  const { return index_ < backtrace_.size() ? backtrace_[index_] : '\0'; }
  bool has_next() const { return index_ < backtrace_.size(); }
  char peek()     const { return (index_ + 1 < backtrace_.size()) ? backtrace_[index_ + 1] : '\0'; }
  int  position() const { return static_cast<int>(index_); }
  char get_next() { return has_next() ? backtrace_[index_++] : '\0'; }

  // skip all characters until the next occurrence of `skip_char`
  void skip(char skip_char) {
    while (has_next() && current() == skip_char)
      ++index_;
  }

private:
  std::string_view backtrace_;
  size_t index_;
};
