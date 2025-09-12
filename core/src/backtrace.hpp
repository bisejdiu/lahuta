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
