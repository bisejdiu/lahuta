#pragma once

#include <ostream>
#include <streambuf>

class TeeBuf : public std::streambuf {
 public:
  TeeBuf(std::streambuf *a, std::streambuf *b) : a_(a), b_(b) {}

 protected:
  int overflow(int ch) override {
    if (ch == EOF) return !EOF;
    const char c = static_cast<char>(ch);
    if (a_->sputc(c) == EOF || b_->sputc(c) == EOF) return EOF;
    return ch;
  }

  int sync() override {
    const int ra = a_->pubsync();
    const int rb = b_->pubsync();
    return (ra == 0 && rb == 0) ? 0 : -1;
  }

 private:
  std::streambuf *a_;
  std::streambuf *b_;
};

class TeeStream : public std::ostream {
 public:
  TeeStream(std::ostream &a, std::ostream &b)
      : std::ostream(&buf_), buf_(a.rdbuf(), b.rdbuf()) {}

 private:
  TeeBuf buf_;
};
