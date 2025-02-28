#ifndef LAHUTA_UTILS_HPP
#define LAHUTA_UTILS_HPP

namespace lahuta {

constexpr static int ITERATIONS = 100;

constexpr bool isspace(const unsigned char c) { //
  return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

constexpr bool isalpha(const unsigned char c) { //
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

constexpr char __toupper__(const unsigned char c) { //
  return (c >= 'a' && c <= 'z') ? (c - 'a' + 'A') : c;
}

template <typename T> constexpr T __pow__(T base, int exp) {
  return (exp == 0) ? 1 : (exp > 0) ? base * __pow__(base, exp - 1) : 1 / __pow__(base, -exp);
}

constexpr double exp_approx(double x, int terms = ITERATIONS) {
  double result = 1.0;
  double term = 1.0;
  for (int n = 1; n < terms; ++n) {
    term *= x / n;
    result += term;
  }
  return result;
}

constexpr double log_approx(double x, double curr, int n, double result) {
  return n > ITERATIONS
             ? result
             : log_approx(x, curr * (x - 1) / (x + 1) * (x - 1) / (x + 1), n + 2, result + curr / n);
}

constexpr double log_approx(double x) {
  return (x > 0) ? 2 * log_approx(x, (x - 1) / (x + 1), 1, 0.0) : -1e10;
}

constexpr double log2_approx(double x) { return log_approx(x) / log_approx(2.0); }

} // namespace lahuta

#endif // LAHUTA_UTILS_HPP
