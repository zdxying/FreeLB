/* Minimal assertion helpers for the FreeLB test suite.
 *
 * Usage:
 *   #include "../../test_util.h"   // from tests/<category>/<name>/
 *   ...
 *   CHECK(cond);
 *   CHECK_EQ(got, expected);
 *   CHECK_NEAR(got, expected, tol);
 *   ...
 *   return flbtest::summary("my_test");
 *
 * summary() prints the check/failure counts and returns 0 iff there were no
 * failures, so `return flbtest::summary(...)` makes the process exit code
 * meaningful to tests/run_tests.sh.
 */
#pragma once

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

namespace flbtest {

inline int& checks() {
  static int n = 0;
  return n;
}
inline int& failures() {
  static int n = 0;
  return n;
}

template <typename A, typename B>
inline void checkEq(const A& a, const B& b, const char* ea, const char* eb,
                    const char* file, int line) {
  ++checks();
  if (!(a == b)) {
    ++failures();
    std::ostringstream os;
    os << "CHECK_EQ FAILED " << file << ":" << line << ": " << ea << " == " << eb
       << "  (got \"" << a << "\" vs \"" << b << "\")";
    std::cerr << os.str() << std::endl;
  }
}

template <typename A, typename B>
inline void checkNear(const A& a, const B& b, double tol, const char* ea,
                      const char* eb, const char* file, int line) {
  ++checks();
  if (!(std::abs(static_cast<double>(a) - static_cast<double>(b)) <= tol)) {
    ++failures();
    std::cerr << "CHECK_NEAR FAILED " << file << ":" << line << ": " << ea
              << " ~= " << eb << " (got " << a << " vs " << b << ", tol " << tol
              << ")" << std::endl;
  }
}

inline int summary(const char* name) {
  std::cout << "[" << name << "] checks: " << checks()
            << ", failures: " << failures() << std::endl;
  std::cout << (failures() == 0 ? "ALL PASSED" : "FAILED") << std::endl;
  return failures() == 0 ? 0 : 1;
}

}  // namespace flbtest

#define CHECK(cond)                                     \
  do {                                                  \
    ++flbtest::checks();                                \
    if (!(cond)) {                                      \
      ++flbtest::failures();                            \
      std::cerr << "CHECK FAILED " << __FILE__ << ":" << __LINE__ << ": " \
                << #cond << std::endl;                  \
    }                                                   \
  } while (0)

#define CHECK_EQ(a, b) \
  flbtest::checkEq((a), (b), #a, #b, __FILE__, __LINE__)

#define CHECK_NEAR(a, b, tol) \
  flbtest::checkNear((a), (b), (tol), #a, #b, __FILE__, __LINE__)
