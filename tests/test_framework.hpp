#pragma once
// Minimal dependency-free test harness.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace octest {

struct Failure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct Case {
  const char* name;
  void (*fn)();
};

inline std::vector<Case>& registry() {
  static std::vector<Case> r;
  return r;
}

struct Registrar {
  Registrar(const char* name, void (*fn)()) { registry().push_back({name, fn}); }
};

inline std::string where(const char* file, int line) { return std::string(file) + ":" + std::to_string(line); }

// Path to a file under data/, valid from any working directory (CTest runs in the build tree).
#ifndef OMNICHEM_TEST_DATA_DIR
#define OMNICHEM_TEST_DATA_DIR "data"
#endif
inline std::string data_path(const std::string& relative) { return std::string(OMNICHEM_TEST_DATA_DIR) + "/" + relative; }

}  // namespace octest

#define OC_CAT2(a, b) a##b
#define OC_CAT(a, b) OC_CAT2(a, b)
#define TEST(name)                                                    \
  static void name();                                                 \
  static ::octest::Registrar OC_CAT(reg_, name)(#name, &name);         \
  static void name()

#define CHECK(cond)                                                                                  \
  do {                                                                                               \
    if (!(cond)) throw ::octest::Failure(::octest::where(__FILE__, __LINE__) + ": CHECK(" #cond ")"); \
  } while (0)

#define CHECK_EQ(a, b)                                                                                          \
  do {                                                                                                          \
    if (!((a) == (b)))                                                                                          \
      throw ::octest::Failure(::octest::where(__FILE__, __LINE__) + ": CHECK_EQ(" #a ", " #b ")");             \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                                                   \
  do {                                                                                                          \
    const double oc_a = (a), oc_b = (b);                                                                        \
    if (!(std::fabs(oc_a - oc_b) <= (tol))) {                                                                   \
      char oc_buf[256];                                                                                         \
      std::snprintf(oc_buf, sizeof oc_buf, ": CHECK_NEAR(" #a ", " #b ") got %.10g vs %.10g (tol %g)", oc_a,  \
                    oc_b, static_cast<double>(tol));                                                            \
      throw ::octest::Failure(::octest::where(__FILE__, __LINE__) + oc_buf);                                    \
    }                                                                                                           \
  } while (0)
