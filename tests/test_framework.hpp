#pragma once
// Minimal dependency-free test harness.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace cltest {

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
#ifndef CHEMLARP_TEST_DATA_DIR
#define CHEMLARP_TEST_DATA_DIR "data"
#endif
inline std::string data_path(const std::string& relative) { return std::string(CHEMLARP_TEST_DATA_DIR) + "/" + relative; }

}  // namespace cltest

#define CL_CAT2(a, b) a##b
#define CL_CAT(a, b) CL_CAT2(a, b)
#define TEST(name)                                                    \
  static void name();                                                 \
  static ::cltest::Registrar CL_CAT(reg_, name)(#name, &name);         \
  static void name()

#define CHECK(cond)                                                                                  \
  do {                                                                                               \
    if (!(cond)) throw ::cltest::Failure(::cltest::where(__FILE__, __LINE__) + ": CHECK(" #cond ")"); \
  } while (0)

#define CHECK_EQ(a, b)                                                                                          \
  do {                                                                                                          \
    if (!((a) == (b)))                                                                                          \
      throw ::cltest::Failure(::cltest::where(__FILE__, __LINE__) + ": CHECK_EQ(" #a ", " #b ")");             \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                                                   \
  do {                                                                                                          \
    const double cl_a = (a), cl_b = (b);                                                                        \
    if (!(std::fabs(cl_a - cl_b) <= (tol))) {                                                                   \
      char cl_buf[256];                                                                                         \
      std::snprintf(cl_buf, sizeof cl_buf, ": CHECK_NEAR(" #a ", " #b ") got %.10g vs %.10g (tol %g)", cl_a,  \
                    cl_b, static_cast<double>(tol));                                                            \
      throw ::cltest::Failure(::cltest::where(__FILE__, __LINE__) + cl_buf);                                    \
    }                                                                                                           \
  } while (0)
