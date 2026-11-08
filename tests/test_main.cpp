#include <cstdio>
#include <exception>

#include "test_framework.hpp"

int main() {
  int failed = 0;
  std::setvbuf(stdout, nullptr, _IOLBF, 0);  // line-buffered: progress survives a hang or crash
  for (const auto& c : octest::registry()) {
    try {
      c.fn();
      std::printf("  ok    %s\n", c.name);
    } catch (const std::exception& e) {
      ++failed;
      std::printf("  FAIL  %s\n        %s\n", c.name, e.what());
    }
  }
  std::printf("\n%zu tests, %d failed\n", octest::registry().size(), failed);
  return failed == 0 ? 0 : 1;
}
