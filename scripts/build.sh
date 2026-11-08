#!/usr/bin/env bash
# Build library, CLI and tests, then run the test suite (macOS, Linux, MSYS2).
# Uses CMake when available, otherwise the Makefile.
set -euo pipefail
cd "$(dirname "$0")/.."
JOBS="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
if command -v cmake >/dev/null 2>&1 && [ "${OMNICHEM_USE_MAKE:-0}" != "1" ]; then
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build --parallel "$JOBS"
  ctest --test-dir build --output-on-failure
else
  make -j"$JOBS" all
  make test
fi
