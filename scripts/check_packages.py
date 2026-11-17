#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""The package version must be the same everywhere it is written (run in CI): CMakeLists.txt project(VERSION),
cl_version() in src/c_api.cpp, pyproject.toml and release/npm/package.json. --version-only prints it."""
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read(p):
    return open(os.path.join(ROOT, p)).read()


found = {
    "CMakeLists.txt": re.search(r"project\(\w+ VERSION ([\d.]+)", read("CMakeLists.txt")).group(1),
    "src/c_api.cpp": re.search(r'cl_version\(void\) \{ return "([\d.]+)"', read("src/c_api.cpp")).group(1),
    "pyproject.toml": re.search(r'^version = "([\d.]+)"', read("pyproject.toml"), re.M).group(1),
    "release/npm/package.json": json.loads(read("release/npm/package.json"))["version"],
}
if sys.argv[1:] == ["--version-only"]:
    print(found["pyproject.toml"])
    sys.exit(0)
for k, v in found.items():
    print("%-28s %s" % (k, v))
if len(set(found.values())) != 1:
    sys.exit("versions differ")
print("all the same")
