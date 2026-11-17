# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""`chemlarp` command for a pip installation: runs the native CLI bundled in the package."""
import os
import platform
import subprocess
import sys
from pathlib import Path


def executable() -> Path:
    name = "chemlarp.exe" if platform.system() == "Windows" else "chemlarp"
    return Path(__file__).resolve().parent / "bin" / name


def main() -> int:
    exe = executable()
    if not exe.exists():
        sys.stderr.write("the chemlarp program is not in this installation (%s)\n" % exe)
        return 1
    args = [str(exe), *sys.argv[1:]]
    if platform.system() == "Windows":
        return subprocess.call(args)
    os.execv(args[0], args)
    return 0  # not reached


if __name__ == "__main__":
    sys.exit(main())
