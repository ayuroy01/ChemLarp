# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Builds the native library and CLI and puts them inside the Python package (see pyproject.toml).

The wheel holds the shared library (loaded with ctypes) and the `chemlarp` command-line program. It does not
depend on the Python version, so it is tagged py3-none-<platform>: one wheel per operating system and
architecture serves every Python from 3.10 on.

macOS and Linux build with the Makefile into build/wheel-<arch> (so the objects never mix with a normal build),
with MACOSX_DEPLOYMENT_TARGET (default 11.0) on macOS and libstdc++ linked statically on Linux. On macOS the
architecture comes from ARCHFLAGS when it is set (cibuildwheel cross-compiles x86_64 wheels on arm64 machines
that way), else from the machine. Windows builds with CMake and the static MSVC runtime, so the bundled
chemlarp.exe needs no Visual C++ redistributable. CHEMLARP_WHEEL_JOBS sets the number of parallel jobs.
"""
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution

try:
    from setuptools.command.bdist_wheel import bdist_wheel
except ImportError:  # setuptools < 70.1
    from wheel.bdist_wheel import bdist_wheel

ROOT = Path(__file__).resolve().parent
MACOS_TARGET = os.environ.get("MACOSX_DEPLOYMENT_TARGET", "11.0")


def macos_arch():
    m = re.search(r"-arch\s+(\S+)", os.environ.get("ARCHFLAGS", ""))
    return m.group(1) if m else platform.machine()


def native_build():
    """Builds the library and CLI; returns (library path, CLI path)."""
    jobs = os.environ.get("CHEMLARP_WHEEL_JOBS", str(os.cpu_count() or 2))
    system = platform.system()
    if system == "Windows":
        out = ROOT / "build" / "wheel-cmake"
        subprocess.check_call(["cmake", "-S", str(ROOT), "-B", str(out), "-DCHEMLARP_BUILD_TESTS=OFF",
                               "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded"])
        subprocess.check_call(["cmake", "--build", str(out), "--config", "Release", "--parallel", jobs,
                               "--target", "chemlarp", "chemlarp_cli"])
        for d in (out / "Release", out):
            if (d / "chemlarp.dll").exists():
                return d / "chemlarp.dll", d / "chemlarp.exe"
        raise RuntimeError("CMake did not produce chemlarp.dll")
    env = dict(os.environ)
    if system == "Darwin":
        arch = macos_arch()
        out = ROOT / "build" / ("wheel-" + arch)
        make = ["make", "-C", str(ROOT), "BUILD=" + str(out), "-j" + jobs, "ARCH_FLAGS=-arch " + arch]
        env["MACOSX_DEPLOYMENT_TARGET"] = MACOS_TARGET
        lib = out / "libchemlarp.dylib"
    else:
        out = ROOT / "build" / ("wheel-" + platform.machine())
        # Self-contained on Linux: no newer libstdc++ needed on the target machine than the one built against.
        make = ["make", "-C", str(ROOT), "BUILD=" + str(out), "-j" + jobs, "EXTRA_LDFLAGS=-static-libstdc++ -static-libgcc"]
        lib = out / "libchemlarp.so"
    subprocess.check_call(make + [str(lib), str(out / "chemlarp")], env=env)
    return lib, out / "chemlarp"


class BuildWithNative(build_py):
    def run(self):
        super().run()
        lib, cli = native_build()
        pkg = Path(self.build_lib) / "chemlarp"
        (pkg / "bin").mkdir(parents=True, exist_ok=True)
        shutil.copy2(lib, pkg / lib.name)
        shutil.copy2(cli, pkg / "bin" / cli.name)


class PlatformWheel(bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, plat = super().get_tag()
        if platform.system() == "Darwin":
            # The native code is built for one architecture, at MACOS_TARGET.
            plat = "macosx_" + MACOS_TARGET.replace(".", "_") + "_" + macos_arch()
        return "py3", "none", plat


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


setup(cmdclass={"build_py": BuildWithNative, "bdist_wheel": PlatformWheel}, distclass=BinaryDistribution)
