#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Checks an installed chemlarp wheel (run outside the source tree, e.g. by cibuildwheel's test step): the
library is the one inside the installed package, the main entry points work, and the `chemlarp` command runs."""
import os
import platform
import shutil
import subprocess
import sys
import sysconfig

import chemlarp

lib = str(chemlarp.library_path())
print("python", sys.version.split()[0], platform.machine(), "library", lib)
assert "site-packages" in lib or "dist-packages" in lib, "the library must come from the installed package"
assert chemlarp.version(), "no version"
assert chemlarp.canonical("C1=CC=CC=C1")["canonical"] == "c1ccccc1"
assert chemlarp.inchi("CCO")["key"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
assert chemlarp.iupac_name("CC(C)Cc1ccc(cc1)C(C)C(=O)O") == "2-[4-(2-methylpropyl)phenyl]propanoic acid"
e = chemlarp.embed("C[C@@H](C(=O)O)N")
assert e["stereo"]["centres_ok"] == 1 and len(e["coordinates"]) == 13
q = chemlarp.qc("O", method="hf", basis="STO-3G")
assert abs(q["energy"] - -74.96) < 0.05, q["energy"]
m = chemlarp.mbar([[0, 0.5, 2, 0.1], [0.2, 1.0, 0.1, 0]], [2, 2])
assert abs(m["f"][1] - -0.2639) < 1e-3
# The console script pip installed: in the interpreter's scripts directory (on PATH only when activated).
name = "chemlarp.exe" if platform.system() == "Windows" else "chemlarp"
exe = os.path.join(sysconfig.get_path("scripts"), name)
if not os.path.exists(exe):
    exe = shutil.which("chemlarp")
assert exe and os.path.exists(exe), "the chemlarp command was not installed"
out = subprocess.run([exe, "canon", "OCC"], capture_output=True, text=True, check=True).stdout.strip()
assert out == "CCO", out
print("wheel smoke test passed:", chemlarp.version(), os.path.basename(lib), exe)
