#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Fills the Homebrew formula and winget manifests for a release.

    python3 scripts/release_manifests.py --version 0.1.0 --dist DIR [--source-tarball FILE] [--out DIR]

DIR holds the release files (chemlarp-<version>-windows-x64.zip, ...). The source tarball is the one GitHub
serves for the tag (https://github.com/ayuroy01/ChemLarp/archive/refs/tags/v<version>.tar.gz): Homebrew checks
its sha256, so it must be that exact file, not a local `git archive`. Placeholders that cannot be filled (no
file given) are left in place and listed, and the exit status is 1, so an incomplete manifest is not mistaken
for a finished one.
"""
import argparse
import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", required=True)
    ap.add_argument("--dist", required=True)
    ap.add_argument("--source-tarball", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    v = a.version
    out = a.out or os.path.join(a.dist, "manifests")
    values = {"VERSION": v}
    win = os.path.join(a.dist, "chemlarp-%s-windows-x64.zip" % v)
    if os.path.exists(win):
        values["SHA256_WINDOWS_X64"] = sha256(win)
    if a.source_tarball:
        values["SHA256"] = sha256(a.source_tarball)
    templates = [("release/homebrew/chemlarp.rb", "homebrew/chemlarp.rb")]
    for f in sorted(os.listdir(os.path.join(ROOT, "release", "winget"))):
        templates.append(("release/winget/" + f, "winget/manifests/a/ayuroy01/ChemLarp/%s/%s" % (v, f)))
    left = []
    for src, dst in templates:
        text = open(os.path.join(ROOT, src)).read()
        for k, val in values.items():
            text = text.replace("@" + k + "@", val)
        left += ["%s: @%s@" % (dst, m) for m in re.findall(r"@([A-Z0-9_]+)@", text)]
        path = os.path.join(out, dst)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").write(text)
        print("wrote", path)
    if left:
        print("not filled:\n  " + "\n  ".join(left))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
