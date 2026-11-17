#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Writes a Windows .ico from PNG files (each stored as a PNG entry, which Windows Vista and later read).

    python3 make_ico.py out.ico in1.png [in2.png ...]   (each PNG at most 256 x 256)
"""
import struct
import sys


def png_size(data):
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit("not a PNG file")
    return struct.unpack(">II", data[16:24])


def main():
    out, pngs = sys.argv[1], sys.argv[2:]
    images = []
    for p in pngs:
        data = open(p, "rb").read()
        w, h = png_size(data)
        if w > 256 or h > 256:
            raise SystemExit("%s is %dx%d; entries can be at most 256 x 256" % (p, w, h))
        images.append((w, h, data))
    header = struct.pack("<HHH", 0, 1, len(images))
    offset = 6 + 16 * len(images)
    entries, blobs = b"", b""
    for w, h, data in images:
        entries += struct.pack("<BBBBHHII", w % 256, h % 256, 0, 0, 1, 32, len(data), offset + len(blobs))
        blobs += data
    open(out, "wb").write(header + entries + blobs)
    print("wrote", out, ", ".join("%dx%d" % (w, h) for w, h, _ in images))


if __name__ == "__main__":
    main()
