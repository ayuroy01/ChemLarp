#!/usr/bin/env python3
"""Pack / inspect ChemLarp Exchange (.clx) containers. See docs/CLX_SPEC.md.

    chemlarp parse 'CCO' --json | scripts/clx_pack.py pack ethanol.clx -
    scripts/clx_pack.py pack mix.clx structure.json --thermo thermo.json
    scripts/clx_pack.py inspect ethanol.clx
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from datetime import datetime, timezone

OCX_FORMAT = "clx"
OCX_VERSION = "0.1"


def _load(path: str) -> dict:
    text = sys.stdin.read() if path == "-" else open(path, encoding="utf-8").read()
    doc = json.loads(text)
    # Accept the CLI's {"ok":..,"structure":..,"diagnostics":..} envelope or a bare payload.
    return doc.get("structure", doc) if isinstance(doc, dict) else doc


def _validate_structure(s: dict) -> None:
    atoms, bonds = s["atoms"], s["bonds"]
    n, m = atoms["count"], bonds["count"]
    for col in ("element", "formal_charge", "implicit_h", "flags", "x", "y", "z", "species"):
        if len(atoms[col]) != n:
            raise ValueError(f"atoms.{col} has {len(atoms[col])} entries, expected {n}")
    for col in ("begin", "end", "order", "flags", "species"):
        if len(bonds[col]) != m:
            raise ValueError(f"bonds.{col} has {len(bonds[col])} entries, expected {m}")
    for k in range(m):
        if not (0 <= bonds["begin"][k] < n and 0 <= bonds["end"][k] < n):
            raise ValueError(f"bond {k} references an atom outside 0..{n - 1}")


def pack(out: str, structure: str, thermo: str | None, process: str | None) -> None:
    payloads = {"structure.json": _load(structure)}
    _validate_structure(payloads["structure.json"])
    if thermo:
        payloads["thermo.json"] = _load(thermo)
    if process:
        payloads["process.json"] = _load(process)
    encoded = {name: json.dumps(body, separators=(",", ":")).encode() for name, body in payloads.items()}
    manifest = {
        "format": OCX_FORMAT,
        "version": OCX_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "chemlarp-core 0.1.0",
        "payloads": [
            {"path": name, "kind": name.split(".")[0], "encoding": "json", "sha256": hashlib.sha256(data).hexdigest()}
            for name, data in encoded.items()
        ],
    }
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("manifest.json", json.dumps(manifest, indent=2))  # first entry, per spec
        for name, data in encoded.items():
            z.writestr(name, data)
    print(f"wrote {out}: {', '.join(encoded)}")


def inspect(path: str) -> None:
    with zipfile.ZipFile(path) as z:
        manifest = json.loads(z.read("manifest.json"))
        if manifest.get("format") != OCX_FORMAT:
            raise ValueError("not an .clx container")
        print(f"{path}: clx v{manifest['version']} ({manifest.get('generator', '?')})")
        for p in manifest["payloads"]:
            data = z.read(p["path"])
            ok = hashlib.sha256(data).hexdigest() == p["sha256"]
            print(f"  {p['path']:<16} {len(data):>8} bytes  sha256 {'ok' if ok else 'MISMATCH'}")
            if p["kind"] == "structure":
                s = json.loads(data)
                for sp in s["species"]:
                    print(f"    {sp['formula']:<10} {sp['smiles']}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("out")
    p.add_argument("structure", help="structure JSON file, or - for stdin")
    p.add_argument("--thermo")
    p.add_argument("--process")
    i = sub.add_parser("inspect")
    i.add_argument("path")
    a = ap.parse_args()
    if a.cmd == "pack":
        pack(a.out, a.structure, a.thermo, a.process)
    else:
        inspect(a.path)


if __name__ == "__main__":
    main()
