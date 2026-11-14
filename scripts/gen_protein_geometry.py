"""Side-chain geometry and rotamer library for src/protein_geometry.inc (needs numpy).

    python3 scripts/gen_protein_geometry.py <dir with PDB files> [--exclude ID,ID,...] > src/protein_geometry.inc

From X-ray structures at 2.2 Å or better (PDB data are CC0) this collects, for every side-chain heavy atom
of the 20 amino acids, the internal coordinates that place it from three atoms already present:

  fixed atoms   mean bond length, bond angle and dihedral (ring atoms, branch partners, CB)
  chi atoms     mean bond length and bond angle; the dihedral is a side-chain torsion chi

and a rotamer library: the chi angles are binned (three 120° sectors for sp3 torsions, 60° sectors for a
terminal sp2 torsion, with the 180° symmetry of Asp, Glu, Phe and Tyr folded), and each bin with at least 5
observations gives one rotamer with its frequency and circular-mean chi angles. The dependence on the
backbone is kept as a factor on chi1: how much more often each chi1 sector occurs in helical, extended and
left-handed backbone conformations than overall.

Only residues with full occupancy, no alternate locations and B factors below 40 Å² are used. The entries
used for the published table are listed in docs/ARCHITECTURE.md; entries named with --exclude are left out
so they can serve as an independent test set.
"""
import math, os, sys
from collections import defaultdict
import numpy as np

# atom: (a, b, c, kind) placing `atom` bonded to c with angle b-c-atom and dihedral a-b-c-atom.
# kind 'x' = chi torsion (free), 'f' = fixed. Order matters: each atom's references come earlier.
TOPO = {
    "ALA": [("CB", "N", "C", "CA", "f")],
    "ARG": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD", "CA", "CB", "CG", "x"), ("NE", "CB", "CG", "CD", "x"),
            ("CZ", "CG", "CD", "NE", "x"), ("NH1", "CD", "NE", "CZ", "f"), ("NH2", "NH1", "NE", "CZ", "f")],
    "ASN": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("OD1", "CA", "CB", "CG", "x"), ("ND2", "OD1", "CB", "CG", "f")],
    "ASP": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("OD1", "CA", "CB", "CG", "x"), ("OD2", "OD1", "CB", "CG", "f")],
    "CYS": [("CB", "N", "C", "CA", "f"), ("SG", "N", "CA", "CB", "x")],
    "GLN": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD", "CA", "CB", "CG", "x"), ("OE1", "CB", "CG", "CD", "x"),
            ("NE2", "OE1", "CG", "CD", "f")],
    "GLU": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD", "CA", "CB", "CG", "x"), ("OE1", "CB", "CG", "CD", "x"),
            ("OE2", "OE1", "CG", "CD", "f")],
    "GLY": [],
    "HIS": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("ND1", "CA", "CB", "CG", "x"), ("CD2", "ND1", "CB", "CG", "f"),
            ("CE1", "CB", "CG", "ND1", "f"), ("NE2", "CB", "CG", "CD2", "f")],
    "ILE": [("CB", "N", "C", "CA", "f"), ("CG1", "N", "CA", "CB", "x"), ("CG2", "CG1", "CA", "CB", "f"), ("CD1", "CA", "CB", "CG1", "x")],
    "LEU": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD1", "CA", "CB", "CG", "x"), ("CD2", "CD1", "CB", "CG", "f")],
    "LYS": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD", "CA", "CB", "CG", "x"), ("CE", "CB", "CG", "CD", "x"),
            ("NZ", "CG", "CD", "CE", "x")],
    "MET": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("SD", "CA", "CB", "CG", "x"), ("CE", "CB", "CG", "SD", "x")],
    "PHE": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD1", "CA", "CB", "CG", "x"), ("CD2", "CD1", "CB", "CG", "f"),
            ("CE1", "CB", "CG", "CD1", "f"), ("CE2", "CB", "CG", "CD2", "f"), ("CZ", "CG", "CD1", "CE1", "f")],
    "PRO": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD", "CA", "CB", "CG", "x")],
    "SER": [("CB", "N", "C", "CA", "f"), ("OG", "N", "CA", "CB", "x")],
    "THR": [("CB", "N", "C", "CA", "f"), ("OG1", "N", "CA", "CB", "x"), ("CG2", "OG1", "CA", "CB", "f")],
    "TRP": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD1", "CA", "CB", "CG", "x"), ("CD2", "CD1", "CB", "CG", "f"),
            ("NE1", "CB", "CG", "CD1", "f"), ("CE2", "CB", "CG", "CD2", "f"), ("CE3", "CE2", "CG", "CD2", "f"),
            ("CZ2", "CG", "CD2", "CE2", "f"), ("CZ3", "CG", "CD2", "CE3", "f"), ("CH2", "CD2", "CE2", "CZ2", "f")],
    "TYR": [("CB", "N", "C", "CA", "f"), ("CG", "N", "CA", "CB", "x"), ("CD1", "CA", "CB", "CG", "x"), ("CD2", "CD1", "CB", "CG", "f"),
            ("CE1", "CB", "CG", "CD1", "f"), ("CE2", "CB", "CG", "CD2", "f"), ("CZ", "CG", "CD1", "CE1", "f"), ("OH", "CD1", "CE1", "CZ", "f")],
    "VAL": [("CB", "N", "C", "CA", "f"), ("CG1", "N", "CA", "CB", "x"), ("CG2", "CG1", "CA", "CB", "f")],
}
# Terminal sp2 torsions: period of the bins (360 = six 60° bins, 180 = symmetric group, three 60° bins).
SP2 = {("ASN", "OD1"): 360, ("ASP", "OD1"): 180, ("GLN", "OE1"): 360, ("GLU", "OE1"): 180, ("HIS", "ND1"): 360,
       ("PHE", "CD1"): 180, ("TYR", "CD1"): 180, ("TRP", "CD1"): 360}


def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0.dot(b1) * b1, b2 - b2.dot(b1) * b1
    return math.degrees(math.atan2(np.cross(b1, v).dot(w), v.dot(w)))


def angle(a, b, c):
    u, v = a - b, c - b
    return math.degrees(math.acos(max(-1.0, min(1.0, u.dot(v) / np.linalg.norm(u) / np.linalg.norm(v)))))


def circ_mean(v):
    r = np.radians(v)
    return math.degrees(math.atan2(np.sin(r).mean(), np.cos(r).mean()))


def bin_of(res, atom, chi):
    period = SP2.get((res, atom))
    if period is None:
        return int((chi % 360.0) // 120.0)
    return int(((chi + 30.0) % period) // 60.0)


def backbone_class(phi, psi):
    """0 helical (phi < 0, -120 <= psi < 50), 1 extended (phi < 0, other psi), 2 left-handed (phi >= 0)."""
    if phi >= 0:
        return 2
    return 0 if -120 <= psi < 50 else 1


def residues(path):
    res, order, bad = {}, [], set()
    resolution = None
    for l in open(path, errors="replace"):
        if l.startswith("REMARK   2 RESOLUTION."):
            try:
                resolution = float(l[23:30])
            except ValueError:
                pass
        if l.startswith("ENDMDL"):
            break
        if not l.startswith("ATOM  "):
            continue
        key = (l[21], l[22:27], l[17:20])
        if l[16] != " " or float(l[54:60] or 1) < 1.0 or float(l[60:66] or 0) > 40.0:
            bad.add(key)
        if key not in res:
            res[key] = {}
            order.append(key)
        res[key][l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    out = []
    for i, k in enumerate(order):
        at, cls = res[k], None
        if 0 < i < len(order) - 1 and all(a in at for a in ("N", "CA", "C")):
            p, n = res[order[i - 1]], res[order[i + 1]]
            if "C" in p and "N" in n and np.linalg.norm(p["C"] - at["N"]) < 1.45 and np.linalg.norm(at["C"] - n["N"]) < 1.45:
                cls = backbone_class(dihedral(p["C"], at["N"], at["CA"], at["C"]), dihedral(at["N"], at["CA"], at["C"], n["N"]))
        if k not in bad:
            out.append((k[2], at, cls))
    return resolution, out


def main():
    d = sys.argv[1]
    exclude = set()
    if "--exclude" in sys.argv:
        exclude = {x.upper() for x in sys.argv[sys.argv.index("--exclude") + 1].split(",")}
    geom = defaultdict(lambda: defaultdict(list))  # (res, atom) -> bond/angle/dihedral lists
    chis = defaultdict(list)
    chi1 = defaultdict(lambda: np.zeros((3, 3)))  # residue -> [backbone class][chi1 sector]
    used = []
    for f in sorted(os.listdir(d)):
        if not f.lower().endswith(".pdb") or f[:4].upper() in exclude:
            continue
        resolution, rs = residues(os.path.join(d, f))
        if resolution is None or resolution > 2.2:
            continue
        used.append(f[:4].upper())
        for name, at, cls in rs:
            topo = TOPO.get(name)
            if not topo or any(a not in at for a in ("N", "CA", "C")) or any(t[0] not in at for t in topo):
                continue
            ok, row = True, []
            for atom, a, b, c, kind in topo:
                bl = float(np.linalg.norm(at[atom] - at[c]))
                if not 1.1 < bl < 2.0:
                    ok = False
                    break
                g = geom[(name, atom)]
                g["b"].append(bl)
                g["a"].append(angle(at[b], at[c], at[atom]))
                dh = dihedral(at[a], at[b], at[c], at[atom])
                g["d"].append(dh)
                if kind == "x":
                    row.append((atom, dh))
            if ok and row:
                chis[name].append(row)
                if cls is not None and name != "PRO":
                    chi1[name][cls][int((row[0][1] % 360.0) // 120.0)] += 1
    out = []
    out.append("// Generated by scripts/gen_protein_geometry.py from %d X-ray structures (resolution <= 2.2 A). Do not edit." % len(used))
    out.append("// {residue, atom, ref a, ref b, ref c, bond c-atom (A), angle b-c-atom (deg), dihedral a-b-c-atom (deg), chi index or -1}")
    out.append("static const SideChainAtom kSideChainAtoms[] = {")
    for name in sorted(TOPO):
        k = 0
        for atom, a, b, c, kind in TOPO[name]:
            g = geom[(name, atom)]
            if not g["b"]:
                raise SystemExit(f"no data for {name} {atom}")
            chi = -1
            if kind == "x":
                chi, k = k, k + 1
            out.append('    {"%s", "%s", "%s", "%s", "%s", %.4f, %.2f, %.2f, %d},' % (
                name, atom, a, b, c, float(np.mean(g["b"])), float(np.mean(g["a"])), circ_mean(g["d"]) if kind == "f" else 0.0, chi))
    out.append("};")
    out.append("// {residue, frequency, number of chi, chi1..chi4 (deg)}; most frequent first within a residue")
    out.append("static const Rotamer kRotamers[] = {")
    total_rot = 0
    for name in sorted(chis):
        rows = chis[name]
        bins = defaultdict(list)
        for row in rows:
            bins[tuple(bin_of(name, atom, chi) for atom, chi in row)].append([chi for _, chi in row])
        kept = [(len(v), v) for v in bins.values() if len(v) >= 5]
        n = sum(c for c, _ in kept)
        for c, v in sorted(kept, key=lambda t: -t[0]):
            arr = np.array(v)
            means = []
            for j, (atom, _) in enumerate(rows[0]):
                col = arr[:, j]
                period = SP2.get((name, atom))
                if period == 180:  # symmetric group: average on the doubled angle
                    means.append(circ_mean(2 * col) / 2)
                else:
                    means.append(circ_mean(col))
            means += [0.0] * (4 - len(means))
            out.append('    {"%s", %.5f, %d, {%.1f, %.1f, %.1f, %.1f}},' % (name, c / n, arr.shape[1], *means))
            total_rot += 1
    out.append("};")
    out.append("// {residue, backbone class (0 helical, 1 extended, 2 left-handed), P(chi1 sector | class) / P(chi1 sector)")
    out.append("//  for the sectors [0,120), [120,240), [240,360) degrees}; one pseudo-count per cell")
    out.append("static const Chi1Bias kChi1Bias[] = {")
    for name in sorted(chi1):
        t = chi1[name] + 1.0
        overall = t.sum(axis=0) / t.sum()
        for cls in range(3):
            r = t[cls] / t[cls].sum() / overall
            out.append('    {"%s", %d, {%.3f, %.3f, %.3f}},' % (name, cls, r[0], r[1], r[2]))
    out.append("};")
    print("\n".join(out))
    print(f"{len(used)} structures, {sum(len(v) for v in chis.values())} side chains, {total_rot} rotamers", file=sys.stderr)
    print("used: " + " ".join(used), file=sys.stderr)


main()
