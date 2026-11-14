"""Validation of protein repair and the hydrogen-bond network (needs Biopython and numpy; build first).

    python3 scripts/validate_protein_build.py sidechains <PDB files...>
    python3 scripts/validate_protein_build.py loops <PDB files...>
    python3 scripts/validate_protein_build.py neutron <mmCIF files of neutron structures...>
    python3 scripts/validate_protein_build.py formats <dir with X.pdb> <dir with X.cif>

sidechains  Side chains are deleted from X-ray structures and rebuilt by `chemlarp protein`: all at once, and
            one residue at a time. Reported: chi1 within 40° of the crystal, chi1 and chi2 both within 40°
            (chi2 of Phe/Tyr modulo 180°; Asn, Asp, His chi2 not counted), for residues with complete side
            chains, full occupancy and B below 40 Å². Use structures that were not used for the rotamer
            library (see scripts/gen_protein_geometry.py --exclude).
loops       Stretches of 3, 4, 6, 8 and 12 residues are deleted, listed in REMARK 465 and rebuilt; backbone
            (N, CA, C) rmsd to the crystal and the lengths of the two closing peptide bonds.
neutron     Structures with experimental hydrogen (deuterium) positions: histidine tautomers, the
            orientation of Asn/Gln amides and His rings (as deposited, and after turning every one of them
            by 180° in the input), and whether each Ser/Thr/Tyr hydroxyl hydrogen donates to the same
            acceptor as in the neutron structure.
formats     The PDB and mmCIF files of the same entries give identical results.
"""
import json, math, os, subprocess, sys, tempfile, warnings
from collections import defaultdict
import numpy as np

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
CL = os.path.join(HERE, "..", "build", "chemlarp")
TMP = tempfile.mkdtemp()
BB = {"N", "CA", "C", "O", "OXT"}
SIDE = {"ALA": "CB", "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2", "ASP": "CB CG OD1 OD2", "CYS": "CB SG",
        "GLN": "CB CG CD OE1 NE2", "GLU": "CB CG CD OE1 OE2", "GLY": "", "HIS": "CB CG ND1 CD2 CE1 NE2", "ILE": "CB CG1 CG2 CD1",
        "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ", "MET": "CB CG SD CE", "PHE": "CB CG CD1 CD2 CE1 CE2 CZ",
        "PRO": "CB CG CD", "SER": "CB OG", "THR": "CB OG1 CG2", "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
        "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "CB CG1 CG2"}


def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0.dot(b1) * b1, b2 - b2.dot(b1) * b1
    return math.degrees(math.atan2(np.cross(b1, v).dot(w), v.dot(w)))


def adiff(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def read_pdb(text):
    """First model: header lines, and residues in order as (key, name, {atom: (xyz, line)})."""
    res, order, head = {}, [], []
    for l in text.splitlines():
        if l.startswith("ENDMDL"):
            break
        if l[:6] in ("ATOM  ", "HETATM"):
            if l[16] not in " A" or l[76:78].strip() in ("H", "D"):
                continue
            key = (l[21], l[22:27])
            if key not in res:
                res[key] = (l[17:20], {}, l[:6] == "HETATM")
                order.append(key)
            res[key][1][l[12:16].strip()] = (np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l)
        elif l[:6] in ("HEADER", "TITLE ", "SSBOND"):
            head.append(l)
    return head, [(k, *res[k]) for k in order]


def write_pdb(head, residues, missing=()):
    out = list(head)
    for name, chain, num in missing:
        out.append(f"REMARK 465     {name} {chain} {num:>5s}")
    for key, name, atoms, het in residues:
        out += [l for _, l in atoms.values()]
    return "\n".join(out) + "\nEND\n"


def run(text, *flags, suffix=".pdb"):
    src, dst = os.path.join(TMP, "in" + suffix), os.path.join(TMP, "out.pdb")
    open(src, "w").write(text)
    subprocess.run([CL, "protein", src, "--out", dst, *flags], capture_output=True, text=True)
    return read_pdb(open(dst).read())[1]


def good(name, atoms, het):
    if het or name not in SIDE or not SIDE[name]:
        return False
    if any(a not in atoms for a in SIDE[name].split()) or any(a not in atoms for a in ("N", "CA", "C")):
        return False
    return all(float(l[54:60]) >= 1.0 and float(l[60:66]) <= 40.0 for _, l in atoms.values())


def chis(name, atoms):
    g = next((a for a in ("CG", "CG1", "OG", "OG1", "SG") if a in atoms), None)
    if not g or name == "PRO":
        return None, None
    x = lambda n: atoms[n][0]
    c1 = dihedral(x("N"), x("CA"), x("CB"), x(g))
    d = next((a for a in ("CD", "CD1", "OD1", "ND1", "SD") if a in atoms), None)
    c2 = dihedral(x("CA"), x("CB"), x(g), x(d)) if d and name not in ("ASN", "ASP", "HIS") else None
    return c1, c2


def strip(residue):
    key, name, atoms, het = residue
    return (key, name, {a: v for a, v in atoms.items() if a in BB}, het) if name in SIDE and not het else residue


def compare(stat, ref, new):
    name, a, b = ref[1], ref[2], new[2]
    c1, c2 = chis(name, a)
    if c1 is None or any(n not in b for n in SIDE[name].split()):
        return
    n1, n2 = chis(name, b)
    ok = adiff(c1, n1) < 40.0
    stat["n1"] += 1
    stat["ok1"] += ok
    if c2 is not None:
        d = adiff(c2, n2)
        if name in ("PHE", "TYR"):
            d = min(d, 180.0 - d)
        stat["n2"] += 1
        stat["ok2"] += ok and d < 40.0
    for n in SIDE[name].split():
        stat["sq"] += float(np.sum((a[n][0] - b[n][0]) ** 2))
        stat["na"] += 1


def sidechains(files):
    tot = {m: defaultdict(float) for m in ("all", "one")}
    for f in files:
        head, residues = read_pdb(open(f).read())
        keep = [i for i, r in enumerate(residues) if good(*r[1:])]
        st = defaultdict(float)
        new = run(write_pdb(head, [strip(r) for r in residues]), "--no-flips")
        index = {r[0]: r for r in new}
        for i in keep:
            compare(st, residues[i], index[residues[i][0]])
        s1 = defaultdict(float)
        for i in keep:
            new = run(write_pdb(head, [strip(r) if j == i else r for j, r in enumerate(residues)]), "--no-flips")
            compare(s1, residues[i], next(r for r in new if r[0] == residues[i][0]))
        for m, s in (("all", st), ("one", s1)):
            for k, v in s.items():
                tot[m][k] += v
        print(f"  {os.path.basename(f):10s} all at once: chi1 {int(st['ok1'])}/{int(st['n1'])}   one at a time: chi1 {int(s1['ok1'])}/{int(s1['n1'])}")
    for m, s in tot.items():
        print(f"{'all side chains removed' if m == 'all' else 'one side chain removed '}: chi1 {int(s['ok1'])}/{int(s['n1'])} = {100 * s['ok1'] / s['n1']:.1f}%, "
              f"chi1+2 {int(s['ok2'])}/{int(s['n2'])} = {100 * s['ok2'] / s['n2']:.1f}%, side-chain rmsd {math.sqrt(s['sq'] / s['na']):.2f} Å")


def loops(files):
    for length in (3, 4, 6, 8, 12):
        rms, bonds, failed, cases = [], [], 0, 0
        for f in files:
            head, residues = read_pdb(open(f).read())
            head = [l for l in head if not l.startswith("SSBOND")]
            i = 3
            while i + length + 3 < len(residues):
                seg = residues[i - 1:i + length + 1]
                ok = all(not r[3] and r[1] in SIDE and all(a in r[2] for a in ("N", "CA", "C")) and r[0][0] == seg[0][0][0] for r in seg)
                ok = ok and all(np.linalg.norm(seg[k][2]["C"][0] - seg[k + 1][2]["N"][0]) < 1.45 for k in range(len(seg) - 1))
                if ok:
                    cases += 1
                    gone = residues[i:i + length]
                    new = run(write_pdb(head, residues[:i] + residues[i + length:], [(r[1], r[0][0], r[0][1].strip()) for r in gone]), "--no-flips")
                    index = {r[0]: r for r in new}
                    if all(r[0] in index for r in gone):
                        d = [index[r[0]][2][a][0] - r[2][a][0] for r in gone for a in ("N", "CA", "C")]
                        rms.append(math.sqrt(np.mean(np.sum(np.array(d) ** 2, axis=1))))
                        chain = [index[r[0]] for r in residues[i - 1:i + length + 1]]
                        bonds += [float(np.linalg.norm(chain[0][2]["C"][0] - chain[1][2]["N"][0])), float(np.linalg.norm(chain[-2][2]["C"][0] - chain[-1][2]["N"][0]))]
                    else:
                        failed += 1
                i += length + 2
        rms.sort()
        print(f"{length:2d} residues: {cases} cases, {failed} not built; backbone rmsd median {rms[len(rms) // 2]:.2f} Å, "
              f"quartiles {rms[len(rms) // 4]:.2f}-{rms[3 * len(rms) // 4]:.2f}; closing C-N bonds {min(bonds):.2f}-{max(bonds):.2f} Å")


def neutron(files):
    from Bio.PDB import MMCIFParser, PDBIO, Select

    class Heavy(Select):
        def accept_atom(self, a):
            return a.element not in ("H", "D") and (not a.is_disordered() or a.get_altloc() in ("A", " "))

        def accept_residue(self, r):
            return True

    T = defaultdict(int)
    for f in files:
        st = MMCIFParser(QUIET=True).get_structure("n", f)
        model = st[0]
        if any(len(ch.id) > 1 for ch in model):
            continue
        hyd = {}  # (chain, resseq, icode) -> {name: xyz} of deposited H/D, D names mapped to H names
        for ch in model:
            for r in ch:
                for a in r:
                    if a.element in ("H", "D"):
                        hyd.setdefault((ch.id, r.id[1], r.id[2]), {})["H" + a.get_id()[1:]] = (a.coord.astype(float), a.get_occupancy() or 0.0)
        io = PDBIO()
        io.set_structure(st)
        path = os.path.join(TMP, "heavy.pdb")
        io.save(path, Heavy())
        text = open(path).read()
        head, residues = read_pdb(text)

        def protonated(flip):
            rs = []
            for key, name, atoms, het in residues:
                atoms = dict(atoms)
                if flip and not het:
                    for a, b in {"ASN": [("OD1", "ND2")], "GLN": [("OE1", "NE2")], "HIS": [("ND1", "CD2"), ("CE1", "NE2")]}.get(name, []):
                        if a in atoms and b in atoms:
                            (xa, la), (xb, lb) = atoms[a], atoms[b]
                            atoms[a] = (xb, la[:30] + lb[30:54] + la[54:])
                            atoms[b] = (xa, lb[:30] + la[30:54] + lb[54:])
                rs.append((key, name, atoms, het))
            src, dst = os.path.join(TMP, "n.pdb"), os.path.join(TMP, "n_out.pdb")
            open(src, "w").write(write_pdb(head, rs))
            p = subprocess.run([CL, "protein", src, "--pH", "7", "--json"], capture_output=True, text=True)
            j = json.loads(p.stdout)
            out = {}
            for l in j["pdb"].splitlines():
                if l[:6] in ("ATOM  ", "HETATM"):
                    out.setdefault((l[21], int(l[22:26]), l[26]), {})[l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            return j, out

        acceptors = [(key, a, v[0]) for key, name, atoms, het in residues for a, v in atoms.items() if v[1][76:78].strip() == "O"]
        for flip in (0, 1):
            j, out = protonated(flip)
            states = {(r["chain"], r["number"]): r["state"] for r in j["residues"]}
            und = {u.split(":")[0] for u in j["undetermined"]}
            for key, name, atoms, het in residues:
                k = (key[0], int(key[1][:4]), key[1][4])
                h, o = hyd.get(k, {}), out.get(k, {})
                label = f"{key[0]} {k[1]} {name}"
                probe = {"ASN": ("ND2", ("HD21", "HD22")), "GLN": ("NE2", ("HE21", "HE22")), "HIS": ("ND1", ("HD1", "HE2"))}.get(name)
                if probe and probe[0] in atoms and probe[0] in o and any(x in h for x in probe[1]):
                    same = np.linalg.norm(o[probe[0]] - atoms[probe[0]][0]) < 0.5
                    if not flip:
                        T["orient_n"] += 1
                        T["orient_wrong"] += not same
                    else:
                        T["flip_n"] += 1
                        T["flip_ok"] += same
                        T["flip_und"] += label in und
                        T["flip_det_ok"] += same and label not in und
                if flip:
                    continue
                if name == "HIS" and ("HD1" in h or "HE2" in h):
                    ref = "HIP" if "HD1" in h and "HE2" in h else "HID" if "HD1" in h else "HIE"
                    ours = states.get((key[0], k[1]))
                    if ref != "HIP" and ours != "HIP":
                        T["taut_n"] += 1
                        T["taut_ok"] += ref == ours
                for hn, tol, tag in (("HG", 45, "oh"), ("HG1", 45, "oh"), ("HH", 90, "tyr")):
                    if name in ("SER", "THR", "TYR") and hn in h and hn in o:
                        O, X, R = {"HG": ("OG", "CB", "CA"), "HG1": ("OG1", "CB", "CA"), "HH": ("OH", "CZ", "CE1")}[hn]
                        if not all(a in atoms for a in (O, X, R)):
                            continue
                        t1 = dihedral(atoms[R][0], atoms[X][0], atoms[O][0], h[hn][0])
                        t2 = dihedral(atoms[R][0], atoms[X][0], atoms[O][0], o[hn])
                        T[tag + "_n"] += 1
                        T[tag + "_ok"] += adiff(t1, t2) <= tol
                        near = min(((np.linalg.norm(h[hn][0] - x), x) for kk, a, x in acceptors if kk != key), key=lambda t: t[0], default=(9, None))
                        if near[0] < 2.5:
                            T["acc_n"] += 1
                            T["acc_ok"] += np.linalg.norm(o[hn] - near[1]) < 2.7
    pct = lambda a, b: f"{T[a]}/{T[b]} = {100 * T[a] / max(1, T[b]):.1f}%"
    print(f"histidine tautomer (both neutral)          {pct('taut_ok', 'taut_n')}")
    print(f"Asn/Gln/His turned away from the neutron orientation   {pct('orient_wrong', 'orient_n')}")
    det = T["flip_n"] - T["flip_und"]
    print(f"input turned by 180°: restored             {pct('flip_ok', 'flip_n')}; reported undetermined {T['flip_und']}; "
          f"of the {det} others restored {T['flip_det_ok']} = {100 * T['flip_det_ok'] / max(1, det):.1f}%")
    print(f"Ser/Thr hydroxyl within 45° of the neutron hydrogen   {pct('oh_ok', 'oh_n')}")
    print(f"Tyr hydroxyl on the same side              {pct('tyr_ok', 'tyr_n')}")
    print(f"hydroxyl donates to the same acceptor      {pct('acc_ok', 'acc_n')}")


def formats(pdb_dir, cif_dir):
    bad = 0
    for f in sorted(os.listdir(pdb_dir)):
        cif = os.path.join(cif_dir, f[:-4] + ".cif")
        if not f.lower().endswith(".pdb") or not os.path.exists(cif):
            continue
        a, b = (json.loads(subprocess.run([CL, "protein", p, "--json"], capture_output=True, text=True).stdout) for p in (os.path.join(pdb_dir, f), cif))
        keys = ("chains", "residues", "disulfides", "chain_breaks", "atoms", "hydrogens_added", "net_charge", "flipped", "undetermined", "repair")
        same = all(a[k] == b[k] for k in keys) and all(abs(x["pKa"] - y["pKa"]) < 1e-6 for x, y in zip(a["pka"], b["pka"]))
        bad += not same
        print(f"  {'ok ' if same else 'BAD'} {f[:-4]}: {a['atoms']} atoms, {len(a['pka'])} titratable groups")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    mode, args = sys.argv[1], sys.argv[2:]
    {"sidechains": lambda: sidechains(args), "loops": lambda: loops(args), "neutron": lambda: neutron(args),
     "formats": lambda: formats(*args)}[mode]()
