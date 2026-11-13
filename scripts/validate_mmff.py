"""MMFF94 validation (run after building build/chemlarp).

1. The official MMFF94 validation suite (Halgren / Merck, CCL archive www.ccl.net/cca/data/MMFF94/,
   MMFF94.tar.gz): atom types, partial charges and every energy term of 759 molecules against Merck's
   OPTIMOL output (MMFF94_opti.log). Pass the directory holding the unpacked suite:
       python3 scripts/validate_mmff.py --suite path/to/MMFF94
   The suite is not redistributed with ChemLarp (its structures derive from the Cambridge Structural
   Database).
2. RDKit's independent MMFF94 implementation on drug-like molecules entered as SMILES (always run; needs
   rdkit): ChemLarp embeds and minimizes each molecule, then both programs evaluate the same coordinates.
"""
import argparse, json, os, re, subprocess, sys, tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CL = os.path.join(ROOT, "build", "chemlarp")
failures = []


def check(label, ok, text):
    print(f"  {'ok ' if ok else 'BAD'} {label:46s} {text}")
    if not ok:
        failures.append(label)


# ------------------------------------------------------------------------------------------- suite
def suite(d):
    def charges(path):
        out, cur = {}, None
        for l in open(path, encoding="latin-1"):
            m = re.match(r"Molecule (\S+)", l)
            if m:
                cur = m.group(1); out[cur] = {}
            m = re.match(r"Atom (\S+), formal charge = (-?\d+)", l)
            if m and cur:
                out[cur][m.group(1)] = int(m.group(2))
        return out

    def to_sdf(mol2, fc):
        recs, out = open(mol2, encoding="latin-1").read().split("@<TRIPOS>MOLECULE")[1:], []
        for r in recs:
            lines = r.strip("\n").split("\n")
            name = lines[0].strip()
            atoms_txt = r.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().split("\n")
            bonds_txt = r.split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0].strip().split("\n")
            atoms = []
            for l in atoms_txt:
                f = l.split()
                el = f[5].split(".")[0]
                atoms.append((f[1], el[0] + el[1:].lower(), float(f[2]), float(f[3]), float(f[4])))
            bonds = [(int(f[1]), int(f[2]), {"1": 1, "2": 2, "3": 3, "am": 1, "ar": 4}[f[3]]) for f in (l.split() for l in bonds_txt)]
            bos = [0] * len(atoms)
            for a, b, o in bonds:
                bos[a - 1] += o; bos[b - 1] += o
            m = [name, "", "", f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
            for k, (_, el, x, y, z) in enumerate(atoms):
                m.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {el:<3} 0  0  0  0  0{(bos[k] or 15):3d}  0  0  0  0  0  0")
            m += [f"{a:3d}{b:3d}{o:3d}  0" for a, b, o in bonds]
            ch = [(i + 1, fc.get(name, {}).get(a[0], 0)) for i, a in enumerate(atoms)]
            ch = [c for c in ch if c[1]]
            for k in range(0, len(ch), 8):
                part = ch[k:k + 8]
                m.append("M  CHG" + f"{len(part):3d}" + "".join(f" {i:3d} {c:3d}" for i, c in part))
            m += ["M  END", "$$$$"]
            out.append("\n".join(m))
        return "\n".join(out) + "\n"

    def precise_coords():
        lines, out, i = open(f"{d}/MMFF94.mmd", encoding="latin-1").read().split("\n"), {}, 0
        while i < len(lines):
            m = re.match(r"\s*(\d+)\s+\[([^,]+),", lines[i])
            if m:
                na = int(m.group(1))
                out[m.group(2)] = [float(v) for l in lines[i + 1:i + 1 + na] for v in l.split()[13:16]]
                i += na + 1
            else:
                i += 1
        return out

    def reference():
        text = open(f"{d}/MMFF94_opti.log", encoding="latin-1").read()
        ref = {}
        for b in re.split(r"Structure Name(?:/Conformational Index)?: ", text)[1:]:
            name = b.split()[0]

            def table(tag, conv):
                m = re.search(r"# " + tag + r"\s*\n(.*?)\n OPTIMOL-LIST>", b, re.S)
                vals = {int(mm.group(2)): conv(mm.group(3)) for mm in re.finditer(r"(\S+) #\s*(\d+)\s+(\S+)", m.group(1))}
                return [vals[i] for i in sorted(vals)]
            e = {}
            for key, lab in [("total", r"Total ENERGY \(Kcal\)"), ("bond", "Bond Stretching"), ("angle", "Angle Bending"),
                             ("out_of_plane", "Out-of-Plane Bending"), ("stretch_bend", "Stretch-Bend"),
                             ("torsion", "Total Torsion"), ("vdw", "Net vdW"), ("electrostatic", "Electrostatic")]:
                m = re.search(lab + r"\s+(-?[\d.]+)", b)
                e[key] = float(m.group(1))
            ref[name] = {"types": table("ty", int), "charges": table("char", float), "energy": e}
        return ref

    ref, coords = reference(), precise_coords()
    print(f"1. MMFF94 validation suite: {len(ref)} molecules (OPTIMOL reference)")
    for form, mol2, fcf in (("hypervalent", "MMFF94_hypervalent.mol2", "MMFF94.fc_hypervalent"),
                            ("dative", "MMFF94_dative.mol2", "MMFF94.fc_dative")):
        sdf = to_sdf(f"{d}/{mol2}", charges(f"{d}/{fcf}"))
        with tempfile.TemporaryDirectory() as tmp:
            open(f"{tmp}/s.sdf", "w").write(sdf)
            open(f"{tmp}/c.txt", "w").write("".join(n + " " + " ".join(map(str, v)) + "\n" for n, v in coords.items()))
            out = json.loads(subprocess.run([CL, "mmff", "--file", f"{tmp}/s.sdf", "--coords", f"{tmp}/c.txt"],
                                            capture_output=True, text=True).stdout)
        n = atoms = tbad = 0
        dq, de = [], {}
        for o in out:
            if "types" not in o:
                continue
            r = ref[o["name"]]
            n += 1
            atoms += len(r["types"])
            tbad += sum(a != b for a, b in zip(o["types"], r["types"]))
            dq.append((max(abs(a - b) for a, b in zip(o["charges"], r["charges"])), o["name"]))
            for k, v in r["energy"].items():
                de.setdefault(k, []).append((abs(o["energy"][k] - v), o["name"]))
        print(f" {form}: {n} molecules read, {atoms} atoms")
        check(f"{form}: atom types", tbad == 0, f"{tbad} of {atoms} differ")
        badq = [x for x in dq if x[0] > 0.0011]
        # FEZPOP: the suite's formal-charge file gives +1 instead of -1 to two oxygens of an anion.
        check(f"{form}: partial charges (to the 3 printed decimals)", [x[1] for x in badq] == ["FEZPOP"],
              f"{len(badq)} molecules differ: {[x[1] for x in badq]}")
        tot = sorted(de["total"])
        exact = sum(1 for x in tot if x[0] <= 1e-4)
        close = [x for x in tot if x[0] > 6e-3 and x[1] != "FEZPOP"]
        check(f"{form}: total energy", not close, f"{exact} of {n} within 1e-4 kcal/mol, worst others "
              f"{[(round(x[0], 4), x[1]) for x in tot[-4:]]}")
        for k in ("bond", "angle", "stretch_bend", "out_of_plane", "torsion", "vdw", "electrostatic"):
            worst = max(x for x in de[k] if x[1] != "FEZPOP")
            check(f"{form}: {k}", worst[0] < 6e-3, f"largest difference {worst[0]:.4f} ({worst[1]})")


# ------------------------------------------------------------------------------------------- RDKit
SMILES = [
    "CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(=O)Nc1ccc(O)cc1",
    "OC(=O)CCC(=O)O", "C[NH3+]", "CC(=O)[O-]", "NC(=[NH2+])N", "c1ccncc1", "c1ccc2[nH]ccc2c1", "c1cnc[nH]1",
    "c1ccoc1", "c1ccsc1", "O=[N+]([O-])c1ccccc1", "CS(=O)(=O)N", "CS(=O)C", "COP(=O)(OC)OC", "C#N", "CC#N", "C=C=C",
    "C1CC1", "C1CCC1", "C1CCCC1", "C1CCCCC1", "c1ccc2ccccc2c1", "OCC(O)CO", "NCC(=O)O", "[NH3+]CC(=O)[O-]",
    "CC[N+](C)(C)C", "Clc1ccccc1", "Brc1ccccc1", "Fc1ccccc1", "Ic1ccccc1", "CC(F)(F)F", "CCOC(=O)C", "CC(=O)N(C)C",
    "O=C1NC(=O)C=C1", "c1ccc(cc1)c1ccccc1", "CCN(CC)CC", "C[S-]", "CC(=S)N", "S=C(N)N", "CN=C=O", "CC=NO",
    "c1ccc[n+](C)c1", "[O-][n+]1ccccc1", "CN1CCC[C@H]1c1cccnc1", "OC(=O)c1ccccc1O", "Oc1ccccc1", "CC(C)(C)O",
    "C[Si](C)(C)C", "CSSC", "CSC", "CS", "O=C(O)C(F)(F)F", "N#Cc1ccccc1", "c1ccc2c(c1)[nH]c1ccccc12", "CC1=CC(=O)CC1",
    "O=c1cc[nH]cc1", "c1cc2ccc3cccc4ccc(c1)c2c34", "CC(C)NCC(O)COc1cccc2ccccc12", "COc1ccc2[nH]cc(CCN(C)C)c2c1",
    "CN1CCN(CC1)C(=O)c1ccccc1", "c1ccc(nc1)N", "Nc1ncnc2[nH]cnc12", "O=C1CCCN1", "OC1CCCCC1", "CC(=O)C(C)=O",
    "c1ccc2nc3ccccc3cc2c1", "C1=CC=CC=C1C=O", "OS(=O)(=O)O", "CS(=O)(=O)[O-]", "C[P+](C)(C)C", "N=[N+]=[N-]",
    "CC(=O)OO", "C(=O)(N)N", "O=C(N)c1cccnc1", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "c1coc2ccccc12", "C1COCCO1",
]


def rdkit_compare():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdForceFieldHelpers as H, rdMolTransforms
    from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties, MMFFGetMoleculeForceField
    RDLogger.DisableLog("rdApp.*")
    print(f"2. RDKit MMFF94 on {len(SMILES)} molecules from SMILES (same coordinates)")
    tbad = qbad = ebad = n = 0
    worst = (0.0, "")
    for smi in SMILES:
        r = json.loads(subprocess.run([CL, "mmff", smi, "--json"], capture_output=True, text=True).stdout)
        if not r["ok"]:
            failures.append(smi)
            print("  BAD could not set up", smi, r["diagnostics"][:1])
            continue
        m = Chem.MolFromMolBlock(r["sdf"], removeHs=False)
        props = MMFFGetMoleculeProperties(m)
        if props is None:
            print("  RDKit has no MMFF setup for", smi)
            continue
        n += 1
        types = [props.GetMMFFAtomType(i) for i in range(m.GetNumAtoms())]
        charges = [props.GetMMFFPartialCharge(i) for i in range(m.GetNumAtoms())]
        e = MMFFGetMoleculeForceField(m, props, nonBondedThresh=1000.0, ignoreInterfragInteractions=False).CalcEnergy()
        if types != r["types"]:
            tbad += 1
            print("   types differ:", smi, [(i, a, b) for i, (a, b) in enumerate(zip(r["types"], types)) if a != b])
        dq = max(abs(a - b) for a, b in zip(r["charges"], charges))
        if smi == "C[P+](C)(C)C":
            # RDKit gives the phosphonium ion no formal charge (every partial charge 0, net 0); ChemLarp
            # keeps the +1 on P. Not covered by the official suite.
            print(f"   known difference: {smi} net charge ChemLarp {sum(r['charges']):+.0f}, RDKit {sum(charges):+.0f}")
            continue
        if dq > 1e-4:
            qbad += 1
            print(f"   charges differ by {dq:.4f}:", smi)
        # The SDF coordinates have four decimals, so energies agree to about 1e-3 kcal/mol.
        de = abs(e - r["energy"]["total"])
        worst = max(worst, (de, smi))
        if de > 0.01:
            ebad += 1
            print(f"   energy differs: ChemLarp {r['energy']['total']:.4f} RDKit {e:.4f}", smi)
    check("RDKit: atom types", tbad == 0, f"{tbad} of {n} molecules differ")
    check("RDKit: partial charges", qbad == 0, f"{qbad} of {n} molecules differ by more than 1e-4")
    check("RDKit: total energy at the same coordinates", ebad == 0, f"largest difference {worst[0]:.4f} kcal/mol ({worst[1]})")


ap = argparse.ArgumentParser()
ap.add_argument("--suite", help="directory with the unpacked MMFF94 validation suite")
a = ap.parse_args()
if a.suite:
    suite(a.suite)
rdkit_compare()
if failures:
    sys.exit("failed: " + "; ".join(failures))
print("\nall checks passed")
