"""Protein handling validation (needs Biopython; run after building build/chemlarp).

    python3 scripts/validate_protein.py <dir with PDB files>

The directory holds PDB-format files downloaded from RCSB (https://files.rcsb.org/download/<ID>.pdb); the
results in docs/ARCHITECTURE.md used 1L2Y, 2JOF, 1VII, 1D3Z and 2K39 (NMR models with hydrogens) and 1CRN,
1UBQ, 1AKI and 4HHB (X-ray structures).

A  Hydrogens: every structure with deposited hydrogens is protonated by `chemlarp protein` (all hydrogens
   removed and rebuilt) and each rebuilt hydrogen is compared with the deposited one of the same name. Methyl
   and NH3+ groups are compared as a set (their rotation is free); hydroxyl and thiol hydrogens, placed by a
   hydrogen-bond search, are reported separately. A methylene pair whose names are swapped would show as a
   distance near 1.8 Å and is counted.
B  Sequences, isoelectric points and net charges against Biopython (PDB parser, IsoelectricPoint, the same
   Bjellqvist pKa set); molecular weights against ProtParam with ChemLarp's residue formulas evaluated at
   Biopython's atomic weights (the two programs use different atomic-weight tables).
"""
import itertools, json, math, os, re, subprocess, sys, warnings
from collections import defaultdict

warnings.simplefilter("ignore")
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

HERE = os.path.dirname(os.path.abspath(__file__))
CL = os.path.join(HERE, "..", "build", "chemlarp")
failures = []


def check(label, ok, text):
    print(f"  {'ok ' if ok else 'BAD'} {label:44s} {text}")
    if not ok:
        failures.append(label)


def first_model_atoms(text):
    res = {}
    for l in text.splitlines():
        if l.startswith("ENDMDL"):
            break
        if l[:6] in ("ATOM  ", "HETATM"):
            key = (l[21], int(l[22:26]), l[26])
            el = l[76:78].strip() or l[12:16].strip()[0]
            res.setdefault(key, {})[l[12:16].strip()] = (float(l[30:38]), float(l[38:46]), float(l[46:54]), el, l[17:20])
    return res


def protonate(path, pH=7.4):
    p = subprocess.run([CL, "protein", path, "--pH", str(pH), "--json"], capture_output=True, text=True)
    return json.loads(p.stdout)


rotors = {("SER", "HG"), ("THR", "HG1"), ("TYR", "HH"), ("CYS", "HG")}
pdb_dir = sys.argv[1]
files = sorted(f for f in os.listdir(pdb_dir) if f.lower().endswith(".pdb"))

print("A  Rebuilt hydrogens against deposited NMR hydrogens")
for f in files:
    text = open(os.path.join(pdb_dir, f)).read()
    ref = first_model_atoms(text)
    if not any(a[3] == "H" for r in ref.values() for a in r.values()):
        continue
    r = protonate(os.path.join(pdb_dir, f))
    ours = first_model_atoms(r["pdb"])
    d = lambda a, b: math.dist(a[:3], b[:3])
    fixed, rot, swaps, pairs, tautomer = [], [], 0, 0, []
    for key, ra in ref.items():
        oa = ours.get(key, {})
        rh = {n: v for n, v in ra.items() if v[3] == "H"}
        oh = {n: v for n, v in oa.items() if v[3] == "H"}
        if set(rh) != set(oh):
            tautomer.append(f"{next(iter(ra.values()))[4]}{key[1]}")
        groups = defaultdict(list)
        for n in rh:
            if n in oh:
                groups[n[:-1] if n[-1].isdigit() else n].append(n)
        for g, ns in groups.items():
            resn = rh[ns[0]][4]
            if len(ns) >= 3:
                best = min(itertools.permutations(ns), key=lambda p: sum(d(rh[a], oh[b]) for a, b in zip(ns, p)))
                fixed += [d(rh[a], oh[b]) for a, b in zip(ns, best)]
            elif len(ns) == 2:
                a, b = ns
                pairs += 1
                swaps += d(rh[a], oh[b]) + d(rh[b], oh[a]) < d(rh[a], oh[a]) + d(rh[b], oh[b]) - 0.5
                fixed += [d(rh[a], oh[a]), d(rh[b], oh[b])]
            else:
                (rot if (resn, ns[0]) in rotors else fixed).append(d(rh[ns[0]], oh[ns[0]]))
    rms = math.sqrt(sum(x * x for x in fixed) / len(fixed))
    within = sum(x < 0.3 for x in fixed) / len(fixed)
    check(f"{f}: {len(fixed)} H, rms", rms < 0.35, f"{rms:.3f} Å, {100 * within:.1f}% within 0.3 Å, max {max(fixed):.2f}")
    check(f"{f}: methylene names", swaps == 0, f"{swaps} of {pairs} pairs swapped")
    if rot:
        print(f"      hydroxyl/thiol rotors: {len(rot)}, mean {sum(rot) / len(rot):.2f} Å (orientation not compared)")
    if tautomer:
        print(f"      different hydrogen sets (histidine tautomer or protonation choice): {', '.join(tautomer)}")

print("B  Sequences, pI, charge and mass against Biopython")
for f in files:
    path = os.path.join(pdb_dir, f)
    st = PDBParser(QUIET=True).get_structure(f, path)
    bio = ["".join(seq1(r.get_resname()) for r in ch if r.id[0] == " ") for ch in st[0]]
    bio = [s for s in bio if s]
    r = protonate(path)
    ours = [c["sequence"] for c in r["chains"]]
    check(f"{f}: sequences", ours == bio, f"{len(ours)} chain(s), {sum(map(len, ours))} residues")
    for c in r["chains"]:
        ip = IsoelectricPoint(c["sequence"])
        ref_pi = ip.pi()
        # Biopython's search is bounded below by pH 4.05; compare only inside its range.
        if ref_pi > 4.06:
            check(f"{f} {c['id']}: pI", abs(ref_pi - c["pI"]) < 2e-3, f"{c['pI']:.4f} vs {ref_pi:.4f}")

seqs = ["ACDEFGHIKLMNPQRSTVWY", "EKKK", "SDEEKHC", "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAENCDKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDSISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"]
worst = 0.0
for s in seqs:
    lines = "".join(f"ATOM  {i + 1:5d}  CA  {r:3s} A{i + 1:4d}    {i * 3.8:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           C\n"
                    for i, r in enumerate({v: k for k, v in __import__('Bio.Data.IUPACData', fromlist=['x']).protein_letters_3to1.items()}[c].upper() for c in s))
    # Net charge from the sequence at several pH values through the titration curve of `chemlarp protein`.
    tmp = os.path.join(pdb_dir, "_seq.pdb.tmp")
    open(tmp, "w").write(lines)
    r = protonate(tmp)
    os.remove(tmp)
    ip = IsoelectricPoint(s)
    for ph, q in r["titration"]:
        worst = max(worst, abs(q - ip.charge_at_pH(ph)))
check("net charge, pH 0-14, 4 sequences", worst < 1e-4, f"largest difference {worst:.1e}")

src = open(os.path.join(HERE, "..", "src", "protein.cpp")).read()
from Bio.Data import IUPACData as I
w = I.atom_weights
bad = []
for m in re.finditer(r"\{\"([A-Z]{3})\", \{'([A-Z])',.*?\{(\d+), (\d+), (\d+), (\d+), (\d+)\}\}\}", src, re.S):
    f = list(map(int, m.groups()[2:]))
    ours = sum(n * w[e] for n, e in zip(f, "CHNOS"))
    if abs(ours - (I.protein_weights[m.group(2)] - 2 * w["H"] - w["O"])) > 1e-3:
        bad.append(m.group(2))
check("residue formulas vs ProtParam weights", not bad, f"20 residues, mismatches: {bad or 'none'}")

print("\nall checks passed" if not failures else f"\n{len(failures)} FAILED: {failures}")
sys.exit(1 if failures else 0)
