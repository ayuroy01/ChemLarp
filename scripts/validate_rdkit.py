#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Structure identity, CIP labels, SMARTS counts and fingerprints against RDKit (needs rdkit; build first).

    python3 scripts/validate_rdkit.py <smiles file> [--limit N] [--examples N]

The file holds one SMILES per line; a leading numeric field (a PubChem CID, say) is skipped. RDKit is an
independent implementation, so the comparisons below test ChemLarp's reading, canonical form, CIP ranking,
substructure matching and fingerprints against another program rather than against values kept here.

A  Identity of the canonical form: ChemLarp's canonical SMILES and the input, both read by RDKit, must give
   the same RDKit canonical SMILES (stereo included). This checks that the canonical form denotes the input
   structure; that it is unique per structure is checked by the random re-spellings in the C++ tests.
   Differences are split into constitution (which atom is bonded to which, and the hydrogen count) and
   stereochemistry only. A constitution difference is always a defect, so the check fails on the first one.
   Stereo-only differences are counted against a rate, because two classes of them are known and documented
   in docs/GAPS.md: stereo markers on three-coordinate P, As, N and S, which ChemLarp does not read, and a
   double-bond configuration that can be dropped when a mancude macrocycle is written in its aromatic form.
B  CIP labels: the multiset of descriptors per molecule (R, S, r, s, E, Z) against RDKit's rdCIPLabeler
   (Hanson et al., J. Chem. Inf. Model. 58 (2018) 1755). Compared as a multiset, not atom by atom, because
   the two programs number atoms differently; a swap between two atoms of the same molecule would be missed.
   Molecules whose input carries a stereo marker on a hetero centre are counted separately, for the same
   reason as in A.
C  SMARTS: unique match counts for 40 patterns (ring, valence, charge, recursive, aromatic, isotope and
   degree primitives) on every molecule.
D  Fingerprints: Morgan radius 2 Tanimoto on 2,000 random pairs. ChemLarp hashes environments with its own
   32-bit function, so individual bits differ from RDKit's by construction; what must agree is the
   similarity itself (mean absolute difference and the ordering of the pairs).
"""
import argparse
import os
import random
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import chemlarp  # noqa: E402

from rdkit import Chem, DataStructs, RDLogger  # noqa: E402
from rdkit.Chem import rdCIPLabeler, rdFingerprintGenerator  # noqa: E402

RDLogger.DisableLog("rdApp.*")

PATTERNS = [
    "[CX4]", "[CX3]=[OX1]", "[#6][CX3](=O)[#6]", "[OX2H]", "[NX3;H2,H1;!$(NC=O)]", "[R2]", "[r5]", "[x3]",
    "[v5]", "[h2]", "[D3;!R]", "[*;R1]@[*;R1]", "[c,n;H1]", "a:a:a", "[!#6;!#1]", "[+,-]", "[N+](=O)[O-]",
    "C=,#C", "[$(C=O)]O", "[R]!@[R]", "*~*~*~*", "[Cl,Br,I]", "[#7;a]", "[CH3]", "[CH2;R]", "c1ccccc1",
    "C1CCCCC1", "[OH]c1ccccc1", "[S;X4](=O)(=O)", "[2H]", "[13C]", "[c;r6;x2]", "[C;!$(C=O)]", "[n;H0;X2]",
    "[a;R2]", "[#6]-,:[#7]", "[!R;D1]", "O=[C,S]", "[NX3;H0;!$(NC=O)]", "[$([CX3]=[CX3])]",
]

# A stereo marker on a three-coordinate hetero centre: ChemLarp does not read these (docs/GAPS.md).
HETERO_STEREO = re.compile(r"\[(?:As|Se|Si|Sn|Te|Ge|[PNSB])(?:H\d?)?@")

failures = []


def check(label, ok, text):
    print("  %s %-44s %s" % ("ok " if ok else "BAD", label, text))
    if not ok:
        failures.append(label)


def read_smiles(path, limit):
    out = []
    for line in open(path):
        f = line.split()
        if not f:
            continue
        s = f[1] if len(f) > 1 and f[0].isdigit() else f[0]
        out.append(s)
        if limit and len(out) >= limit:
            break
    return out


def cip_multiset(labels):
    d = {}
    for x in labels:
        if x:
            d[x] = d.get(x, 0) + 1
    return d


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            mean = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = mean
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den if den else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles_file")
    ap.add_argument("--limit", type=int, default=1500)
    ap.add_argument("--examples", type=int, default=5)
    a = ap.parse_args()

    smis = read_smiles(a.smiles_file, a.limit)
    mols, keep = {}, []
    for i, s in enumerate(smis):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            mols[i] = m
            keep.append(i)
    print("%d SMILES, %d read by RDKit (%d not, and not compared)" % (len(smis), len(keep), len(smis) - len(keep)))

    # A  identity of the canonical form
    print("\nA  canonical SMILES denotes the input structure (judged by RDKit)")
    same = 0
    ours_failed = []
    constitution, stereo_hetero, stereo_other = [], [], []
    for i in keep:
        try:
            r = chemlarp.canonical(smis[i])
        except Exception:
            ours_failed.append(smis[i])
            continue
        c = Chem.MolFromSmiles(r["canonical"])
        if c is None:
            constitution.append((smis[i], r["canonical"], "RDKit cannot read our canonical SMILES"))
            continue
        if Chem.MolToSmiles(c) == Chem.MolToSmiles(mols[i]):
            same += 1
        elif Chem.MolToSmiles(c, isomericSmiles=False) != Chem.MolToSmiles(mols[i], isomericSmiles=False):
            constitution.append((smis[i], r["canonical"], Chem.MolToSmiles(c)))
        elif HETERO_STEREO.search(smis[i]):
            stereo_hetero.append((smis[i], r["canonical"], Chem.MolToSmiles(c)))
        else:
            stereo_other.append((smis[i], r["canonical"], Chem.MolToSmiles(c)))
    n = len(keep) - len(ours_failed)
    nstereo = len(stereo_hetero) + len(stereo_other)
    check("constitution preserved", not constitution,
          "%d/%d molecules; %d with a different constitution" % (n - len(constitution), n, len(constitution)))
    check("stereochemistry preserved", n and nstereo <= 0.02 * n,
          "%d/%d (%.2f%%); %d stereo on hetero centres, %d other" %
          (n - len(constitution) - nstereo, n, 100.0 * (n - len(constitution) - nstereo) / max(n, 1),
           len(stereo_hetero), len(stereo_other)))
    if ours_failed:
        print("      %d not read by ChemLarp: %s" % (len(ours_failed), ", ".join(ours_failed[: a.examples])))
    for name, group in (("constitution", constitution), ("stereo, hetero centre", stereo_hetero),
                        ("stereo, other", stereo_other)):
        for s, c, got in group[: a.examples]:
            print("      [%s] input %s\n        ours %s\n        RDKit reads ours as %s" % (name, s, c, got))

    # B  CIP labels
    print("\nB  CIP descriptors against rdCIPLabeler (multiset per molecule)")
    with_stereo = agree = 0
    units_ours = units_both = 0
    bad = []
    for i in keep:
        try:
            cip = chemlarp.canonical(smis[i])["cip"]
        except Exception:
            continue
        ours = [c["label"] for c in cip["centres"]] + [c["label"] for c in cip["double_bonds"]]
        ours = [x for x in ours if x]
        m = Chem.Mol(mols[i])
        try:
            rdCIPLabeler.AssignCIPLabels(m)
        except Exception:
            continue
        ref = [at.GetProp("_CIPCode") for at in m.GetAtoms() if at.HasProp("_CIPCode")]
        for b in m.GetBonds():
            s = str(b.GetStereo())
            if s in ("STEREOE", "STEREOTRANS"):
                ref.append("E")
            elif s in ("STEREOZ", "STEREOCIS"):
                ref.append("Z")
        if not ours and not ref:
            continue
        with_stereo += 1
        a_, b_ = cip_multiset(ours), cip_multiset(ref)
        units_ours += sum(a_.values())
        units_both += sum(min(v, b_.get(k, 0)) for k, v in a_.items())
        if a_ == b_:
            agree += 1
        else:
            bad.append((smis[i], a_, b_))
    hetero = [x for x in bad if HETERO_STEREO.search(x[0])]
    other = [x for x in bad if not HETERO_STEREO.search(x[0])]
    if with_stereo == 0:
        print("      no molecule in this sample carries stereochemistry: nothing compared")
    else:
        check("same descriptors", len(other) <= 0.01 * with_stereo,
              "%d/%d molecules with stereo (%.2f%%), %d/%d units; %d differ on a hetero centre, %d otherwise" %
              (agree, with_stereo, 100.0 * agree / with_stereo, units_both, units_ours, len(hetero), len(other)))
    for name, group in (("hetero centre", hetero), ("other", other)):
        for s, x, y in group[: a.examples]:
            print("      [%s] %s\n        ours %s  RDKit %s" % (name, s, x, y))

    # C  SMARTS
    print("\nC  SMARTS unique match counts, %d patterns" % len(PATTERNS))
    queries = [Chem.MolFromSmarts(p) for p in PATTERNS]
    for p, q in zip(PATTERNS, queries):
        assert q is not None, p
    pairs = agree_pairs = 0
    worst = {}
    for i in keep[: min(len(keep), 400)]:
        for p, q in zip(PATTERNS, queries):
            ref = len(mols[i].GetSubstructMatches(q, uniquify=True, maxMatches=5000))
            try:
                got = len(chemlarp.smarts_match(smis[i], p, max_matches=5000))
            except Exception:
                got = -1
            pairs += 1
            if got == ref:
                agree_pairs += 1
            else:
                worst[p] = worst.get(p, 0) + 1
    check("same counts", pairs and agree_pairs >= 0.999 * pairs,
          "%d/%d pattern-molecule pairs (%.3f%%)" % (agree_pairs, pairs, 100.0 * agree_pairs / max(pairs, 1)))
    if worst:
        print("      patterns that differ:", dict(sorted(worst.items(), key=lambda kv: -kv[1])[: a.examples]))

    # D  fingerprints
    print("\nD  Morgan radius 2 Tanimoto on random pairs")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = {i: gen.GetFingerprint(mols[i]) for i in keep}
    random.seed(5)
    ours, theirs = [], []
    for _ in range(2000):
        i, j = random.sample(keep, 2) if len(keep) > 1 else (keep[0], keep[0])
        try:
            o = chemlarp.similarity(smis[i], smis[j], radius=2, bits=2048)["tanimoto"]
        except Exception:
            continue
        ours.append(o)
        theirs.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    if ours:
        mean_abs = sum(abs(x - y) for x, y in zip(ours, theirs)) / len(ours)
        worst = max(abs(x - y) for x, y in zip(ours, theirs))
        rho = spearman(ours, theirs)
        # The rank correlation is reported but not a criterion: a sample of enumerated stereoisomers holds
        # many pairs with identical fingerprints, and those ties pull Spearman down whatever the agreement.
        check("similarity agrees", mean_abs <= 0.02 and worst <= 0.25,
              "%d pairs, mean |difference| %.4f, largest %.3f, Spearman %.4f (ties included)" %
              (len(ours), mean_abs, worst, rho))
    else:
        check("similarity agrees", False, "no pairs scored")

    if failures:
        sys.exit("outside tolerance: " + "; ".join(failures))
    print("\nall checks passed")


if __name__ == "__main__":
    main()
