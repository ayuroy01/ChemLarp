"""Weights of the structure-based pKa model in src/protein_pka.inc (needs numpy; run after building).

    python3 scripts/fit_protein_pka.py <pkad.json> <dir with PDB files> [--write] [--propka propka3] [--cache dir]

pkad.json holds the measured pKa values of wild-type proteins from PKAD-2 (Ancona, Bastola & Alexov,
Database 2023, baad033; http://compbio.clemson.edu/PKAD-2) as a list of {"pdb", "res", "chain", "num",
"pka"}; the directory holds the PDB files of those entries from RCSB. For every measured group the features
come from `chemlarp protein <file> --json` ("pka"[i]["features"]).

Model: pKa − reference = offset(group) + w(class) · features, classes carboxyl (ASP, GLU, C-terminus),
HIS and amine (LYS, N-terminus); ridge regression. TYR and CYS keep their reference values: the 25 and 15
measured values are too few, and the cysteines are mostly reactive active-site residues (mean pKa 6.5),
so an offset fitted to them would wrongly deprotonate ordinary cysteines. The error is estimated by 10-fold
cross-validation in which all structures of one protein (same sequence of the first chain) stay together,
and compared with the null model (reference values only) and, with --propka, with PROPKA 3 on the same
groups. --write replaces src/protein_pka.inc with weights fitted to all data.
"""
import json, math, os, subprocess, sys, tempfile
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CL = os.path.join(HERE, "..", "build", "chemlarp")
GROUPS = ["ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG", "N+", "C-"]
CLASS = {"ASP": 0, "GLU": 0, "C-": 0, "HIS": 1, "LYS": 2, "N+": 2, "TYR": 3, "CYS": 3}
NF = 10
RIDGE = 3.0


def sites(path):
    p = subprocess.run([CL, "protein", path, "--json"], capture_output=True, text=True)
    r = json.loads(p.stdout)
    if not r.get("ok"):
        return None, None
    return r["pka"], r["chains"][0]["sequence"] if r["chains"] else ""


def propka(exe, path, cache):
    out = os.path.join(cache, os.path.basename(path)[:-4] + ".pka")
    if not os.path.exists(out):
        subprocess.run([exe, "-q", os.path.abspath(path)], cwd=cache, capture_output=True)
    res = {}
    if not os.path.exists(out):
        return res
    on = False
    for l in open(out):
        if l.startswith("SUMMARY OF THIS PREDICTION"):
            on = True
            continue
        if on and len(l) > 20 and l[:6].strip() in ("ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG", "N+", "C-"):
            t = l.split()
            try:
                if float(t[3]) < 50.0:  # 99.99 marks a group PROPKA does not titrate (e.g. a disulfide)
                    res[(t[0], t[2], int(t[1]))] = float(t[3])
            except ValueError:
                pass
    return res


def fit(rows, idx):
    """Per-class ridge fit on rows[idx]; returns offsets by group and weights by class."""
    off, w = {}, {}
    for c in range(3):  # class 3 (TYR, CYS) is not fitted: see the module text
        sel = [i for i in idx if CLASS[rows[i]["group"]] == c]
        groups = sorted({rows[i]["group"] for i in sel})
        if not sel:
            continue
        X = np.array([rows[i]["f"] + [1.0 if rows[i]["group"] == g else 0.0 for g in groups] for i in sel])
        y = np.array([rows[i]["pka"] - rows[i]["ref"] for i in sel])
        reg = np.diag([RIDGE] * NF + [1e-6] * len(groups))
        beta = np.linalg.solve(X.T @ X + reg, X.T @ y)
        w[c] = beta[:NF]
        for g, b in zip(groups, beta[NF:]):
            off[g] = b
    return off, w


def predict(row, off, w):
    c = CLASS[row["group"]]
    if c not in w:
        return row["ref"]
    return row["ref"] + off.get(row["group"], 0.0) + float(np.dot(w[c], row["f"]))


def main():
    data = json.load(open(sys.argv[1]))
    pdb_dir = sys.argv[2]
    exe = sys.argv[sys.argv.index("--propka") + 1] if "--propka" in sys.argv else None
    by_pdb = defaultdict(list)
    for d in data:
        by_pdb[d["pdb"].upper()].append(d)
    rows, skipped = [], 0
    cache = sys.argv[sys.argv.index("--cache") + 1] if "--cache" in sys.argv else tempfile.mkdtemp()
    os.makedirs(cache, exist_ok=True)
    for pdb in sorted(by_pdb):
        path = os.path.join(pdb_dir, pdb + ".pdb")
        if not os.path.exists(path):
            skipped += len(by_pdb[pdb])
            continue
        ss, seq = sites(path)
        if ss is None:
            skipped += len(by_pdb[pdb])
            continue
        pk = propka(exe, path, cache) if exe else {}
        index = defaultdict(list)
        for s in ss:
            index[(s["chain"], s["number"], s["group"])].append(s)
        for d in by_pdb[pdb]:
            g = {"N-TERM": "N+", "C-TERM": "C-"}.get(d["res"], d["res"])
            try:
                num = int(d["num"])
            except ValueError:
                skipped += 1
                continue
            hit = index.get((d["chain"], num, g))
            if not hit and g in ("N+", "C-"):  # termini are listed with the residue name or number of the chain end
                hit = [s for s in ss if s["chain"] == d["chain"] and s["group"] == g]
            if not hit or g not in CLASS:
                skipped += 1
                if os.environ.get("SHOW"): print("unmatched", pdb, d)
                continue
            s = hit[0]
            rows.append({"pdb": pdb, "seq": seq, "group": g, "pka": float(d["pka"]), "ref": s["reference"] if not WRITTEN else REFS[g],
                         "f": s["features"], "propka": pk.get((g, d["chain"], s["number"]))})
    print(f"{len(rows)} measured groups in {len({r['pdb'] for r in rows})} structures ({skipped} not matched)")
    proteins = sorted({r["seq"] for r in rows})
    fold_of = {s: i % 10 for i, s in enumerate(proteins)}
    cv = np.zeros(len(rows))
    for k in range(10):
        train = [i for i, r in enumerate(rows) if fold_of[r["seq"]] != k]
        off, w = fit(rows, train)
        for i, r in enumerate(rows):
            if fold_of[r["seq"]] == k:
                cv[i] = predict(r, off, w)
    y = np.array([r["pka"] for r in rows])
    null = np.array([r["ref"] for r in rows])
    pk = np.array([r["propka"] if r["propka"] is not None else np.nan for r in rows])
    rm = lambda a, m: math.sqrt(np.nanmean((a[m] - y[m]) ** 2)) if m.any() else float("nan")
    ma = lambda a, m: float(np.nanmean(np.abs(a[m] - y[m]))) if m.any() else float("nan")
    print(f"{len(proteins)} distinct proteins; 10-fold cross-validation by protein")
    print(f"{'group':6s} {'n':>5s}  {'null rmse':>9s} {'cv rmse':>8s} {'cv mae':>7s}  {'PROPKA n':>8s} {'rmse':>6s} {'mae':>6s} {'ours same n':>11s}")
    for g in GROUPS + ["all"]:
        m = np.array([g == "all" or r["group"] == g for r in rows])
        if not m.any():
            continue
        mp = m & ~np.isnan(pk)
        print(f"{g:6s} {m.sum():5d}  {rm(null, m):9.3f} {rm(cv, m):8.3f} {ma(cv, m):7.3f}  {mp.sum():8d} {rm(pk, mp):6.3f} {ma(pk, mp):6.3f} {rm(cv, mp):11.3f}")
    big = np.abs(y - null) > 1.0
    print(f"groups shifted by more than 1 unit from the reference: {big.sum()}; cv rmse {rm(cv, big):.3f}, null {rm(null, big):.3f}, PROPKA {rm(pk, big & ~np.isnan(pk)):.3f}")
    off, w = fit(rows, list(range(len(rows))))
    lines = ["// Generated by scripts/fit_protein_pka.py from %d measured pKa values (PKAD-2). Do not edit." % len(rows),
             "// Weights by class (carboxyl, HIS, amine, TYR/CYS) for the features: buried fraction, local density,",
             "// backbone NH, backbone O, side-chain NH, hydroxyl, amide O, Coulomb exposed, Coulomb buried, charged contacts.",
             "static const double kPkaWeights[4][kPkaFeatures] = {"]
    for c in range(4):
        lines.append("    {" + ", ".join("%.4f" % v for v in w.get(c, np.zeros(NF))) + "},")
    lines += ["};", "// Offset of the reference pKa by group: " + " ".join(GROUPS),
              "static const double kPkaOffset[9] = {" + ", ".join("%.4f" % off.get(g, 0.0) for g in GROUPS) + "};"]
    if "--write" in sys.argv:
        open(os.path.join(HERE, "..", "src", "protein_pka.inc"), "w").write("\n".join(lines) + "\n")
        print("wrote src/protein_pka.inc")
    else:
        print("\n".join(lines))


# Reference values in water (PROPKA 3 model values); used instead of the JSON "reference", which already
# includes the fitted offset once src/protein_pka.inc has been written.
REFS = {"ASP": 3.80, "GLU": 4.50, "HIS": 6.50, "CYS": 9.00, "TYR": 10.00, "LYS": 10.50, "ARG": 12.50, "N+": 8.00, "C-": 3.20}
WRITTEN = True
main()
