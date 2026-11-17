#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Explicit-water simulation and hydration free energies (build first; needs FreeSolv's database.txt).

    python3 scripts/validate_hydration.py FREESOLV_DATABASE_TXT [--ps P] [--molecules a,b,...] [--rerun]

Every simulation is slow (bulk water: a few minutes; one hydration free energy: about a quarter of an hour on
8 cores), so results are kept as JSON next to the FreeSolv file (runs/) and reused unless --rerun is given.

A  Bulk TIP3P at 298.15 K and 1 atm, two seeds, against the published TIP3P liquid with Ewald and a
   long-range Lennard-Jones correction: density 0.982 g/cm3 and potential energy -9.66 kcal/mol at 300 K
   (Vega, Abascal, Conde & Aragones, Faraday Discuss. 141 (2009) 251, table 6). The model, not water, is the
   reference here: TIP3P is known to be 1.5 % too light.
B  A TIP3P molecule decoupled from TIP3P water (the model's own hydration free energy): the estimators must
   agree and neighbouring windows overlap. The value is printed beside the published TIP3P results, -5.4 to
   -5.6 kcal/mol (-22.6 and -23.43 kJ/mol, quoted in arXiv 0809.1628), but those used cut-off electrostatics
   and 216 molecules, so they are not a pass/fail reference; experiment is -6.3.
C  Small molecules from FreeSolv (Mobley & Guthrie, J. Comput.-Aided Mol. Des. 28 (2014) 711; database.txt,
   v0.52) against experiment, with FreeSolv's own GAFF/AM1-BCC/TIP3P values alongside. ChemLarp's solute model
   (MMFF94 charges and intramolecular terms, UFF Lennard-Jones) is not GAFF, so it is judged against
   experiment with a tolerance at the level of standard fixed-charge force fields (GAFF's error on FreeSolv
   is about 1.5 kcal/mol RMS), not against GAFF's numbers.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import chemlarp  # noqa: E402

failures = []
DEFAULT_MOLECULES = ["methane", "ethane", "methanol", "ethanol", "acetone", "benzene", "phenol", "acetamide"]


def check(label, ok, text):
    print("  %s %-44s %s" % ("ok " if ok else "BAD", label, text), flush=True)
    if not ok:
        failures.append(label)


def cached(path, rerun, compute):
    if not rerun and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    r = compute()
    with open(path, "w") as f:
        json.dump(r, f)
    return r


def read_freesolv(path):
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = [x.strip() for x in line.split(";")]
            rows[p[2]] = dict(id=p[0], smiles=p[1], name=p[2], expt=float(p[3]), expt_err=float(p[4]),
                              gaff=float(p[5]), gaff_err=float(p[6]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("freesolv", help="FreeSolv database.txt")
    ap.add_argument("--ps", type=float, default=50.0, help="production per window (ps)")
    ap.add_argument("--molecules", default=",".join(DEFAULT_MOLECULES), help="FreeSolv names, comma-separated")
    ap.add_argument("--rerun", action="store_true", help="ignore stored results")
    a = ap.parse_args()
    runs = os.path.join(os.path.dirname(os.path.abspath(a.freesolv)), "runs")
    os.makedirs(runs, exist_ok=True)
    print("threads", chemlarp.threads(), "; results kept in", runs)

    print("\nA  bulk TIP3P, 400 molecules, 298.15 K, 1 atm, 20 + 100 ps")
    boxes = []
    for seed in (1, 2):
        r = cached(os.path.join(runs, "water_box_%d.json" % seed), a.rerun,
                   lambda: chemlarp.water_box(400, 298.15, 1.0, 20.0, 100.0, seed))
        boxes.append(r)
        print("  seed %d: density %.4f +- %.4f g/cm3, U %.3f +- %.3f kcal/mol, kinetic T %.1f K" % (
            seed, r["density"], r["density_error"], r["potential_per_molecule"], r["potential_error"], r["temperature"]))
    rho = sum(b["density"] for b in boxes) / len(boxes)
    u = sum(b["potential_per_molecule"] for b in boxes) / len(boxes)
    check("density vs TIP3P 0.982 g/cm3", abs(rho - 0.982) < 0.01, "%.4f" % rho)
    check("potential energy vs TIP3P -9.66 kcal/mol", abs(u + 9.66) < 0.15, "%.3f" % u)
    peak = max(boxes[0]["g_oo"], key=lambda p: p[1])
    check("first O-O peak 2.75 +- 0.1 A", abs(peak[0] - 2.75) < 0.1, "%.3f A, g = %.2f" % (peak[0], peak[1]))

    def hydrate(name, smiles, seed=1):
        key = "%s_%gps_seed%d" % (name.replace(" ", "_"), a.ps, seed)
        return cached(os.path.join(runs, "hydration_%s.json" % key), a.rerun,
                      lambda: chemlarp.hydration_free_energy(smiles, 400, 298.15, 20.0, a.ps, seed))

    def show(name, r):
        d = r["dg_hyd"]
        print("  %-12s MBAR %6.2f +- %.2f  BAR %6.2f  TI %6.2f  (charging %6.2f, vdW %5.2f)  overlap %.2f" % (
            name, d["mbar"]["dg"], d["mbar"]["error"], d["bar"]["dg"], d["ti"]["dg"], r["dg_elec"], r["dg_vdw"],
            r["min_neighbour_overlap"]))

    def estimators_agree(name, r):
        d = r["dg_hyd"]
        spread = max(abs(d["mbar"]["dg"] - d["bar"]["dg"]), abs(d["mbar"]["dg"] - d["ti"]["dg"]))
        check(name + ": MBAR, BAR and TI agree", spread < 0.5, "max difference %.2f kcal/mol" % spread)
        check(name + ": neighbouring windows overlap", r["min_neighbour_overlap"] > 0.03,
              "smallest %.3f" % r["min_neighbour_overlap"])

    print("\nB  TIP3P water in TIP3P water")
    w = hydrate("water", "water", seed=2)
    show("water", w)
    estimators_agree("water", w)
    g = w["dg_hyd"]["mbar"]["dg"]
    # Not a pass/fail check: the published TIP3P values were computed without Ewald summation (cut-off
    # electrostatics, 216 molecules), and no PME value was found to compare with (docs/GAPS.md).
    print("  --  water vs cut-off TIP3P results -5.4 to -5.6, experiment -6.3: %.2f (not comparable setups)" % g)

    print("\nC  FreeSolv small molecules (kcal/mol)")
    db = read_freesolv(a.freesolv)
    names = [n for n in a.molecules.split(",") if n]
    errs, gaff_errs = [], []
    for n in names:
        if n not in db:
            check(n + ": in FreeSolv", False, "not found")
            continue
        row = db[n]
        r = hydrate(n, row["smiles"])
        show(n, r)
        estimators_agree(n, r)
        g = r["dg_hyd"]["mbar"]["dg"]
        errs.append(g - row["expt"])
        gaff_errs.append(row["gaff"] - row["expt"])
        print("  %-12s experiment %6.2f +- %.2f   FreeSolv GAFF %6.2f   ChemLarp - experiment %+.2f" % (
            "", row["expt"], row["expt_err"], row["gaff"], g - row["expt"]))
        check(n + ": within 2.5 kcal/mol of experiment", abs(g - row["expt"]) < 2.5, "%+.2f" % (g - row["expt"]))
    if errs:
        rms = math.sqrt(sum(e * e for e in errs) / len(errs))
        mean = sum(errs) / len(errs)
        grms = math.sqrt(sum(e * e for e in gaff_errs) / len(gaff_errs))
        print("  %d molecules: ChemLarp RMS error %.2f (mean %+.2f); FreeSolv GAFF on the same set %.2f" % (
            len(errs), rms, mean, grms))
        check("RMS error vs experiment below 2 kcal/mol", rms < 2.0, "%.2f" % rms)

    print("\n%s" % ("all checks passed" if not failures else "%d failed: %s" % (len(failures), ", ".join(failures))))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
