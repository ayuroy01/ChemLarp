#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Redocking of public crystal complexes (build first; needs the PDB entries and their ligands).

    python3 scripts/validate_docking.py COMPLEX_DIR [--runs N] [--complexes 1STP,3PTB,...] [--rerun]

COMPLEX_DIR holds, for each entry, the deposited structure `<ID>.pdb` (or `.cif`) and its ligand as
`<ID>_<CODE>.sdf` with the crystal coordinates. Both come from RCSB: the structure from
files.rcsb.org/download/<ID>.pdb and the ligand from the model server,
models.rcsb.org/v1/<id>/ligand?auth_seq_id=<n>&encoding=sdf. Nothing is downloaded here.

The test is the standard one for a docking program. The crystal ligand gives only its chemical identity: the
SMILES is taken from the deposited ligand and a conformation is generated from scratch, so the crystal
geometry is not an input. That molecule is then docked into the receptor with the crystal ligand's own site
cleared, and the top-ranked pose is compared with the deposited pose by heavy-atom RMSD with the atoms
matched by graph isomorphism (symmetric groups do not inflate it). Under 2 A counts as reproducing the
crystal pose, the usual threshold.

A failure is split into its two causes, which need different work:

  search    the crystal pose, minimised where it stands, scores better than the best pose the search found,
            so the function knows the right answer and the search did not reach it;
  scoring   the pose found scores better than the crystal pose, so the search did its job and the function
            prefers somewhere else.

Each complex is docked twice, into the structure as deposited and again with its waters taken out, because
the difference is large and is worth reporting: a resolved water is part of the site, and taking it away
leaves a cavity that the diffuse part of the score sinks a ligand into. The pass/fail criterion is the
structure as deposited, which is also the engine's default.

There is no comparison against AutoDock Vina's own numbers here: that program has no macOS build on PyPI and
its binaries were not downloaded. The reference is the crystal pose, which is an experiment.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import chemlarp  # noqa: E402

failures = []

# PDB entry -> the ligand's chemical component code. A spread of sizes and binding sites: a vitamin in
# streptavidin, small inhibitors in trypsin and thrombin, a cofactor-bound reductase, a cytochrome, a sugar
# transporter, two HIV protease inhibitors, an acetylcholinesterase inhibitor and a kinase inhibitor.
COMPLEXES = [
    ("1STP", "BTN"),  # streptavidin + biotin
    ("3PTB", "BEN"),  # beta-trypsin + benzamidine
    ("4DFR", "MTX"),  # dihydrofolate reductase + methotrexate
    ("2CPP", "CAM"),  # cytochrome P450cam + camphor
    ("1ABE", "ARA"),  # L-arabinose binding protein + arabinose
    ("1HSG", "MK1"),  # HIV-1 protease + indinavir
    ("1ACJ", "THA"),  # acetylcholinesterase + tacrine
    ("1HVR", "XK2"),  # HIV-1 protease + XK263
    ("1DWB", "BEN"),  # thrombin + benzamidine
    ("2BRB", "PFQ"),  # PDE4 + a furanone inhibitor
    ("1PPH", "0ZG"),  # trypsin + a peptidyl inhibitor
    ("1ETS", "MID"),  # thrombin + NAPAP
]

GOOD_RMSD = 2.0         # A: the pose counts as the crystal pose
SUCCESS_FRACTION = 0.5  # of the complexes that ran

# The four this engine is documented not to reproduce (docs/ARCHITECTURE.md, "Docking: validation"), with
# the RMSD recorded there. They are checked against that number rather than against 2 A, so the suite is a
# guard against losing ground without pretending the module gets them right.
EXPECTED_MISS = {"4DFR": 3.38, "1HVR": 3.80, "1DWB": 2.35, "1PPH": 4.26}
MISS_SLACK = 1.0  # A of room above the recorded value, for the scatter between seeds


def check(label, ok, text):
    print("  %s %-52s %s" % ("ok " if ok else "BAD", label, text), flush=True)
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


def receptor_path(directory, pdb_id):
    for ext in (".pdb", ".cif", ".ent"):
        p = os.path.join(directory, pdb_id + ext)
        if os.path.exists(p):
            return p
    return None


def redock(directory, pdb_id, code, runs, seed, remove_waters):
    """Dock the deposited ligand's molecule back into its own structure."""
    structure = receptor_path(directory, pdb_id)
    ligand_file = os.path.join(directory, "%s_%s.sdf" % (pdb_id, code))
    if not structure or not os.path.exists(ligand_file):
        return None
    with open(structure) as f:
        receptor = f.read()
    with open(ligand_file) as f:
        crystal = f.read()
    # The identity of the ligand, and nothing else, comes from the deposited coordinates: its canonical
    # SMILES, from which a conformation is then generated.
    records = chemlarp.read_molecules(crystal)
    if not records or not records[0].get("ok"):
        return None
    smiles = records[0]["smiles"]
    r = chemlarp.dock(receptor, smiles, reference=crystal, runs=runs, seed=seed, margin=4.0,
                      remove_waters=remove_waters)
    best = r["poses"][0]
    return dict(pdb=pdb_id, code=code, smiles=smiles, heavy_atoms=r["heavy_atoms"],
                rotatable_bonds=r["rotatable_bonds"], receptor_atoms=r["receptor_atoms"], seconds=r["seconds"],
                affinity=best["affinity"], rmsd=best["rmsd"], internal=best["internal"],
                rmsd_any=min([p["rmsd"] for p in r["poses"] if p["rmsd"] is not None] or [None]),
                poses=len(r["poses"]), reference=r.get("reference", {}),
                waters_removed=r["receptor"]["waters_removed"],
                box=[round(v, 1) for v in r["box"]["size"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory", help="directory with <ID>.pdb and <ID>_<CODE>.sdf for each complex")
    ap.add_argument("--runs", type=int, default=16, help="independent Monte Carlo runs per complex")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--complexes", default="", help="comma-separated PDB ids, default all of them")
    ap.add_argument("--rerun", action="store_true", help="ignore stored results")
    a = ap.parse_args()
    runs_dir = os.path.join(a.directory, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    wanted = [c.strip().upper() for c in a.complexes.split(",") if c.strip()]
    chosen = [(i, c) for i, c in COMPLEXES if not wanted or i in wanted]
    print("threads", chemlarp.threads(), "; %d runs per complex; results kept in %s" % (a.runs, runs_dir))
    print("\nA  Redocking into the structure as deposited: the top-ranked pose against the deposited pose")
    print("   (the last column is the same run with the structure's waters taken out)\n")

    rows, stripped, missing = [], {}, []
    for pdb_id, code in chosen:
        results = {}
        for no_water in (False, True):
            key = os.path.join(runs_dir, "%s_%s_r%d_s%d%s.json" % (
                pdb_id, code, a.runs, a.seed, "_nowat" if no_water else ""))
            try:
                results[no_water] = cached(
                    key, a.rerun, lambda nw=no_water: redock(a.directory, pdb_id, code, a.runs, a.seed, nw))
            except chemlarp.ChemLarpError as e:
                check("%s %s: docked" % (pdb_id, code), False, str(e)[:70])
                results[no_water] = None
        r, dry = results.get(False), results.get(True)
        if r is None:
            if dry is None:
                missing.append(pdb_id)
            continue
        rows.append(r)
        if dry is not None:
            stripped[pdb_id] = dry
        # The crystal pose minimised where it stands: which half of the job failed, when one did.
        ref = r.get("reference") or {}
        ref_affinity = ref.get("affinity")
        cause = ""
        if r["rmsd"] is not None and r["rmsd"] > GOOD_RMSD and ref_affinity is not None:
            cause = " (search)" if r["affinity"] > ref_affinity else " (scoring)"
        print("  %-5s %-4s %2d heavy %2d rot  %4d waters  affinity %7.2f  crystal %7.2f  RMSD %6.2f A%-10s"
              "  no water %6.2f" % (
                  r["pdb"], r["code"], r["heavy_atoms"], r["rotatable_bonds"],
                  (dry or {}).get("waters_removed", 0), r["affinity"],
                  ref_affinity if ref_affinity is not None else float("nan"), r["rmsd"], cause,
                  (dry or {}).get("rmsd", float("nan"))))
        if pdb_id in EXPECTED_MISS:
            limit = EXPECTED_MISS[pdb_id] + MISS_SLACK
            check("%s %s: no worse than the %.2f A recorded for it" % (pdb_id, code, EXPECTED_MISS[pdb_id]),
                  r["rmsd"] is not None and r["rmsd"] <= limit,
                  "%.2f A in %.1f s (this one is a known miss)" % (r["rmsd"], r["seconds"]))
        else:
            check("%s %s: top pose within %.1f A of the crystal pose" % (pdb_id, code, GOOD_RMSD),
                  r["rmsd"] is not None and r["rmsd"] <= GOOD_RMSD,
                  "%.2f A in %.1f s" % (r["rmsd"], r["seconds"]))

    if missing:
        print("\n  no files for: %s" % ", ".join(missing))
    if not rows:
        print("\nnothing ran: the directory holds no <ID>.pdb with its <ID>_<CODE>.sdf")
        return 1

    good = [r for r in rows if r["rmsd"] is not None and r["rmsd"] <= GOOD_RMSD]
    rmsds = sorted(r["rmsd"] for r in rows if r["rmsd"] is not None)
    median = rmsds[len(rmsds) // 2] if len(rmsds) % 2 else 0.5 * (rmsds[len(rmsds) // 2 - 1] + rmsds[len(rmsds) // 2])
    print("\nB  The set as a whole\n")
    print("  %d of %d complexes within %.1f A; median RMSD %.2f A; %.0f s in all" % (
        len(good), len(rows), GOOD_RMSD, median, sum(r["seconds"] for r in rows)))
    searched = [r for r in rows if r["rmsd"] > GOOD_RMSD and (r.get("reference") or {}).get("affinity") is not None
                and r["affinity"] > (r["reference"] or {})["affinity"]]
    scored = [r for r in rows if r["rmsd"] > GOOD_RMSD and (r.get("reference") or {}).get("affinity") is not None
              and r["affinity"] <= (r["reference"] or {})["affinity"]]
    if searched or scored:
        print("  of the %d that missed: %d the search did not reach (%s), %d the function prefers elsewhere (%s)" % (
            len(searched) + len(scored), len(searched), ", ".join(r["pdb"] for r in searched) or "-",
            len(scored), ", ".join(r["pdb"] for r in scored) or "-"))
    check("at least %.0f%% of the complexes within %.1f A" % (100 * SUCCESS_FRACTION, GOOD_RMSD),
          len(good) >= math.ceil(SUCCESS_FRACTION * len(rows)), "%d of %d" % (len(good), len(rows)))
    if stripped:
        dry_good = set(p for p, d in stripped.items() if d["rmsd"] is not None and d["rmsd"] <= GOOD_RMSD)
        dry_rmsds = sorted(d["rmsd"] for d in stripped.values() if d["rmsd"] is not None)
        dry_median = (dry_rmsds[len(dry_rmsds) // 2] if len(dry_rmsds) % 2 else
                      0.5 * (dry_rmsds[len(dry_rmsds) // 2 - 1] + dry_rmsds[len(dry_rmsds) // 2]))
        print("  with the waters taken out: %d of %d within %.1f A, median %.2f A" % (
            len(dry_good), len(stripped), GOOD_RMSD, dry_median))
        lost = sorted(set(r["pdb"] for r in good) - dry_good)
        gained = sorted(dry_good - set(r["pdb"] for r in good))
        if lost:
            print("  taking the waters out loses the crystal pose for: %s" % ", ".join(lost))
        if gained:
            print("  and finds it only without them for: %s" % ", ".join(gained))
    # Any pose in the list, not only the top one: a weaker claim about the function, a stronger one about
    # whether the right basin was visited at all.
    any_good = [r for r in rows if r["rmsd_any"] is not None and r["rmsd_any"] <= GOOD_RMSD]
    print("  the crystal pose is somewhere in the list for %d of %d" % (len(any_good), len(rows)))
    # Every crystal pose must survive being minimised where it stands: a few tenths of an angstrom, not a
    # slide into something else. That checks the local minimiser and the receptor preparation, not the search.
    moved = [r for r in rows if (r.get("reference") or {}).get("rmsd") is not None and r["reference"]["rmsd"] > 1.0]
    check("every crystal pose stays put when minimised", not moved,
          "worst %.2f A" % max([(r.get("reference") or {}).get("rmsd") or 0.0 for r in rows]))

    print("\n%s" % ("all checks passed" if not failures else "%d failed: %s" % (len(failures), ", ".join(failures))))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
