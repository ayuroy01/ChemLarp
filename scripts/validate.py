#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""One runner for every validation suite in the project, with a report.

    python3 scripts/validate.py --list
    python3 scripts/validate.py [--only id,id] [--skip id,id] [--data DIR] [--config FILE]
                                [--python PATH] [--report out.md] [--json out.json]

Each suite is one of the scripts under scripts/ (or the C++ test binary) run as a subprocess; this file only
decides what can run, runs it, and writes the report. Nothing is downloaded: a suite whose reference data,
Python packages or external tools are missing is skipped and the report says what was missing and where to
get it. Exit status is 0 only if every suite that ran passed.

Two kinds of suite, kept apart in the report because they carry different weight:

  independent  compared with another implementation of the same model (PHREEQC, RDKit, CDK, Biopython,
               thermo, scipy), with published measurements, or with a definition computed a second way.
  internal     in-repo expectations: the C++ test suite. Its reference values come from outside (JANAF
               tables, Merck's MMFF94 output, crystal structures), but the comparison is frozen in the
               repository rather than recomputed, so a wrong value would stay wrong.

The report also lists the roadmap modules with no independent suite at all, and the validations that
docs/ARCHITECTURE.md describes but that no script here reproduces (see NOT_SCRIPTED).

Reference data is not redistributed. Point the runner at a directory with this layout (--data, or the
environment variable CHEMLARP_VALIDATION_DATA):

    pdb/        PDB files with deposited hydrogens (1L2Y 2JOF 1VII 1D3Z 2K39) and X-ray (1CRN 1UBQ 1AKI 4HHB)
    cif/        the mmCIF file of each entry in pdb/, same names
    heldout/    PDB files not used to fit anything: 1BPI 2TRX 1IGD 1POH 2CPL 1HRC 4ICB 1LZ1 7RSA 1XNB
    neutron/    neutron/X-ray structures as mmCIF: 1L2K 1VCX 2ZOI 2ZWB 3A1R 3INS 3KYY 3X2P 4AR3 4FC1 4NY6
                5CE4 5MOQ 5PTI 5ZN0 6EYM 6F2Q 6L46 7AZ9 7VOS
    mmff94/     the unpacked MMFF94 validation suite from www.ccl.net/cca/data/MMFF94/ (MMFF94.tar.gz)
    jars/       cdk-*.jar, centres-core-*.jar, centres-cdk-*.jar
    freesolv/   database.txt from github.com/MobleyLab/FreeSolv (runs/ is created beside it)
    smiles.txt  one SMILES per line (first whitespace-separated field), e.g. a PubChem sample

or write the paths explicitly in a JSON file (--config, or validation.json in the data directory):

    {"python": "/path/to/venv/bin/python", "pdb_dir": "...", "mmff94_dir": "...", "cdk_jar": "..."}

The Python packages the suites need (phreeqpython, thermo, scipy, rdkit, biopython, pymbar) usually live in a
virtual environment outside the repository: --python points at that interpreter.
"""
import argparse
import glob
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Paths each suite may ask for: key -> (what it is, default location under --data)
PATHS = {
    "pdb_dir": ("PDB files with deposited hydrogens", "pdb"),
    "cif_dir": ("mmCIF files of the same entries", "cif"),
    "heldout_dir": ("PDB files held out of every fit", "heldout"),
    "neutron_dir": ("neutron structures (mmCIF)", "neutron"),
    "mmff94_dir": ("unpacked MMFF94 validation suite", "mmff94"),
    "cdk_jar": ("CDK jar", "jars/cdk*.jar"),
    "centres_core_jar": ("centres-core jar", "jars/centres-core*.jar"),
    "centres_cdk_jar": ("centres-cdk jar", "jars/centres-cdk*.jar"),
    "smiles_file": ("SMILES list, one per line", "smiles.txt"),
    "freesolv_file": ("FreeSolv database.txt (github.com/MobleyLab/FreeSolv)", "freesolv/database.txt"),
    "dock_dir": ("crystal complexes: <ID>.pdb with its ligand as <ID>_<CODE>.sdf", "dock"),
}

# Validations docs/ARCHITECTURE.md reports that no script here reproduces.
NOT_SCRIPTED = [
    ("InChI and InChIKey", "1483/1500 identical to PubChem", "needs the PubChem sample and its InChIs"),
    ("IUPAC naming", "14,188/14,193 names parsed back by OPSIN", "needs OPSIN (Java) and the PubChem sample"),
    ("Quantum chemistry", "HF/DFT energies against PySCF on 19 molecules x 4 basis sets", "needs pyscf"),
    ("2D layout", "every stereocentre's wedge read back from the drawing", "done inside the C++ tests only"),
    ("Structure editor", "scripted browser checks", "needs a browser harness"),
    ("Readers", "same structure as RDKit on 999 PubChem 2D and 865 3D files", "needs those file sets"),
    ("Pure-component properties", "constants and correlations against DIPPR/thermo", "not kept as a script"),
    ("Equations of state, NRTL/UNIFAC", "against thermo and DDBst data", "not kept as a script"),
]


def suites():
    """The registry. argv(cfg) returns the command; cfg holds `python` and the resolved paths."""
    return [
        dict(
            id="tests",
            title="C++ test suite",
            kind="internal",
            modules=["*"],
            reference="in-repo expectations and bundled reference data (JANAF, MMFF reference SDF, COD structures)",
            needs=dict(build=["chemlarp_tests"]),
            argv=lambda c: [os.path.join(ROOT, "build", "chemlarp_tests")],
            parse="tests",
        ),
        dict(
            id="aqueous",
            title="Aqueous speciation, saturation indices, solubility",
            kind="independent",
            modules=["Solution chemistry (pH, ions, solubility)"],
            reference="PHREEQC (phreeqpython) on the same database; measured pKw, calcite and gypsum solubility",
            needs=dict(python=["phreeqpython"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_aqueous.py"],
        ),
        dict(
            id="column",
            title="Distillation columns",
            kind="independent",
            modules=["Distillation columns"],
            reference="the same MESH equations solved by scipy with thermo's correlations and UNIFAC",
            needs=dict(python=["numpy", "scipy", "thermo"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_column.py"],
        ),
        dict(
            id="flowsheet",
            title="Flowsheet solver",
            kind="independent",
            modules=["Flowsheet solver", "Unit operations (mixers, exchangers, reactors)"],
            reference="closed-form balances (Felder & Rousseau) and scipy/thermo equation-oriented solutions",
            needs=dict(python=["numpy", "scipy", "thermo"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_flowsheet.py"],
        ),
        dict(
            id="md",
            title="Molecular dynamics",
            kind="independent",
            modules=["Molecular dynamics"],
            reference="an independent Monte Carlo sampler of the same ensemble; normal modes; equipartition",
            needs=dict(python=["numpy"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_md.py"],
        ),
        dict(
            id="mmff-rdkit",
            title="MMFF94 against RDKit",
            kind="independent",
            modules=["MMFF94, partial charges, electrostatics"],
            reference="RDKit's MMFF94 on the same coordinates",
            needs=dict(python=["rdkit"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_mmff.py"],
        ),
        dict(
            id="mmff-suite",
            title="MMFF94 against Merck's OPTIMOL output (759 molecules)",
            kind="independent",
            modules=["MMFF94, partial charges, electrostatics"],
            reference="the official MMFF94 validation suite (Halgren/Merck, CCL archive)",
            needs=dict(python=["rdkit"], build=["chemlarp"], paths=["mmff94_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_mmff.py", "--suite", c["mmff94_dir"]],
        ),
        dict(
            id="protein-h",
            title="Rebuilt hydrogens, sequences, pI and mass",
            kind="independent",
            modules=["Protein handling (PDB, residues, protonation)"],
            reference="deposited NMR hydrogens; Biopython (PDB parser, IsoelectricPoint, ProtParam)",
            needs=dict(python=["Bio"], build=["chemlarp"], paths=["pdb_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_protein.py", c["pdb_dir"]],
        ),
        dict(
            id="protein-sidechains",
            title="Side-chain rebuilding (chi1/chi2 against the crystal)",
            kind="independent",
            modules=["Protein handling (PDB, residues, protonation)"],
            reference="X-ray side chains deleted and rebuilt; structures held out of the rotamer library",
            needs=dict(python=["Bio", "numpy"], build=["chemlarp"], paths=["heldout_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_protein_build.py", "sidechains"]
            + sorted(glob.glob(os.path.join(c["heldout_dir"], "*.pdb"))),
        ),
        dict(
            id="protein-loops",
            title="Loop rebuilding (backbone rmsd, closing bonds)",
            kind="independent",
            modules=["Protein handling (PDB, residues, protonation)"],
            reference="stretches of 3-12 residues deleted from held-out X-ray structures and rebuilt",
            needs=dict(python=["Bio", "numpy"], build=["chemlarp"], paths=["heldout_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_protein_build.py", "loops"]
            + sorted(glob.glob(os.path.join(c["heldout_dir"], "*.pdb")))[:4],
        ),
        dict(
            id="protein-neutron",
            title="Protonation against neutron structures",
            kind="independent",
            modules=["Protein handling (PDB, residues, protonation)"],
            reference="experimental hydrogen (deuterium) positions: His tautomers, amide flips, OH donors",
            needs=dict(python=["Bio", "numpy"], build=["chemlarp"], paths=["neutron_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_protein_build.py", "neutron"]
            + sorted(glob.glob(os.path.join(c["neutron_dir"], "*.cif"))),
        ),
        dict(
            id="protein-formats",
            title="PDB and mmCIF of the same entry give the same result",
            kind="independent",
            modules=["Protein handling (PDB, residues, protonation)", "Readers for MOL/SDF, PDB, CML"],
            reference="the two deposited formats of each entry",
            needs=dict(python=["Bio"], build=["chemlarp"], paths=["pdb_dir", "cif_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_protein_build.py", "formats", c["pdb_dir"], c["cif_dir"]],
        ),
        dict(
            id="fep",
            title="Free energies: MBAR, BAR, TI, alchemical LJ particle",
            kind="independent",
            modules=["Free-energy perturbation (MBAR, BAR, TI; alchemical LJ)"],
            reference="exact oscillator free energies; pymbar on the same matrices; Widom insertion in the pure fluid",
            needs=dict(python=["numpy", "pymbar"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_fep.py"],
        ),
        dict(
            id="hydration",
            title="Explicit TIP3P water and hydration free energies (hours)",
            kind="independent",
            modules=["Hydration free energies in explicit water (FreeSolv)"],
            reference="published TIP3P liquid (Vega et al. 2009) and TIP3P self-hydration; FreeSolv experiment",
            needs=dict(build=["chemlarp"], paths=["freesolv_file"]),
            argv=lambda c: [c["python"], "scripts/validate_hydration.py", c["freesolv_file"]],
        ),
        dict(
            id="docking",
            title="Redocking crystal complexes (minutes)",
            kind="independent",
            modules=["Docking"],
            reference="the deposited poses of twelve RCSB PDB complexes (heavy-atom RMSD under 2 A)",
            needs=dict(build=["chemlarp"], paths=["dock_dir"]),
            argv=lambda c: [c["python"], "scripts/validate_docking.py", c["dock_dir"]],
        ),
        dict(
            id="lab",
            title="Virtual lab: mixtures, titrations, heat of neutralisation",
            kind="independent",
            modules=["Interactive lab simulations and safety rules"],
            reference="PHREEQC (MIX + EQUILIBRIUM_PHASES) on the same solutions; CODATA heat of neutralisation",
            needs=dict(python=["phreeqpython"], build=["chemlarp"]),
            argv=lambda c: [c["python"], "scripts/validate_lab.py"],
        ),
        dict(
            id="axial",
            title="Axial chirality (allenes, atropisomers)",
            kind="independent",
            modules=["Axial chirality (allenes, atropisomers)"],
            reference="CDK with the centres CIP labeller (M/P), on SMILES and on 3D files",
            needs=dict(tools=["java", "c++"], build=["libchemlarp.a"],
                       paths=["cdk_jar", "centres_core_jar", "centres_cdk_jar"]),
            argv=lambda c: ["sh", "scripts/validate_axial/run.sh", c["cdk_jar"], c["centres_core_jar"], c["centres_cdk_jar"]],
            parse="axial",
        ),
        dict(
            id="rdkit-structure",
            title="Structure identity, CIP labels, SMARTS counts, fingerprints",
            kind="independent",
            # Not "InChI and InChIKey; CIP rules 1b, 3-5 (r/s)": the CIP part of that module is exercised
            # here, but nothing reproduces the InChI comparison with PubChem, so it counts as uncovered.
            modules=["Canonical SMILES and CIP R/S, E/Z", "Substructure search (SMARTS) and fingerprints",
                     "SMILES parser, valence checks, error vectors"],
            reference="RDKit (canonical SMILES, rdCIPLabeler, substructure counts, Morgan fingerprints)",
            needs=dict(python=["rdkit"], build=["chemlarp"], paths=["smiles_file"]),
            argv=lambda c: [c["python"], "scripts/validate_rdkit.py", c["smiles_file"]],
        ),
    ]


def resolve_paths(cfg, data_dir):
    """Fill missing path keys from the documented layout under `data_dir` (globs allowed)."""
    if not data_dir:
        return
    for key, (_, rel) in PATHS.items():
        if cfg.get(key):
            continue
        pattern = os.path.join(data_dir, rel)
        hits = sorted(glob.glob(pattern))
        if hits:
            cfg[key] = hits[0]
        elif "*" not in rel and os.path.exists(pattern):
            cfg[key] = pattern


def have_module(python, name):
    r = subprocess.run([python, "-c", "import " + name], capture_output=True)
    return r.returncode == 0


def missing(suite, cfg, modcache):
    """Everything the suite needs and does not have, as readable phrases."""
    out = []
    needs = suite["needs"]
    for b in needs.get("build", []):
        if not os.path.exists(os.path.join(ROOT, "build", b)):
            out.append("build/" + b + " (run make)")
    for t in needs.get("tools", []):
        if not shutil.which(t):
            out.append("the " + t + " command")
    for m in needs.get("python", []):
        if m not in modcache:
            modcache[m] = have_module(cfg["python"], m)
        if not modcache[m]:
            out.append("the " + m + " package for " + cfg["python"])
    for k in needs.get("paths", []):
        p = cfg.get(k)
        if not p or not os.path.exists(p):
            out.append(PATHS[k][0] + " (" + k + ")")
    return out


def parse_output(text, how):
    """(checks ok, checks total, headline). total 0 means the suite prints no per-check lines."""
    ok = len(re.findall(r"^\s*ok\s", text, re.M))
    bad = len(re.findall(r"^\s*BAD\s", text, re.M))
    if how == "tests":
        m = re.search(r"(\d+) tests, (\d+) failed", text)
        if m:
            return int(m.group(1)) - int(m.group(2)), int(m.group(1)), m.group(0)
    lines = [l.strip() for l in text.strip().splitlines() if l.strip()]
    # Headline: the script's own closing summary when it has one, else its last line.
    for l in reversed(lines):
        if l.startswith(("all checks passed", "largest differences", "worst ")):
            return ok, ok + bad, l
    return ok, ok + bad, lines[-1] if lines else ""


def run_suite(suite, cfg, timeout):
    argv = suite["argv"](cfg)
    t0 = time.time()
    try:
        p = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        out, code = p.stdout + p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or "") + (e.stderr or "") if isinstance(e.stdout, str) else ""
        out += "\n[timed out after %d s]" % timeout
        code = 124
    ok, total, headline = parse_output(out, suite.get("parse"))
    return dict(status="passed" if code == 0 else "failed", exit_code=code, seconds=round(time.time() - t0, 1),
                checks_ok=ok, checks_total=total, headline=headline, output=out, command=" ".join(argv))


def roadmap_coverage(results):
    """module name -> list of independent suite ids that cover it."""
    try:
        d = json.load(open(os.path.join(ROOT, "data", "roadmap.json")))
    except OSError:
        return []
    cover = {}
    for s in suites():
        if s["kind"] != "independent":
            continue
        for m in s["modules"]:
            cover.setdefault(m, []).append(s["id"])
    rows = []
    for area in d["areas"]:
        for m in area["modules"]:
            rows.append(dict(area=area["name"], module=m["name"], status=m["status"],
                             suites=cover.get(m["name"], [])))
    return rows


def write_report(path, meta, results, skipped, rows):
    L = []
    w = L.append
    w("# ChemLarp validation report")
    w("")
    w("Written by `scripts/validate.py`. Every line below is the result of running a suite now, not a")
    w("recorded claim. Suites whose reference data or packages were missing are listed as skipped; they are")
    w("not evidence of anything.")
    w("")
    w("| | |")
    w("|---|---|")
    for k, v in meta.items():
        w("| %s | %s |" % (k, v))
    w("")
    passed = [r for r in results if r["status"] == "passed"]
    failed = [r for r in results if r["status"] == "failed"]
    w("**%d suites ran: %d passed, %d failed. %d skipped.**" % (len(results), len(passed), len(failed), len(skipped)))
    w("")
    for kind, head in (("independent", "Independent references"), ("internal", "In-repo expectations")):
        part = [r for r in results if r["kind"] == kind]
        if not part:
            continue
        w("## " + head)
        w("")
        if kind == "independent":
            w("Another implementation of the same model, published measurements, or a definition computed a")
            w("second way.")
        else:
            w("Values frozen in the repository. A wrong expectation would stay wrong, so these do not")
            w("establish correctness on their own.")
        w("")
        w("| Suite | Checks | Time | Result | Compared with |")
        w("|---|---|---|---|---|")
        for r in part:
            checks = "%d/%d" % (r["checks_ok"], r["checks_total"]) if r["checks_total"] else "—"
            mark = "passed" if r["status"] == "passed" else "**FAILED** (exit %d)" % r["exit_code"]
            w("| %s | %s | %.0f s | %s | %s |" % (r["title"], checks, r["seconds"], mark, r["reference"]))
        w("")
        for r in part:
            if r["headline"]:
                w("- %s: %s" % (r["id"], r["headline"]))
        w("")
    if failed:
        w("## Failures")
        w("")
        for r in failed:
            w("### " + r["title"])
            w("")
            w("`" + r["command"] + "` exited %d." % r["exit_code"])
            w("")
            w("```")
            tail = r["output"].strip().splitlines()[-40:]
            L.extend(tail)
            w("```")
            w("")
    if skipped:
        w("## Skipped")
        w("")
        w("| Suite | Missing | Compared with |")
        w("|---|---|---|")
        for s in skipped:
            w("| %s | %s | %s |" % (s["title"], "; ".join(s["missing"]), s["reference"]))
        w("")
        w("See the header of `scripts/validate.py` for where the reference data comes from.")
        w("")
    w("## Roadmap modules without an independent suite")
    w("")
    w("Modules marked done in `data/roadmap.json` that no suite here compares with anything outside the")
    w("repository while it runs. Several are checked against outside values that were brought into the")
    w("repository instead and frozen in the C++ test suite (JANAF tables for thermochemistry, crystal")
    w("structures for planar chirality, Merck's output for MMFF94); others were compared once by hand and")
    w("are listed in the last table. Neither is the same as a suite that recomputes the comparison.")
    w("")
    w("| Area | Module |")
    w("|---|---|")
    for r in rows:
        if r["status"] == "done" and not r["suites"]:
            w("| %s | %s |" % (r["area"], r["module"]))
    w("")
    w("## Described in the architecture notes but not scripted here")
    w("")
    w("| Module | Reported | Why it is not re-run |")
    w("|---|---|---|")
    for mod, claim, why in NOT_SCRIPTED:
        w("| %s | %s | %s |" % (mod, claim, why))
    w("")
    open(path, "w").write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Run every validation suite and write a report.")
    ap.add_argument("--list", action="store_true", help="show the suites and whether they can run")
    ap.add_argument("--only", default="", help="comma-separated suite ids")
    ap.add_argument("--skip", default="", help="comma-separated suite ids")
    ap.add_argument("--data", default=os.environ.get("CHEMLARP_VALIDATION_DATA", ""), help="reference-data directory")
    ap.add_argument("--config", default="", help="JSON file with explicit paths")
    ap.add_argument("--python", default="", help="interpreter that has the validation packages")
    ap.add_argument("--report", default="", help="write a Markdown report here")
    ap.add_argument("--json", default="", help="write the results as JSON here")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds per suite (default 3600)")
    a = ap.parse_args()

    cfg = {}
    config = a.config or (os.path.join(a.data, "validation.json") if a.data else "")
    if config and os.path.exists(config):
        cfg.update(json.load(open(config)))
    resolve_paths(cfg, a.data)
    cfg["python"] = a.python or cfg.get("python") or sys.executable

    only = {s for s in a.only.split(",") if s}
    skip = {s for s in a.skip.split(",") if s}
    chosen = [s for s in suites() if (not only or s["id"] in only) and s["id"] not in skip]
    modcache = {}
    plan = [(s, missing(s, cfg, modcache)) for s in chosen]

    if a.list:
        print("interpreter for the suites:", cfg["python"])
        print("reference data:", a.data or "(none given: --data or CHEMLARP_VALIDATION_DATA)")
        for s, miss in plan:
            print("  %-20s %-8s %s" % (s["id"], s["kind"], "ready" if not miss else "needs " + "; ".join(miss)))
        return 0

    results, skipped = [], []
    for s, miss in plan:
        if miss:
            print("skip %-20s needs %s" % (s["id"], "; ".join(miss)))
            skipped.append(dict(id=s["id"], title=s["title"], kind=s["kind"], reference=s["reference"], missing=miss))
            continue
        print("run  %-20s %s" % (s["id"], s["title"]), flush=True)
        r = run_suite(s, cfg, a.timeout)
        r.update(id=s["id"], title=s["title"], kind=s["kind"], reference=s["reference"], modules=s["modules"])
        print("     %s in %.0f s%s" % (r["status"], r["seconds"],
                                       "  " + r["headline"] if r["headline"] else ""), flush=True)
        results.append(r)

    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        except OSError:
            return "?"

    meta = {
        "run at": time.strftime("%Y-%m-%d %H:%M %Z"),
        "commit": git("rev-parse", "--short", "HEAD") + (" (modified working tree)" if git("status", "--porcelain") else ""),
        "platform": platform.platform(),
        "interpreter": cfg["python"],
        "reference data": a.data or "(paths given individually)",
    }
    rows = roadmap_coverage(results)
    if a.report:
        write_report(a.report, meta, results, skipped, rows)
        print("report:", a.report)
    if a.json:
        json.dump(dict(meta=meta, results=[{k: v for k, v in r.items() if k != "output"} for r in results],
                       skipped=skipped, coverage=rows), open(a.json, "w"), indent=1)
        print("json:  ", a.json)
    failed = [r for r in results if r["status"] == "failed"]
    print("\n%d ran, %d passed, %d failed, %d skipped" % (len(results), len(results) - len(failed), len(failed), len(skipped)))
    for r in failed:
        print("  FAILED", r["id"], "exit", r["exit_code"])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
