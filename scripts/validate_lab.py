#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Virtual laboratory against PHREEQC and reference data (needs phreeqpython; build first).

    python3 scripts/validate_lab.py

A  Mixtures: each experiment is posed to PHREEQC (IPhreeqc through phreeqpython, the same phreeqc.dat) as two
   SOLUTIONs at the same molalities ChemLarp uses, a MIX of them, and EQUILIBRIUM_PHASES for the precipitate
   or gas at saturation index 0 with nothing present initially (so it can only form). Compared: pH, moles
   precipitated or given off. PHREEQC lets the mass of water change with the reactions; ChemLarp holds it.
B  Titration curves at several volumes, the same way.
C  Heat of neutralisation (H+ + OH- -> H2O) per mole, at 0.001 to 0.1 mol/L: against the CODATA value
   (-55.815 kJ/mol from the standard enthalpies of formation of H2O(l) and OH-(aq), Cox, Wagman & Medvedev,
   CODATA Key Values for Thermodynamics, 1989) and against what phreeqc.dat's log K(T) for water implies.
D  The classic calorimetry experiment: 50 mL of 1 mol/L HCl with 50 mL of 1 mol/L NaOH warms by about 6.7 K
   (0.05 mol x 55.8 kJ/mol over 100 g x 4.184 J/(g K)).
"""
import json
import math
import os
import subprocess
import sys

import phreeqpython

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CL = os.path.join(ROOT, "build", "chemlarp")
RHO = 0.99705  # kg water per litre, as in src/lab.cpp
pp = phreeqpython.PhreeqPython()
failures = []


def check(label, ok, text):
    print("  %s %-46s %s" % ("ok " if ok else "BAD", label, text))
    if not ok:
        failures.append(label)


def chemlarp(kind, text):
    p = subprocess.run([CL, kind, "-"], input=text, capture_output=True, text=True)
    r = json.loads(p.stdout)
    assert r["ok"], r["diagnostics"]
    return r


# PHREEQC SOLUTION lines for catalogue reagents (element totals per mol of reagent; pH from the charge balance).
ELEMENTS = {
    "hcl": {"Cl": 1}, "naoh": {"Na": 1}, "h3po4": {"P": 1}, "na2co3": {"Na": 2, "C(4)": 1}, "nahco3": {"Na": 1, "C(4)": 1},
    "bacl2": {"Ba": 1, "Cl": 2}, "na2so4": {"Na": 2, "S(6)": 1}, "na2s": {"Na": 2, "S(-2)": 1}, "cacl2": {"Ca": 1, "Cl": 2},
    "nh3": {"N(-3)": 1}, "nh4cl": {"N(-3)": 1, "Cl": 1}, "pbno3": {"Pb": 1, "N(5)": 2},
}


def phreeqc_mix(portions, phases, T=25.0):
    """portions: [(volume_mL, {reagent: mol/L})]; returns pH and the moles formed of each phase at T (°C)."""
    lines = []
    for i, (v, sol) in enumerate(portions, 1):
        kgw = v / 1000.0 * RHO
        tot = {}
        for reagent, c in sol.items():
            for el, k in ELEMENTS[reagent].items():
                tot[el] = tot.get(el, 0.0) + k * c * (v / 1000.0) / kgw
        lines += ["SOLUTION %d" % i, " units mol/kgw", " temp %.10g" % T, " -water %.10g" % kgw, " pH 7 charge"]
        lines += [" %s %.12g" % (el, m) for el, m in tot.items()]
    lines += ["MIX 1"] + [" %d 1" % i for i in range(1, len(portions) + 1)]
    if phases:
        lines += ["EQUILIBRIUM_PHASES 1"] + [" %s 0 0" % p for p in phases]
    lines += ["SELECTED_OUTPUT", " -reset false", " -pH true"]
    if phases:
        lines += [" -equilibrium_phases " + " ".join(phases)]
    lines += ["END"]
    pp.ip.run_string("\n".join(lines))
    out = pp.ip.get_selected_output_array()
    head, row = out[0], out[-1]
    res = {"pH": row[head.index("pH")]}
    for p in phases:
        res[p] = row[head.index("d_" + p)] if ("d_" + p) in head else 0.0
    return res


def mix_text(portions):
    t = ""
    for i, (v, sol) in enumerate(portions):
        t += "portion p%d %g\n" % (i, v)
        t += "".join("%s %g\n" % (r, c) for r, c in sol.items())
    return t


def main():
    print("A  mixtures against PHREEQC (MIX + EQUILIBRIUM_PHASES)")
    cases = [
        ("HCl + NaOH, 1 mol/L, 50 + 50 mL", [(50, {"hcl": 1}), (50, {"naoh": 1})], []),
        ("BaCl2 + Na2SO4, 0.1 mol/L: barite", [(50, {"bacl2": 0.1}), (50, {"na2so4": 0.1})], ["Barite"]),
        ("NaHCO3 + HCl, 1 mol/L: CO2 bubbles", [(50, {"nahco3": 1}), (50, {"hcl": 1})], ["CO2(g)"]),
        ("Na2S 0.5 + HCl 1 mol/L: H2S bubbles", [(50, {"na2s": 0.5}), (50, {"hcl": 1})], ["H2S(g)"]),
        ("CaCl2 + Na2CO3, 0.05 mol/L: calcite", [(40, {"cacl2": 0.05}), (60, {"na2co3": 0.05})], ["Calcite"]),
    ]
    for label, portions, phases in cases:
        r = chemlarp("lab", mix_text(portions))
        # ChemLarp reports the mixture at the temperature its reaction heat brings it to; PHREEQC mixes the same
        # solutions at that temperature (so the energy balance is ChemLarp's, the chemistry compared).
        ref = phreeqc_mix(portions, phases, r["T_final_C"])
        check(label + ": pH at %.1f °C" % r["T_final_C"], abs(r["pH"] - ref["pH"]) < 0.02,
              "ChemLarp %.4f  PHREEQC %.4f" % (r["pH"], ref["pH"]))
        for p in phases:
            ours = sum(x["mol"] for x in r["precipitates"] + r["gases"] if x["phase"] == p)
            theirs = abs(ref[p])
            rel = abs(ours - theirs) / max(theirs, 1e-12)
            check(label + ": " + p, rel < 0.02, "ChemLarp %.6g mol  PHREEQC %.6g mol  rel %.1e" % (ours, theirs, rel))

    print("\nB  titration curves against PHREEQC")
    for label, analyte, titrant, vols in [
        ("25 mL 0.1 M H3PO4 with 0.1 M NaOH", (25, {"h3po4": 0.1}), ("naoh", 0.1), [0, 5, 12.5, 24, 26, 37.5, 49, 51, 60]),
        ("25 mL 0.1 M Na2CO3 with 0.1 M HCl", (25, {"na2co3": 0.1}), ("hcl", 0.1), [0, 10, 24, 26, 40, 49, 51]),
    ]:
        text = "analyte %g\n" % analyte[0] + "".join("%s %g\n" % kv for kv in analyte[1].items())
        text += "titrant %s %g\nmax %g\npoints %d\n" % (titrant[0], titrant[1], vols[-1], int(vols[-1] * 4))
        r = chemlarp("titrate", text)
        curve = {round(v, 6): ph for v, ph in r["curve"]}
        worst = 0.0
        for v in vols:
            ours = curve[round(v, 6)]
            # Titrations take no heat into account: PHREEQC at 25 °C throughout.
            portions = [analyte] + ([(v, {titrant[0]: titrant[1]})] if v > 0 else [])
            ref = phreeqc_mix(portions, ["CO2(g)"] if "na2co3" in analyte[1] else [])
            worst = max(worst, abs(ours - ref["pH"]))
        check(label, worst < 0.03, "largest pH difference over %d points %.4f" % (len(vols), worst))
        print("      equivalence points:", ", ".join("%.2f mL (pH %.2f)" % (e["volume_mL"], e["pH"]) for e in r["equivalence"]))

    print("\nC  heat of neutralisation per mole")
    a = [293.29227, 0.1360833, -10576.913, -123.73158, 0.0, -6.996455e-05]  # phreeqc.dat, H2O = OH- + H+
    lk = lambda T: a[0] + a[1] * T + a[2] / T + a[3] * math.log10(T) + a[4] / T ** 2 + a[5] * T * T
    T, h, R = 298.15, 0.01, 8.314462618e-3
    db = -R * T * T * math.log(10) * (lk(T + h) - lk(T - h)) / (2 * h)
    for c in (0.1, 0.01, 0.001):
        r = chemlarp("lab", mix_text([(100, {"hcl": c}), (100, {"naoh": c})]))
        dh = -r["heat_released_kJ"] / (0.1 * c)
        check("%g mol/L: against phreeqc.dat's log K(T)" % c, abs(dh - db) < 0.05, "%.3f kJ/mol  database %.3f" % (dh, db))
        check("%g mol/L: against CODATA -55.815" % c, abs(dh + 55.815) < 0.02 * 55.815,
              "%.3f kJ/mol  (%.1f %% from CODATA; the database's own value is %.3f)" % (dh, 100 * abs(dh + 55.815) / 55.815, db))

    print("\nD  calorimetry: 50 mL 1 M HCl + 50 mL 1 M NaOH")
    r = chemlarp("lab", mix_text([(50, {"hcl": 1}), (50, {"naoh": 1})]))
    dT = r["T_final_C"] - r["T_initial_C"]
    expect = 0.05 * 55.815e3 / (100 * 4.184)
    check("temperature rise", abs(dT - expect) < 0.15, "%.2f K  (from CODATA: %.2f K)" % (dT, expect))

    if failures:
        sys.exit("outside tolerance: " + "; ".join(failures))
    print("\nall checks passed")


if __name__ == "__main__":
    main()
