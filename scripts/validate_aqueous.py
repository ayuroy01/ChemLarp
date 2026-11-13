"""Aqueous speciation validation against PHREEQC (needs phreeqpython, which bundles IPhreeqc and the same
phreeqc.dat; run after building build/chemlarp: python3 scripts/validate_aqueous.py).

PHREEQC is an independent implementation of the same ion-association model and database; differences come
from model details ChemLarp leaves out (PHREEQC lets the mass of water change, computes the Debye-Hueckel
A and B from its own water properties, and adds ionic-strength terms from -Vm parameters) and from the
solution itself. Measured data checked at the end: pKw (IUPAC), calcite-CO2 equilibrium (Plummer &
Busenberg 1982) and gypsum solubility (Marshall & Slusher 1966).
"""
import json, os, subprocess, sys
import phreeqpython

CL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "chemlarp")
pp = phreeqpython.PhreeqPython()


def ref(comp, t=25):
    """PHREEQC solution from component totals (mmol/kgw) with the pH from the charge balance, as ChemLarp's
    default. (add_solution_simple starts from water fixed at pH 7.000, which is not charge-balanced and
    leaves a 1e-9 eq/kg excess in every solution made from it.)"""
    d = {"temp": t, "units": "mmol/kgw", "pH": "7 charge"}
    d.update(comp)
    return pp.add_solution(d)


def chemlarp(text):
    p = subprocess.run([CL, "solution", "-", "--json"], input=text, capture_output=True, text=True)
    r = json.loads(p.stdout)
    assert r["ok"], r["diagnostics"]
    return r


def sp(r, name):
    for s in r["species"]:
        if s["name"] == name:
            return s["molality"]
    return 0.0


def si(r, phase):
    return next(s["si"] for s in r["saturation"] if s["phase"] == phase)


worst = {"pH": 0.0, "I": 0.0, "species": 0.0, "SI": 0.0}
failures = []


def cmp(case, label, got, ref, tol, kind, rel=False):
    e = abs(got - ref) / max(abs(ref), 1e-300) if rel else abs(got - ref)
    worst[kind] = max(worst[kind], e)
    ok = e <= tol
    if not ok:
        failures.append(f"{case}: {label}")
    print(f"  {'ok ' if ok else 'BAD'} {label:28s} ChemLarp {got:13.6g}  PHREEQC {ref:13.6g}  {'rel ' if rel else ''}{e:.2g}")


def compare(case, text, sol, species=(), phases=(), tol_ph=2e-3, tol_sp=5e-3, tol_si=3e-3):
    print(case)
    r = chemlarp(text)
    cmp(case, "pH", r["pH"], sol.pH, tol_ph, "pH")
    cmp(case, "ionic strength", r["ionic_strength"], sol.I, tol_sp, "I", rel=True)
    for s in species:
        cmp(case, s, sp(r, s), sol.species_molalities.get(s, 0.0), tol_sp, "species", rel=True)
    for p in phases:
        cmp(case, "SI " + p, si(r, p), sol.si(p), tol_si, "SI")
    return r


# Pure water at several temperatures.
for t in (0, 25, 50, 100):
    compare(f"pure water, {t} C", f"temperature {t}\n", ref({}, t))

compare("NaCl 100 mmol/kg", "add NaCl 100\n", ref({"Na": 100, "Cl": 100}), ["Na+", "Cl-", "NaCl"])
compare("CaCl2 1 + NaHCO3 2", "add CaCl2 1\nadd NaHCO3 2\n", ref({"Ca": 1, "Cl": 2, "Na": 2, "C(4)": 2}),
        ["HCO3-", "CO3-2", "CO2", "CaHCO3+", "CaCO3", "Ca+2"], ["Calcite", "Aragonite"])
compare("same at 60 C", "temperature 60\nadd CaCl2 1\nadd NaHCO3 2\n",
        ref({"Ca": 1, "Cl": 2, "Na": 2, "C(4)": 2}, 60), ["HCO3-", "CO3-2", "CO2", "CaCO3"], ["Calcite"])
compare("HCl 1 mmol/kg", "add HCl 1\n", ref({"Cl": 1}), ["Cl-"])
compare("NH4Cl 10", "add NH4Cl 10\n", ref({"N(-3)": 10, "Cl": 10}), ["NH4+", "NH3"])
compare("NH3 10 (component N(-3) with charge balance)", "N(-3) 10\n",
        ref({"N(-3)": 10}), ["NH4+", "NH3", "OH-"])
compare("phosphate buffer NaH2PO4 10 + Na2HPO4 10", "add NaH2PO4 10\nadd Na2HPO4 10\n",
        ref({"Na": 30, "P": 20}), ["H2PO4-", "HPO4-2", "NaHPO4-"])
compare("AlCl3 0.1 (hydrolysis)", "add AlCl3 0.1\n", ref({"Al": 0.1, "Cl": 0.3}),
        ["Al+3", "AlOH+2", "Al(OH)2+"], ["Gibbsite"])
compare("FeCl3 0.1", "add FeCl3 0.1\n", ref({"Fe(3)": 0.1, "Cl": 0.3}), ["Fe+3", "FeOH+2", "Fe(OH)2+"],
        ["Fe(OH)3(a)"])
compare("MgSO4 5 + K2SO4 2", "add MgSO4 5\nadd K2SO4 2\n", ref({"Mg": 5, "K": 4, "S(6)": 7}),
        ["MgSO4", "KSO4-", "SO4-2"])

# Equilibrium with phases.
for t in (10, 25, 50, 80):
    s = ref({}, t)
    s.equalize(["Calcite", "CO2(g)"], [0, -3.5])
    # PHREEQC raises the CO2(g) solubility constant by 0.002 (10 C) to 0.017 (80 C) log units with a correction
    # that is not in phreeqc.dat; ChemLarp uses the database's constant (Plummer & Busenberg 1982). Measured
    # Henry constants at 25 C span about 0.008 log units, so the tolerance here covers that difference.
    r = compare(f"calcite + CO2(g) 10^-3.5 atm, {t} C", f"temperature {t}\nequilibrium Calcite 0\nequilibrium CO2(g) -3.5\n",
                s, ["Ca+2", "HCO3-", "CO3-2", "CaHCO3+"], ["Calcite", "CO2(g)"], tol_ph=0.012, tol_sp=0.03)
    cmp("calcite", "dissolved Ca (mol/kg)", next(t2["molality"] for t2 in r["totals"] if t2["component"] == "Ca"),
        s.total_element("Ca", units="mol") / s.mass, 0.012, "species", rel=True)
for gas, aq, target in (("H2S(g)", "H2S", -3), ("NH3(g)", "NH3", -5), ("N2(g)", "N2", -1)):
    s = ref({})
    s.equalize([gas], [target])
    compare(f"{gas} at 10^{target} atm", f"equilibrium {gas} {target}\n", s, [aq], [gas])
for phase in ("Gypsum", "Fluorite", "Anhydrite"):
    s = ref({})
    s.saturate(phase, 0)
    sat = next(p for p in compare(f"{phase} in pure water", f"equilibrium {phase} 0\n", s, [], [phase])["phases"])
    want = s.total_element({"Gypsum": "Ca", "Fluorite": "Ca", "Anhydrite": "Ca", "Barite": "Ba"}[phase], units="mol") / s.mass
    cmp(phase, "dissolved (mol/kg)", sat["dissolved_mol"], want, 5e-3, "species", rel=True)

# Barite: PHREEQC's solubility constant differs from the analytic expression phreeqc.dat gives for it
# (-9.97 against -9.84 at 25 C, 0.02-0.37 log units over 0-100 C) for a reason not visible in the file;
# ChemLarp evaluates the expression as written, so its barite solubility is 16 % higher. Reported, not tested.
s = ref({})
s.saturate("Barite", 0)
r = chemlarp("equilibrium Barite 0\n")
print(f"Barite in pure water (known difference): ChemLarp {r['phases'][0]['dissolved_mol']:.4g} mol/kg, "
      f"PHREEQC {s.total_element('Ba', units='mol') / s.mass:.4g} mol/kg")

# Every phase's saturation index in a solution with traces of every component, at 25 and 60 C.
comps = {"Na": 10, "Cl": 10, "Ca": 0.5, "Mg": 0.3, "K": 0.2, "Fe(2)": 0.01, "Fe(3)": 0.001, "Mn(2)": 0.01, "Al": 0.001,
         "Ba": 0.01, "Sr": 0.01, "Si": 0.1, "C(4)": 1.0, "S(6)": 0.5, "N(5)": 0.1, "N(-3)": 0.05, "B": 0.01, "P": 0.01,
         "F": 0.05, "Li": 0.01, "Br": 0.01, "Zn": 0.001, "Cd": 0.001, "Pb": 0.001, "Cu(2)": 0.001, "S(-2)": 0.01}
for t in (25, 60):
    d = {"temp": t, "units": "mmol/kgw", "pH": 7.5}
    d.update(comps)
    s = pp.add_solution(d)
    r = chemlarp(f"temperature {t}\npH 7.5\n" + "".join(f"{k} {v}\n" for k, v in comps.items()))
    n = 0
    for x in r["saturation"]:
        if x["phase"] in s.phases and x["phase"] not in ("Barite", "CO2(g)", "CH4(g)"):
            n += 1
            e = abs(x["si"] - s.phases[x["phase"]])
            worst["SI"] = max(worst["SI"], e)
            if e > 5e-3:
                failures.append(f"SI {x['phase']} at {t} C")
                print(f"  BAD SI {x['phase']} at {t} C: {x['si']:.4f} vs {s.phases[x['phase']]:.4f}")
    print(f"all-phase check at {t} C: {n} saturation indices compared")

# Seawater (PHREEQC example 1 composition, ppm converted to mmol/kgw by PHREEQC's own SOLUTION).
sea = pp.add_solution({"temp": 25, "pH": "8.22 charge", "units": "mmol/kgw", "Ca": 10.6, "Mg": 54.7, "Na": 485.8,
                       "K": 10.6, "Cl": 565.8, "S(6)": 29.3, "C(4)": 2.4})
compare("seawater (I = 0.7)", "units mmol/kgw\npH 8.22 charge\nCa 10.6\nMg 54.7\nNa 485.8\nK 10.6\nCl 565.8\nS(6) 29.3\nC(4) 2.4\n",
        sea, ["Na+", "Mg+2", "SO4-2", "MgSO4", "NaSO4-", "HCO3-", "CO3-2", "MgCO3"], ["Calcite", "Gypsum"],
        tol_ph=0.05, tol_sp=0.05, tol_si=0.05)

print("\nMeasured data")
r = chemlarp("temperature 25\n")
cmp("pKw", "pKw at 25 C (IUPAC 13.995)", 2 * r["pH"], 13.995, 0.01, "pH")
r = chemlarp("equilibrium Calcite 0\nequilibrium CO2(g) -3.5\n")
cmp("calcite", "Ca, calcite + 10^-3.5 atm CO2", next(t["molality"] for t in r["totals"] if t["component"] == "Ca"),
    0.00049, 0.05, "species", rel=True)  # about 0.49 mmol/kg (Plummer & Busenberg 1982 constants)
r = chemlarp("equilibrium Gypsum 0\n")
cmp("gypsum", "gypsum solubility 25 C", r["phases"][0]["dissolved_mol"], 0.01545, 0.03, "species", rel=True)

print("\nlargest differences:", {k: f"{v:.2g}" for k, v in worst.items()})
if failures:
    sys.exit("outside tolerance: " + "; ".join(failures))
