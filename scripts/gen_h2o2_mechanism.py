#!/usr/bin/env python3
"""Generate the example H2/O2 mechanism in every supported format from one source.

    python3 scripts/gen_h2o2_mechanism.py

Writes:
    data/mechanisms/h2o2.ocm               OmniChem native (SMILES-bound species)
    data/mechanisms/h2o2/chem.inp          CHEMKIN-II reactions
    data/mechanisms/h2o2/therm.dat         CHEMKIN-II NASA-7 thermo cards (fixed columns)
    data/mechanisms/h2o2/species.smi       NAME -> SMILES map for the CHEMKIN / YAML imports
    data/mechanisms/h2o2.yaml              Cantera YAML

All four describe the same mechanism; the test suite checks that the three
importers produce identical rates.

ILLUSTRATIVE PARAMETERS: rate constants follow the structure and approximate
values of the Li, Zhao, Kazakov & Dryer (2004) H2/O2 mechanism; thermo is
GRI-Mech 3.0. Replace with a validated mechanism for quantitative work.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# name, SMILES, composition, (Tlow, Tmid, Thigh), low-T coeffs, high-T coeffs   (GRI-Mech 3.0)
SPECIES = [
    ("H2", "[H][H]", {"H": 2}, (200, 1000, 3500),
     [2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01],
     [3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00]),
    ("O2", "O=O", {"O": 2}, (200, 1000, 3500),
     [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00],
     [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00]),
    ("H", "[H]", {"H": 1}, (200, 1000, 3500),
     [2.50000000E+00, 7.05332819E-13, -1.99591964E-15, 2.30081632E-18, -9.27732332E-22, 2.54736599E+04, -4.46682853E-01],
     [2.50000001E+00, -2.30842973E-11, 1.61561948E-14, -4.73515235E-18, 4.98197357E-22, 2.54736599E+04, -4.46682914E-01]),
    ("O", "[O]", {"O": 1}, (200, 1000, 3500),
     [3.16826710E+00, -3.27931884E-03, 6.64306396E-06, -6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+00],
     [2.56942078E+00, -8.59741137E-05, 4.19484589E-08, -1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00]),
    ("OH", "[OH]", {"O": 1, "H": 1}, (200, 1000, 3500),
     [3.99201543E+00, -2.40131752E-03, 4.61793841E-06, -3.88113333E-09, 1.36411470E-12, 3.61508056E+03, -1.03925458E-01],
     [3.09288767E+00, 5.48429716E-04, 1.26505228E-07, -8.79461556E-11, 1.17412376E-14, 3.85865700E+03, 4.47669610E+00]),
    ("H2O", "O", {"H": 2, "O": 1}, (200, 1000, 3500),
     [4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01],
     [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]),
    ("HO2", "[O]O", {"H": 1, "O": 2}, (200, 1000, 3500),
     [4.30179801E+00, -4.74912051E-03, 2.11582891E-05, -2.42763894E-08, 9.29225124E-12, 2.94808040E+02, 3.71666245E+00],
     [4.01721090E+00, 2.23982013E-03, -6.33658150E-07, 1.14246370E-10, -1.07908535E-14, 1.11856713E+02, 3.78510215E+00]),
    ("H2O2", "OO", {"H": 2, "O": 2}, (200, 1000, 3500),
     [4.27611269E+00, -5.42822417E-04, 1.67335701E-05, -2.15770813E-08, 8.62454363E-12, -1.77025821E+04, 3.43505074E+00],
     [4.16500285E+00, 4.90831694E-03, -1.90139225E-06, 3.71185986E-10, -2.87908305E-14, -1.78617877E+04, 2.91615662E+00]),
    ("AR", "[Ar]", {"Ar": 1}, (300, 1000, 5000),
     [2.5, 0.0, 0.0, 0.0, 0.0, -7.45375E+02, 4.366],
     [2.5, 0.0, 0.0, 0.0, 0.0, -7.45375E+02, 4.366]),
]

# Reactions: (equation, A [cm,mol,s], b, Ea [cal/mol], extras)
# extras: eff={...}, low=(A,b,Ea), troe=(a,T3,T1[,T2]), dup=True
M_EFF = {"H2": 2.5, "H2O": 12.0}
REACTIONS = [
    ("H + O2 <=> O + OH", 3.547e15, -0.406, 16599, {}),
    ("O + H2 <=> H + OH", 5.080e4, 2.67, 6290, {}),
    ("H2 + OH <=> H2O + H", 2.160e8, 1.51, 3430, {}),
    ("O + H2O <=> OH + OH", 2.970e6, 2.02, 13400, {}),
    ("H2 + M <=> H + H + M", 4.577e19, -1.40, 104380, {"eff": {**M_EFF, "AR": 0.0}}),
    ("H2 + AR <=> H + H + AR", 5.840e18, -1.10, 104380, {}),
    ("O + O + M <=> O2 + M", 6.165e15, -0.50, 0, {"eff": {**M_EFF, "AR": 0.0}}),
    ("O + O + AR <=> O2 + AR", 1.886e13, 0.0, -1788, {}),
    ("O + H + M <=> OH + M", 4.714e18, -1.00, 0, {"eff": {**M_EFF, "AR": 0.75}}),
    ("H + OH + M <=> H2O + M", 3.800e22, -2.00, 0, {"eff": {**M_EFF, "AR": 0.38}}),
    ("H + O2 (+M) <=> HO2 (+M)", 1.475e12, 0.60, 0,
     {"low": (6.366e20, -1.72, 524.8), "troe": (0.8, 1e-30, 1e30), "eff": {"H2": 2.0, "H2O": 11.0, "O2": 0.78, "AR": 0.67}}),
    ("HO2 + H <=> H2 + O2", 1.660e13, 0.0, 823, {}),
    ("HO2 + H <=> OH + OH", 7.079e13, 0.0, 295, {}),
    ("HO2 + O <=> O2 + OH", 3.250e13, 0.0, 0, {}),
    ("HO2 + OH <=> H2O + O2", 2.890e13, 0.0, -497, {}),
    ("HO2 + HO2 <=> H2O2 + O2", 4.200e14, 0.0, 11982, {"dup": True}),
    ("HO2 + HO2 <=> H2O2 + O2", 1.300e11, 0.0, -1629.3, {"dup": True}),
    ("H2O2 (+M) <=> OH + OH (+M)", 2.951e14, 0.0, 48430,
     {"low": (1.202e17, 0.0, 45500), "troe": (0.5, 1e-30, 1e30), "eff": {**M_EFF, "AR": 0.64}}),
    ("H2O2 + H <=> H2O + OH", 2.410e13, 0.0, 3970, {}),
    ("H2O2 + H <=> HO2 + H2", 4.820e13, 0.0, 7950, {}),
    ("H2O2 + O <=> OH + HO2", 9.550e6, 2.00, 3970, {}),
    ("H2O2 + OH <=> H2O + HO2", 1.000e12, 0.0, 0, {"dup": True}),
    ("H2O2 + OH <=> H2O + HO2", 5.800e14, 0.0, 9557, {"dup": True}),
]

HEADER = """ILLUSTRATIVE H2/O2 mechanism (9 species, 23 reactions) generated by
scripts/gen_h2o2_mechanism.py. Rate constants follow the structure and
approximate values of Li, Zhao, Kazakov & Dryer (2004); thermo is GRI-Mech 3.0.
Not validated: replace with a published mechanism for quantitative work."""


def aux_lines(extra, indent):
    out = []
    if "low" in extra:
        out.append(f"{indent}LOW / {extra['low'][0]:.4e} {extra['low'][1]} {extra['low'][2]} /")
    if "troe" in extra:
        out.append(f"{indent}TROE / " + " ".join(f"{v:g}" for v in extra["troe"]) + " /")
    if "eff" in extra:
        out.append(indent + " ".join(f"{k}/{v}/" for k, v in extra["eff"].items()))
    if extra.get("dup"):
        out.append(f"{indent}DUPLICATE")
    return out


def reaction_block(indent):
    lines = []
    for eq, A, b, Ea, extra in REACTIONS:
        lines.append(f"{indent}{eq:<30} {A:.4e} {b:>7} {Ea:>9}")
        lines += aux_lines(extra, indent + "  ")
    return lines


def write_ocm():
    lines = ["# " + l for l in HEADER.splitlines()]
    lines += ["#", "# Reverse rates come from detailed balance with the species thermo.", "", "species"]
    lines += [f"  {n:<6} {smi}" for n, smi, *_ in SPECIES]
    lines += ["end", "", "reactions  cal/mol  cm3"] + reaction_block("  ") + ["end", ""]
    (ROOT / "data/mechanisms/h2o2.ocm").write_text("\n".join(lines))


def card(name, comp, T, low, high):
    elems = "".join(f"{sym.upper():<2}{n:>3d}" for sym, n in comp.items())
    elems = elems.ljust(20)
    l1 = f"{name:<18}{'OMNI':<6}{elems}G{T[0]:>10.3f}{T[2]:>10.3f}{T[1]:>8.2f}".ljust(79) + "1"
    c = high[:5], high[5:] + low[:3], low[3:]
    l2 = "".join(f"{v:15.8E}" for v in c[0]).ljust(79) + "2"
    l3 = "".join(f"{v:15.8E}" for v in c[1]).ljust(79) + "3"
    l4 = "".join(f"{v:15.8E}" for v in c[2]).ljust(79) + "4"
    return [l1, l2, l3, l4]


def write_chemkin():
    d = ROOT / "data/mechanisms/h2o2"
    d.mkdir(parents=True, exist_ok=True)
    chem = ["! " + l for l in HEADER.splitlines()]
    chem += ["ELEMENTS", "H O AR", "END", "SPECIES", " ".join(n for n, *_ in SPECIES), "END",
             "REACTIONS  CAL/MOLE  MOLES"]
    chem += reaction_block("") + ["END", ""]
    (d / "chem.inp").write_text("\n".join(chem))
    therm = ["THERMO", "   300.000  1000.000  5000.000"]
    for n, _, comp, T, low, high in SPECIES:
        therm += card(n, comp, T, low, high)
    therm += ["END", ""]
    (d / "therm.dat").write_text("\n".join(therm))
    (d / "species.smi").write_text("# NAME SMILES map for the CHEMKIN and YAML versions\n" +
                                   "\n".join(f"{n} {smi}" for n, smi, *_ in SPECIES) + "\n")


def yrate(A, b, Ea):
    return f"{{A: {A:.6g}, b: {b:g}, Ea: {Ea:g}}}"


def write_yaml():
    y = ["description: |-"] + ["  " + l for l in HEADER.splitlines()]
    y += ["", "units: {length: cm, time: s, quantity: mol, activation-energy: cal/mol}", "",
          "phases:", "- name: gas", "  thermo: ideal-gas", "  elements: [H, O, Ar]",
          "  species: [" + ", ".join(n for n, *_ in SPECIES) + "]", "  kinetics: gas",
          "  state: {T: 300.0, P: 1 atm}", "", "species:"]
    for n, _, comp, T, low, high in SPECIES:
        y += [f"- name: {n}", "  composition: {" + ", ".join(f"{k}: {v}" for k, v in comp.items()) + "}",
              "  thermo:", "    model: NASA7", f"    temperature-ranges: [{T[0]:.1f}, {T[1]:.1f}, {T[2]:.1f}]",
              "    data:"]
        for arr in (low, high):
            vals = [f"{v:.9e}" for v in arr]
            y += ["    - [" + ", ".join(vals[:4]) + ",", "      " + ", ".join(vals[4:]) + "]"]
    y += ["", "reactions:"]
    for eq, A, b, Ea, extra in REACTIONS:
        y.append(f"- equation: {eq}")
        if "low" in extra:
            y += ["  type: falloff", f"  low-P-rate-constant: {yrate(*extra['low'])}",
                  f"  high-P-rate-constant: {yrate(A, b, Ea)}"]
            t = extra["troe"]
            y.append("  Troe: {" + ", ".join(f"{k}: {v:g}" for k, v in zip(["A", "T3", "T1", "T2"], t)) + "}")
        else:
            if " M " in f" {eq} ":
                y.append("  type: three-body")
            y.append(f"  rate-constant: {yrate(A, b, Ea)}")
        if "eff" in extra:
            y.append("  efficiencies: {" + ", ".join(f"{k}: {v}" for k, v in extra["eff"].items()) + "}")
        if extra.get("dup"):
            y.append("  duplicate: true")
    (ROOT / "data/mechanisms/h2o2.yaml").write_text("\n".join(y) + "\n")


if __name__ == "__main__":
    write_ocm()
    write_chemkin()
    write_yaml()
    print("wrote h2o2.ocm, h2o2/chem.inp, h2o2/therm.dat, h2o2/species.smi, h2o2.yaml")
