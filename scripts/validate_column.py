"""Distillation column validation against an independent solver (needs numpy, scipy and thermo; run after
building build/chemlarp: python3 scripts/validate_column.py).

The same MESH model (equilibrium stages, total condenser, kettle reboiler, reflux ratio and distillate flow
specified) is written out here as one system of equations and solved with scipy's hybrid Powell method,
with thermo's correlations instead of ChemLarp's data: vapour pressures (VaporPressure), ideal-gas heat
capacities (HeatCapacityGas), enthalpies of vaporization (EnthalpyVaporization; liquid enthalpy
H_V(T) - dHvap(T)) and, for the non-ideal case, thermo's original UNIFAC. Differences therefore combine
the property data and the solution; the stage-by-stage check in tests/test_column.cpp isolates the solver.
"""
import json, os, subprocess, sys
import numpy as np
from scipy.optimize import root
from thermo import VaporPressure, HeatCapacityGas, EnthalpyVaporization
from thermo.unifac import UNIFAC

CL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "chemlarp")

def chemlarp(text):
    p = subprocess.run([CL, "flowsheet", "-", "--json"], input=text, capture_output=True, text=True)
    r = json.loads(p.stdout)
    assert r["ok"], [d for d in r["diagnostics"] if d["severity"] == "error"]
    return r

def solve_mesh(cas, gammas, N, fstage, feed, Tfeed, R, D, P, start):
    C = len(cas)
    ps = [VaporPressure(CASRN=c) for c in cas]
    cg = [HeatCapacityGas(CASRN=c) for c in cas]
    hv = [EnthalpyVaporization(CASRN=c) for c in cas]
    HV = lambda i, T: cg[i].T_dependent_property_integral(298.15, T)
    HL = lambda i, T: HV(i, T) - hv[i].T_dependent_property(T)
    F = sum(feed)
    Hf = sum(feed[i] * HL(i, Tfeed) for i in range(C))
    B = F - D
    m = 2 * C + 1
    def K(T, x):
        g = gammas(T, x)
        return [g[i] * ps[i].T_dependent_property(T) / P for i in range(C)]
    def res(X):
        X = X.reshape(N, m)
        r = np.zeros((N, m))
        for j in range(N):
            T, l, v = X[j, 0], X[j, 1:1 + C], X[j, 1 + C:]
            L, V = l.sum(), v.sum()
            x = l / L
            for i in range(C):
                inn = (feed[i] if j == fstage - 1 else 0) + (X[j - 1, 1 + i] if j > 0 else 0) + (X[j + 1, 1 + C + i] if j < N - 1 else 0)
                r[j, 1 + i] = (l[i] + v[i] - inn) / F
            k = K(T, x)
            if j == 0:  # total condenser: v holds the distillate
                for i in range(C - 1):
                    r[j, 1 + C + i] = v[i] / V - x[i]
                r[j, 1 + C + C - 1] = sum(k[i] * x[i] for i in range(C)) - 1
                r[j, 0] = (L - R * V) / F
                continue
            for i in range(C):
                r[j, 1 + C + i] = v[i] / V - k[i] * x[i]
            if j == N - 1:
                r[j, 0] = (L - B) / F
            else:
                H = sum(l[i] * HL(i, T) + v[i] * HV(i, T) - X[j - 1, 1 + i] * HL(i, X[j - 1, 0]) - X[j + 1, 1 + C + i] * HV(i, X[j + 1, 0]) for i in range(C))
                r[j, 0] = (H - (Hf if j == fstage - 1 else 0)) / (F * 1e4)
        return r.ravel()
    sol = root(res, start.ravel(), method="hybr", options={"xtol": 1e-12})
    X = sol.x.reshape(N, m)
    assert np.abs(res(sol.x)).max() < 1e-8, (sol.message, np.abs(res(sol.x)).max())
    Qc = sum((X[0, 1 + i] + X[0, 1 + C + i]) * HL(i, X[0, 0]) - X[1, 1 + C + i] * HV(i, X[1, 0]) for i in range(C))
    Qr = sum(X[-1, 1 + i] * HL(i, X[-1, 0]) + X[-1, 1 + C + i] * HV(i, X[-1, 0]) - X[-2, 1 + i] * HL(i, X[-2, 0]) for i in range(C))
    return X, Qc, Qr

def start_from(r, C):
    col = next(u for u in r["units"] if u["type"] == "column")
    st = col["profile"]
    D = next(s for s in r["streams"] if s["from"] == col["name"])
    X = []
    for j, s in enumerate(st):
        l = [s["L"] * x for x in s["x"]]
        v = D["flows"] if j == 0 else [s["V"] * y for y in s["y"]]
        X.append([s["T_K"], *l, *v])
    return np.array(X), col

worst = {}
def check(case, label, got, ref, tol, absolute=False):
    e = abs(got - ref) if absolute else abs(got - ref) / abs(ref)
    worst[case] = max(worst.get(case, 0), e / tol)
    print(f"  {'ok ' if e <= tol else 'BAD'} {label:34s} ChemLarp {got:12.6g}  independent {ref:12.6g}  {'diff' if absolute else 'rel'} {e:.2g}")
    if e > tol:
        sys.exit(f"{case}: {label} outside tolerance {tol}")

ideal = lambda T, x: [1.0] * len(x)

print("A  benzene/toluene, 15 stages, feed on 8 (360 K liquid), R = 2, D = 50 mol/s, 1 atm, ideal solution")
txt = """components c1ccccc1 Cc1ccccc1
feed F T=360 P=1atm flows c1ccccc1=50 Cc1ccccc1=50
column T1 F -> D B stages=15 feed=8 R=2 D=50 P=1atm
"""
r = chemlarp(txt)
X0, col = start_from(r, 2)
X, Qc, Qr = solve_mesh(["71-43-2", "108-88-3"], ideal, 15, 8, [50, 50], 360.0, 2.0, 50.0, 101325.0, X0)
d = next(s for s in r["streams"] if s["name"] == "D")
check("A", "distillate benzene (mol/s)", d["flows"][0], X[0, 3], 0.01)
check("A", "bottoms benzene (mol/s)", 50 - d["flows"][0], 50 - X[0, 3], 0.05)
check("A", "top temperature (K)", X0[0, 0], X[0, 0], 0.3, True)
check("A", "bottom temperature (K)", X0[-1, 0], X[-1, 0], 0.5, True)
check("A", "feed-stage temperature (K)", X0[7, 0], X[7, 0], 0.5, True)
check("A", "condenser duty (kW)", col["values"]["condenser_duty_W"] / 1e3, Qc / 1e3, 0.02)
check("A", "reboiler duty (kW)", col["values"]["reboiler_duty_W"] / 1e3, Qr / 1e3, 0.02)

print("B  ethanol/water, 20 stages, feed on 12 (350 K), R = 3, D = 11 mol/s, 1 atm, UNIFAC")
txt = """components CCO O
model unifac
feed F T=350 P=1atm flows CCO=10 O=90
column T1 F -> D B stages=20 feed=12 R=3 D=11 P=1atm
"""
r = chemlarp(txt)
X0, col = start_from(r, 2)
groups = [{1: 1, 2: 1, 14: 1}, {16: 1}]   # CH3 CH2 OH; H2O (original UNIFAC subgroups)
gam = lambda T, x: UNIFAC.from_subgroups(T=T, xs=list(x), chemgroups=groups, version=0).gammas()
X, Qc, Qr = solve_mesh(["64-17-5", "7732-18-5"], gam, 20, 12, [10, 90], 350.0, 3.0, 11.0, 101325.0, X0)
d = next(s for s in r["streams"] if s["name"] == "D")
check("B", "distillate ethanol mole fraction", d["flows"][0] / 11, X[0, 3] / 11, 0.01, True)
check("B", "bottoms ethanol (mol/s)", 10 - d["flows"][0], 10 - X[0, 3], 0.15)
check("B", "top temperature (K)", X0[0, 0], X[0, 0], 0.3, True)
check("B", "bottom temperature (K)", X0[-1, 0], X[-1, 0], 0.5, True)
check("B", "condenser duty (kW)", col["values"]["condenser_duty_W"] / 1e3, Qc / 1e3, 0.03)
check("B", "reboiler duty (kW)", col["values"]["reboiler_duty_W"] / 1e3, Qr / 1e3, 0.03)
print("\nworst deviation as a fraction of its tolerance:", {k: round(v, 2) for k, v in worst.items()})
