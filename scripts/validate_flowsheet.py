"""Flowsheet solver validation against independent solutions (needs numpy, scipy and thermo; run after building
build/chemlarp: python3 scripts/validate_flowsheet.py).

A  Felder & Rousseau Example 4.7-2 (propane dehydrogenation, recycle): closed-form balance + published figures.
B  CO2 + 3 H2 -> CH3OH + H2O loop with condenser and purge (Felder & Rousseau Example 4.7-3 statement):
   the problem is solved by hand algebra for the reactor-feed basis, then posed to ChemLarp in sequential form.
C  benzene/toluene flash with liquid recycle: an equation-oriented Newton solve (scipy) with thermo's vapour
   pressures (independent correlations) and Raoult's law.
D  N2 gas loop with a heater and a recycle splitter: steady-state energy balance with thermo's ideal-gas Cp.
"""
import json, os, subprocess, sys
import numpy as np
from scipy.optimize import brentq, fsolve
from thermo import VaporPressure, HeatCapacityGas

CL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "chemlarp")

def run(text):
    p = subprocess.run([CL, "flowsheet", "-", "--json"], input=text, capture_output=True, text=True)
    r = json.loads(p.stdout)
    errs = [d for d in r["diagnostics"] if d["severity"] == "error"]
    assert r["ok"], errs
    return r

def stream(r, name):
    return next(s for s in r["streams"] if s["name"] == name)

def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)

worst = {}
def check(case, label, got, ref, tol):
    e = rel(got, ref)
    worst[case] = max(worst.get(case, 0), e)
    flag = "ok " if e <= tol else "BAD"
    print(f"  {flag} {label:38s} ChemLarp {got:14.8g}  reference {ref:14.8g}  rel {e:.1e}")
    if e > tol:
        sys.exit(f"{case}: {label} outside tolerance {tol}")

# ---------------------------------------------------------------- A
print("A  propane dehydrogenation (Felder & Rousseau Ex. 4.7-2)")
a, b = 1 - 0.00555, 0.05               # fraction of reactor-outlet propane / propylene recycled
X = 95 / (95 + 5 / 0.00555)            # single-pass conversion for 95 % overall conversion
F = 100.0
for method in ("wegstein", "direct"):
    txt = f"""components CCC C=CC [HH]
feed F T=300 P=1bar flows CCC={F}
mixer M F R -> Q1
reactor RX Q1 -> Q2 key=CCC X={X:.17g} T=800 nu CCC=-1 C=CC=1 [HH]=1
separator S Q2 -> P R T=300 split CCC=0.00555 C=CC=0.95 [HH]=1
method {method}
tolerance 1e-10
max-iterations 1000
"""
    r = run(txt)
    Pin = F / (1 - a * (1 - X))                    # closed form
    Yout = X * Pin / (1 - b)
    q1, rr, pp = stream(r, "Q1")["flows"], stream(r, "R")["flows"], stream(r, "P")["flows"]
    print(f" {method}: {r['iterations']} iterations")
    check("A", "reactor feed propane", q1[0], Pin, 1e-7)
    check("A", "reactor feed propylene", q1[1], b * Yout, 1e-7)
    check("A", "product propane", pp[0], (1 - a) * (1 - X) * Pin, 1e-6)
    check("A", "product propylene", pp[1], X * Pin, 1e-7)
    check("A", "recycle propane", rr[0], a * (1 - X) * Pin, 1e-7)
# Published figures (rounded): reactor feed 995.9 mol C3H8 and 5 mol C3H6 per 100 mol fresh feed,
# recycle 895.9 + 5 mol, single-pass conversion 9.54 %, recycle ratio 9.0.
check("A", "published reactor-feed propane 995.9", q1[0], 995.9, 1e-4)
check("A", "published recycle ratio 9.0", (rr[0] + rr[1]) / F, 9.0, 2e-3)

# ---------------------------------------------------------------- B
print("B  CO2 hydrogenation loop with condenser and purge (Felder & Rousseau Ex. 4.7-3 statement)")
# Reactor feed 28.0 % CO2, 70.0 % H2, 2.00 % inerts; 60 % single-pass H2 conversion; fresh feed 0.400 % inerts;
# the condenser removes all methanol and water. Basis 100 mol/s reactor feed, solved by hand algebra:
xr = np.array([28.0, 70.0, 2.0])                     # CO2, H2, I in the reactor feed
ext = 0.6 * 70.0 / 3
gas = xr - ext * np.array([1, 3, 0])                 # condenser gas: 14, 28, 2
# n0 + r = 100 and 0.004 n0 + r * gas_I/sum(gas) = 2
yI = gas[2] / gas.sum()
n0 = (2 - 100 * yI) / (0.004 - yI)
r_ = 100 - n0
fresh = xr - r_ * gas / gas.sum()
purge_frac = 1 - r_ / gas.sum()
txt = f"""components O=C=O [HH] [Ar] CO O
feed F T=300 P=50bar flows O=C=O={fresh[0]:.17g} [HH]={fresh[1]:.17g} [Ar]={fresh[2]:.17g}
mixer M F R -> Q1
reactor RX Q1 -> Q2 key=[HH] X=0.6 T=500 nu O=C=O=-1 [HH]=-3 CO=1 O=1
separator C Q2 -> G L T=300 split O=C=O=1 [HH]=1 [Ar]=1
splitter SP G -> R PURGE fractions={1 - purge_frac:.17g},{purge_frac:.17g}
tolerance 1e-10
max-iterations 500
"""
r = run(txt)
q1 = stream(r, "Q1")["flows"]
print(f"  {r['iterations']} iterations; fresh feed {fresh.sum():.4f} mol/s, x_I = {fresh[2]/fresh.sum():.5f}")
check("B", "reactor feed CO2", q1[0], 28.0, 1e-7)
check("B", "reactor feed H2", q1[1], 70.0, 1e-7)
check("B", "reactor feed inerts", q1[2], 2.0, 1e-7)
check("B", "methanol product", stream(r, "L")["flows"][3], ext, 1e-7)
check("B", "purge", sum(stream(r, "PURGE")["flows"]), gas.sum() - r_, 1e-6)

# ---------------------------------------------------------------- C
print("C  benzene/toluene flash with liquid recycle (equation-oriented Newton, thermo vapour pressures)")
T, P, frec = 368.0, 101325.0, 0.6
txt = f"""components c1ccccc1 Cc1ccccc1
feed F T=300 P=1atm flows c1ccccc1=50 Cc1ccccc1=50
mixer M F R -> S1
heater H S1 -> S2 T={T}
flash D S2 -> V L T={T} P={P}
splitter SP L -> R B fractions={frec},{1-frec}
"""
r = run(txt)
psat = [VaporPressure(CASRN="71-43-2").T_dependent_property(T), VaporPressure(CASRN="108-88-3").T_dependent_property(T)]
K = np.array(psat) / P
f = np.array([50.0, 50.0])
def eqs(u):
    l = u[:2]                                         # liquid from the flash
    z = f + frec * l                                  # flash feed
    beta = brentq(lambda b: np.sum(z / z.sum() * (K - 1) / (1 + b * (K - 1))), 0, 1)
    x = z / z.sum() / (1 + beta * (K - 1))
    return l - (1 - beta) * z.sum() * x
l = fsolve(eqs, [40.0, 60.0], xtol=1e-13)
v = f + frec * l - l
Vc = stream(r, "V")["flows"]
# The same equations with ChemLarp's own K-values isolate the recycle solution from the vapour-pressure data.
K_thermo = K
K = np.array(stream(r, "V")["y"]) / np.array(stream(r, "L")["x"])
l2 = fsolve(eqs, [40.0, 60.0], xtol=1e-13)
check("C", "vapour benzene (ChemLarp K)", Vc[0], (f + frec * l2 - l2)[0], 1e-7)
check("C", "vapour toluene (ChemLarp K)", Vc[1], (f + frec * l2 - l2)[1], 1e-7)
K = K_thermo
print(f"  {r['iterations']} iterations; K ChemLarp/thermo: "
      f"{stream(r,'V')['y'][0]/stream(r,'L')['x'][0]/K[0]:.4f}, {stream(r,'V')['y'][1]/stream(r,'L')['x'][1]/K[1]:.4f}")
check("C", "vapour benzene", Vc[0], v[0], 0.02)
check("C", "vapour toluene", Vc[1], v[1], 0.02)
check("C", "recycle total", sum(stream(r, "R")["flows"]), frec * l.sum(), 0.01)

# ---------------------------------------------------------------- D
print("D  N2 loop: heater and 50 % recycle; steady state F (h(T) - h(300 K)) = Q")
Q = 500e3
txt = f"""components N#N
feed F T=300 P=1bar flows N#N=100
mixer M F R -> S1
heater H S1 -> S2 Q={Q}
splitter SP S2 -> R OUT fractions=0.5,0.5
"""
r = run(txt)
cp = HeatCapacityGas(CASRN="7727-37-9")
Tref = brentq(lambda t: 100 * cp.T_dependent_property_integral(300, t) - Q, 300, 2000)
print(f"  {r['iterations']} iterations")
check("D", "product temperature (K)", stream(r, "OUT")["T_K"], Tref, 1e-3)
check("D", "energy balance: H_out - H_in = Q", stream(r, "OUT")["H_W"] - stream(r, "F")["H_W"], Q, 1e-7)

print("\nworst relative deviations:", {k: f"{v:.1e}" for k, v in worst.items()})
