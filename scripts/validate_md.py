"""Molecular dynamics validation (needs numpy; run after building build/chemlarp: python3 scripts/validate_md.py).

A  Lennard-Jones fluid: ChemLarp's MD (Bussi and Langevin thermostats) against an independent Metropolis Monte
   Carlo simulation written here (same truncated-and-shifted potential, same N, NVT). Different algorithms
   sampling the same canonical ensemble must give the same potential energy and pressure within statistics.
B  Harmonic limit: at 1 K the dynamics of a molecule is harmonic, so the peaks of the velocity spectrum
   must sit at the normal-mode frequencies from the finite-difference Hessian.
C  Equipartition: for a classical harmonic system the mean potential energy above the minimum is
   (f/2) kB T, with f the number of vibrational degrees of freedom (3N - 6).
D  Integrator: the energy fluctuation of velocity Verlet scales as dt^2.
"""
import json, os, subprocess, sys
import numpy as np

CL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "chemlarp")
failures = []


def run(args):
    p = subprocess.run([CL, *args], capture_output=True, text=True)
    return json.loads(p.stdout)


def check(label, ok, text):
    print(f"  {'ok ' if ok else 'BAD'} {label:44s} {text}")
    if not ok:
        failures.append(label)


# ---------------------------------------------------------------- A
def lj_mc(cells, rho, T, rc, sweeps, equil, seed=7):
    rng = np.random.default_rng(seed)
    N = 4 * cells ** 3
    L = (N / rho) ** (1 / 3)
    a = L / cells
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    x = np.array([(np.array([i, j, k]) + b + 0.25) * a for i in range(cells) for j in range(cells) for k in range(cells) for b in basis])
    rc2 = rc * rc
    ushift = 4 * (rc ** -12 - rc ** -6)

    def pair_u(d2):
        m = d2 < rc2
        ir6 = np.where(m, 1.0 / np.where(m, d2, 1.0) ** 3, 0.0)
        return np.where(m, 4 * ir6 * (ir6 - 1) - ushift, 0.0)

    def particle_energy(i, pos):
        d = x - pos
        d -= L * np.round(d / L)
        d2 = (d * d).sum(1)
        d2[i] = np.inf
        return pair_u(d2).sum()

    def totals():
        d = x[:, None, :] - x[None, :, :]
        d -= L * np.round(d / L)
        d2 = (d * d).sum(2)
        iu = np.triu_indices(N, 1)
        r2 = d2[iu]
        m = r2 < rc2
        ir6 = 1.0 / r2[m] ** 3
        U = (4 * ir6 * (ir6 - 1) - ushift).sum()
        W = (48 * ir6 * (ir6 - 0.5)).sum()  # sum of r * F(r)
        return U, W

    step = 0.12
    us, ps = [], []
    acc = 0
    for sweep in range(equil + sweeps):
        for i in rng.permutation(N):
            old = x[i].copy()
            e0 = particle_energy(i, old)
            new = old + step * (rng.random(3) - 0.5) * 2
            e1 = particle_energy(i, new)
            if e1 <= e0 or rng.random() < np.exp(-(e1 - e0) / T):
                x[i] = new % L
                acc += 1
        if sweep < equil and sweep % 50 == 49:  # tune the step for about 40 % acceptance
            rate = acc / (50 * N)
            step *= 1.1 if rate > 0.4 else 0.9
            acc = 0
        if sweep >= equil:
            U, W = totals()
            us.append(U / N)
            ps.append((N * T + W / 3) / L ** 3)

    def block(y):
        y = np.array(y)
        b = y[: len(y) // 20 * 20].reshape(20, -1).mean(1)
        return y.mean(), b.std(ddof=1) / np.sqrt(20)

    return block(us), block(ps)


print("A  Lennard-Jones fluid (N = 108, rc = 2.5, truncated and shifted): MD against Monte Carlo")
for rho, T in ((0.8, 1.0), (0.4, 1.5)):
    (umc, uerr), (pmc, perr) = lj_mc(3, rho, T, 2.5, sweeps=6000, equil=1000)
    for th in ("bussi", "langevin"):
        r = run(["lj", "--cells", "3", "--density", str(rho), "-T", str(T), "--equil", "4000", "--steps", "60000",
                 "--thermostat", th, "--seed", "3"])
        su = np.hypot(uerr, r["energy_error"])
        sp = np.hypot(perr, r["pressure_error"])
        check(f"rho {rho} T {T} {th}: U/N", abs(r["energy_per_particle"] - umc) < 4 * su,
              f"MD {r['energy_per_particle']:.4f} ± {r['energy_error']:.4f}   MC {umc:.4f} ± {uerr:.4f}")
        check(f"rho {rho} T {T} {th}: pressure", abs(r["pressure"] - pmc) < 4 * sp,
              f"MD {r['pressure']:.4f} ± {r['pressure_error']:.4f}   MC {pmc:.4f} ± {perr:.4f}")
        check(f"rho {rho} T {T} {th}: kinetic temperature", abs(r["temperature"] - T) < 0.02 * T, f"{r['temperature']:.4f}")
r = run(["lj", "--cells", "3", "--density", "0.8", "-T", "1.0", "--equil", "0", "--steps", "20000", "--thermostat", "none", "--dt", "0.002"])
check("NVE drift per particle per unit time", abs(r["drift"]) < 1e-4, f"{r['drift']:.2e}")

# ---------------------------------------------------------------- B
print("B  harmonic limit: velocity-spectrum peaks against normal modes (1 K, NVE)")
for smiles in ("O", "C=O", "C#N"):
    r = run(["md", smiles, "-T", "1", "--thermostat", "none", "--steps", "80000", "--dt", "0.25", "--every", "4", "--modes", "--json"])
    spec = np.array(r["spectrum"]["intensity"])
    grid = np.linspace(0, r["spectrum"]["max_cm1"], len(spec))
    peaks = [grid[i] for i in range(2, len(spec) - 2) if spec[i] > 0.01 and spec[i] == spec[i - 2:i + 3].max()]
    worst = 0.0
    for f in r["modes"]:
        near = min(peaks, key=lambda p: abs(p - f))
        worst = max(worst, abs(near - f))
    check(f"{smiles}: {len(r['modes'])} modes found in the spectrum", worst <= 6.0,
          f"modes {[round(f, 1) for f in r['modes']]}; largest offset {worst:.1f} cm-1 (grid 5 cm-1)")

# ---------------------------------------------------------------- C
print("C  equipartition at 10 K: <U - Umin> = (3N - 6) kB T / 2")
kB = 0.0019872043
for smiles, th in (("CCO", "bussi"), ("CCO", "langevin"), ("c1ccccc1", "bussi")):
    r = run(["md", smiles, "-T", "10", "--thermostat", th, "--tau", "50", "--steps", "200000", "--dt", "0.5", "--every", "20", "--json"])
    n = len(r["elements"])
    f = 3 * n - 6
    excess = r["mean_potential"] - r["minimized_energy"]
    want = 0.5 * f * kB * 10
    check(f"{smiles} {th}: mean potential above minimum", abs(excess - want) < 0.05 * want,
          f"{excess:.5f} kcal/mol, expected {want:.5f} (mean T {r['mean_temperature_K']:.2f} K)")

# ---------------------------------------------------------------- D
print("D  velocity Verlet: energy fluctuation scales as dt^2 (ethanol, 300 K, NVE)")
sd = []
for dt in (0.4, 0.2, 0.1):
    r = run(["md", "CCO", "-T", "300", "--thermostat", "none", "--steps", str(int(2000 / dt)), "--dt", str(dt), "--every", "5", "--json"])
    sd.append(r["fluctuation"])
check("fluctuation ratio dt 0.4 / 0.2", 3.5 < sd[0] / sd[1] < 4.5, f"{sd[0] / sd[1]:.2f}")
check("fluctuation ratio dt 0.2 / 0.1", 3.5 < sd[1] / sd[2] < 4.5, f"{sd[1] / sd[2]:.2f}")

if failures:
    sys.exit("failed: " + "; ".join(failures))
print("\nall checks passed")
