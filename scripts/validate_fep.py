#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""Free-energy estimators and the alchemical Lennard-Jones particle (needs numpy and pymbar; build first).

    python3 scripts/validate_fep.py [--seeds N]

A  Exact answers: five harmonic oscillators u_k(x) = K_k (x - x0_k)^2 / 2 (kT units) sampled exactly. Their
   reduced free energies are f_k = ln(K_k / 2 pi) / 2, so f_k - f_0 is known in closed form. ChemLarp's MBAR
   must be unbiased (mean over replicates within three standard errors of the exact value) and its reported
   uncertainty must match the scatter between replicates.
B  pymbar (Shirts & Chodera's own implementation, an independent code) on the same matrices: free energies,
   their uncertainties and the overlap matrix, including a state that is evaluated but not sampled.
C  The alchemical Lennard-Jones particle: ChemLarp's MBAR on the matrix its own run produced, against pymbar
   on that matrix.
D  The same excess chemical potential by two different routes: MBAR over soft-core lambda windows, and Widom
   test-particle insertion into a separate simulation of the pure fluid. They share no samples, so their
   agreement over independent seeds checks the simulation and the coupling as well as the estimator.
"""
import argparse
import math
import os
import sys

import numpy as np
import pymbar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import chemlarp  # noqa: E402

failures = []


def check(label, ok, text):
    print("  %s %-48s %s" % ("ok " if ok else "BAD", label, text))
    if not ok:
        failures.append(label)


def oscillators(rng, n, Kk, x0, sampled=None):
    sampled = sampled if sampled is not None else [n] * len(Kk)
    xs = np.concatenate([x0[k] + rng.standard_normal(sampled[k]) / math.sqrt(Kk[k]) for k in range(len(Kk))])
    u = np.array([0.5 * Kk[k] * (xs - x0[k]) ** 2 for k in range(len(Kk))])
    return u, np.array(sampled)


def pymbar_result(u, nk):
    m = pymbar.MBAR(u, nk, verbose=False)
    r = m.compute_free_energy_differences()
    ov = m.compute_overlap()["matrix"]
    return np.asarray(r["Delta_f"]), np.asarray(r["dDelta_f"]), np.asarray(ov)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=6, help="independent seeds for part D (each run takes a few s)")
    a = ap.parse_args()
    print("pymbar", pymbar.__version__)

    # A  exact answers
    print("\nA  harmonic oscillators against the exact free energies")
    Kk = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    x0 = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    exact = 0.5 * np.log(Kk / Kk[0])
    rng = np.random.default_rng(11)
    est, err = [], []
    for rep in range(100):
        u, nk = oscillators(rng, 500, Kk, x0)
        r = chemlarp.mbar(u.tolist(), nk.tolist())
        est.append(r["f"])
        err.append([r["ddf"][0][k] if k else 0.0 for k in range(len(Kk))])
    est, err = np.array(est), np.array(err)
    mean, spread = est.mean(axis=0), est.std(axis=0, ddof=1)
    sem = spread / math.sqrt(len(est))
    bias_ok = all(abs(mean[k] - exact[k]) <= 3 * sem[k] for k in range(1, len(Kk)))
    check("unbiased (within 3 standard errors)", bias_ok,
          "f_4 - f_0 = %.4f (exact %.4f), standard error %.4f" % (mean[-1], exact[-1], sem[-1]))
    ratio = err[:, 1:].mean(axis=0) / spread[1:]
    check("reported error matches the scatter", all(0.8 <= x <= 1.25 for x in ratio),
          "reported / observed: " + " ".join("%.3f" % x for x in ratio))

    # B  pymbar on the same matrices
    print("\nB  pymbar %s on the same matrices" % pymbar.__version__)
    worst_f = worst_e = worst_o = 0.0
    for rep in range(10):
        u, nk = oscillators(rng, 400, Kk, x0)
        ours = chemlarp.mbar(u.tolist(), nk.tolist())
        df, ddf, ov = pymbar_result(u, nk)
        worst_f = max(worst_f, np.abs(np.array(ours["df"]) - df).max())
        worst_e = max(worst_e, np.abs(np.array(ours["ddf"]) - ddf).max())
        worst_o = max(worst_o, np.abs(np.array(ours["overlap"]) - ov).max())
    check("free energies", worst_f < 1e-6, "largest |difference| %.2e kT" % worst_f)
    check("uncertainties", worst_e < 1e-5, "largest |difference| %.2e kT" % worst_e)
    check("overlap matrix", worst_o < 1e-6, "largest |difference| %.2e" % worst_o)
    # A state with no samples (the middle one): free energy only, its uncertainty is not reported here.
    u, nk = oscillators(rng, 400, Kk, x0, sampled=[400, 400, 0, 400, 400])
    ours = chemlarp.mbar(u.tolist(), nk.tolist())
    df, ddf, _ = pymbar_result(u, nk)
    d = np.abs(np.array(ours["df"]) - df).max()
    check("unsampled state", d < 1e-6, "largest |difference| %.2e kT" % d)

    # C  pymbar on the Lennard-Jones matrices
    print("\nC  alchemical Lennard-Jones matrices, ChemLarp's MBAR against pymbar")
    for rho, T in ((0.8, 1.2), (0.5, 1.5)):
        r = chemlarp.lj_solvation(density=rho, temperature=T, production=20000, insertions=0, keep_matrix=True)
        u = np.array(r["u_kn"])
        nk = np.array(r["n_k"])
        df, ddf, ov = pymbar_result(u, nk)
        ours_mu, ours_err = r["mu_ex"]["mbar"]["kT"]["df"], r["mu_ex"]["mbar"]["kT"]["error"]
        check("rho* %.1f T* %.1f: mu_ex" % (rho, T), abs(ours_mu - df[0, -1]) < 1e-6 and abs(ours_err - ddf[0, -1]) < 1e-5,
              "ChemLarp %.6f +- %.6f, pymbar %.6f +- %.6f (kT; %d windows, %d samples)" %
              (ours_mu, ours_err, df[0, -1], ddf[0, -1], len(nk), int(nk.sum())))

    # D  MBAR against Widom insertion, independent seeds
    print("\nD  MBAR over lambda windows against Widom insertion in the pure fluid (%d seeds each)" % a.seeds)
    for rho, T in ((0.8, 1.2), (0.5, 1.5)):
        mb, wi = [], []
        for seed in range(1, a.seeds + 1):
            r = chemlarp.lj_solvation(density=rho, temperature=T, seed=seed)
            mb.append(r["mu_ex"]["mbar"]["kT"]["df"])
            wi.append(r["mu_ex"]["widom"]["kT"]["df"])
        mb, wi = np.array(mb), np.array(wi)
        diff = wi - mb
        se = diff.std(ddof=1) / math.sqrt(len(diff))
        check("rho* %.1f T* %.1f: Widom - MBAR" % (rho, T), abs(diff.mean()) <= 3 * se + 0.02,
              "mu_ex/kT MBAR %.4f, Widom %.4f, difference %.4f +- %.4f" % (mb.mean(), wi.mean(), diff.mean(), se))

    if failures:
        sys.exit("outside tolerance: " + "; ".join(failures))
    print("\nall checks passed")


if __name__ == "__main__":
    main()
