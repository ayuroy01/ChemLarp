# OmniChem Core: architecture

## Provenance: clean-room, literature-based

This project was scoped against ChemDraw, Labster, Schrödinger (Glide/Jaguar), and Aspen Plus/HYSYS. None of those products were decompiled or reverse-engineered.

- No installed copies were available to analyze. Labster is hosted SaaS, and Aspen is Windows-only.
- Every capability listed is already documented in open literature, so nothing needs to be recovered from a binary.
- Decompiling commercial chemistry suites to build a competitor would break their EULAs. It would also put the open-source license at risk.

Each algorithm cites its public source in its header. Use the same rule for every new module: **cite the paper or standard, never the product.**

| Capability | Public basis |
|---|---|
| SMILES parsing, aromaticity | OpenSMILES spec (opensmiles.org); Weininger, *JCICS* 28 (1988) 31 |
| InChI | IUPAC/InChI Trust technical manual + open-source reference library (MIT-style license) |
| R/S, E/Z | IUPAC 2013 recommendations P-92 (CIP rules); Hanson et al., *JCIM* 58 (2018) 1755 (CIP algorithm validation) |
| 2D→3D embedding | Crippen & Havel, *Distance Geometry and Molecular Conformation* (1988); Riniker & Landrum ETKDG, *JCIM* 55 (2015) 2562 |
| Force fields | MMFF94: Halgren, *J. Comput. Chem.* 17 (1996), 5 parts; UFF: Rappé et al., *JACS* 114 (1992) 10024; OPLS-AA: Jorgensen et al., *JACS* 118 (1996) 11225 |
| MD integrators | Verlet (1967); Swope et al. velocity Verlet, *JCP* 76 (1982) 637 |
| Docking scoring | AutoDock Vina, Trott & Olson, *J. Comput. Chem.* 31 (2010) 455 (open, Apache-2) |
| FEP / MBAR | Zwanzig (1954); Shirts & Chodera, *JCP* 129 (2008) 124105 |
| Thermo data | GRI-Mech 3.0; NASA Glenn (McBride, Zehe & Gordon, NASA/TP-2002-211556); Burcat database |
| Gibbs minimisation | White, Johnson & Dantzig, *JCP* 28 (1958) 751; Reynolds STANJAN (1986); Smith & Missen, *Chemical Reaction Equilibrium Analysis* (1982) |
| Kinetics | modified Arrhenius / mass action; CHEMKIN-II format (Kee et al., Sandia SAND89-8009) |
| Activity models / EoS | NRTL: Renon & Prausnitz (1968); UNIFAC: Fredenslund et al. (1975), group tables published by DDBST; Peng–Robinson (1976) |
| Flash / distillation | Rachford & Rice (1952); Michelsen, *Fluid Phase Equil.* 9 (1982) (stability); Naphtali–Sandholm (1971) |
| Flowsheet convergence | Westerberg et al., *Process Flowsheeting* (1979); Biegler, Grossmann & Westerberg (1997) |

## Data flow

```
 SMILES / PDB ──▶ Topology graph ──▶ Thermo state ──▶ Force-field minimiser ──▶ Continuity solver
  smiles.cpp      system_state.hpp     thermo.cpp        (planned: ff/)           (planned: process/)
                  formula.cpp          joback.cpp
                                       gibbs.cpp
                                       reaction.cpp
            └──────────────── one ChemicalSystemState, no copies between stages ────────────────┘
```

Every stage reads and writes the same `ChemicalSystemState`, so data passes between stages without copying:

- **SoA columns** (`AtomTable`, `BondTable`) are `AlignedVector<T>` with 64-byte alignment. Per-atom kernels such as `x[i] += dt*v[i]` auto-vectorize with `-O2 -march=native`; there are no hand-written intrinsics yet.
- **Contiguous species slices.** Species *k* owns `[atom_begin, atom_begin+atom_count)`, so a worker receives a `std::span` and never needs a copy.
- **CSR adjacency** (`Adjacency`) gives O(1) neighbour iteration with no per-atom heap allocation. It is rebuilt in O(N+M) after topology edits.
- **`revision`** is bumped on every mutation. Downstream caches (Newton warm starts, neighbour lists, FF parameter assignment) key on it.
- **C ABI** (`omnichem_c.h`) exposes raw column pointers. Python wraps them with `numpy.ctypeslib.as_array`, which is verified zero-copy: a write from NumPy is visible in C++. Node can do the same with `ArrayBuffer` via N-API.

## Identity: canonical SMILES and CIP

`canonical.cpp` produces one string per structure, stereo included. Every species carries `canonical_smiles`, set at parse time, and the thermo database is keyed by it. The pipeline:

1. **Normalised view.** A plain `[H]` with a single heavy neighbour is folded into that atom's H count. Aromaticity is re-perceived from the Kekulé structure, so input notation never leaks into identity.
2. **Ranking.** Initial atom invariants are (degree, Z, isotope, charge, H count, ring, aromatic, CIP label, atom class). They are refined iteratively by sorted (neighbour rank, bond code) multisets until stable. While ties remain, the lowest tied class is split and refinement runs again (Weininger 1989).
3. **Writing.** DFS starts from the lowest rank, neighbours are visited in rank order, and ring closures take the lowest free digit. Tetrahedral parity is recomputed as the permutation parity between the stored reference order and the output neighbour order. `/` `\` marks come from a small constraint search, so conjugated systems (shared single bonds) stay consistent; the choice is anchored on output order only.
4. **Stereo is part of identity.** CIP labels are rank invariants, so meso and chiral diastereomers, and the members of meso pairs, are distinguished deterministically. Stereo marks that CIP rules 1–2 show are not stereogenic are dropped, with a warning.

How it is verified: the randomized writer (`SmilesWriteOptions{canonical=false}`) emits arbitrary atom orders in aromatic or Kekulé form. The test suite requires `canonical(random(x)) == canonical(x)` for 33 molecules × 40 seeds; an ad-hoc run of 14 harder molecules × 1000 seeds also passed.

Stereo storage (`TetrahedralStereo`, `DoubleBondStereo`) keeps references to explicit neighbours plus implicit-H and lone-pair markers. It never depends on write order, so the same structure can be re-emitted in any order.

## Platforms

Windows, macOS and Linux are first-class targets, and the code base follows a few rules to stay that way:

- **ISO C++20 and the standard library only.** No compiler builtins (`std::countr_zero` instead of `__builtin_ctzll`), no POSIX calls, no `<windows.h>` outside `src/cli/main.cpp`, where it is used only to switch the console to UTF-8.
- **One exported C ABI.** `OMNICHEM_API` expands to `__declspec(dllexport/dllimport)` on Windows and to default visibility elsewhere. `OMNICHEM_BUILD_DLL` and `OMNICHEM_STATIC` select the mode, and CMake sets them per target.
- **Builds.** CMake is canonical (MSVC: `/utf-8 /permissive- /W4`). The Makefile covers macOS, Linux and MSYS2 MinGW. MinGW DLLs are linked with static libstdc++, so the DLL is self-contained for Python's `ctypes`.
- **Files and text.** Mechanism files are read in text mode and accept CRLF. `.gitattributes` pins line endings for scripts. Tests resolve `data/` through `OMNICHEM_TEST_DATA_DIR`, so CTest can run from any build directory.

## Error vector

No stage throws for bad chemistry. Every stage appends a `Diagnostic{code, severity, position, atom, species, element, delta, message}` to a shared `Diagnostics` object.

- Parsing is **transactional**. On error the state is untouched and *all* errors are reported, both syntax errors (1xx) and independent chemistry errors (2xx).
- Codes are stable integers: 1xx syntax, 2xx valence/aromaticity, 3xx stoichiometry, 4xx thermo/equilibrium. Front ends can switch on them.
- Positions are character offsets in the original input, including SMILES embedded in a reaction equation, so the CLI can draw a caret under the offending character.

## Threading model (designed, not yet implemented)

The current code is single-threaded and its state is not shared across threads. The planned design:

1. **Reaction dependency graph.** Nodes are reactions; an edge joins two reactions that share a species. Connected components, and independent sets inside a component for split-operator kinetics, are the units of parallel work.
2. **Work-stealing pool.** There is one Chase–Lev deque per worker. Tasks are `(kernel, species span)` pairs. Chase–Lev (SPAA 2005) is lock-free, with a CAS only on steal.
3. **Ownership rather than locks.** A task gets exclusive `std::span<double>` access to the SoA columns of its species slice. A reduction across tasks, such as total energy, writes to per-worker accumulators that are combined at the end. No mutexes are used on the hot path.
4. **Determinism.** Reductions run in a fixed tree order, so results do not depend on thread count.

Design rule: build a correct serial kernel first, then parallelize it behind the same API. Lock-free code is held to a higher testing bar (TSan plus stress tests) before it merges.

## Status

| Module | Status |
|---|---|
| OpenSMILES parser (organic subset, brackets, rings, branches, charges, isotopes, `@`/`@@`, `/` `\`) | **done**, tested |
| Kekulization, implicit H, charge-aware valence validation | **done**, tested |
| Ring perception: bridge test, SSSR (GF(2) basis over shortest cycles) | **done** |
| Aromaticity perception (Hückel 4n+2 over SSSR rings and fused-pair envelopes) | **done**; Kekulé and aromatic inputs give identical state |
| Canonical SMILES (invariant refinement + tie-breaking, stereo-aware) | **done**; 14,000-case randomized round-trip fuzz with zero mismatches |
| Randomized SMILES writer (aromatic or Kekulé) | **done** (used as the test oracle) |
| CIP R/S and E/Z, rules 1a + 2 (hierarchical digraph, BFS by sphere) | **done**; rules 1b, 3-5 (r/s, like/unlike, seqcis) not implemented and reported as `CipUndetermined` |
| Stereo error vector: conflicting `/` `\` marks, small-ring E/Z, non-tetrahedral `@` | **done** |
| Thermo database keyed by canonical SMILES (formula fallback warns `ThermoFormulaMatch`) | **done** |
| Composition, Hill formula, molar mass | **done** |
| Reaction parser, element/charge balance error vector, ΔrH°/ΔrS°/ΔrG°/K | **done**, tested |
| NASA-7 thermo (10 GRI-Mech species) | **done**, verified against JANAF to <0.6 kJ/mol |
| Joback ΔfH°/ΔfG° from graph | **done**; ethanol, propane, acetone, benzene, and acetic acid fall within the tolerances shown in tests |
| Ideal-gas Gibbs minimisation (element-potential Newton) | **done**, tested (K consistency, conservation, Le Chatelier, rank deficiency) |
| C ABI + Python ctypes binding; `.ocx` v0.1 container | **done** |
| InChI | planned (canonical SMILES now covers structure identity inside OmniChem) |
| Allene / atropisomer / square-planar stereo (@AL, @SP, @TB, @OH) | planned |
| Explicit hydrogens with stereo remapping | **done** |
| Distance-geometry 3D embedding (bounds + triangle smoothing, metric matrix in 4D, chiral/planar refinement) | **done**; stereo read back from 3D and enforced, best-of-N conformers |
| UFF force field (bond, angle, torsion, inversion, vdW; analytic gradient) + L-BFGS | **done**; gradients exact vs finite differences, geometries vs experiment in tests |
| SDF (V2000) and XYZ output | **done** |
| Multiphase / condensed-phase equilibrium, activity models | planned |
| Mechanism format `.ocm` (SMILES-bound species, CHEMKIN-style reactions, units, third-body efficiencies) | **done**, with error vector (5xx) |
| Mass-action rates, modified Arrhenius, third bodies, detailed-balance reverse rates, analytic Jacobian | **done**; relaxes exactly to the Gibbs equilibrium (tested) |
| Batch reactors: isothermal/adiabatic × constant V/P, ignition delay | **done** |
| Variable-order (1–5) variable-step BDF with dense output | **done**; Robertson, Prothero–Robinson, ignition benchmarks |
| Fall-off (Lindemann, Troe), chemically activated (HIGH), specific colliders `(+AR)`, PLOG, explicit REV, DUPLICATE | **done**, with analytic Jacobians (Troe dF/dPr, PLOG dk/dP) |
| CHEMKIN-II import (chem.inp + fixed-column NASA-7 therm.dat, inline THERMO, optional NAME→SMILES map) | **done**; composition-only species when no SMILES is given |
| Cantera YAML import (NASA7; elementary, three-body, falloff/Troe, chemically-activated, pressure-dependent-Arrhenius) | **done** (built-in YAML subset reader) |
| Cross-format check: one mechanism as .ocm / CHEMKIN / YAML gives identical rates | **done** (`scripts/gen_h2o2_mechanism.py`, test `three_formats_give_identical_rates`) |
| SRI fall-off, Chebyshev, FORD/RORD, NASA9 thermo, ionic species | planned (`UnsupportedReactionType` today) |
| Sparse LU for large mechanisms (>~300 species) | planned; dense LU is O(S³) per step |
| MMFF94 / OPLS force fields; partial charges & electrostatics; velocity-Verlet MD | planned (UFF done) |
| ETKDG-style torsion preferences; conformer ensembles with RMSD pruning | planned |
| Docking scoring, FEP/MBAR | planned |
| Peng–Robinson and SRK cubic EoS, mixture fugacities, 30-compound critical table | **done**; vapour pressures within 0.1–0.7% of NIST for non-polar compounds |
| Michelsen stability test + Rachford–Rice two-phase flash; pure-component vapour pressure | **done** |
| NRTL / UNIFAC activity models; column & flowsheet solvers | planned |
| Installable web app (WebAssembly, PWA, GitHub Pages) | **done**: https://ayuroy01.github.io/ChemLarp/ |
| Worker pool | planned (design above) |
| Node.js addon | planned (C ABI is ready for it) |
| IUPAC naming | planned, and the hardest item on this list. It needs the full P-rules engine; there is no shortcut "hash map". |
