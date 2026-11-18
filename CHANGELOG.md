# Changelog

All notable changes to ChemLarp. The format follows [Keep a Changelog](https://keepachangelog.com/1.1.0/),
and versions follow [semantic versioning](https://semver.org/).

## [Unreleased]

### Added

- **Protein–ligand docking** (`include/chemlarp/docking.hpp`): a flexible ligand placed into a rigid
  receptor. The scoring function has the form published for AutoDock Vina and is written from the paper —
  gauss1, gauss2, repulsion, hydrophobic and hydrogen-bond terms in the surface distance over heavy-atom
  pairs within 8 Å, X-Score atom types and radii, the affinity divided by the rotatable-bond factor. The
  receptor's share is precomputed on 0.375 Å grid maps and read back by trilinear interpolation; the ligand
  comes from SMILES (a conformation is generated) or from MOL/SDF (it is kept) and is varied over its
  position, orientation and the torsions of a torsion tree; the search is Monte Carlo with a local L-BFGS
  after every step on a landscape whose repulsion is squashed, Metropolis acceptance and independent runs in
  parallel, and the poses are clustered by heavy-atom RMSD with the atoms matched by graph isomorphism.
  New `chemlarp dock` command, `cl_dock_json` in the C API and `dock()` in Python. Redocking of twelve
  crystal complexes is reported in `docs/ARCHITECTURE.md`, including the complexes it gets wrong.
- Validation suite `docking` (`scripts/validate_docking.py`): the deposited pose of each complex against the
  top-ranked docked pose, with each failure attributed to the search or to the scoring function. One complex
  (PDB 3PTB, trypsin with benzamidine) is committed under `data/dock` and redocked by the test suite, so
  every CI job checks a real complex end to end without downloading anything.
- Two more entry points in the fuzz target (docking receptors and ligands) and docking cases in
  `tests/test_robustness.cpp`.

- **Explicit water and hydration free energies** (`include/chemlarp/solvation.hpp`,
  `include/chemlarp/ewald.hpp`): rigid TIP3P water with SETTLE constraints, Ewald electrostatics with smooth
  particle-mesh Ewald over a mixed-radix FFT written for the project, Langevin BAOAB dynamics at 2 fs, a Monte
  Carlo barostat and an analytic long-range dispersion correction. A neutral solute from SMILES is decoupled
  over 19 λ windows (charges first, then soft-core van der Waals) and analysed with MBAR, with BAR and
  thermodynamic integration as checks. New `chemlarp water` and `chemlarp hydrate` commands,
  `cl_water_box_json` and `cl_hydration_json` in the C API, and `water_box()` and `hydration_free_energy()` in
  Python.
- Validation suite `hydration` (`scripts/validate_hydration.py`): bulk TIP3P against published values, the
  model's own water-in-water free energy, and FreeSolv small molecules against experiment.
- `docs/AUDITS.md`: the correctness, security, integrity, git and licensing audits, with findings.
- `SECURITY.md`, `CONTRIBUTING.md`, `LICENSING.md`, `NOTICE`, this changelog.
- A Content-Security-Policy and a "Source code" link in the web app.
- README rewritten with screenshots and animations of every tool, an accuracy table and a comparison with the
  usual commercial packages; `INSTALL.md`, a step-by-step installation guide for non-technical users.

### Changed

- **Licence: AGPL-3.0-only** from this version, with commercial licences available from the copyright holder
  (see `LICENSING.md`). Version 0.1.0 and everything before it remain available under Apache-2.0.
- Every source file carries an SPDX licence header.
- GitHub Actions in the workflows are pinned to commit SHAs instead of tags.

### Fixed

- Aqueous solver: mass-balance rows are scaled by the largest term, so a phase that takes nearly all of a
  component (calcite near 25 °C) no longer reports false non-convergence.

## [0.1.0] - 2026-11-15

First release, under the Apache License 2.0.

### Added

- Structures and identifiers: OpenSMILES parser with valence checks and an error vector, canonical and random
  SMILES, relevant-cycle ring perception and Hückel aromaticity, CIP R/S and E/Z including axial and planar
  chirality, standard InChI and InChIKey, systematic IUPAC names, SMARTS substructure search, Morgan/ECFP
  fingerprints, readers and writers for MOL V2000/V3000, SDF, PDB, mmCIF and CML, 2D depiction with wedges.
- Molecular modelling: distance-geometry 3D embedding with stereo checked in 3D, UFF and MMFF94/MMFF94s force
  fields with analytic gradients, molecular dynamics (velocity Verlet, Bussi, Langevin), normal modes and
  vibrational spectra, Hartree–Fock and DFT (SVWN, PBE, B3LYP) single points.
- Proteins: PDB and mmCIF readers, side-chain and loop repair, pKa prediction from structure, protonation at a
  pH and hydrogen placement with hydrogen-bond optimisation.
- Thermodynamics and kinetics: NASA-7 thermochemistry, reaction balancing and ΔrG°/K, Gibbs-minimisation
  equilibrium, stiff kinetics with `.ocm`, CHEMKIN-II and Cantera YAML mechanisms and a variable-order BDF
  integrator.
- Process and solution chemistry: cubic equations of state, NRTL and UNIFAC activity models, flash and unit
  operations, flowsheets with recycle convergence, rigorous distillation columns, aqueous speciation with the
  PHREEQC database, and a virtual lab with safety rules and GHS data.
- Free energies: EXP, BAR, MBAR with uncertainties and overlap, thermodynamic integration, and an alchemical
  Lennard-Jones particle with Widom insertion as a second route.
- Platforms: C++ library with a C ABI, `chemlarp` command-line program, Python bindings, a WebAssembly web app
  at https://ayuroy01.github.io/ChemLarp/, and macOS, Windows and Android applications.

[Unreleased]: https://github.com/ayuroy01/ChemLarp/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/ayuroy01/ChemLarp/releases/tag/v0.1.0
