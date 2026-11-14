# Known gaps

Everything that is not 100% done, in one place, so it can be picked up later without re-reading the code.
Each entry: what is missing, where the code is, how to resume, and what data it needs. The measured accuracy
of finished work is in docs/ARCHITECTURE.md; this file only lists what is short of complete.

Update this file whenever a module is finished with known shortfalls, and strike entries (delete them) when
they are fixed. Roadmap modules not yet started are in data/roadmap.json (status "planned") and are listed at
the end for convenience.

Downloaded reference data is never committed. Each entry says where it came from so it can be fetched again.

---

## Proteins (`src/protein*.cpp`, `include/chemlarp/protein.hpp`)

| Gap | Now | To resume |
|---|---|---|
| Nucleic acids and ligands get no hydrogens or pKa | read and kept, listed under `unknown_residues` | Add DNA/RNA templates (DA, DC, DG, DT, A, C, G, U) to `residue_defs()` in `src/protein_internal.hpp` with H definitions; for ligands, build a molecule from HETATM coordinates (bond perception exists in the PDB molecule reader) and protonate with the aqueous/pKa machinery. Validate against neutron structures with DNA. |
| Arg, Tyr, Cys pKa = reference values | not fitted: Tyr 25, Cys 15 (mostly reactive), Arg 0 measurements in PKAD-2 | Needs more measured values (PKAD-2 mutants, literature Cys pKa sets). `scripts/fit_protein_pka.py` class 3 is skipped in `fit()`; re-enable when n ≥ ~80. |
| Large pKa shifts (> 1 unit) | cross-validated error 1.9 units (PROPKA 1.8) | Features are linear; try coupled-site titration (iterate neighbour charges at the pH) and a buried-pair desolvation term. Refit with `fit_protein_pka.py --write`; keep the CV-by-protein protocol. |
| pKa of a group ignores coupled titration | neighbours at reference charge | Iterate predict_pka with neighbour charges at pH (Tanford–Kirkwood style), or Monte Carlo over states. |
| Loops > ~6 residues, chain ends | built as placeholders (median 3.6 Å at 8 residues, 6.4 Å at 12) | Would need a knowledge-based loop score (DFIRE/statistical potential) or fragment library; current scoring is overlap only (`src/protein_build.cpp`, LoopResult scoring). |
| Side-chain packing, all removed | chi1 84.8% (one at a time 90.6%) | Backbone-dependent library (only a chi1 factor now, `kChi1Bias`), better energy (H-bond directionality, electrostatics), dead-end elimination instead of iterated improvement. Regenerate with `scripts/gen_protein_geometry.py`. |
| Undetermined Asn/Gln/His (~1/3 of them) | left as read, listed in `undetermined` | Fundamentally undecidable from structure alone; could use B-factor/electron density or crystal pH. Not a bug. |
| Ser/Thr OH within 45° of neutron H | 63.7% (same acceptor 91.5%) | Waters are treated as free acceptors/donors without orientation; optimise water hydrogens too. |

Reference data (re-download into a scratch folder; never into the repo):
- PKAD-2 measured pKa values: http://compbio.clemson.edu/PKAD-2 (page has `var wildData = [[...]]`; parse to `{pdb,res,chain,num,pka}`), plus the 161 PDB files from `https://files.rcsb.org/download/<ID>.pdb`.
- Neutron structures (20): 1L2K 1VCX 2ZOI 2ZWB 3A1R 3INS 3KYY 3X2P 4AR3 4FC1 4NY6 5CE4 5MOQ 5PTI 5ZN0 6EYM 6F2Q 6L46 7AZ9 7VOS as `.cif`.
- Side-chain/loop test set (held out): 1BPI 2TRX 1IGD 1POH 2CPL 1HRC 4ICB 1LZ1 7RSA 1XNB; weight-setting set 1CRN 1UBQ 1AKI 2LZT 1A2P 4PTI 3RN3 1PGA 2CI2.
- PROPKA 3.5.1 (`pip install propka`) for comparison.
- Scripts: `scripts/validate_protein.py`, `scripts/validate_protein_build.py`, `scripts/fit_protein_pka.py`, `scripts/gen_protein_geometry.py`.

## Stereo and identity

| Gap | Where | To resume |
|---|---|---|
| Planar and helical chirality (cyclophanes, metallocenes, helicenes) | roadmap module | See roadmap list below. |
| CIP: `[P@H]` phosphines and quinoid `c(=C)` input not parsed; rule-5 double bonds print E/Z, not seqCis/seqTrans | `src/cip*.cpp`, SMILES parser | 47/6312 enumerated stereoisomers differ from RDKit; list them with the CIP validation script and fix one class at a time. |
| InChI: isotopic stereo, allene stereo, /f and /r layers, polymers, some charge-migration/salt rules (diazonium, N-oxides, macrocycle metal complexes) | `src/inchi*.cpp` | Compare against the reference InChI on PubChem subsets that contain these features. |
| 2D layout: bridged polycycles legible but not template-quality | `src/depict*.cpp` | Add cage templates (adamantane, morphinan, bicyclo[2.2.2]). |
| Readers: query features, R-groups, S-group expansion, enhanced stereo not read | MOL/SDF reader | Add on demand. |

## Axial chirality (`AxialStereo`; `src/cip.cpp`, `src/smiles.cpp`, `src/canonical.cpp`, `src/molio.cpp`, `src/conformer.cpp`, `src/uff.cpp`)

Validated (docs/ARCHITECTURE.md "Axial chirality: validation"): 367/367 allene labels and 106/106 atropisomer labels agree with centres 1.5, with no wrong atropisomer.

| Gap | Now | To resume |
|---|---|---|
| InChI layer for allenes/atropisomers | InChI omits it, with a StereoIgnored warning (`src/inchi.cpp`, `inchi()`) | Standard InChI encodes allenes (even cumulenes) as a /t centre on the middle atom. Map AxialStereo to a Centre whose ligands are the terminal substituents; the parity code needs non-neighbour refs (`cut` check). Reference: RDKit cannot generate allene stereo, so get reference InChIs from the InChI Trust binary or PubChem records. |
| IUPAC names with Ra/Sa | NamingUnsupported (`src/iupac.cpp`, "axial chirality ... not part of the names") | Add the "(Ra)-"/"(Sa)-" prefix (P-93.5.1.1) and check that OPSIN accepts it; OPSIN may not parse axial descriptors, in which case use a hand-checked set. |
| Axial units from 2D wedges (MOL files) | only 3D coordinates are read (`Builder::add_stereo`, `dims == 3`) | Allenes: a wedge on a terminal substituent. Atropisomers: the wedge convention RDKit 2024+ and CDK use (a wedge on a bond next to the axis). Validate with CDK `StereoElementFactory.using2DCoordinates`. |
| 2D drawings of axial units | Ra/Sa shown as a label; no wedges drawn (`src/depict.cpp`) | Wedge one substituent on each allene terminal (MOL convention) once 2D input is read the same way, so drawings round-trip. |
| Atropisomers in SMILES | not written (there is no OpenSMILES notation), with a warning | Optional CXSMILES `wU`/`wD` extension, which RDKit and CDK read. |
| Atropisomer test is structural | 3 or more substituted ortho positions, aromatic ring atoms, non-ring single bond | Ring biaryls with sp2 non-aromatic ends, N–C axes (anilides) and 2-ortho biaryls with bulky groups are missed. A computed rotation barrier (UFF scan of the axis) would be the general test. |
| Cyclic allenes with a ring-closure digit at a terminal or the centre | read as for a tetrahedral centre; CDK 2.13 reads them differently and fails to round-trip one of them | Recheck against a fixed CDK release or OpenSMILES clarification; 3 of 734 strings affected (`scripts/validate_axial/run.sh`). |
| Axial units in CIP rules 4–5 | they give no auxiliary descriptors to other centres | Add axial nodes to `Digraph::compute_aux` (`src/cip.cpp`); needs test molecules with an allene inside a pseudo-asymmetric arrangement. |
| Distance-geometry bounds for allenes | none: the chiral volume keeps the sign and the UFF cumulene term makes the end groups perpendicular; `minimize=false` gives flat allenes | Add 1-5 distance bounds across the cumulene in `build_bounds` (`src/conformer.cpp`). |
| UFF energies of cumulenes | include the added cumulene torsion, so they differ from published UFF and RDKit | Intentional (geometry is right only with it); keep documented. |

Reference data (re-download): CDK 2.13 jar (`gh release download cdk-2.13 -R cdk/cdk -p cdk-2.13.jar`, sha256 5a42221d…2411); centres 1.5 from Maven Central (`com/simolecule/centres/centres-core/1.5`, `centres-cdk/1.5`). Run `scripts/validate_axial/run.sh <cdk.jar> <centres-core.jar> <centres-cdk.jar>`.

## Quantum chemistry

| Gap | To resume |
|---|---|
| UHF/UKS of NO in larger bases does not converge reliably (several nearby solutions) | Add level shifting / second-order SCF; see the QC row in ARCHITECTURE.md for the full list of unconverged cases. |

## Force fields and dynamics

| Gap | To resume |
|---|---|
| MMFF: no OPLS-AA, no distance-dependent dielectric, no periodic/solvated systems, no cutoffs (all pairs summed) | `src/mmff*.cpp`; cutoffs need a neighbour list (see MD). |
| MD: no periodic boxes for molecules, solvent, SHAKE/RATTLE, NPT barostat, electrostatics with UFF, neighbour lists, multiple time steps | `src/md*.cpp`. Neighbour lists are prerequisite for "Multithreading and sparse solvers" and FEP. |

## Threads and sparse solvers (`src/parallel.cpp`, `src/sparse.cpp`)

| Gap | Now | To resume |
|---|---|---|
| Serial kernels | only QC (integrals, J/K, DFT grid) is parallel | Next candidates by cost: conformer embedding attempts (`embed_species` loop over attempts; each attempt independent, keep the lowest-energy in attempt order), MD/MMFF non-bonded pair loops (partition by atom, private force arrays added in order), protein repair (loop trials), SMARTS screening over many molecules (`cl_search_json`), side-chain scoring. Use `parallel_for` / `parallel_sum` with a fixed partition count so results stay bit-identical. |
| QC scaling | 3.1× on 4 threads, 3.7× on 10 threads (72.5 s to 19.7 s) on an Apple M2 with 4 performance cores; not measured on a larger machine | The serial part is the SCF linear algebra (`symmetric_eigen`, `mul`, O(n³)) and the grid setup; parallelise `mul` by row blocks; the eigensolver needs a blocked Householder or divide-and-conquer. |
| No work stealing | one atomic counter per region | Enough for coarse loops; fine-grained nested parallelism would need the Chase–Lev design in docs/ARCHITECTURE.md history (git log). |
| WebAssembly is single-threaded | `__EMSCRIPTEN_PTHREADS__` not set | Needs `-pthread` builds and cross-origin isolation (COOP/COEP headers), which GitHub Pages cannot send; a service-worker header shim (coi-serviceworker) is the usual workaround. |
| Sparse LU only in BDF | kinetics with ≥ 60 species and a sparse Newton matrix | Distillation (block-tridiagonal, already O(N) by block Thomas) does not need it; an equation-oriented flowsheet solver would. The Newton matrix is still assembled dense (n² memory): pass a sparse Jacobian through `OdeJacobian` for very large mechanisms. |
| Ordering cost | exact minimum degree on an explicit elimination graph (std::set adjacency) | Fine to a few thousand unknowns; for larger systems implement AMD (Amestoy, Davis & Duff 1996) with quotient graphs. |
| Finite-difference Jacobian is dense | n right-hand-side evaluations when no analytic Jacobian is given | Curtis–Powell–Reid column colouring once a sparsity pattern is known (from the caller, since a pattern detected at one state can miss entries that are zero there). |
| ThreadSanitizer only on Linux CI | crashes on this macOS (Darwin 25) even for a trivial program | Re-check with a newer Xcode clang. |

## Thermodynamics and process

| Gap | To resume |
|---|---|
| NRTL: no parameter database; no modified UNIFAC (Dortmund); no vapour non-ideality; no liquid-liquid splits | Dortmund parameters are published (DDBST public set); LLE needs a stability test (tangent-plane) in the flash. |
| Pure components: no estimation methods (Chung, Lucas, Brock–Bird), no pressure corrections for vapour transport, no solids, no web UI | `src/pure*.cpp`. |
| Unit operations: no equilibrium or kinetic (CSTR/PFR) reactors, no real-gas vapour, no LLE | Reactors can reuse the kinetics module. |
| Flowsheet: no design specs/controllers, no equation-oriented solve, no Broyden acceleration, no minimal tear sets, no graphical editor | `src/flowsheet*.cpp`. |
| Aqueous: no redox (pe), no Pitzer (valid below I ≈ 0.5–1 mol/kg), water mass fixed, no pressure effects, kinetics, surface complexation, ion exchange, organic acids | phreeqc.dat / pitzer.dat (public domain, USGS). |
| Distillation: only the current specifications; no purity/recovery/boilup/duty specs, side draws, pumparounds, stage efficiencies, two liquids, non-condensables, shortcut FUG design | `src/distillation*.cpp`. |

## Roadmap modules not started (data/roadmap.json)

| Module | Size | Notes for starting |
|---|---|---|
| Planar and helical chirality (cyclophanes, metallocenes, helicenes) | S | Split from the axial module on 2026-11-14. Rp/Sp with a pilot atom (P-93.5.1.2) and M/P helicity; from 3D only (no SMILES notation). centres 1.5 labels neither, so the reference has to be literature descriptors of known compounds ([2.2]paracyclophanes, planar-chiral ferrocenes, hexahelicene). The AxialStereo machinery (3D parity, 3D perception in `Builder::add_stereo`) can be reused. |
| Docking | M | Needs protein module (done) + MMFF/UFF ligand; Vina-like empirical score, grid maps, Monte Carlo + BFGS. Validate redocking RMSD on PDBbind core set or CASF (licence check before download). |
| Free-energy perturbation (FEP/MBAR) | M | Needs MD with solvent, periodic boxes and electrostatics (see MD gaps). Validate MBAR on analytic harmonic oscillators, then hydration free energies vs FreeSolv. |
| Interactive lab simulations and safety rules | M | Web feature on top of kinetics/thermo; GHS hazard rules need a data source (public PubChem GHS). |
| Native installers (.apk, .exe, .dmg) | S | CI packaging; signing needs the user's certificates. |
| Packages (pip, Homebrew, winget, Node) | S | Publishing needs the user's accounts; prepare manifests only. |
| Validation against reference data | M | Collect the per-module validation scripts into one runner and a report. |
