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
| Canonical SMILES of a mancude macrocycle can drop one double-bond configuration: an exocyclic C=C whose ring neighbour is written aromatic loses its `/` mark (InChI b layer shows `?` where the input has a configuration). 8 of 1,500 molecules in the validation sample, all stereoisomers of one chlorophyll-like macrocycle | `src/canonical.cpp` (double-bond marks), `src/rings.cpp` (aromatic perception) | Reproduce with `scripts/validate_rdkit.py`, class "stereo, other": `C=CC1=C(C)C2=NC/1=C/C1=NC(=C(CC)C/1=C/O)/C=C1N=C3C(=C\1C)C(=O)[C@@H](C(=O)OC)/C3=C1/N/C(=C\2)[C@H](C)[C@H]1CCC(=O)OC`. The reference atom of the double bond is a ring atom that becomes aromatic, so the writer has no non-aromatic neighbour to hang the mark on; keep the Kekule neighbour for the mark, or write such rings in Kekule form. |
| CIP: `[P@H]` phosphines and quinoid `c(=C)` input not parsed; rule-5 double bonds print E/Z, not seqCis/seqTrans | `src/cip*.cpp`, SMILES parser | 47/6312 enumerated stereoisomers differ from RDKit; list them with the CIP validation script and fix one class at a time. |
| InChI: isotopic stereo, allene stereo, /f and /r layers, polymers, some charge-migration/salt rules (diazonium, N-oxides, macrocycle metal complexes) | `src/inchi*.cpp` | Compare against the reference InChI on PubChem subsets that contain these features. |
| 2D layout: bridged polycycles legible but not template-quality | `src/depict*.cpp` | Add cage templates (adamantane, morphinan, bicyclo[2.2.2]). |
| Readers: query features, R-groups, S-group expansion, enhanced stereo not read | MOL/SDF reader | Add on demand. |

## Planar and helical chirality (`PlanarStereo`; `src/planar.cpp`, `src/cip.cpp`, `src/molio.cpp`)

Done 2026-11-14: cyclophanes (pilot-atom Rp/Sp), metallocenes (Schlögl Rp/Sp), helicenes (P/M) from 3D files; checked against 3 COD crystal structures and mirror images (docs/ARCHITECTURE.md, "Planar and helical chirality: validation").

| Gap | Where | To resume |
|---|---|---|
| Validation covers 3 compounds ([2.2]paracyclophane, two ferrocenes); helicenes only against the torsion definition | `tests/test_planar.cpp`, `data/planar_chirality_cod.sdf` | More COD entries with descriptors in their names (download permission needed): ferrocenes 4073056-4073059 (pR), 2220286 (Sp), 2220288 (Rp), 7000979 (Sp,Sp), 7020919 (Sp), 2015329 (SFc); search `https://www.crystallography.net/cod/result?format=json&text=ferrocen&year=YYYY` and filter names for Rp/Sp/pR/SFc. Convert with `scripts/validate_planar/cif2sdf.py` (caps each Cp ring with an H away from the metal, RDKit DetermineBondOrders charge 0, removes the cap as C−). Need helicene and arene-Cr(CO)3 entries with stated P/M or Rp/Sp. |
| No notation: SMILES, InChI and IUPAC names cannot carry the units (warning StereoIgnored / NamingUnsupported); both enantiomers have the same canonical SMILES | `src/canonical.cpp`, `src/inchi.cpp`, `src/iupac.cpp` | None exists in SMILES or standard InChI; names would need Rp/Sp or P/M with locants (P-93.5.1.4). |
| Embedding from SMILES gives a random enantiomer of a planar or helical unit; units carried into a new species are not enforced as constraints | `src/conformer.cpp` | Add a chirality-plane constraint (sign of a pilot–ring volume) to the distance-geometry refinement, and an API to request a descriptor. |
| A plane counts only if it carries a substituent or ring heteroatom: an unsubstituted deck in a chiral molecule (e.g. the second deck of 4-bromo[2.2]paracyclophane) is not labelled, as in usual practice, but CIP would distinguish it | `perceive_planar_stereo` | Decide from the full CIP digraph plus a symmetry test of the deck. |
| Configurational stability not judged: large ansa rings that let the ring rotate through the loop, long-bridged [n]paracyclophanes, low-barrier helicenes are all labelled; [4]helicene is excluded by the five-ring minimum | `src/planar.cpp` | Use bridge length / ring size thresholds from the literature (e.g. [n]paracyclophane racemisation barriers). |
| Cyclophane bridges must start with a non-aromatic atom at least 0.25 Å off the ring plane; only aromatic 5/6-rings are planes | `perceive_planar_stereo` | Extend when a real case needs it (biphenylophanes, heteroaromatic decks larger than 6). |
| Metallocenes: Schlögl's convention only (the IUPAC 2013 alternative treats each ring carbon as a centre and can differ); substituents ranked on the ionic form (metal bonds cut); arene complexes use the same rule without a reference check | `src/cip.cpp` (`metallocene_rp`) | Add the P-93.5.1.4.2 per-atom description as an option; check against an arene-Cr(CO)3 crystal structure. |
| Hapto bonds cut on reading: a carbon left with three bonds and no multiple bond becomes C− and the metal gains +1 (Cp− / Fe2+). Other organometallic bonding patterns (allyl, η4-dienes, carbonyl clusters) get that rule too and may not give the intended charge split | `build_species` in `src/molio.cpp` | Compare against InChI's metal disconnection on a set of organometallic SDFs. |
| `cip_compare_atoms` copies the whole state per comparison (adds a dummy root) | `src/cip.cpp` | Give the Digraph a two-root mode. |
| Helicenes: aromatic 5/6-rings only; a double helicene is reported as separate chains; the unit is not tied to helicene numbering | `perceive_planar_stereo` | |
| Output: the CLI `read` command does not print the units; the web app lists them without atom numbers (it re-parses the SMILES, whose atom order differs) and the drawing has no Rp/Sp/P/M notes there | `src/cli/main.cpp`, `web/app.js` | Map file atom order to the re-parsed SMILES, or draw from the file's state. |

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

## Free energies (`src/free_energy.cpp`, `include/chemlarp/free_energy.hpp`)

Done 2026-11-14: EXP, BAR, MBAR (uncertainties, overlap), TI, statistical inefficiency; soft-core LJ particle
in the LJ fluid with Widom insertion as a second route. MBAR matches pymbar 4.0.3 to 5e-10 kT.

| Gap | Where | To resume |
|---|---|---|
| No binding free energies (protein–ligand) | — | Needs docking poses and restraints on top of `solvation.hpp`; see the next section for the solvent side. |
| BAR summed over neighbouring windows adds the pair errors as if independent; neighbouring pairs share a window's samples, so the error was 0.49 of the scatter between seeds | `lj_solvation` | Report MBAR's error for the BAR route too, or compute the covariance of adjacent pair estimates. |
| TI error 0.65 of the scatter: the statistical inefficiency of the slowest windows (λ ≈ 0.25, g ≈ 20) is underestimated from 4,000 frames | `statistical_inefficiency` | Longer runs, or block averaging as a second estimate of g; pymbar's warning applies equally. |
| Widom insertion error is about twice the scatter (too conservative) | `lj_solvation` | Compute the error from blocks of frames instead of the delta method on the per-frame means. |
| TI uses the trapezoidal rule only | `thermodynamic_integration` | Add Simpson or cubic-spline integration for unevenly spaced λ. |
| μ_ex of the truncated-and-shifted LJ fluid not compared with a literature equation of state | `scripts/validate_fep.py` | Implement the Thol et al. (Int. J. Thermophys. 36 (2015) 25) LJTS equation of state (coefficients from the paper) and compare a few state points. |
| The LJ simulation is O(N²) without neighbour lists; one solute; no replica exchange between λ windows | `src/free_energy.cpp` | Cell lists for N > 500; Hamiltonian replica exchange for windows with poor overlap. |
| Uncertainties of unsampled states are reported as NaN | `mbar` | Use the full N × N form of Θ restricted to the needed columns (pymbar does). |
| No web UI | `web/` | A panel that runs `lj_solvation` and plots ∂u/∂λ and the overlap matrix. |

## Explicit water and hydration free energies (`src/solvation.cpp`, `src/ewald.cpp`, `include/chemlarp/solvation.hpp`)

Done 2026-11-15: TIP3P with SETTLE, PME, Langevin BAOAB, Monte Carlo barostat; neutral solutes from SMILES
decoupled in 19 windows; MBAR/BAR/TI. Validation: `scripts/validate_hydration.py` (suite `hydration`, needs
FreeSolv's `database.txt` from https://github.com/MobleyLab/FreeSolv, v0.52, under `freesolv/` in the data
directory; results are cached in `freesolv/runs/`).

| Gap | Where | To resume |
|---|---|---|
| TIP3P water in TIP3P water gives −6.07 ± 0.18 and −6.11 ± 0.22 kcal/mol (two seeds); the published TIP3P values found, −5.4 and −5.6, used cut-off electrostatics and 216 molecules, so the 0.5 kcal/mol difference is not resolved either way | `scripts/validate_hydration.py` part B | Find a PME value for TIP3P self-hydration (water-model papers' supporting information), or run the same decoupling in OpenMM/GROMACS (TIP3P, PME, 9 Å, LJ tail correction) and compare like for like. |
| **The solute's Lennard-Jones parameters are wrong for solvation, so hydration free energies of aromatics and carbonyls are 2.5-3 kcal/mol too negative.** UFF's well depths were fitted to geometries, not condensed phases: against GAFF's, UFF's epsilon is 2.9x larger for hydrogen (0.044 vs 0.015 kcal/mol), 1.2x for aromatic carbon and 0.29x for carbonyl oxygen, so hydrogen-rich and aromatic solutes attract water far too strongly (benzene -4.00 against -0.90 measured; acetone -6.35 against -3.80). Alkanes and alcohols are within 1.3 kcal/mol. The free-energy machinery is not implicated: see "Hydration free energies: validation" in docs/ARCHITECTURE.md | `add_solute` in `src/solvation.cpp` (`assign_uff_types`, `uff_parameters`) | Implement GAFF atom types with AM1-BCC charges (AM1 is not in the QC module, so semiempirical support would come first), or OPLS-AA, and revalidate against FreeSolv. A cheaper interim step: use a published Lennard-Jones set for the common elements instead of UFF's, and state plainly that charges remain MMFF94. Do not tune parameters against FreeSolv: that would fit the test set. |
| FreeSolv's calculated column is GAFF/AM1-BCC, so ChemLarp's numbers can only be compared with experiment, never used to reproduce FreeSolv's own calculations | `scripts/validate_hydration.py` part C | Follows from the row above. |
| Charged solutes are rejected: no finite-size or net-charge corrections, no ion parameters | `add_solute` | Add the Rocklin et al. (J. Chem. Phys. 139 (2013) 184103) corrections and Joung–Cheatham ions. |
| Kinetic temperature at 2 fs is 296–297 K for a 298.15 K thermostat (BAOAB's known O(Δt²) bias in the momenta; 299 ± 2 K at 1 fs). The configurational temperature, which the free energies depend on, was not measured | `Dynamics::temperature` | Compute the configurational temperature ⟨|∇U|²⟩/⟨∇²U⟩ for the water box, or report the kinetic temperature at the O step. |
| Speed: about 5.6 ms per step on one thread for 1,200 atoms; the neighbour list is O(N²); windows run one per thread, so 19 windows on 8 threads leave threads idle in the last round | `Engine::compute`, `build_list`, `hydration_free_energy` | Cell lists, molecule-based water–water kernels (9 site pairs with one O–O cutoff test), SIMD; let a window use idle threads or pick a window count that fills the pool. |
| Cubic boxes only, one water model (TIP3P), no constraints on the solute's bonds, no replica exchange between windows | `src/solvation.cpp` | Rectangular boxes are a small change (PME already takes three edges); TIP4P-Ew needs a virtual site; SHAKE for X–H bonds would allow 2 fs without hydrogen-mass repartitioning. |
| No web UI and not in the web build | `web/` | Too slow for the browser at this size; a panel could show stored results or run a short water box. |

## Packages (`pyproject.toml`, `setup.py`, `release/`, `.github/workflows/packages.yml`, docs/PACKAGING.md)

Done 2026-11-14: wheels for five platforms, npm (WebAssembly), Homebrew formula, winget manifests, CLI archives,
draft release on a tag. All built and installed in CI; macOS wheels also tested locally.

| Gap | Where | To resume |
|---|---|---|
| Nothing is published: PyPI, npm, the Homebrew tap and winget-pkgs need the owner's accounts | docs/PACKAGING.md, "Releasing" | Tag a version, publish the draft release, then the four manual steps there (or add a PyPI trusted-publisher job). |
| The winget manifests and the Homebrew stable URL have never been used: they need a published release (the formula is only tested as `--HEAD`) | `release/winget`, `release/homebrew` | After the first release: `wingetcreate validate` on the filled manifests, `brew install` (not `--HEAD`) from the tap. |
| The release job (draft release, checksums) runs only on a tag, so it has not run yet | `.github/workflows/packages.yml`, job `release` | Check it on the first tag; push a throwaway tag such as `v0.0.0-test` and delete the draft if needed. |
| No code signing: macOS CLI archive not notarized, Windows exe not Authenticode-signed | CI | Needs an Apple Developer ID and a Windows code-signing certificate from the owner. |
| The Windows exe is linked with the static MSVC runtime so that it needs no redistributable, but this is only exercised on CI runners, which have the runtime installed anyway | `cli` job | Run `dumpbin /dependents chemlarp.exe` in a Developer shell and check that no VCRUNTIME/MSVCP DLL is listed. |
| No Windows on Arm, musl Linux, 32-bit or Python < 3.10 wheels; no conda-forge recipe | `packages.yml` | Add matrix entries / a feedstock if users ask. |
| The npm package exposes 16 of the C API's functions (the web worker's set plus MBAR) | `release/npm/index.js` | Generate the wrapper from `chemlarp_c.h` instead of listing calls by hand. |
| Homebrew HEAD builds the main branch, not the commit under test, so a pull request's formula change is tested against main | `homebrew` job | Point `head` at the checkout (`brew install --HEAD` with a local `url`) in CI. |

## Native apps (`release/macos`, `release/windows`, `release/android`, scripts/build_macos_app.sh, scripts/build_android.sh)

Done 2026-11-14: macOS .dmg, Windows installer, Android .apk, each the offline web app in the platform's web
view, each run in CI with an engine call through its window.

| Gap | Where | To resume |
|---|---|---|
| No code signing or notarization: macOS ad hoc only (Gatekeeper asks), Windows unsigned (SmartScreen), Android signed with a throwaway key per build (a later APK will not install over an earlier one) | `scripts/build_macos_app.sh`, `installer-windows`, `scripts/build_android.sh` | Needs the owner's Apple Developer ID (codesign + notarytool), a Windows code-signing certificate (signtool), and an Android release keystore (CHEMLARP_KEYSTORE as a CI secret). |
| Android: exports (SVG, PNG, MOL) do not save: the WebView ignores the page's blob downloads | `release/android/.../MainActivity.java` | Add a JavaScript interface that receives the file as base64 and writes it through the Storage Access Framework (ACTION_CREATE_DOCUMENT); the web app would call it when present. |
| Android: Google Play needs an app bundle (.aab) signed with the owner's key; only an APK is built | `release/android` | `gradle bundleRelease` with the release keystore. |
| Windows: the WebView2 Runtime is assumed (Windows 11, current Windows 10); the installer does not bundle the bootstrapper | `release/windows/chemlarp.iss` | Bundle MicrosoftEdgeWebview2Setup.exe and run it from [Run] when the runtime is missing. |
| Windows: x64 only (no Arm64 build) | `installer-windows` | cl.exe /arm64 cross-compiler and the arm64 WebView2 loader. |
| The apps do not update themselves (the web app on GitHub Pages does) | all three | Sparkle (macOS), winget upgrade (Windows), Play Store (Android) once published. |
| No iOS app (needs an Apple developer account); no Linux desktop package (the installable web app covers Linux) | — | WKWebView shell like macOS for iOS; WebKitGTK or an AppImage for Linux if wanted. |
| Accessibility and keyboard behaviour inside the native windows were not reviewed separately from the web app | all three | Check VoiceOver/Narrator/TalkBack on each shell. |

## Virtual lab and safety rules (`src/lab.cpp`, `include/chemlarp/lab.hpp`, `src/ghs_data.inc`, web Virtual Lab)

Done 2026-11-14: mixing, titrations, calorimetry, precipitates and gases on the aqueous model; computed hazards; NRC
incompatibility rules; PubChem GHS for 207 chemicals. Validated against PHREEQC and CODATA.

| Gap | Where | To resume |
|---|---|---|
| No redox between valence states, so bleach + acid (chlorine), permanganate or dichromate oxidations, metals dissolving in acid are not simulated | `src/aqueous.cpp` | Couple valence states through electron activity (pe) as PHREEQC does; phreeqc.dat has the half-reactions. |
| Bleach + acid has no incompatibility rule: *Prudent Practices* (2011) does not state it in the chapters used | `kRules` in `src/lab.cpp` | Add with a second citable source (e.g. the CAMEO Chemicals reactivity profile in the PubChem record of sodium hypochlorite; needs download approval). |
| Reagents outside phreeqc.dat (hypochlorite, permanganate, peroxides, cyanide, silver, organics, metals) take part in the safety rules only | reagent catalogue | Extend the database (llnl.dat has more species) or give such reagents their own reaction rules. |
| Molality from molarity assumes the water mass equals the volume times 0.99705 kg/L: a few per cent high at 1 mol/L | `src/lab.cpp` | Use solution densities (e.g. from Pitzer apparent molar volumes) when concentrated. |
| Heat: no excess (non-ideal) enthalpies or heats of dilution; heat capacity 4.184 J/(g K) and 1 g/mL for the whole volume; no heat loss | `lab_mix` | Apparent molar heat capacities and enthalpies of dilution for concentrated solutions. |
| The GHS data are the pure substances' (hydrochloric acid as hydrogen chloride gas, aqueous ammonia as ammonia gas); concentration limits for solutions are not applied | `ghs_hazard` | Apply GHS concentration limits (generic cut-offs or the ECHA specific limits) to the solution's concentration. |
| PubChem's name lookup resolved "white phosphorus" to the red phosphorus record; that entry is dropped as a duplicate | `scripts/gen_ghs_data.py` | Look such compounds up by CID instead of name. |
| No computed warning thresholds for temperature short of boiling (burn risk), gas volumes in a closed vessel (pressure), or toxic-gas exposure limits | `computed_hazards` | Need sourced thresholds (e.g. ISO 13732-1 for hot surfaces, occupational exposure limits). |
| No kinetics: everything reaches equilibrium at once (no slow precipitation, no dissolution rates, no gas transfer rate) | `lab_mix` | The kinetics module could drive time-resolved lab experiments. |

## Validation runner (`scripts/validate.py`, `scripts/validate_rdkit.py`)

Done 2026-11-14: one runner for every validation suite, with a Markdown/JSON report that separates
independent references from in-repo expectations, lists what was skipped and why, and lists the roadmap
modules with no independent suite. `make validate DATA=... PY=...`.

| Gap | Where | To resume |
|---|---|---|
| Validations that docs/ARCHITECTURE.md reports but no script reproduces: InChI against PubChem (1483/1500), IUPAC names through OPSIN (14,188/14,193), quantum chemistry against PySCF, readers against RDKit on PubChem 2D/3D files, pure-component properties and EoS/activity models against thermo and DDBst | the `NOT_SCRIPTED` list in `scripts/validate.py` is the authoritative list | Each needs its data source re-fetched (PubChem sample with InChIs, OPSIN jar, pyscf, DDBst extracts) and the ad-hoc script rewritten against the current API. The old throwaway versions are gone, so these are the largest remaining pieces of this module. |
| CIP labels are compared as a multiset per molecule, not atom by atom, because the two programs number atoms differently; a swap between two centres of one molecule would pass | `scripts/validate_rdkit.py` part B | Map atoms through a canonical ranking (or match our canonical SMILES to RDKit's output atom order) and compare per atom. |
| SMARTS comparison uses 40 patterns on the first 400 molecules, not RDKit's full fragment set on the whole sample (the documented 187,375 pairs) | `scripts/validate_rdkit.py` part C | Read the patterns from RDKit's `FragmentDescriptors.csv` as the earlier run did, and drop the 400-molecule cap. |
| Fingerprints are compared through similarity values (mean difference, Spearman), not bit by bit, because ChemLarp uses its own 32-bit hash | `scripts/validate_rdkit.py` part D | Compare the sets of environment identifiers instead of the hashed bits. |
| The runner does not run in CI: the reference data is not redistributable and the Python packages are not installed there | `.github/workflows/ci.yml` | A job could install thermo/rdkit/biopython and run the suites that need no downloads (md, rdkit-structure with a committed SMILES sample, mmff-rdkit). |
| Tolerances live inside each script, so the report shows pass/fail but not how close each number was to its limit | the suite scripts | Have each script print a machine-readable line (name, value, tolerance) and collect those in the report. |
| No history: each run overwrites the report, so a slow drift in an agreement rate is invisible | `scripts/validate.py` | Append a row per run to a CSV beside the report. |
| Docking has no suite because the module is not built | — | Add when it lands (redocking RMSD on a public complex set). |

## Licensing and bundled data (`LICENSE`, `NOTICE`, `LICENSING.md`, `src/*_data.inc`)

Done 2026-11-17: relicensed to AGPL-3.0-only with commercial licences offered separately; NOTICE, LICENSING.md
(name and logo terms, earlier Apache-2.0 releases), CONTRIBUTING.md with a sign-off that keeps relicensing
possible, SPDX headers in every source file. Full audit in docs/AUDITS.md.

| Gap | Where | To resume |
|---|---|---|
| **`src/transport_data.inc` stores DIPPR correlation coefficients (as printed in Perry's Handbook) and VDI Heat Atlas PPDS coefficients.** Both are commercially licensed databases. The coefficients are fitted constants, which US law treats as facts, and they reach ChemLarp through the MIT-licensed `chemicals` package, but the EU database right protects substantial extraction regardless, and selling a commercial licence to software containing them is more exposed than giving it away | `src/transport_data.inc`, `scripts/gen_transport_data.py` | Before commercial licences are sold, either have a lawyer look at it, or replace those columns with openly licensed correlations (NIST/TDE, or fits of published measurements made here), or move the file out of the repository and load it from a user-supplied data directory. The Yaws Collection is already safe: it only votes between candidates and none of its values are stored. |
| Commits and tags are not signed, so GitHub does not show "Verified" and a forged commit in someone's fork could be passed off as the author's | the repository | Add an SSH or GPG signing key to the owner's GitHub account, then `git config commit.gpgsign true` and `git config tag.gpgsign true`, and sign future tags with `git tag -s`. |
| Releases are not signed or notarised: macOS warns on first launch, Windows SmartScreen may warn, and a download cannot be proved to come from the author beyond `SHA256SUMS` published on the same page | `.github/workflows/packages.yml` | Sigstore/cosign signing of the release assets is free and needs no certificate; an Apple Developer ID (paid) for notarisation and an Authenticode certificate (paid) remove the warnings. |
| Builds are not reproducible, so a release cannot be rebuilt bit-for-bit from its tag to prove it matches the source | the build | Set `SOURCE_DATE_EPOCH`, strip build paths (`-ffile-prefix-map`), and pin the toolchain versions in the workflow. |

## Roadmap modules not started (data/roadmap.json)

| Module | Size | Notes for starting |
|---|---|---|
| Docking | M | Needs protein module (done) + MMFF/UFF ligand; Vina-like empirical score, grid maps, Monte Carlo + BFGS. Validate redocking RMSD on PDBbind core set or CASF (licence check before download). |
