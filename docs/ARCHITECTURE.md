# ChemLarp: architecture

## Provenance: clean-room, literature-based

This project was scoped against ChemDraw, Labster, Schrödinger (Glide/Jaguar), and Aspen Plus/HYSYS. None of those products were decompiled or reverse-engineered.

- No installed copies were available to analyze. Labster is hosted SaaS, and Aspen is Windows-only.
- Every capability listed is already documented in open literature, so nothing needs to be recovered from a binary.
- Decompiling commercial chemistry suites to build a competitor would break their EULAs. It would also put the open-source license at risk.

Each algorithm cites its public source in its header. Use the same rule for every new module: **cite the paper or standard, never the product.**

| Capability | Public basis |
|---|---|
| SMILES parsing, aromaticity | OpenSMILES spec (opensmiles.org); Weininger, *JCICS* 28 (1988) 31 |
| MOL/SDF, PDB, CML | BIOVIA CTfile Formats (2020); wwPDB Atomic Coordinate Entry Format 3.3; Murray-Rust & Rzepa, *JCICS* 39 (1999) 928 (CML); Cordero et al., *Dalton Trans.* (2008) 2832 (covalent radii) |
| SMARTS, fingerprints | Daylight Theory Manual ch. 4 (SMARTS); Rogers & Hahn, *JCIM* 50 (2010) 742 (ECFP); Vismara, *Electron. J. Combin.* 4 (1997) R9 (relevant cycles); Willett, Barnard & Downs, *JCICS* 38 (1998) 983 (similarity coefficients) |
| Quantum chemistry | Helgaker, Jorgensen & Olsen, *Molecular Electronic-Structure Theory* (2000) ch. 9 (McMurchie-Davidson integrals); Pulay, *CPL* 73 (1980) 393 (DIIS); Becke, *JCP* 88 (1988) 2547 (grid partition); Treutler & Ahlrichs, *JCP* 102 (1995) 346 (radial grid); Vosko, Wilk & Nusair, *CJP* 58 (1980) 1200; Perdew, Burke & Ernzerhof, *PRL* 77 (1996) 3865; Perdew & Wang, *PRB* 45 (1992) 13244; Becke, *PRA* 38 (1988) 3098; Lee, Yang & Parr, *PRB* 37 (1988) 785 with Miehlich et al., *CPL* 157 (1989) 200; Stephens et al., *JPC* 98 (1994) 11623; basis sets from the original papers listed in src/qc_basis_data.inc, digitised by the Basis Set Exchange (Pritchard et al., *JCIM* 59 (2019) 4814) |
| IUPAC names | Favre & Powell, *Nomenclature of Organic Chemistry: IUPAC Recommendations and Preferred Names 2013* (RSC 2014): P-14 (locants, multiplying prefixes, alphanumerical order), P-16.5 (enclosing marks), P-21/P-22/P-25 (parent hydrides, retained ring names, fused-ring numbering), P-23.2 (von Baeyer), P-31.1.4 (numbering, indicated and hydro prefixes), P-41 (seniority of classes), P-44/P-45 (parent selection), P-6 (characteristic groups), P-93 (stereodescriptors). Validation by parsing the names with OPSIN (Lowe et al., *JCIM* 51 (2011) 739) |
| InChI | Heller et al., *J. Cheminform.* 7 (2015) 23; InChI Technical Manual. Written from the published description; the reference software is used only to produce test strings |
| R/S, E/Z | IUPAC 2013 recommendations P-92 (CIP rules); Hanson et al., *JCIM* 58 (2018) 1755 (CIP algorithm validation) |
| 2D depiction | Helson, *Rev. Comput. Chem.* 13 (1999) 313; Clark, Labute & Santavy, *JCIM* 46 (2006) 1107; IUPAC graphical representation of stereochemical configuration, Brecher et al., *Pure Appl. Chem.* 78 (2006) 1897; ACS 1996 drawing settings |
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
- **C ABI** (`chemlarp_c.h`) exposes raw column pointers. Python wraps them with `numpy.ctypeslib.as_array`, which is verified zero-copy: a write from NumPy is visible in C++. Node can do the same with `ArrayBuffer` via N-API.

## Identity: canonical SMILES and CIP

`canonical.cpp` produces one string per structure, stereo included. Every species carries `canonical_smiles`, set at parse time, and the thermo database is keyed by it. The pipeline:

1. **Normalised view.** A plain `[H]` with a single heavy neighbour is folded into that atom's H count. Aromaticity is re-perceived from the Kekulé structure, so input notation never leaks into identity.
2. **Ranking.** Initial atom invariants are (degree, Z, isotope, charge, H count, ring, aromatic, CIP label, atom class). They are refined iteratively by sorted (neighbour rank, bond code) multisets until stable. While ties remain, the lowest tied class is split and refinement runs again (Weininger 1989). With stereo present, each stereo element's parity relative to the current ranks is refined in once its deciding neighbours have distinct ranks; while any element is still undecided, every member of the tied class is tried and the split giving the largest rank-ordered parity vector wins. Without this, mirror-image spellings of an achiral molecule with ring cis/trans stereo (C[C@H]1CC[C@@H](O)CC1) gave two strings.
3. **Writing.** DFS starts from the lowest rank, neighbours are visited in rank order, and ring closures take the lowest free digit. Tetrahedral parity is recomputed as the permutation parity between the stored reference order and the output neighbour order. `/` `\` marks come from a small constraint search, so conjugated systems (shared single bonds) stay consistent; the choice is anchored on output order only.
4. **Stereo is part of identity.** CIP labels are rank invariants, so meso and chiral diastereomers, and the members of meso pairs, are distinguished deterministically. Stereo marks on centres that are neither CIP-stereogenic nor ring para-stereocentres are dropped, with a warning.

How it is verified: the randomized writer (`SmilesWriteOptions{canonical=false}`) emits arbitrary atom orders in aromatic or Kekulé form. The test suite requires `canonical(random(x)) == canonical(x)` for 33 molecules × 40 seeds; an ad-hoc run of 14 harder molecules × 1000 seeds also passed. Spelling invariance across 5 random spellings each (2026-11-12): 0 of 617 hard symmetric stereo cases, 9 of 6322 enumerated stereoisomers and 0 of 1500 PubChem compounds give more than one string (were 109, 181 and 1). Ring para-stereocentres without a CIP descriptor (outer carbons of 1,3,5-trisubstituted cyclohexanes) are kept; isomer counts for ten ring families equal RDKit's. E/Z next to aromatic ring bonds (quinone imines, dyes) is now written on the Kekulé-single ring bond instead of being dropped; that also removed 10 of 12 cases where different RDKit InChIKeys shared one string (the 2 left are an orthoester cage that InChI gives bridgehead stereo). The 9 left: E/Z on meso double bonds of chlorins, whose 18-pi macrocycle the ring-by-ring aromaticity perception does not see, and one phosphazane cage.

Stereo storage (`TetrahedralStereo`, `DoubleBondStereo`) keeps references to explicit neighbours plus implicit-H and lone-pair markers. It never depends on write order, so the same structure can be re-emitted in any order.

## Platforms

Windows, macOS and Linux are first-class targets, and the code base follows a few rules to stay that way:

- **ISO C++20 and the standard library only.** No compiler builtins (`std::countr_zero` instead of `__builtin_ctzll`), no POSIX calls, no `<windows.h>` outside `src/cli/main.cpp`, where it is used only to switch the console to UTF-8.
- **One exported C ABI.** `CHEMLARP_API` expands to `__declspec(dllexport/dllimport)` on Windows and to default visibility elsewhere. `CHEMLARP_BUILD_DLL` and `CHEMLARP_STATIC` select the mode, and CMake sets them per target.
- **Builds.** CMake is canonical (MSVC: `/utf-8 /permissive- /W4`). The Makefile covers macOS, Linux and MSYS2 MinGW. MinGW DLLs are linked with static libstdc++, so the DLL is self-contained for Python's `ctypes`.
- **Files and text.** Mechanism files are read in text mode and accept CRLF. `.gitattributes` pins line endings for scripts. Tests resolve `data/` through `CHEMLARP_TEST_DATA_DIR`, so CTest can run from any build directory.

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
| Canonical SMILES (invariant refinement + tie-breaking, stereo-aware) | **done**; 14,000-case randomized round-trip fuzz with zero mismatches; mirror spellings of achiral ring-stereo molecules unified; ring para-stereocentres kept (2026-11-12). Open: E/Z in chlorin-type macrocycles (needs macrocycle aromaticity) |
| Randomized SMILES writer (aromatic or Kekulé) | **done** (used as the test oracle) |
| CIP R/S, r/s and E/Z, rules 1a-5 (hierarchical digraph, mancude duplicates, auxiliary descriptors, like/unlike) | **done**; agrees with RDKit's CIP labeller on 816/816 PubChem stereo structures and 6265/6312 enumerated stereoisomers. The differences: ring-closure digits on three-coordinate N/P centres (the two SMILES readers disagree), `[P@H]` phosphines and quinoid `c(=C)` input (not parsed). Rule-5 double bonds print E/Z, not seqCis/seqTrans; no axial or planar units |
| Stereo error vector: conflicting `/` `\` marks, small-ring E/Z, non-tetrahedral `@` | **done** |
| Thermo database keyed by canonical SMILES (formula fallback warns `ThermoFormulaMatch`) | **done** |
| Composition, Hill formula, molar mass | **done** |
| Reaction parser, element/charge balance error vector, ΔrH°/ΔrS°/ΔrG°/K | **done**, tested |
| NASA-7 thermo (10 GRI-Mech species) | **done**, verified against JANAF to <0.6 kJ/mol |
| Joback ΔfH°/ΔfG° from graph | **done**; ethanol, propane, acetone, benzene, and acetic acid fall within the tolerances shown in tests |
| Ideal-gas Gibbs minimisation (element-potential Newton) | **done**, tested (K consistency, conservation, Le Chatelier, rank deficiency) |
| C ABI + Python ctypes binding; `.clx` v0.1 container | **done** |
| 2D layout (ring templates by arc closure, 120° chains, E/Z by mirroring, overlap relief by subtree reflection) and SVG drawing with IUPAC-style wedges, ring-fusion H, CIP labels; 2D molfile export | **done**; every stereocentre's wedge is checked by reading the chirality back from the drawing. Bridged polycycles (morphine, adamantane) are legible but not template-quality |
| Standard InChI and InChIKey (normalisation, mobile H, canonical numbering, /c /h /q /p /b /t /m /s /i layers) | **done** for organic structures; 1483/1500 (98.9%) identical to PubChem on a random sample never used for tuning. Not covered: isotopic stereo, allene stereo, /f and /r layers, polymers, several charge-migration and salt rules (diazonium, N-oxides, metal complexes of macrocycles) |
| Allene / atropisomer / square-planar stereo (@AL, @SP, @TB, @OH) | planned |
| Explicit hydrogens with stereo remapping | **done** |
| Distance-geometry 3D embedding (bounds + triangle smoothing, metric matrix in 4D, chiral/planar refinement) | **done**; stereo read back from 3D and enforced, best-of-N conformers; inverted structures (an sp3 carbon past 150 degrees, from geometrically impossible stereo such as a wrong norbornane bridgehead pair) are rejected and reported instead of returned |
| UFF force field (bond, angle, torsion, inversion, vdW; analytic gradient) + L-BFGS | **done**; gradients exact vs finite differences, geometries vs experiment in tests |
| SDF (V2000) and XYZ output | **done**; writers mark open double bonds "either" and give the valence field where MDL implicit-H rules would differ |
| Readers: MOL V2000/V3000, SDF (data items), PDB (CONECT, LINK, amino-acid templates, first model/altloc), CML | **done**; stereo from wedges (2D) or coordinates (3D). Same structure as RDKit reading the same files: 999/999 PubChem 2D, 865/866 PubChem 3D and 1465/1468 of our own 3D files by canonical SMILES (the rest: cage/imine stereo outside the CIP model, or a distorted 3D input); crambin and ubiquitin identical InChIKey. Not read: query features, R-groups, S-group expansion, enhanced stereo, nucleic-acid bond orders, XYZ |
| Structure editor (web): bonds, wedge/hash, chains, rings 3-8 and benzene (fused or substituent), element and charge tools and hover hotkeys, select/move/merge, undo/redo, zoom/pan, clean-up via the 2D layout | **done**; the drawing is exchanged as a V2000 molfile, so validation, stereo from wedges and the canonical SMILES come from the engine's reader; scripted browser checks (zig-zag, phenyl, naphthalene, chain, eraser, valence errors, charges, (S)-alanine round trip) |
| Quantum chemistry: RHF/UHF and RKS/UKS (SVWN5, PBE, B3LYP) with STO-3G, 3-21G, 6-31G, 6-31G* (H-Ar); McMurchie-Davidson integrals, SAD guess, DIIS, level shift; energy, orbitals, Mulliken charges, dipole, <S^2> | **done**; against PySCF on the same geometries: HF energies within 2e-7 hartree for 19 molecules x 4 basis sets (closed shells, ions, radicals, triplet O2), DFT within 2e-6 for all three functionals at 6-31G* on the fine grid (benzoic acid B3LYP: 1.6e-8), dipoles and Mulliken charges within 1e-5. Default grid error about 4e-6 hartree. Not converged: UHF/UKS of NO in larger bases (several nearby solutions; PySCF also struggles). No geometry optimisation or gradients yet; elements beyond Ar not covered |
| IUPAC names (structure to name): chains, monocyclic carbocycles, benzene, about 140 retained, Hantzsch-Widman and fused ring parents (pyridine ... imidazo[1,2-a]pyridine, thieno[2,3-d]pyrimidine, 1,4-benzodiazepine, cyclopenta[a]phenanthrene, adamantane) with indicated hydrogen and hydro prefixes, ortho- and peri-fused systems of any number of rings by fusion nomenclature (base component by FR-3.1, attached components to second order, numbering from the FR-5.2 orientation rules on an eight-direction grid: thieno[2,3-b]pyridine, 9H-pyrido[2,3-b]indole, 4H-[1,2,4]triazolo[4,3-a][1,4]benzodiazepine, cyclopenta[4,5]thieno[2,3-d]pyrimidine, chrysene), Hantzsch-Widman names for any 3-10 membered heteromonocycle, von Baeyer names for bridged systems of any number of rings (cubane as pentacyclo[4.2.0.0^{2,5}.0^{3,8}.0^{4,7}]octane, porphyrins, morphinans; compound locants such as 7(18) for double bonds) and for fused systems without a fusion name here, monospiro systems, large heteromonocycles by replacement names; acids (incl. sulfinic), anhydrides, esters (carboxylic, thio, sulfonic, sulfuric, phosphoric, phosphonic, di- and triphosphoric monoesters), acid halides, silanes and silyl ethers, acyl ureas and thioureas, N,N'-diacylhydrazines, azo compounds (diazenes, diazenyl), exocyclic imines on rings and N-acyl imines, semicarbazones and thiosemicarbazones, phosphanes, phosphine oxides (λ5-phosphanones), phosphinic acids, phosphorothioate and phosphorodithioate triesters, disulfides, sulfonylureas and N-acyl sulfonamides, addition compounds of neutral components (A—B (1/2)), N-oxides (oxido...ium), N-amino heterocycles, amides, imides, thioamides, hydrazides, N-hydroxy and N-alkoxy amides, amidoximes, hydrazines and hydrazones, hydroxylamines, ring sulfoxides and sulfones (lambda convention), nitriles, aldehydes, ketones and thiones (incl. lactams/lactones), alcohols, thiols, amines, imines/oximes, sulfonic acids, sulfonamides, sulfones/sulfoxides, phosphate monoesters, phosphonic acids, carbonates, carbamates, ureas, thioureas, amidines, guanidines, nitro, halogens; anions (-oate, -sulfonate, -olate), cations (-aminium, ring -ium, including -ium beside another suffix), zwitterions, salts, hydrochlorides and hydrates; CIP descriptors | **done**; every name parsed back with OPSIN gives the input structure (same InChIKey, stereo and charge included): 13,795 of 13,800 named structures from 15,000 PubChem compounds; the two others are names OPSIN cannot read (a 53-fold multiplying prefix, and an indicated-hydrogen tautomer whose double bonds OPSIN fails to place). Coverage: 86.3% of CIDs 1-2000 (development set), 89.9% of CIDs 2001-4000, 93.2% of 3,000 random CIDs and 91.0% of a further 8,000 random CIDs (the last three never used for tuning). Anything else gets no name and a NamingUnsupported diagnostic saying why; the largest gaps are uncommon ions and metal complexes, isothioureas and formamidines, spiro systems with fused components, di- and triphosphate diesters (NAD, CoA), mixed and cyclic phosphates, phosphoramides, isothioureas, semicarbazides, boron and metals, substituted or cyclic amidines, hydrazines and azo compounds. Not always the preferred IUPAC name: bridged fused systems (morphine) and fused systems with several peri-fused components get von Baeyer names instead of bridged-fusion or fusion names, ring assemblies are named substitutively (phenylbenzene, not 1,1'-biphenyl) and saturation is written with hydro prefixes instead of added hydrogen (uracil: 1,2,3,4-tetrahydropyrimidine-2,4-dione, PIN pyrimidine-2,4(1H,3H)-dione). Agreement with PubChem's own (non-PIN) names is 17-48% character for character |
| SMARTS substructure search (full Daylight atom/bond logic, recursive $(), ring/valence primitives, optional chirality) | **done**; unique match counts identical to RDKit on 187,338 of 187,375 pattern-molecule pairs (99.98%: 125 patterns incl. RDKit's fragment set × 1,499 PubChem compounds); the rest are aromaticity-model differences in 5 molecules. 1,500 fuzzed patterns: no crash or hang under UBSan |
| Relevant cycles (unique ring set) | **done**; used by SMARTS R/r so ring counts do not depend on atom order (bicyclo[2.2.2]octane: 3 rings) |
| Morgan / ECFP fingerprints, Tanimoto, Dice | **done**; total environment counts identical to RDKit's Morgan for 1,499/1,499 compounds at radius 1-3 and distinct counts for 99.7% (rest: aromaticity); own 32-bit hash, so bit positions differ from RDKit |
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
