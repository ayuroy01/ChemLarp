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
| Force fields | MMFF94: Halgren, *J. Comput. Chem.* 17 (1996) 490, 520, 553, 587, 616 and MMFF94s 20 (1999) 720, 730; parameter tables (Copyright Merck & Co., 1994–1996, published as the papers' supplementary material) converted by `scripts/gen_mmff_data.py` from the verbatim copies in RDKit's `Params.cpp` (BSD; only the data are used, the typing and energy code is written from the papers and checked against the official validation suite); UFF: Rappé et al., *JACS* 114 (1992) 10024; OPLS-AA: Jorgensen et al., *JACS* 118 (1996) 11225 |
| MD integrators | Verlet (1967); Swope et al. velocity Verlet, *JCP* 76 (1982) 637; Bussi, Donadio & Parrinello, *JCP* 126 (2007) 014101 (stochastic velocity rescaling); Leimkuhler & Matthews, *Appl. Math. Res. Express* 2013, 34 (BAOAB Langevin); Marsaglia & Tsang, *ACM TOMS* 26 (2000) 363 (gamma deviates) |
| Docking scoring | AutoDock Vina, Trott & Olson, *J. Comput. Chem.* 31 (2010) 455 (open, Apache-2) |
| FEP / MBAR | Zwanzig (1954); Shirts & Chodera, *JCP* 129 (2008) 124105 |
| Thermo data | GRI-Mech 3.0; NASA Glenn (McBride, Zehe & Gordon, NASA/TP-2002-211556); Burcat database |
| Gibbs minimisation | White, Johnson & Dantzig, *JCP* 28 (1958) 751; Reynolds STANJAN (1986); Smith & Missen, *Chemical Reaction Equilibrium Analysis* (1982) |
| Kinetics | modified Arrhenius / mass action; CHEMKIN-II format (Kee et al., Sandia SAND89-8009) |
| Process streams | Rachford & Rice (1952); liquid enthalpy referenced at the normal boiling point; validation against IAPWS-95 (Wagner & Pruß, *J. Phys. Chem. Ref. Data* 31 (2002) 387) |
| Transport properties | Perry's Chemical Engineers' Handbook, 8th ed. (2008), Tables 2-312 to 2-315 (DIPPR); Mulero & Cachadiña, *J. Phys. Chem. Ref. Data* 41 (2012) 043105; Jasper, *J. Phys. Chem. Ref. Data* 1 (1972) 841; validation against VDI Heat Atlas (PPDS) and REFPROP |
| Pure-component data | Critical constants: IUPAC critical-property series (Ambrose, Tsonopoulos, Young et al., *J. Chem. Eng. Data* 1995–2015), CRC Handbook, PSRK appendix (Horstmann et al., *Fluid Phase Equilib.* 227 (2005) 157); vapour-pressure equations of Poling et al. and McGarry (*Ind. Eng. Chem. Process Des. Dev.* 22 (1983) 313); Watson (1943); Rackett (1970), Spencer & Danner (1972); validation against NIST-JANAF (1998), TRC and VDI Heat Atlas (PPDS) data |
| Activity models / EoS | NRTL: Renon & Prausnitz, *AIChE J.* 14 (1968) 135; UNIFAC: Fredenslund, Jones & Prausnitz, *AIChE J.* 21 (1975) 1086, subgroup and interaction tables published by DDBST (Hansen et al., *Ind. Eng. Chem. Res.* 30 (1991) 2352, and revisions); Antoine coefficients: Poling, Prausnitz & O'Connell, *The Properties of Gases and Liquids*, 5th ed. (2001); Peng–Robinson (1976) |
| Flash / distillation | Rachford & Rice (1952); Michelsen, *Fluid Phase Equil.* 9 (1982) (stability); Naphtali–Sandholm (1971) |
| Flowsheet convergence | Westerberg et al., *Process Flowsheeting* (1979); Biegler, Grossmann & Westerberg (1997) |
| Aqueous chemistry | Database: PHREEQC's phreeqc.dat (U.S. Geological Survey, public domain; Parkhurst & Appelo, USGS Techniques and Methods 6-A43, 2013), converted by `scripts/gen_aqueous_data.py`; ion-association model with WATEQ Debye–Hückel (Truesdell & Jones 1974) and Davies activity coefficients; dielectric constant of water: Malmberg & Maryott, *J. Res. NBS* 56 (1956) 1; density: Kell, *J. Chem. Eng. Data* 20 (1975) 97. PHREEQC itself is used only as a validation reference |

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
| IUPAC names (structure to name): chains, monocyclic carbocycles, benzene, about 140 retained, Hantzsch-Widman and fused ring parents (pyridine ... imidazo[1,2-a]pyridine, thieno[2,3-d]pyrimidine, 1,4-benzodiazepine, cyclopenta[a]phenanthrene, adamantane) with indicated hydrogen and hydro prefixes, ortho- and peri-fused systems of any number of rings by fusion nomenclature (base component by FR-3.1, attached components to second order, numbering from the FR-5.2 orientation rules on an eight-direction grid: thieno[2,3-b]pyridine, 9H-pyrido[2,3-b]indole, 4H-[1,2,4]triazolo[4,3-a][1,4]benzodiazepine, cyclopenta[4,5]thieno[2,3-d]pyrimidine, chrysene), Hantzsch-Widman names for any 3-10 membered heteromonocycle (O, S, N, P), ring phosphorus with the λ5 convention in monocycles and fused systems (N,N-bis(2-chloroethyl)-2-oxo-1,3,2λ5-oxazaphosphinan-2-amine, 4H-2λ5-furo[3,2-d][1,3,2]dioxaphosphinine for cyclic nucleotides), von Baeyer names for bridged systems of any number of rings (cubane as pentacyclo[4.2.0.0^{2,5}.0^{3,8}.0^{4,7}]octane, porphyrins, morphinans; compound locants such as 7(18) for double bonds) and for fused systems without a fusion name here, monospiro systems including fused components (1',2'-dihydrospiro[cyclohexane-1,3'-indole]-2'-one, 3',4'-dihydro-2'H-spiro[imidazolidine-4,1'-naphthalene]-2,5-dione), large heteromonocycles by replacement names, ring assemblies of two or three identical benzene or aromatic heteromonocycle units ([1,1'-biphenyl]-4-ol, 2-([1,1'-biphenyl]-4-yl)propanoic acid, [1,1':4',1''-terphenyl]-4-ol, 2,3'-bipyridine, 2,2':6',2''-terpyridine, [2,2'-bipyridin]-3-ol, 2,2'-bi-1,3-thiazole; junction locants lowest first), multiplying prefixes that cannot fuse with a numeral (tris(decyl), not tridecyl); acids (incl. sulfinic), anhydrides, esters (carboxylic, thio, sulfonic, sulfuric, phosphoric, phosphonic, di- and triphosphoric monoesters, and their diesters such as nucleotide sugars as phosphoryl phosphates: [hydroxy(methoxy)phosphoryl] methyl hydrogen phosphate), acid halides, silanes and silyl ethers, acyl ureas and thioureas, N,N'-diacylhydrazines, azo compounds (diazenes, diazenyl), exocyclic imines on rings and N-acyl imines, diazo compounds (2-amino-6-diazo-5-oxohexanoic acid), amidinium cations (-carboximidamidium), diazeniumdiolates as N-hydroxy-N-nitroso hydrazines, mancude parents with up to four indicated hydrogens (2H,6H,10H-dipyrano[2,3-f:2',3'-h][1]benzopyran), semicarbazides, semicarbazones and their thio analogues as hydrazinecarboxamides (2-(propan-2-ylidene)hydrazine-1-carboxamide, 2-acetylhydrazine-1-carbothioamide), phosphanes, phosphine oxides (λ5-phosphanones), phosphinic acids, phosphorothioate and phosphorodithioate triesters and acid esters (O,O-diethyl O-hydrogen phosphorothioate; identical organyls multiplied across locant letters: O,S-dimethyl), disulfides, sulfonylureas and N-acyl sulfonamides, addition compounds of neutral components (A—B (1/2)), N-oxides (oxido...ium), N-amino heterocycles and amides on ring nitrogen (N-(piperidin-1-yl)benzamide), N-hydroxy and N-alkoxy ring nitrogen (1-hydroxypiperidine), amides, imides, thioamides, hydrazides, N-hydroxy and N-alkoxy amides, amidoximes, hydrazines and hydrazones (incl. 2-acylhydrazinylidene prefixes), hydroxylamines (O-(2-methylpropyl)hydroxylamine), several N-substituted amide or amine suffixes with N1/N4 locants (N1-methyl-N4-phenylbenzene-1,4-dicarboxamide), ring sulfoxides and sulfones (lambda convention), nitriles, aldehydes, ketones and thiones (incl. lactams/lactones), alcohols, thiols, peroxols and peroxides (propane-2-peroxol, (methylperoxy)methane), amines, imines/oximes, sulfonic acids, sulfonamides (incl. on ring nitrogen: piperidine-1-sulfonamide, N-chloro), sulfonohydrazides, sulfuric diamides, sulfamate esters, thionitrites (nitrososulfanyl), sulfones/sulfoxides, phosphate monoesters, phosphonic acids, phosphonic and H-phosphinic acids (methylphosphinic acid, hydroxyphosphinoyl), boronic acids (phenylboronic acid, borono), carbamoyl phosphate, phosphate-sulfate anhydrides ([hydroxy(methoxy)phosphoryl] hydrogen sulfate), di- and triphosphate anions (methyl diphosphate, trihydrogen diphosphate), phosphoric acid, phosphate and phosphonate anions (hexadecyl 2-(trimethylazaniumyl)ethyl phosphate), acyl phosphates (acetyl dihydrogen phosphate), thiosulfate S-esters (S-ethyl hydrogen sulfurothioate, sulfosulfanyl), hydrodisulfides (disulfanyl), small inorganic molecules and ions by exact match against a table of accepted names (water, ammonia, sulfurous acid, thiosulfate, hydroxylammonium; formonitrile keeps its organic name), sulfuric and sulfamic acids, sulfates and sulfamates (sodium ethyl sulfate, N-methylsulfamic acid), carbonates, carbamates and carbamic acids, hydrogen carbonates, dithiocarbamates, xanthates and thiocarbamates (O-/S- locants), ureas, thioureas, amidines (incl. ring-N carboximidamides), guanidines and acylguanidines (N-(diaminomethylidene)acetamide), isothioureas as carbamimidothioates (2-aminoethyl carbamimidothioate), nitro, nitroso (C- and N-nitroso), azido, cyanato, thiocyanato, isocyanato, isothiocyanato, halogens; anions (-oate, -sulfonate, -olate, ring -ide: 1H-indol-1-ide, piperidin-1-ide, 9H-fluoren-9-ide), cations (-aminium, ring -ium, including -ium beside another suffix; sulfanium, phosphanium, -diazonium), borates, sulfonyl azanides, nitric acid and nitrate, zwitterions, salts (multiplied ions as bis(...) where di/tri would be ambiguous or malformed: bis(azanium) sulfate, calcium bis(N,N-dimethylcarbamodithioate); calcium diacetate), hydrochlorides and hydrates; CIP descriptors | **done**; every name parsed back with OPSIN gives the input structure (same InChIKey, stereo and charge included): 14,188 of 14,193 named structures from 14,996 PubChem compounds; the five others are names OPSIN cannot read (a 53-fold multiplying prefix, an indicated-hydrogen tautomer whose double bonds OPSIN fails to place, and three addition compounds with ethanedioic acid). Coverage: 95.9% of CIDs 1-2000 (development set), 95.3% of CIDs 2001-4000, 95.4% of 3,000 random CIDs and 93.9% of a further 8,000 random CIDs. The OPSIN round trip shows that each name denotes the input structure; it does not show that the name follows the IUPAC rules or is the preferred IUPAC name. Rule conformance (seniority, numbering, name construction) is checked only by the hand-reviewed cases in tests/test_iupac.cpp and by the refusals described below. The last three sets were never used to tune rules, but their failure lists were consulted to choose which groups to implement in several batches up to November 2026, so they are not fully independent; failure classification now uses the development set only. Anything else gets no name and a NamingUnsupported diagnostic saying why, including a principal characteristic group that could only be cited as a prefix; the largest gaps are uncommon ions and metal complexes, formamidines, spiro systems with several spiro atoms, di- and triphosphates carrying anions, mixed phosphates, ring phosphorus with a ring double bond (phosphinines, phospholes), phosphoramides, boron and metals, substituted or cyclic amidines, hydrazines and azo compounds. Not always the preferred IUPAC name: bridged fused systems (morphine) and fused systems with several peri-fused components get von Baeyer names instead of bridged-fusion or fusion names, ring assemblies of four or more rings, of different rings, of fused units or of rings needing indicated hydrogen are named substitutively (1-phenyl-4-(4-phenylphenyl)benzene for quaterphenyl; 2-(thiophen-2-yl)pyridine), polyol esters get multiplicative names when every ester group is on one parent and all acid parts are identical (ethane-1,2-diyl diacetate, propane-1,2,3-triyl tris(decanoate), cyclohexane-1,2,3,4,5,6-hexayl hexakis(dihydrogen phosphate), 1,4-phenylene bis(dihydrogen phosphate)); esters on different parents (ring and CH2, as in sugar bisphosphates) or of different acids cite the further ester groups as prefixes and saturation is written with hydro prefixes instead of added hydrogen (uracil: 1,2,3,4-tetrahydropyrimidine-2,4-dione, PIN pyrimidine-2,4(1H,3H)-dione). Agreement with PubChem's own (non-PIN) names is 10-36% character for character; it measures closeness to PubChem's naming style, not correctness, and drops wherever the namer moves closer to the preferred name (PubChem writes 4-phenylphenyl, not [1,1'-biphenyl]-4-yl) |
| SMARTS substructure search (full Daylight atom/bond logic, recursive $(), ring/valence primitives, optional chirality) | **done**; unique match counts identical to RDKit on 187,338 of 187,375 pattern-molecule pairs (99.98%: 125 patterns incl. RDKit's fragment set × 1,499 PubChem compounds); the rest are aromaticity-model differences in 5 molecules. 1,500 fuzzed patterns: no crash or hang under UBSan |
| Relevant cycles (unique ring set) | **done**; used by SMARTS R/r so ring counts do not depend on atom order (bicyclo[2.2.2]octane: 3 rings) |
| Morgan / ECFP fingerprints, Tanimoto, Dice | **done**; total environment counts identical to RDKit's Morgan for 1,499/1,499 compounds at radius 1-3 and distinct counts for 99.7% (rest: aromaticity); own 32-bit hash, so bit positions differ from RDKit |
| Multiphase / condensed-phase equilibrium, activity models | **in progress**: cubic-EoS vapour-liquid flash and NRTL / UNIFAC activity-coefficient VLE done (below); liquid-liquid and solid-liquid equilibrium planned |
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
| MMFF94 and MMFF94s (`include/chemlarp/mmff.hpp`): the 99 MMFF atom types from structure (hypervalent S=O / P=O or dative S+–O− drawings, charges the drawing leaves implicit, MMFF's own aromaticity for 5- and 6-rings including fused, cationic and anionic rings), formal charges shared over resonance groups, bond-charge-increment partial charges, all seven energy terms with Halgren's step-down parameter search and the empirical rules for missing bond, angle, stretch-bend and torsion parameters, analytic gradient (forward-mode automatic differentiation), minimization, MMFF dynamics; `chemlarp mmff` (SMILES, or `--file` for the terms of every record of an SDF), `cl_mmff_json`, Python `mmff()`, force-field choice in `chemlarp md`, `cl_md_json` and `molecular_dynamics()`; web 3D Model view (UFF / MMFF94 / MMFF94s, energy terms, types and charges) and Molecular Dynamics view | **done** (`src/mmff.cpp`); validation in "MMFF94: validation" below. Not covered: OPLS-AA, MMFF's distance-dependent dielectric option, periodic or solvated systems, non-bonded cutoffs (all pairs are summed) |
| Molecular dynamics (`include/chemlarp/md.hpp`): velocity Verlet (NVE), Bussi stochastic velocity rescaling with its conserved effective energy, and BAOAB Langevin dynamics; Maxwell–Boltzmann start with the centre-of-mass motion and overall rotation removed; platform-independent random numbers; harmonic normal modes (finite-difference Hessian, translations and rotations projected out, linear molecules detected from the moments of inertia); vibrational spectrum from the velocity autocorrelation; UFF molecules from SMILES; periodic Lennard-Jones fluid with block-averaged energy and pressure; `chemlarp md`, `chemlarp lj`, `cl_md_json`, `cl_lj_json`, Python `molecular_dynamics()` and `lennard_jones()`, web Molecular Dynamics panel with an animated trajectory, energy, temperature and spectrum charts and XYZ export | **done** (`src/md.cpp`); validation in "Molecular dynamics: validation" below. Not covered: periodic boxes for molecules, solvent, constraints (SHAKE/RATTLE), barostats (NPT), electrostatics (UFF has none here), neighbour lists for large systems, multiple time steps |
| ETKDG-style torsion preferences; conformer ensembles with RMSD pruning | planned |
| Docking scoring, FEP/MBAR | planned |
| Peng–Robinson and SRK cubic EoS, mixture fugacities, 30-compound critical table, then the pure-component database (319 compounds with Tc, Pc and ω) | **done**; vapour pressures within 0.1–0.7% of NIST for non-polar compounds; with database constants, median 1.3% (Tr = 0.6) and 0.9% (Tr = 0.9) from McGarry's Wagner equations over about 160 compounds (below) |
| Michelsen stability test + Rachford–Rice two-phase flash; pure-component vapour pressure | **done** |
| NRTL and original UNIFAC activity coefficients (any number of components), UNIFAC subgroups from structure, Antoine vapour pressures (320 compounds), bubble and dew T and P by modified Raoult's law; `chemlarp vle`, `cl_vle_text_json`, Python `vle()`; web Phase Behaviour panel with T–x–y / P–x–y diagrams and azeotrope marking for binaries | **done** (`src/activity.cpp`); validation in "Activity models: validation" below. Not covered: NRTL parameter database (parameters are supplied by the caller), modified UNIFAC (Dortmund), vapour non-ideality, liquid-liquid splits (dew points that need one fail with VleNotConverged) |
| Pure-component property database (365 compounds keyed by structure: Tc, Pc, Vc, ω, Tm, Tb, ΔHvap(Tb), ideal-gas Cp polynomial, liquid Cp(298), Rackett Z_RA, ideal-gas ΔfH°, ΔfG°, S°) with ideal-gas Cp, H, S, Watson ΔHvap(T), saturated liquid volume (VDI Heat Atlas PPDS density for 168 compounds, else Rackett); transport properties from published correlations with their fitted ranges: liquid and low-pressure vapour viscosity and thermal conductivity (Perry 8th ed., DIPPR equations, about 180 compounds each) surface tension (Mulero–Cachadiña, else Jasper; 225 compounds) and liquid heat capacity (Perry Table 2-153, DIPPR 100; 167 compounds); `chemlarp props`, `cl_pure_component_json`, Python `pure_component()` | **done** (`src/properties.cpp`); validation in "Pure-component properties: validation" below. Scope: published data and the correlations above; a missing value is reported (NoPropertyData), never estimated, and use outside a fitted range is a warning. Not covered: estimation methods for compounds without data (Chung, Lucas, Brock–Bird), pressure corrections to vapour transport properties, solid properties, web UI |
| Material streams and unit operations: TP and PH flash (modified Raoult, ideal or UNIFAC liquid, non-condensable components without a vapour-pressure equation; ideal-gas vapour; liquid enthalpy from ΔHvap at Tb, liquid Cp and v_L(P − P_sat)), mixer, splitter, heater/cooler (to a temperature or by duty), valve, flash drum (isothermal or adiabatic), counter-current heat exchanger (zone profile, temperature-cross check, UA), pump, compressor/expander (isentropic ideal-gas path with efficiency), conversion reactor (element-balanced stoichiometry, heat of reaction from ΔfH°, isothermal or adiabatic), component separator (split fraction per component) | **done** (`src/process.cpp`); validation in "Unit operations: validation" below; used from the CLI, C API, Python and the web through the flowsheet solver. Not covered: equilibrium and kinetic (CSTR/PFR) reactors, real-gas vapour, liquid-liquid splits |
| Flowsheet solver: units connected by named streams in a text format (`include/chemlarp/flowsheet.hpp`), sequential-modular calculation, tear streams from the back edges of a depth-first search (or chosen by the user, with initial guesses), Wegstein acceleration (q bounded to [-5, 0]) or direct substitution on component flows and enthalpy, overall mass and energy balance; `chemlarp flowsheet`, `cl_flowsheet_json`, Python `flowsheet()`, web Flowsheet panel with an automatic process flow diagram | **done** (`src/flowsheet.cpp`); validation in "Flowsheet solver: validation" below. Not covered: design specifications and controllers (adjusting a unit parameter to meet a target), simultaneous (equation-oriented) solution, Broyden acceleration, minimal tear sets for nested loops, a graphical editor |
| Aqueous solution chemistry: speciation (36 components including valence states, 188 aqueous species, 64 minerals and gases from phreeqc.dat), pH from the charge balance or fixed, ionic strength, activity coefficients (WATEQ Debye–Hückel or Davies; 0.1 I for neutral species; water activity), log K(T) from analytic expressions or van 't Hoff over 0–100 °C, saturation indices, equilibrium with minerals and gases (dissolution limited by the amount present, precipitation), alkalinity; input by component totals or by formula (`add NaHCO3 2`, valence states from the formula); `chemlarp solution`, `cl_aqueous_json`, Python `solution()`, web Aqueous Chemistry panel with species, saturation indices and a distribution-with-pH chart | **done** (`src/aqueous.cpp`); validation in "Aqueous chemistry: validation" below. Not covered: redox equilibria between valence states (pe), Pitzer equations for brines (the model is for ionic strength below about 0.5–1 mol/kg), changes in the mass of water, pressure effects, kinetics, surface complexation and ion exchange, organic acids (not in phreeqc.dat) |
| Distillation columns: rigorous equilibrium-stage model (MESH equations with the stream model of the unit operations: modified Raoult with ideal or UNIFAC liquid, liquid and ideal-gas enthalpies), total or partial condenser, kettle reboiler, several feeds, linear pressure profile, reflux ratio with distillate or bottoms flow; Naphtali–Sandholm Newton solution with a block-tridiagonal finite-difference Jacobian and line search, started from constant molar overflow and bubble-point (Wang–Henke) passes (`include/chemlarp/column.hpp`); `column` unit in flowsheets (CLI, C API, Python), web stage-profile chart | **done** (`src/column.cpp`); validation in "Distillation columns: validation" below. Not covered: other specifications (purity, recovery, boilup or duty), side draws and pumparounds, stage efficiencies, two liquid phases, non-condensable components in the column, shortcut (Fenske–Underwood–Gilliland) design |
| Installable web app (WebAssembly, PWA, GitHub Pages) | **done**: https://ayuroy01.github.io/ChemLarp/ |
| Worker pool | planned (design above) |
| Node.js addon | planned (C ABI is ready for it) |
| IUPAC naming | **done** for structure to name (`src/iupac.cpp`; scope, validation and limits in the "IUPAC names" row above). Done criteria, all met on 2026-11-13: no name that OPSIN parses to a different structure on any of the four validation sets (14,996 PubChem compounds); at least 93% of each set named; every structure outside the scope gets a NamingUnsupported diagnostic instead of a guessed name, including a principal group that could only be cited as a prefix; preferred-name rules for seniority, numbering, ring assemblies and multiplicative esters, with the remaining non-preferred forms listed above. Not part of this module: name-to-structure parsing, and coordination and organometallic nomenclature |

## Activity models: validation

All checks below use data not produced by ChemLarp. Tests: `tests/test_activity.cpp`.

**Equations against an independent implementation.** ln γ from `nrtl_ln_gamma` and `unifac_ln_gamma` were compared with thermo 0.6.1 (C. Bell, MIT licence) on 300 random 2-4 component NRTL systems (250-450 K, including a, b, e, f and α terms and infinitely dilute components) and 400 random 2-4 component UNIFAC systems built from subgroups with published interaction parameters. Largest difference: 5.0e-14 in ln γ. thermo writes α = c + d T, ChemLarp α = c + d (T − 273.15 K) (the Aspen convention); the comparison shifts c accordingly. Five of these cases are in the tests.

**Limits (tests).** γ = 1 for a pure component and for a mixture of identical components; NRTL with τ = 0 gives γ = 1 and Raoult's law for bubble and dew pressures (to 1e-12); NRTL infinite dilution matches ln γ∞ = τ21 + τ12 exp(−α τ12); UNIFAC satisfies Gibbs–Duhem (Σ x_i d ln γ_i = 0 to 1e-7 by central differences); the dew point of a bubble-point vapour returns the bubble temperature and liquid (1e-6). Invalid input (wrong length, negative or zero mole fractions, non-positive T or P, unknown subgroup, missing NRTL α) gives InvalidMixtureInput; a pair of main groups without a published a_mn (water + iodomethane) gives NoInteractionParameter and is never taken as zero; a compound without Antoine data gives NoVaporPressureData; use outside an Antoine range sets `extrapolated` and warns.

**Group assignment against DDBST.** DDBST's published original-UNIFAC assignments for 28,846 compounds (InChIKeys, as distributed with thermo/chemicals; SMILES from chemicals' PubChem table) were split at random into a development half (14,418, used while writing the SMARTS rules) and a held-out half (14,405, evaluated once, on 2026-11-13, after the rules were fixed).

| Set | Assigned | Identical to DDBST | Different | Refused |
|---|---|---|---|---|
| development | 12,814 (88.9%) | 12,778 (99.7% of assigned) | 36 | 1,604 |
| held-out | 12,749 (88.5%) | 12,718 (99.8% of assigned) | 31 | 1,656 |

The differences are mostly cases where DDBST's own choice is not determined by the group definitions: CH2 next to both an aromatic ring and an ether O or N (ACCH2 vs CH2O / CH2NH), non-benzenoid aromatics (azulene, cyclophanes) that DDBST splits into olefinic and aromatic groups, and two radicals. Refusals are structures with atoms or groups outside the supported subgroups (silicon, amides, carbonates, most heteroaromatics other than pyridine and thiophene); they get NoUnifacGroups, never a partial assignment.

**Vapour pressures.** Three Antoine rows that miss a reference normal boiling point lying in (or within 5 K of) their fitted range by more than 2 K are dropped by the generator (trans-1,4-dimethylcyclohexane: 429.7 K against 392.5 K; nitrogen trifluoride; dinitrogen tetroxide). The Antoine fits reproduce CRC normal boiling points of the 12 compounds used below within 0.42 K (acetone, −0.42 K; all others within 0.25 K); 11 are in the tests with a 0.5 K tolerance.

**Azeotropes at 101.325 kPa against experiment.** Reference: CRC Handbook of Chemistry and Physics (44th ed.) and Lange's Handbook of Chemistry (10th ed.), as compiled in the Wikipedia article "Azeotrope tables" (ethanol + cyclohexane from the Edinburgh Azeotrope Databank; methanol + hexane from CRC 101st ed.); composition converted from weight % to mole fraction. Homogeneous azeotropes only. The azeotrope is located where y − x changes sign along the UNIFAC bubble-point curve. x_B is the mole fraction of the second component.

| System (A + B) | T exp / K | T UNIFAC / K | ΔT / K | x_B exp | x_B UNIFAC | Δx_B |
|---|---|---|---|---|---|---|
| water + ethanol | 351.25 | 351.16 | −0.09 | 0.892 | 0.897 | +0.005 |
| ethanol + cyclohexane | 338.05 | 337.84 | −0.21 | 0.555 | 0.561 | +0.006 |
| ethanol + hexane | 331.85 | 331.69 | −0.16 | 0.668 | 0.663 | −0.005 |
| ethanol + toluene | 349.85 | 349.94 | +0.09 | 0.190 | 0.191 | 0.000 |
| ethanol + chloroform | 332.55 | 332.69 | +0.14 | 0.837 | 0.858 | +0.022 |
| methanol + benzene | 331.45 | 331.02 | −0.43 | 0.386 | 0.406 | +0.021 |
| methanol + toluene | 336.95 | 337.00 | +0.05 | 0.135 | 0.111 | −0.024 |
| methanol + acetone | 328.85 | 328.38 | −0.47 | 0.800 | 0.781 | −0.019 |
| methanol + chloroform | 326.65 | 326.80 | +0.15 | 0.651 | 0.657 | +0.006 |
| methanol + tetrachloromethane | 328.85 | 329.01 | +0.16 | 0.445 | 0.454 | +0.009 |
| acetone + chloroform (maximum-boiling) | 337.85 | 337.58 | −0.27 | 0.661 | 0.628 | −0.033 |
| acetone + hexane | 322.95 | 322.89 | −0.06 | 0.319 | 0.369 | +0.050 |
| acetone + cyclohexane | 326.15 | 326.74 | +0.59 | 0.254 | 0.244 | −0.010 |
| benzene + cyclohexane | 350.95 | 350.47 | −0.48 | 0.432 | 0.446 | +0.015 |
| water + 1-propanol | 360.85 | 361.42 | +0.57 | 0.432 | 0.422 | −0.009 |
| water + 2-propanol | 353.55 | 353.37 | −0.18 | 0.681 | 0.690 | +0.009 |
| methanol + hexane | 333.15 | 322.42 | −10.73 | 0.484 | 0.492 | +0.008 |

Sixteen of the 17 systems agree within 0.6 K and 0.05 in mole fraction, including the maximum-boiling acetone + chloroform azeotrope, and methanol + water correctly shows none. This is UNIFAC's predictive accuracy for these well-studied systems, not a numerical tolerance, and it will be worse for systems far from UNIFAC's fitting data. Methanol + hexane is unresolved: the tabulated 60 °C lies above the methanol + heptane azeotrope in the same table (59.1 °C), against the usual trend in a homologous series, and UNIFAC gives 49.3 °C. The source value has not been checked against a second compilation, so the system is reported here and left out of the tests. Five systems (water + ethanol, methanol + benzene, acetone + chloroform, benzene + cyclohexane, water + 2-propanol) are in the tests with tolerances of 0.6 K and 0.05.

## Pure-component properties: validation

Data come from published compilations (`scripts/gen_property_data.py` names the source of each field); each correlation below is compared with a source that did not supply its inputs. Scratch scripts compute the comparisons through `cl_pure_component_json`; spot values are in `tests/test_properties.cpp`.

**Errors found in the compilations, and how the generator handles them.** These are confirmed by a disagreement between independent sources, not by our own estimates.

- IUPAC critical table as distributed: rows of cyclohexane and methylcyclopentane exchanged; decane Tc = 611.7 K against 617.6-618.1 K elsewhere. Rule: IUPAC first, then CRC, then PSRK; an IUPAC row that disagrees with CRC (Tc 0.5%, Pc 5%) is replaced when another source (PSRK, Yaws, or the Tc and Pc published with a Wagner equation) confirms CRC. This replaces 10 IUPAC rows. Older compilations never outvote an evaluated one on their own, since they often share an origin.
- Poling ideal-gas Cp polynomials: a polynomial that misses the databank's own Cp(298.15 K) by more than 2% is kept only if a third source (NIST-JANAF, else TRC) confirms it. Kept: acetylene, propadiene, chlorine, two perfluoroalkanes, butylamine, 2-methylnaphthalene (the tabulated Cp was the wrong value). Dropped: 3-methylpentane (coefficients mis-scaled, 33.5 against 140.1 J/(mol K)), cis- and trans-2-butene (polynomials apparently exchanged), sulfur dioxide, trideuteroammonia.
- Acentric factor: a tabulated ω is tied to the Tc and Pc of its own table (PSRK gives 2-methyl-1-pentanol ω = 0.726 with Tc = 582 K; IUPAC has 604.4 K). ω is therefore computed from its definition, −log10(Psat(0.7 Tc)/Pc) − 1, with the chosen Tc and Pc and a Wagner (Poling, McGarry) or Antoine equation valid at 0.7 Tc; a tabulated ω is used only when its table's Tc and Pc agree with ours, and otherwise ω is left missing. Against Yaws's ω (311 compounds), median |Δω| fell from 0.005 to 0.002 and the largest from 0.216 to 0.137.

**Results** (median and 90th percentile of the absolute deviation):

| Property | Reference | n | Median | 90th pct. | Notes |
|---|---|---|---|---|---|
| Ideal-gas Cp, 300-900 K | TRC polynomial fits | 1,475 points | 0.19% | 0.65% | worst: formic acid (Poling 53.4 against TRC 45.7 J/(mol K) at 298 K) and 3,5-lutidine, where the two sources disagree |
| Ideal-gas Cp, 300-900 K | NIST-JANAF tables | 225 points | 0.30% | 1.7% | hydrazine and SF4 3-7% below JANAF (Poling and TRC agree with each other) |
| ΔHvap(T), Watson, 0.55-0.85 Tc | VDI PPDS correlations | 552 points | 1.0% | 3.9% | carboxylic acids up to 28% (vapour association; Watson does not apply); ethyl propyl ether 25% at high Tr |
| ΔHvap(298 K), Watson from Tb | CRC ΔHvap(298) | 180 | 0.75% | 4.7% | long extrapolations (eicosane −11%) and acids |
| Saturated liquid density, Rackett with Z_RA, 0.5-0.8 Tc | VDI PPDS | 316 points | 0.51% | 3.4% | N2O4 (equilibrium mixture with NO2). Since the unit-operations work, the VDI equation itself is used where tabulated (Rackett gives water 1065 kg/m3 against 997, which made pump work 6% low); Rackett remains the fallback |
| Saturated liquid volume, VDI PPDS, Tm to 0.9 Tc | Perry DIPPR-105 | 414 points | 0.08% | 0.47% | 1-heptanol −4.4% near 0.9 Tc; water within 0.25% of IAPWS-95 at 300, 350, 450 K |
| Saturated liquid density, Rackett with Zc (no Z_RA) | VDI PPDS | 356 points | 2.4% | 12.7% | reported as a warning; HF and formic acid up to 38% |
| Tc | second compilation | 287 | 0.10 K | 1.0 K | the largest are the IUPAC errors above |
| Pc | second compilation | 260 | 0.30% | 2.4% | pyrrole (5.7 against 8.0 MPa) and propadiene unresolved |
| Peng–Robinson Psat with database Tc, Pc, ω, Tr = 0.6 | McGarry Wagner equations | 156 | 1.3% | 8.2% | alcohols worst (2-methyl-2-butanol 36%), a known limit of the PR α function |
| Peng–Robinson Psat, Tr = 0.9 | McGarry Wagner equations | 163 | 0.92% | 4.7% | acetic anhydride −32% |

**Unresolved.** Critical pressures of alkylphenols: IUPAC and CRC give about 4.1 MPa for p-cresol, PSRK, Yaws and McGarry about 5.15 MPa; the evaluated value is used and ω is computed consistently with it, but which Pc is right is not settled here. Ideal-gas Cp of formic acid and 3,5-lutidine (Poling against TRC), and ΔHvap of ethyl propyl ether at high reduced temperature (VDI against Watson from the CRC value), likewise need a fourth source.

**Transport properties.** Two checks. (1) Implementation: ChemLarp's evaluation of the Perry, Mulero–Cachadiña and Jasper coefficients agrees with thermo 0.6.1's evaluation of the same coefficients to rounding at 2,826 points (three temperatures across each fitted range). (2) Data: comparison with the VDI Heat Atlas (PPDS) correlations, an independent compilation, at three temperatures inside both fitted ranges (below 0.9 Tc for liquids, at most 1000 K for vapours), and with thermo's fits to REFPROP for about 100 fluids.

| Property | Against VDI: n, median, 90th pct. | Against REFPROP fits: n, median, 90th pct. | Notes |
|---|---|---|---|
| Liquid viscosity | 423, 0.63%, 8.9% | 199, 2.8%, 25% | VDI's p-cresol entry is broken (range 150.7-155.7 K, constant 1e-5 Pa s); Perry agrees with Viswanath & Natarajan and Dutt & Prasad for the cresols and phenol; light gases near their critical point (CF4, 1-butene, dimethyl ether) differ by 50-77% from REFPROP |
| Vapour viscosity, low pressure | 429, 0.12%, 0.91% | 116, 1.3%, 6.6% | |
| Liquid thermal conductivity | 420, 0.03%, 0.16% | 194, 3.1%, 13% | the near-identity with VDI suggests both were fitted to the same data, so this is a weak independence test |
| Vapour thermal conductivity, low pressure | 432, 0.48%, 1.8% | 109, 1.6%, 6.5% | chloroform and CCl4 diverge from VDI above about 400 K (up to 109%); carboxylic acids (vapour dimerization) up to 39% |
| Surface tension | 435, 0.95%, 7.3% | 267, 0.00%, 0.64% | REFPROP uses the same Mulero–Cachadiña fits, so that column is not independent; Jasper's dimethylamine is 27-47% above VDI, which itself differs from a third fit |

Liquid heat capacity (Perry Table 2-153) against the evaluated data of Zabransky et al. at three temperatures in each Zabransky range below 0.9 Tc: 179 points, median 0.60%, 90th percentile 6.3%. The largest deviation, 3-methyl-1-butanol at −88%, is on the reference side (Perry's row gives 210.8 J/(mol K) at 298 K, matching Poling's 211). The generator drops a Perry row that agrees with neither Poling's nor CRC's Cp(298.15 K) within 5% (2-methylpropanoic acid, 1-butyne, acetic anhydride). Water against IAPWS-95: 75.29 against 75.32 J/(mol K) at 300 K.

For water, against REFPROP (IAPWS formulations) at 300 K: liquid viscosity +2.7% (876.5 against 853.7 µPa s), liquid thermal conductivity −0.1%, surface tension +0.1%. None of these disagreements has been settled by a third, evaluated source, so no transport data are dropped; the reference values in the tests carry the tolerances above.

## Unit operations: validation

Model (`include/chemlarp/process.hpp`): ideal-gas vapour, H_V = ΔfH° + ∫Cp,ig; liquid referenced at the normal boiling point, H_L = H_V(Tb) − ΔHvap(Tb) + ∫Cp,L, with no excess enthalpy or pressure effect; phase split by modified Raoult's law with ideal or UNIFAC activity coefficients and Antoine vapour pressures. The liquid reference was chosen over H_V(T) − ΔHvap(T) (Watson) because the latter gives the sensible heat of liquid water about 25% too high (10 bar → 1 atm throttling would give a vapour fraction near 0.19 instead of 0.152).

Water against IAPWS-95 (computed with the `chemicals` implementation; tests/test_process.cpp):

| Case | ChemLarp | IAPWS-95 |
|---|---|---|
| Saturated liquid throttled to 1 atm from 2 / 5 / 10 / 20 bar: vapour fraction | 0.0380 / 0.0979 / 0.1519 / 0.2159 | 0.0380 / 0.0980 / 0.1522 / 0.2169 |
| Saturation temperature at 2 / 10 / 20 bar (Antoine, extrapolated above 373 K) | 393.48 / 453.01 / 485.37 K | 393.36 / 453.03 / 485.53 K |
| Liquid heated 300 → 360 K at 1 atm | 4522.1 J/mol | 4524.8 J/mol |
| Adiabatic mixing of equal flows at 290 K and 350 K | 320.004 K | 320.020 K |
| Steam heated 400 → 600 K at 1 atm | 7075 J/mol | 7178 J/mol (−1.4%: ideal gas and the Cp polynomial) |
| 20 kJ/mol added to saturated liquid at 1 atm: vapour fraction | 0.4921 | 0.4920 |

Mixtures: the benzene + toluene split at 368 K and 1 atm (ideal solution) agrees with an independent Rachford–Rice solve to 1e-10; below the bubble point and above the dew point the stream is one phase. For ethanol + water (UNIFAC) the units conserve every component flow to 1e-12 and the energy balance to 1e-7 W per mol/s; a heater by duty inverts the heater to temperature to 1e-6 K, and remixing the two phases of an adiabatic flash returns its state.

Further units (tests/test_process.cpp), against IAPWS-95 unless stated:

| Case | ChemLarp | Reference |
|---|---|---|
| Pump, water 300 K, 1 → 100 bar: work | 179.0 J/mol | 178.6 J/mol |
| Same at 75% efficiency: outlet temperature | 300.79 K | 301.00 K (the model omits the 0.2 K isentropic rise) |
| Compressor, argon 300 K, 1 → 5 bar, isentropic | 571.096 K | 571.096 K (ideal gas, T1 (P2/P1)^(2/5)) |
| Compressor, steam 450 K, 0.1 → 0.5 bar: isentropic work, outlet | 7312.9 J/mol, 653.58 K | 7313.3 J/mol, 654.40 K |
| Same at 80% efficiency: outlet | 702.69 K | 703.42 K |
| Heat exchanger, liquid water 360 → 320 K against 290 K: duty, cold outlet | 3015.8 J/mol, 330.048 K | 3018.9 J/mol, 330.079 K |
| Methane combustion at 298.15 K, all gas | −802.5 kJ/mol | −802.3 kJ/mol (CODATA ΔfH°) |
| Same with liquid water (10 bar), corrected for the 0.3% that stays as vapour | −890.0 kJ/mol | −890.56 kJ/mol (NIST) |

The exchanger follows both temperature profiles through 20 zones, so a pinch inside a phase change (condensing steam against water) is found, and a temperature cross is an error. The conversion reactor checks element balance from the structures and requires ΔfH° for every reacting component.

Not validated yet: mixtures against measured enthalpy data (heats of mixing are neglected by the model), and pressures where the ideal-gas vapour fails (the steam case above shows the size of the error at 1 atm; it grows with pressure).

## Flowsheet solver: validation

Each case is posed in the flowsheet format, solved by the sequential-modular solver, and compared with a solution that does not use it (`scripts/validate_flowsheet.py`, which needs numpy, scipy and thermo; tests/test_flowsheet.cpp keeps the closed-form cases).

| Case | Reference | Deviation |
|---|---|---|
| Propane dehydrogenation with recycle (Felder & Rousseau, *Elementary Principles of Chemical Processes*, 3rd ed., Example 4.7-2): 95% overall conversion, 0.555% of the reactor-outlet propane and 95% of the propylene to product | Closed-form recycle balance; published figures: 995.9 mol propane to the reactor per 100 mol fresh feed, recycle ratio 9.0 | 1e-9 (Wegstein, 26 iterations; direct substitution, 197) against the closed form; the published figures to their rounding |
| CO2 + 3 H2 → CH3OH + H2O loop with condenser and purge (Felder & Rousseau Example 4.7-3 statement: reactor feed 28/70/2% CO2/H2/inerts, 60% single-pass H2 conversion, 0.4% inerts in the fresh feed) | Hand solution on a 100 mol reactor-feed basis, posed back as fresh feed and purge fraction | 1e-10 on the reactor feed, methanol and purge (20 iterations) |
| Benzene + toluene flash at 368 K, 1 atm, with 60% of the liquid recycled | Equation-oriented Newton solve (scipy) of the loop's balances with Raoult's law: (a) with ChemLarp's K-values, (b) with thermo's vapour pressures | (a) 2e-12; (b) 0.8–1.0% in the vapour flows, from vapour pressures 0.14% apart (the recycle amplifies them) |
| Nitrogen loop: heater with 0.5 MW duty, 50% of its outlet recycled | Steady state F (h(T) − h(300 K)) = Q with thermo's ideal-gas Cp | 0.09 K (2e-4) in the product temperature; energy balance 8e-14 |

Other checks: re-running the loop's units from the converged state reproduces the recycle; Wegstein and direct substitution reach the same point; a feed–effluent exchanger loop (toluene hydrodealkylation) closes its energy balance to 1e-3 W; a loop whose specifications fix the product enthalpy independently of the heater duty (ill-posed) is reported as not converged rather than solved; a closed loop that no feed reaches is flagged. The convergence test is on the relative change of the tear streams (default 1e-7); with direct substitution the remaining error can be larger by 1/(1 − loop gain), so slow loops should use Wegstein or a tighter tolerance.

Starting values matter for units with outlet specifications: a recycle is started empty unless it has a `guess`, and a unit that cannot meet its specification on the first passes (an exchanger outlet temperature above its hot inlet) stops the solve with a hint to add a guess.

## Distillation columns: validation

`scripts/validate_column.py` writes the same stage model as one system of equations and solves it with scipy's hybrid Powell method using thermo's correlations instead of ChemLarp's data: vapour pressures, ideal-gas Cp, enthalpies of vaporization (liquid H = H_V − ΔHvap(T)) and, for ethanol–water, thermo's original UNIFAC. The comparison therefore covers both the solution method and the property data.

| Case | Quantity | ChemLarp | Independent |
|---|---|---|---|
| Benzene/toluene, 15 stages, feed (360 K liquid) on 8, R = 2, D = 50 mol/s, 1 atm, ideal | benzene in distillate | 48.734 mol/s | 48.721 mol/s |
| | top / feed stage / bottom T | 353.67 / 366.24 / 382.57 K | 353.73 / 366.26 / 382.55 K |
| | condenser / reboiler duty | −4631.5 / 4785.8 kW | −4636.5 / 4798.5 kW |
| Ethanol/water, 20 stages, feed (350 K) on 12, R = 3, D = 11 mol/s, 1 atm, UNIFAC | distillate ethanol fraction | 0.8402 | 0.8364 |
| | top / bottom T | 351.22 / 370.27 K | 351.35 / 370.07 K |
| | condenser / reboiler duty | −1718 / 1858 kW | −1741 / 1887 kW |
| | ethanol left in the bottoms | 0.757 mol/s | 0.800 mol/s (a small trace, sensitive to γ at infinite dilution) |

Solver checks (tests/test_column.cpp): the inflows of every intermediate stage, flashed with the separate PH flash at the stage's pressure and enthalpy, reproduce the stage temperature to 1e-5 K and its vapour flow and compositions to 1e-7; component and energy balances close to 1e-8 mol/s and 1e-3 W; a partial condenser's vapour distillate satisfies y = K x with its reflux; more reflux gives a purer distillate at a higher reboiler duty; a distillate or a bottoms specification gives the same column; the ethanol–water distillate stays below the azeotrope. Newton converges in 3 iterations for benzene/toluene and about 25 (with line search) for ethanol–water, in milliseconds. A column also converges inside a flowsheet recycle.

## Aqueous chemistry: validation

Reference: PHREEQC 3 (IPhreeqc through phreeqpython) with the same phreeqc.dat, an independent implementation of the same model; `scripts/validate_aqueous.py` runs the comparison. Solutions are defined in PHREEQC with the pH from the charge balance.

| Case | Largest difference |
|---|---|
| Pure water 0, 25, 50, 100 °C | pH 4e-9 |
| NaCl 100 mmol/kg; HCl; NH4Cl; NH3 (charge balance); phosphate buffer; AlCl3 and FeCl3 hydrolysis; MgSO4 + K2SO4 | pH 0.0016, species 0.1% |
| CaCl2 1 + NaHCO3 2 mmol/kg at 25 and 60 °C | pH 2e-5, species 0.06%, calcite SI 0.0003 |
| Seawater (I = 0.67) | pH 0.0003, species 0.25%, gypsum SI 0.018 |
| Gypsum, fluorite, anhydrite in pure water | dissolved amount 0.1% |
| H2S(g), NH3(g), N2(g) at fixed partial pressure | pH 2e-5, dissolved gas 5e-5 |
| Every phase's saturation index in a solution with traces of all 26 components, 25 and 60 °C (54 phases) | 0.0025 |
| Calcite with CO2(g) at 10^-3.5 atm, 10–80 °C | pH 0.002–0.010, dissolved Ca 0.3–1.0% (see below) |

The remaining differences in dilute solutions come from details ChemLarp leaves out: PHREEQC computes the Debye–Hückel A and B from its own water properties (log γ differ by about 2e-4 at I = 0.1), lets the mass of water change, and adds ionic-strength terms from the -Vm parameters at high ionic strength.

Two known differences in the constants: PHREEQC raises the CO2(g) solubility constant by 0.002 (10 °C) to 0.017 (80 °C) log units and lowers CH4(g)'s by 0.04 with corrections not given in phreeqc.dat; ChemLarp uses the database constants (for CO2, Plummer & Busenberg 1982: log K_H = −1.468 at 25 °C, within the spread of measured values). For barite, PHREEQC's constant (−9.970 at 25 °C) differs from the analytic expression phreeqc.dat gives (−9.844), which ChemLarp evaluates as written; barite solubility is therefore 16% higher than PHREEQC's (unresolved).

Against measurements: pKw at 25 °C 13.9948 (IUPAC 13.995); gypsum solubility 15.10 mmol/kg (measured 15.45, Marshall & Slusher 1966; PHREEQC 15.09); rain water at 10^-3.5 atm CO2 pH 5.66; calcite with atmospheric CO2 0.493 mmol/kg Ca, pH 8.28.

## Molecular dynamics: validation

`scripts/validate_md.py` (numpy) checks the dynamics against results that do not depend on it:

| Check | Reference | Result |
|---|---|---|
| Lennard-Jones fluid, N = 108, r_c = 2.5 (truncated and shifted), ρ* = 0.8, T* = 1.0 | Metropolis Monte Carlo of the same model written in the script (6000 sweeps) | U/N −4.689 ± 0.005 (Bussi), −4.691 ± 0.004 (Langevin) against −4.702 ± 0.005; P 1.65 ± 0.02 and 1.65 ± 0.02 against 1.60 ± 0.03 (within 2σ) |
| Same at ρ* = 0.4, T* = 1.5 | Monte Carlo | U/N −2.238 / −2.247 against −2.230 ± 0.004; P 0.451 / 0.443 against 0.442 ± 0.005 |
| Kinetic temperature with both thermostats | target T | within 1.5% |
| NVE energy drift, Lennard-Jones, dt = 0.002 | energy conservation | 1e-7 per particle per unit time |
| Harmonic limit (1 K, NVE): peaks of the velocity spectrum of water, formaldehyde and HCN | normal-mode frequencies from the finite-difference Hessian | every mode found within the 5 cm⁻¹ spectral grid (HCN: 3N − 5 = 4 modes with the degenerate bend) |
| Equipartition at 10 K: mean potential energy above the minimum | (3N − 6) kB T / 2 (classical harmonic limit) | ethanol 0.2110 / 0.2099 kcal/mol (Bussi / Langevin) against 0.2087; benzene 0.3017 against 0.2981 (the excess is anharmonicity) |
| Velocity Verlet energy fluctuation, ethanol at 300 K, dt = 0.4, 0.2, 0.1 fs | O(dt²) | ratios 3.92 and 4.00 |

Tests (tests/test_md.cpp) also check time reversibility (forward 2000 steps, reverse the velocities, back to the start within 1e-8 Å), reproducibility for a seed, and equipartition for water with both thermostats. The UFF frequencies themselves are those of the force field, which is not fitted to vibrational spectra (water: 1757, 3706, 3760 cm⁻¹ against 1595, 3657, 3756 measured).

## MMFF94: validation

`scripts/validate_mmff.py` runs two independent comparisons.

**Official MMFF94 validation suite** (Halgren / Merck, 761 molecules, www.ccl.net/cca/data/MMFF94/; not redistributed because the structures derive from the Cambridge Structural Database): the hypervalent and dative mol2 files with their formal-charge lists, evaluated at the 6-decimal coordinates of `MMFF94.mmd`, against Merck's OPTIMOL output (`MMFF94_opti.log`).

| Quantity | Result (both representations) |
|---|---|
| Molecules read | 759 of 761 (COVXIU and DOTVEN give a three-valent N− in both files, which the structure reader rejects) |
| Atom types | 17,231 of 17,231 atoms identical |
| Partial charges (3 printed decimals) | identical except FEZPOP, whose formal-charge file gives +1 instead of −1 to two oxygens of an anion |
| Total energy | 752 of 759 within 1e-4 kcal/mol; ERULE_01, 04, 06 within 4e-4; ERULE_03 0.005 (empirical Si–P bond length 0.0013 Å shorter); FEZPOP off because of its charges |
| Each term (bond, angle, stretch-bend, out-of-plane, torsion, van der Waals, electrostatic) | largest difference 0.0044 kcal/mol (ERULE_03 bond); out-of-plane, van der Waals and electrostatics below 5e-5 |
| Torsion parameters, interaction by interaction | 41,127 of 41,133 from the tables as OPTIMOL assigns them; the other six use the empirical rule |

Rules settled against the suite rather than the papers alone: the torsion step-down tries the exact types, then one outer atom at equivalence level 1 with the other a wildcard, then both wildcards (no level-2 step); the bond-type torsion classes take precedence over the five-ring class, and class 2 needs a single central bond; an N bonded to N=C (hydrazone), N–NO2 or N–N=O stays an amine, while N–N=N is type 10; charges of nitro, N-oxide, azide, isonitrile and phosphoryl groups are counted by the group's net charge; angle constants use 143.9325·(π/180)² exactly (the rounded 0.043844 leaves 1e-4–7e-4 kcal/mol errors).

**RDKit's MMFF94** (independent implementation, BSD) on 80 molecules entered as SMILES (drugs, ions, zwitterions, heteroaromatics, sulfonamides, phosphates, nitro and N-oxide groups), at the same coordinates: atom types identical in all 80, partial charges within 1e-4 in 79, total energies within 1e-4 kcal/mol. The exception, tetramethylphosphonium, is a difference in RDKit: it gives the ion no formal charge (net charge 0), while ChemLarp keeps +1. `data/mmff_reference.sdf` stores twelve of these molecules with RDKit's MMFF94 and MMFF94s types, charges and energies for tests/test_mmff.cpp, which also checks the gradient against finite differences, minimization and energy conservation in MMFF dynamics.

