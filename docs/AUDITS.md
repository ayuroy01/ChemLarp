# Audits

The audits to run before ChemLarp is called production-ready, what each covers, how to run it, and where it
stands. Findings are recorded in the "Findings" column; open items move to docs/GAPS.md once they are
understood. Status: **open** (not started), **partial**, **done** (date).

Started 2026-11-17.

## 1. Correctness (bugs)

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 1.1 | Unit tests on every platform | `make && build/chemlarp_tests`; CI on Windows (MSVC, MinGW), macOS, Linux | done (CI) | 241 tests pass |
| 1.2 | Undefined behaviour | UBSan build of the test suite (command in docs/ARCHITECTURE.md; ASan hangs on macOS) | done 2026-11-15 | clean |
| 1.3 | Memory errors | ASan and LeakSanitizer on Linux (CI job), Valgrind on the CLI | open | ASan hangs locally on macOS; needs a Linux run |
| 1.4 | Threading | ThreadSanitizer on the parallel kernels (QC, PME, solvation, sparse LU) | open | |
| 1.5 | Fuzzing of every text input | libFuzzer targets for SMILES, SMARTS, MOL/SDF, PDB/mmCIF, CML, .ocm, CHEMKIN, YAML, flowsheet, lab and aqueous inputs; at least 1 h each, crashes and hangs fixed | partial | SMILES round trip (14,000 cases) and 1,500 SMARTS fuzzed; no coverage-guided fuzzing yet |
| 1.6 | Scientific validation | `make validate DATA=… PY=…` (`scripts/validate.py`): all suites against independent references | done 2026-11-14, hydration suite running | see docs/ARCHITECTURE.md |
| 1.7 | Numerical robustness | Extreme inputs: 0 K, very high T and P, zero amounts, huge molecules, empty files; must fail with a diagnostic, never crash or loop | open | |
| 1.8 | Static analysis | clang-tidy (bugprone, cert, performance), cppcheck, MSVC /analyze; compiler warnings at -Wall -Wextra -Wpedantic | open | |
| 1.9 | Web app | Every view with every example in Chrome, Safari, Firefox, Android and iOS; console errors; offline mode; very large files | partial | headless Chrome run of all 15 views on 2026-11-15, no errors seen |
| 1.10 | Packages | Install every release asset on a clean machine (Windows, macOS arm64/x86_64, Linux, Android) and run its smoke test | partial | CI smoke-tests all; not yet on clean physical machines |

## 2. Security

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 2.1 | Secrets in the repository | Scan all history for keys, tokens, passwords (gitleaks or trufflehog) | partial 2026-11-17 | Pattern scan of all history: no keys, tokens or passwords. Run gitleaks for a thorough pass |
| 2.2 | Personal data in the repository | Emails, home paths, machine names in history | partial 2026-11-17 | Commits use the GitHub noreply address; no personal email or home path found in history |
| 2.3 | Web app: cross-site scripting | `innerHTML` uses (5 in `web/app.js`) take engine-generated SVG; check that every string from user files (SDF titles, PDB headers, names) is escaped in the SVG and HTML; add a Content-Security-Policy | open | No CSP in `web/index.html` |
| 2.4 | Web app: network | The app should only contact its own origin and RCSB (protein download); confirm in the network log; no analytics | open | |
| 2.5 | Native apps | macOS/Android local servers bind 127.0.0.1 only, refuse path traversal (macOS `--selftest` checks this), no other open ports; Windows virtual host mapping | partial | macOS self-test in CI |
| 2.6 | Input handling in the C API | Null pointers, very long strings, malformed JSON; results freed with `cl_free` | partial | `guarded()` wraps every call |
| 2.7 | Supply chain: CI | Third-party GitHub Actions are pinned by tag (e.g. `mymindstorm/setup-emsdk@v14`, `reactivecircus/android-emulator-runner@v2`); pin to commit SHAs, keep `permissions:` minimal, enable Dependabot for actions | open | 14 actions, all pinned by tag only |
| 2.8 | Supply chain: builds | Reproducible builds; SHA256SUMS published (done); sign releases (Sigstore/cosign or GPG); code-sign and notarise the macOS app, Authenticode for Windows | partial | SHA256SUMS in v0.1.0; no signatures |
| 2.9 | Repository settings | Branch protection on `main` (no force push, require CI), 2FA on the owner account, secret scanning and push protection, Dependabot alerts, private vulnerability reporting | open | `main` is not protected |
| 2.10 | Security policy | `SECURITY.md` with how to report a vulnerability | open | |

## 3. Integrity of builds and of this machine

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 3.1 | Release assets match the source | Rebuild a release from its tag and compare checksums; check SHA256SUMS against the downloads | open | |
| 3.2 | Development machine | macOS up to date, FileVault on, firewall on, Gatekeeper/SIP enabled (`csrutil status`, `spctl --status`, `fdesetup status`), no unknown login items or launch agents, Homebrew packages from known taps | open | needs the owner's permission to inspect system settings |
| 3.3 | Toolchain | Compilers, Emscripten SDK, Python venvs from official sources and versions recorded | open | |
| 3.4 | Downloaded reference data | Every validation download kept outside the repository, in its own folder, with its source URL (docs/GAPS.md and `scripts/validate.py` record them) | done | |
| 3.5 | Git object integrity | `git fsck --full`; signed commits and tags (SSH or GPG signing, "Verified" on GitHub) | open | commits are not signed |

## 4. Git hygiene

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 4.1 | Authorship | All commits by the owner; no third-party or tool attribution lines | done 2026-11-17 | 83 commits, all `ayuroy01`; no `Co-Authored-By` lines |
| 4.2 | Large files | Nothing over a few MB in history | done 2026-11-17 | largest: two README GIFs (3.1 and 2.6 MB) |
| 4.3 | Ignored files | Build output, venvs, `.DS_Store`, local agent files (`CLAUDE.md`, `AGENTS.md`) never committed | partial | `CLAUDE.md` and `AGENTS.md` are untracked but not in `.gitignore` |
| 4.4 | Tags and releases | Annotated, signed tags; release notes per version; CHANGELOG | partial | v0.1.0 lightweight tag; no CHANGELOG |

## 5. Licensing and ownership

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 5.1 | Project licence chosen deliberately | Decide between permissive (Apache-2.0, current), copyleft (GPL-3.0 / AGPL-3.0) and source-available (PolyForm Noncommercial, BSL) | open, owner's decision | Apache-2.0 lets anyone reuse, modify and sell the code, including in closed products, as long as they keep the licence and notices. v0.1.0 is already out under Apache-2.0 and stays so |
| 5.2 | Copyright and attribution notices | `NOTICE` file (Apache-2.0 §4 requires derivative works to carry it), copyright line in `LICENSE`, SPDX header in each source file | open | No NOTICE file; no per-file headers |
| 5.3 | Name and logo | "ChemLarp" as a trademark: trademark policy in the repo; registration is optional and outside the code | open | |
| 5.4 | Contributions | Contributor licence agreement or Developer Certificate of Origin, so that the owner can relicense later | open | |
| 5.5 | Bundled data: licences | Check each bundled data set can be redistributed under the project licence, and credit it | partial | See the table below |

### Bundled data provenance

| Data | Where | Source | Redistribution |
|---|---|---|---|
| Aqueous species, minerals (phreeqc.dat) | `src/aqueous_data.inc` | USGS PHREEQC | US government work, public domain |
| MMFF94 parameters | `src/mmff_data.inc` | Halgren, J. Comput. Chem. 17 (1996), © Merck | Published supplementary material, redistributed by RDKit and Open Babel; keep the copyright line |
| Pure-component constants, transport, Antoine | `src/property_data.inc`, `transport_data.inc`, `antoine_data.inc` | via the `chemicals` library (MIT), which compiles DIPPR, Perry's, Poling–Prausnitz–O'Connell, VDI, Yaws | Facts are not copyrightable in the US, but some compilations (DIPPR, VDI, Yaws) are commercial; review each column's source and remove or replace any that are not openly licensed. **Highest licensing risk** |
| UNIFAC tables | `src/unifac_data.inc` | DDBST published tables via `thermo` (MIT) | published parameters; credit DDBST |
| GHS classifications | `src/ghs_data.inc` | PubChem (mostly ECHA C&L Inventory) | PubChem: public; ECHA requires attribution and no implied endorsement |
| Thermo (GRI-Mech 3.0), H₂/O₂ mechanism | `data/mechanisms/` | GRI-Mech (free with citation), Li et al. 2004 | credit given in the files |
| Protein pKa fit | `src/protein_pka.inc` | PKAD-2 measured values | check the PKAD-2 licence; fitted weights only |
| Test structures | `data/*.pdb`, `data/planar_chirality_cod.sdf` | RCSB PDB (CC0), COD (CC0) | free |
| Lab incompatibility rules | `src/lab.cpp` | NRC *Prudent Practices in the Laboratory* (2011) | rules paraphrased with section citations, no text copied |
