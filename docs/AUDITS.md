# Audits

The audits to run before ChemLarp is called production-ready, what each covers, how to run it, and where it
stands. Findings are recorded in the "Findings" column; open items move to docs/GAPS.md once they are
understood. Status: **open** (not started), **partial**, **done** (date).

Started 2026-11-17. Last run 2026-11-17: of 34 audits, 23 are done and 11 partial; none is untouched; four defects were found and fixed (audits 1.1, 1.7 and 1.8).

## 1. Correctness (bugs)

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 1.1 | Unit tests on every platform | `make && build/chemlarp_tests`; CI on Windows (MSVC, MinGW), macOS, Linux | done 2026-11-17 | 251 tests pass. CI caught a Windows build break the macOS build could not: `tests/test_solvation.cpp` used `M_PI`, which MSVC does not define, so both Windows jobs failed from the hydration commit until it was replaced with `std::numbers::pi`. An MSVC-only warning (C4390, an empty `if` in `src/iupac.cpp` whose condition was pure and therefore a no-op) was cleared at the same time |
| 1.2 | Undefined behaviour | UBSan build of the test suite (command in docs/ARCHITECTURE.md; ASan hangs on macOS) | done 2026-11-15 | clean |
| 1.3 | Memory errors | AddressSanitizer + LeakSanitizer in CI (`asan` job, Linux) | done 2026-11-17 | CI job added. On macOS ASan spins inside its own scheduler with this thread pool (confirmed again 2026-11-17), so it runs on Linux. Valgrind on the CLI: still open |
| 1.4 | Threading | ThreadSanitizer in CI (`tsan` job, Linux, 4 threads) | done | Runs on every push; was already in CI when this list was written |
| 1.5 | Fuzzing of every text input | `tests/fuzz_parsers.cpp`: one libFuzzer target covering 14 entry points (SMILES, SMARTS, MOL/SDF, PDB/mmCIF, mechanisms, flowsheet, lab, aqueous, MBAR, reactions, equilibrium, VLE), run in CI (`fuzz` job) with ASan and UBSan and a corpus seeded from the repository's own inputs | partial 2026-11-17 | Target written and wired into CI for 4 minutes a push. Locally, 10.5 million mutated inputs under UBSan in 120 s found no crash (Apple's clang ships no libFuzzer, so the local run uses a random-mutation driver). Still to do: a long soak (hours) with a kept corpus, and crashes triaged from a stored corpus rather than from scratch each run |
| 1.6 | Scientific validation | `make validate DATA=… PY=…` (`scripts/validate.py`): all suites against independent references | done 2026-11-14, hydration suite running | see docs/ARCHITECTURE.md |
| 1.7 | Numerical robustness | `tests/test_robustness.cpp`: ~60 hostile strings (empty, 20,000 atoms, 5,000 unclosed brackets, invalid UTF-8, NaN/inf text), malformed MOL/PDB/mmCIF/mechanisms, and out-of-range temperatures, pressures, amounts and step counts through every C API entry point | done 2026-11-17 | 10 tests, all passing. **Two defects found and fixed**: (a) `cl_water_box_json` and `cl_hydration_json` treated a negative argument as "use the default", so `-1` silently started a 400-water, 120 ps simulation instead of returning an error — they now reject negative, non-finite and absurd values at once; (b) the header did not say `mechanism_text` and `model` of `cl_kinetics_json` may not be NULL, which is now documented and pinned by a test |
| 1.8 | Static analysis | clang static analyzer (`clang --analyze`) over every file in `src/` | partial 2026-11-17 | 5 warnings in 3 files. **Fixed**: a possible division by zero in `src/iupac.cpp` (ring size used as a modulus without a guard; unreachable today because ring perception never returns a ring smaller than 3, now asserted explicitly). **Reviewed, no change**: two dead stores (cannot affect behaviour) and three `cplusplus.Move` reports on `std::sort` with a custom comparator, which are a known false positive of that checker. clang-tidy, cppcheck and MSVC /analyze still to run |
| 1.9 | Web app | Every view with every example in Chrome, Safari, Firefox, Android and iOS; console errors; offline mode; very large files | partial | headless Chrome run of all 15 views on 2026-11-15, no errors seen |
| 1.10 | Packages | Install every release asset on a clean machine (Windows, macOS arm64/x86_64, Linux, Android) and run its smoke test | partial | CI smoke-tests all; not yet on clean physical machines |

## 2. Security

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 2.1 | Secrets in the repository | 12 credential patterns plus a Shannon-entropy sweep over every blob in history, and GitHub secret scanning with push protection | done 2026-11-17 | 1,670 blobs scanned: no credential patterns and no high-entropy strings. GitHub's own secret scanning reports 0 alerts. Push protection is on, so a future commit containing a key is refused |
| 2.2 | Personal data in the repository | Emails, home paths, machine names in history | done 2026-11-17 | All 83 commits use the GitHub noreply address; no personal email, home path or machine name anywhere in history |
| 2.3 | Web app: cross-site scripting | A Content-Security-Policy, and tracing every string from a user's file to where it reaches the DOM | done 2026-11-17 | CSP added and tested on the live build: external scripts, external connections and `eval` are all refused, while the engine, WebAssembly, the worker and the RCSB download still work, with no violation in any of the 15 views. The five `innerHTML` uses all take engine-generated SVG, whose only variable text is the canonical SMILES (escaped for `& < >` in `src/depict.cpp`) and element labels; the one `html:` in the DOM helper builds an orbital diagram from numbers and two hard-coded Greek letters. Every string from a file (SDF titles, PDB headers) goes through `text:`, i.e. `textContent`. Tested with an SDF whose title is `<img src=x onerror=...>"><script>...`: it is escaped to text, no element is created and nothing runs |
| 2.4 | Web app: network | The app should only contact its own origin and RCSB; confirmed by the CSP's `connect-src` and by grepping for external URLs | done 2026-11-17 | `files.rcsb.org` is the only outside host in the source, and the CSP now enforces exactly that. No analytics, no telemetry |
| 2.5 | Native apps | macOS/Android local servers bind 127.0.0.1 only, refuse path traversal (macOS `--selftest` checks this), no other open ports; Windows virtual host mapping | partial | macOS self-test in CI |
| 2.6 | Input handling in the C API | Null pointers, very long strings, malformed input, out-of-range numbers | done 2026-11-17 | `guarded()` wraps every call; `tests/test_robustness.cpp` drives ~60 hostile strings and every out-of-range number through the API, and the fuzz target hammers the parsers. Two defects came out of it and are fixed (audit 1.7) |
| 2.7 | Supply chain: CI | Pin every action to a commit SHA, minimal `permissions:`, Dependabot for actions | done 2026-11-17 | All 14 actions pinned to full commit SHAs with the tag in a comment; `permissions: contents: read` added to `ci.yml` (the other two already had it); `.github/dependabot.yml` added so the pins are kept current |
| 2.8 | Supply chain: builds | Reproducible builds; SHA256SUMS published (done); sign releases (Sigstore/cosign or GPG); code-sign and notarise the macOS app, Authenticode for Windows | partial | SHA256SUMS in v0.1.0; no signatures |
| 2.9 | Repository settings | Branch protection, secret scanning, push protection, Dependabot alerts, private vulnerability reporting | done 2026-11-17 | `main` now refuses force pushes and deletion and requires linear history; secret scanning, push protection, Dependabot alerts and automated security updates, and private vulnerability reporting are all enabled. 2FA on the owner account is the owner's to check |
| 2.10 | Security policy | `SECURITY.md` with how to report a vulnerability | done 2026-11-17 | Private reporting through GitHub Security Advisories, what counts as a vulnerability for a local program, and how to verify a release against `SHA256SUMS` |

## 3. Integrity of builds and of this machine

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 3.1 | Release assets match the source | Check the published `SHA256SUMS` against a download, and that the release was built by CI from the tagged commit | partial 2026-11-17 | v0.1.0 was built by the `packages` workflow (run 43664060281) from commit `c9b698d`, which is exactly what tag `v0.1.0` points at, and a downloaded asset matches its published checksum. A bit-for-bit rebuild from source was not attempted: the builds are not reproducible yet (timestamps, paths) |
| 3.2 | Development machine | Read-only checks of the security posture and the launch agents | done 2026-11-17 | macOS 26.6.2; System Integrity Protection enabled, Gatekeeper assessments enabled, FileVault on, firewall on. Every user launch agent belongs to recognisable software (Google, Grammarly, Steam, Homebrew services, VirtualBox, Adobe, Discord, Teams, OpenAI, Anthropic, ssh-agent) — the owner should confirm they installed each. Nothing was changed |
| 3.3 | Toolchain | Record the versions used, and confirm they come from official sources | done 2026-11-17 | Local: Apple clang 17.0.0 (Command Line Tools), Python 3.12.3 (python.org), Node 22.23.4; CMake is not installed locally, so `scripts/build.sh` uses the Makefile. CI pins its toolchains through the SHA-pinned actions, and the web build uses Emscripten 3.1.56 (pinned in `.github/workflows/web.yml`) |
| 3.4 | Downloaded reference data | Every validation download kept outside the repository, in its own folder, with its source URL (docs/GAPS.md and `scripts/validate.py` record them) | done | |
| 3.5 | Git object integrity | `git fsck --full --strict`; signed commits and tags | partial 2026-11-17 | `git fsck --full --strict` is clean (only dangling objects from an amended commit, which is normal). Commits and tags are still unsigned: that needs a signing key on the owner's account |

## 4. Git hygiene

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 4.1 | Authorship | All commits by the owner; no third-party or tool attribution lines | done 2026-11-17 | 83 commits, all `ayuroy01`; no `Co-Authored-By` lines |
| 4.2 | Large files | Nothing over a few MB in history | done 2026-11-17 | largest: two README GIFs (3.1 and 2.6 MB) |
| 4.3 | Ignored files | Build output, venvs, `.DS_Store`, local agent files never committed | done 2026-11-17 | `CLAUDE.md`, `AGENTS.md` and `build*/` added to `.gitignore` (a scratch build directory was caught by this before it was pushed) |
| 4.4 | Tags and releases | Annotated, signed tags; release notes per version; CHANGELOG | partial 2026-11-17 | `CHANGELOG.md` added (Keep a Changelog, semantic versioning); v0.1.0 has written release notes. Tags are still lightweight and unsigned |

## 5. Licensing and ownership

| # | Audit | How | Status | Findings |
|---|---|---|---|---|
| 5.1 | Project licence chosen deliberately | Decide between permissive, copyleft and source-available | done 2026-11-17 | **AGPL-3.0-only**, with commercial licences offered separately (`LICENSING.md`). Anyone who distributes ChemLarp, or runs a modified copy as a network service, must publish their source under the same licence. v0.1.0 and everything before the change stay available under Apache-2.0 |
| 5.2 | Copyright and attribution notices | `NOTICE` file, copyright line, SPDX header in each source file | done 2026-11-17 | `NOTICE` holds the copyright, the licence summary, the earlier Apache-2.0 releases and the third-party data credits; an SPDX header is in all 180 source files; `NOTICE` ships in the Python and npm packages |
| 5.3 | Name and logo | Trademark policy in the repo | partial 2026-11-17 | `LICENSING.md` states that the licence covers the code and not the name: others may say "based on ChemLarp", but may not name a modified version ChemLarp or imply endorsement, and forks must rename. Registering the mark is a separate, optional step for the owner |
| 5.4 | Contributions | Contributor licence agreement or Developer Certificate of Origin | done 2026-11-17 | `CONTRIBUTING.md` requires a `Signed-off-by` line and grants the copyright holder the right to relicense, which is what keeps commercial licensing possible |
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
