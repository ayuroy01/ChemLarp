# ChemLarp

ChemLarp is a terminal-native chemistry engine written in C++20. It covers the path from SMILES to molecular graph to thermochemistry to Gibbs-minimized equilibrium. Every stage shares one zero-copy structure-of-arrays state, and every error is collected into a single vector.

## Use it in your browser (any OS, installable)

**https://ayuroy01.github.io/ChemLarp/**: the full engine compiled to WebAssembly with a point-and-click interface (molecule analysis, opening MOL/SDF/PDB/CML files, substructure and similarity search, 2D structure drawings, 3D models, reactions, equilibrium, kinetics, progress). It runs entirely on your device and works offline after the first visit.

- **Windows / macOS / Linux (Chrome, Edge):** click **Install app** in the page header or the install icon in the address bar.
- **Android (Chrome):** menu ⋮ → **Install app** / **Add to Home screen**.
- **iPhone / iPad (Safari):** Share → **Add to Home Screen**.

Build it yourself with the Emscripten SDK: `scripts/build_web.sh`, then serve `build/web` with any static file server.

## Build

Supported platforms: **Windows** (Visual Studio 2022 / MSVC, clang-cl, or MSYS2 MinGW-w64), **macOS** (Apple clang, Xcode 15+), and **Linux** (GCC 11+ / Clang 15+). The only requirements are a C++20 compiler and CMake 3.20+. There are no other dependencies.

**macOS / Linux**

```bash
scripts/build.sh            # CMake if installed, otherwise the Makefile; runs the tests
```

**Windows (PowerShell, from a "Developer PowerShell for VS 2022" or with CMake on PATH)**

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1
```

**Any platform, by hand**

```bash
cmake -S . -B build
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

| Output | Windows | macOS | Linux |
|---|---|---|---|
| CLI | `build\Release\chemlarp.exe` (VS) / `build\chemlarp.exe` (Ninja, MinGW) | `build/chemlarp` | `build/chemlarp` |
| C ABI shared library | `chemlarp.dll` + `chemlarp.lib` | `libchemlarp.dylib` | `libchemlarp.so` |

The Python binding finds the library automatically in any of these locations. Set `CHEMLARP_LIB` to point it elsewhere.

**Quoting SMILES on the command line:** use single quotes in bash, zsh and PowerShell (`'[C-]#[O+]'`), and double quotes in `cmd.exe`. In zsh, unquoted `[...]` is treated as a glob.

CI (`.github/workflows/ci.yml`) builds and tests on Windows (MSVC and MinGW), macOS and Linux, including CLI and Python smoke tests.

## Try it

```bash
./build/chemlarp parse 'c1ccncc1'
```

```bash
./build/chemlarp parse 'C(C)(C)(C)(C)C.c1ccnc1'      # shows the error vector
```

```bash
./build/chemlarp canon 'C1=CC=C(C=C1)/C=C/C(=O)O' --random 3   # canonical + equivalent random forms
```

```bash
./build/chemlarp kinetics data/mechanisms/h2o2.ocm -m adiabatic-v -T 1000 -P 1 -t 1e-3 H2=2 O2=1 AR=7
```

The same mechanism in CHEMKIN-II and Cantera YAML form (detected from content):

```bash
./build/chemlarp kinetics data/mechanisms/h2o2/chem.inp --thermo data/mechanisms/h2o2/therm.dat --smiles data/mechanisms/h2o2/species.smi -T 1000 -t 1e-3 H2=2 O2=1 AR=7
```

```bash
./build/chemlarp kinetics data/mechanisms/h2o2.yaml -T 1000 -t 1e-3 H2=2 O2=1 AR=7
```

```bash
./build/chemlarp inchi 'C[C@@H](C(=O)O)N'                                              # InChI and InChIKey
./build/chemlarp read compounds.sdf --inchi             # MOL/SDF, PDB, CML -> canonical SMILES, title, InChIKey
./build/chemlarp match 'CC(=O)Oc1ccccc1C(=O)O' '[CX3](=O)[OX2H1]'                     # SMARTS matches
./build/chemlarp search compounds.sdf --smarts 'c1ccccc1C(=O)O' --like 'OC(=O)c1ccccc1O' --min 0.3   # screen a file
./build/chemlarp depict 'C[C@@H](C(=O)O)N' -o alanine.svg --mol alanine.mol --cip   # 2D drawing with wedges
```

```bash
./build/chemlarp embed 'C[C@@H](C(=O)O)N' --sdf alanine.sdf    # 3D: distance geometry + UFF, stereo checked
```

```bash
./build/chemlarp reaction '1 C + 2 O=O -> 1 O=C=O + 2 O' -T 298.15
```

```bash
./build/chemlarp equilibrate -T 2500 C 1 O=O 2 O=C=O 0 O 0 '[C-]#[O+]' 0 '[H][H]' 0 '[OH]' 0 '[H]' 0 '[O]' 0
```

```bash
./build/chemlarp parse CCO --json | scripts/clx_pack.py pack ethanol.clx -
```

Python (no compiled extension, just ctypes):

```python
import sys; sys.path.insert(0, "python")
import chemlarp
chemlarp.equilibrate({"[C-]#[O+]": 1, "O": 1, "O=C=O": 0, "[H][H]": 0}, temperature_K=1000)
with chemlarp.System() as s:
    s.add("c1ccccc1O")
    xyz = s.atom_x()        # numpy view straight into the C++ SoA column
chemlarp.kinetics("data/mechanisms/h2o2.ocm", {"H2": 2, "O2": 1, "AR": 7}, 1000, 1e5, t_end=1e-3)["ignition_time"]
chemlarp.inchi("N[C@@H](C)C(=O)O")["key"]   # 'QNAYBMKLOCPYGJ-REOHCLBHSA-N'
chemlarp.read_molecules("compounds.sdf")     # [{'title', 'smiles', 'formula', 'properties': {...}}, ...]
chemlarp.smarts_match("CC(=O)Oc1ccccc1C(=O)O", "[CX3](=O)[OX2H1]")   # [(10, 11, 12)]
chemlarp.similarity("CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O")["tanimoto"]   # ECFP4 Tanimoto
chemlarp.depict("N[C@@H](C)C(=O)O")["svg"]  # 2D drawing (SVG + molfile with wedges)
chemlarp.embed("N[C@@H](C)C(=O)O")["sdf"]   # 3D molfile with explicit hydrogens
chemlarp.canonical("N[C@@H](C)C(=O)O")   # {'canonical': 'C[C@H](N)C(=O)O', 'cip': {... 'S' ...}}
```

## Layout

```
include/chemlarp/   public headers (C++ API + chemlarp_c.h C ABI)
  system_state.hpp    ChemicalSystemState, SoA Atom/Bond tables, CSR adjacency
  diagnostics.hpp     error-vector codes (1xx syntax, 2xx valence, 3xx stoichiometry, 4xx thermo)
  smiles.hpp          OpenSMILES parser, kekulizer, valence validator, stereo references
  rings.hpp           SSSR and Hückel aromaticity perception
  canonical.hpp       canonical and randomised SMILES writer
  cip.hpp             CIP priorities (rules 1a-5), R/S, r/s, E/Z
  inchi.hpp           standard InChI and InChIKey
  molio.hpp           readers for MOL V2000/V3000, SDF, PDB, CML
  smarts.hpp          SMARTS parser and substructure matcher
  fingerprint.hpp     Morgan (ECFP) fingerprints, Tanimoto and Dice
  reaction.hpp        equation parser, element/charge balance, ΔrG°, K
  thermo.hpp          NASA-7 polynomials, built-in GRI-Mech data, thermo binding
  joback.hpp          group-additivity ΔfH°/ΔfG° from the graph
  gibbs.hpp           element-potential Gibbs minimiser
  kinetics.hpp        mechanisms (.ocm / CHEMKIN-II / Cantera YAML), fall-off, PLOG, batch reactors
  ode.hpp             variable-order BDF stiff integrator, LU
  depict.hpp          2D layout, wedges, SVG drawing, 2D molfile
  conformer.hpp       explicit H, distance-geometry 3D embedding, stereo from 3D, SDF/XYZ
  uff.hpp             Universal Force Field (energy, analytic gradient, minimiser)
  optimize.hpp        L-BFGS
src/                implementations; src/cli/main.cpp is the `chemlarp` terminal tool
tests/              dependency-free test suite (93 tests: random-SMILES round trips, Robertson, Troe, cross-format kinetics, 2D drawing, 3D/UFF)
data/mechanisms/    example H2/O2 mechanism in .ocm, CHEMKIN-II and YAML form (illustrative parameters; regenerate with scripts/gen_h2o2_mechanism.py)
python/chemlarp/    ctypes binding
scripts/            build.sh, clx_pack.py
docs/               ARCHITECTURE.md (design, provenance, status), MATH.md (SMILES→Gibbs math), CLX_SPEC.md
```

Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for what is built versus planned, and [docs/MATH.md](docs/MATH.md) for the full SMILES → Gibbs matrix derivation.

## License

Apache License 2.0 — see [LICENSE](LICENSE).
