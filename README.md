# OmniChem Core

OmniChem Core is a terminal-native chemistry engine written in C++20. It covers the path from SMILES to molecular graph to thermochemistry to Gibbs-minimized equilibrium. Every stage shares one zero-copy structure-of-arrays state, and every error is collected into a single vector.

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
| CLI | `build\Release\omnichem.exe` (VS) / `build\omnichem.exe` (Ninja, MinGW) | `build/omnichem` | `build/omnichem` |
| C ABI shared library | `omnichem.dll` + `omnichem.lib` | `libomnichem.dylib` | `libomnichem.so` |

The Python binding finds the library automatically in any of these locations. Set `OMNICHEM_LIB` to point it elsewhere.

**Quoting SMILES on the command line:** use single quotes in bash, zsh and PowerShell (`'[C-]#[O+]'`), and double quotes in `cmd.exe`. In zsh, unquoted `[...]` is treated as a glob.

CI (`.github/workflows/ci.yml`) builds and tests on Windows (MSVC and MinGW), macOS and Linux, including CLI and Python smoke tests.

## Try it

```bash
./build/omnichem parse 'c1ccncc1'
```

```bash
./build/omnichem parse 'C(C)(C)(C)(C)C.c1ccnc1'      # shows the error vector
```

```bash
./build/omnichem canon 'C1=CC=C(C=C1)/C=C/C(=O)O' --random 3   # canonical + equivalent random forms
```

```bash
./build/omnichem kinetics data/mechanisms/h2o2.ocm -m adiabatic-v -T 1000 -P 1 -t 1e-3 H2=2 O2=1 AR=7
```

```bash
./build/omnichem reaction '1 C + 2 O=O -> 1 O=C=O + 2 O' -T 298.15
```

```bash
./build/omnichem equilibrate -T 2500 C 1 O=O 2 O=C=O 0 O 0 '[C-]#[O+]' 0 '[H][H]' 0 '[OH]' 0 '[H]' 0 '[O]' 0
```

```bash
./build/omnichem parse CCO --json | scripts/ocx_pack.py pack ethanol.ocx -
```

Python (no compiled extension, just ctypes):

```python
import sys; sys.path.insert(0, "python")
import omnichem
omnichem.equilibrate({"[C-]#[O+]": 1, "O": 1, "O=C=O": 0, "[H][H]": 0}, temperature_K=1000)
with omnichem.System() as s:
    s.add("c1ccccc1O")
    xyz = s.atom_x()        # numpy view straight into the C++ SoA column
omnichem.kinetics("data/mechanisms/h2o2.ocm", {"H2": 2, "O2": 1, "AR": 7}, 1000, 1e5, t_end=1e-3)["ignition_time"]
omnichem.canonical("N[C@@H](C)C(=O)O")   # {'canonical': 'C[C@H](N)C(=O)O', 'cip': {... 'S' ...}}
```

## Layout

```
include/omnichem/   public headers (C++ API + omnichem_c.h C ABI)
  system_state.hpp    ChemicalSystemState, SoA Atom/Bond tables, CSR adjacency
  diagnostics.hpp     error-vector codes (1xx syntax, 2xx valence, 3xx stoichiometry, 4xx thermo)
  smiles.hpp          OpenSMILES parser, kekulizer, valence validator, stereo references
  rings.hpp           SSSR and Hückel aromaticity perception
  canonical.hpp       canonical and randomised SMILES writer
  cip.hpp             CIP priorities, R/S, E/Z
  reaction.hpp        equation parser, element/charge balance, ΔrG°, K
  thermo.hpp          NASA-7 polynomials, built-in GRI-Mech data, thermo binding
  joback.hpp          group-additivity ΔfH°/ΔfG° from the graph
  gibbs.hpp           element-potential Gibbs minimiser
  kinetics.hpp        .ocm mechanisms, mass-action rates, batch reactors
  ode.hpp             variable-order BDF stiff integrator, LU
src/                implementations; src/cli/main.cpp is the `omnichem` terminal tool
tests/              dependency-free test suite (55 tests: random-SMILES round trips, Robertson, ignition invariants)
data/mechanisms/    example mechanisms (h2o2.ocm: illustrative H2/O2 parameters)
python/omnichem/    ctypes binding
scripts/            build.sh, ocx_pack.py
docs/               ARCHITECTURE.md (design, provenance, status), MATH.md (SMILES→Gibbs math), OCX_SPEC.md
```

Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for what is built versus planned, and [docs/MATH.md](docs/MATH.md) for the full SMILES → Gibbs matrix derivation.
