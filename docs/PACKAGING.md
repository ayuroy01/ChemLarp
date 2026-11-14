# Packages

ChemLarp is packaged five ways. Every one of them is built and installed by the `packages` workflow
(`.github/workflows/packages.yml`) on each push, so a release attaches files that were already tested.
Nothing is published automatically: publishing needs accounts that belong to the project owner, and the
steps are at the end of this file.

| Package | What it holds | Built by | Tested by |
|---|---|---|---|
| Python wheel `chemlarp` | the shared library (loaded with ctypes), the Python module, and the `chemlarp` command | `setup.py` (make on macOS/Linux, CMake on Windows) via cibuildwheel | `scripts/wheel_smoke_test.py` in a clean environment |
| npm package `chemlarp` | the engine compiled to WebAssembly for Node, with a JavaScript wrapper | `scripts/build_npm.sh` (Emscripten) | `release/npm/test.js`, then an install of the packed tarball |
| Homebrew formula | the CLI, shared library, headers and data under the Homebrew prefix | `make install` | the formula's `test do` block (`brew test`) |
| winget manifests | the portable `chemlarp.exe` from the release zip | `release/winget/*.yaml`, filled by `scripts/release_manifests.py` | YAML parse; installation needs a published release |
| Command-line archives | `chemlarp` / `chemlarp.exe` with the licence, for Linux x64, macOS arm64, Windows x64 | make; CMake with the static MSVC runtime on Windows | runs `chemlarp canon` on the build machine |

The version is written in four places (CMakeLists.txt, `cl_version()`, pyproject.toml, the npm package.json);
`scripts/check_packages.py` fails if they differ.

## Building locally

```bash
# Python wheel (no download: uses the setuptools already installed)
python3 -m pip wheel . --no-build-isolation --no-deps -w dist
# x86_64 wheel on an Apple-silicon Mac
ARCHFLAGS="-arch x86_64" python3 -m pip wheel . --no-build-isolation --no-deps -w dist

# Install under a prefix (what the Homebrew formula runs)
make install PREFIX=$HOME/.local

# npm package (needs the Emscripten SDK)
scripts/build_npm.sh
```

The wheel is tagged `py3-none-<platform>`: the binding uses ctypes, not a compiled extension, so one wheel per
operating system and architecture serves every Python from 3.10. macOS wheels target macOS 11; Linux wheels
link libstdc++ statically and are repaired to `manylinux_2_28` by auditwheel in CI.

## Releasing

1. Set the new version in the four places above and push; check that the `packages` workflow passes.
2. Tag and push the tag: `git tag v0.1.0 && git push origin v0.1.0`. The workflow then makes a **draft** GitHub
   release with the wheels, the npm tarball, the command-line archives, `manifests.tar.gz` (the Homebrew
   formula and winget manifests with their checksums filled in) and `SHA256SUMS`.
3. Look the draft over and publish it on GitHub.

Then, with the project owner's accounts:

- **PyPI**: `python3 -m twine upload dist/*.whl` (or set up a trusted publisher for this repository on
  pypi.org and add an upload job to the workflow).
- **npm**: `npm publish chemlarp-<version>.tgz`.
- **Homebrew**: create the tap repository `ayuroy01/homebrew-chemlarp` and copy `homebrew/chemlarp.rb` from
  `manifests.tar.gz` to `Formula/chemlarp.rb` there. Users then run `brew install ayuroy01/chemlarp/chemlarp`.
  Before the first release, `brew install --HEAD` builds from the main branch.
- **winget**: open a pull request to `microsoft/winget-pkgs` with the three files under
  `winget/manifests/a/ayuroy01/ChemLarp/<version>/` from `manifests.tar.gz` (`wingetcreate submit` can do it).

## Not covered

No code signing: the macOS command-line archive is not notarized (Gatekeeper asks before running a
downloaded copy; the wheel and Homebrew builds are not affected), and the Windows executable is not
Authenticode-signed (SmartScreen may warn). No Windows on Arm, musl Linux or 32-bit wheels, no conda-forge
recipe, and the npm package exposes 16 of the C API's functions (those of the web app's worker plus MBAR).
