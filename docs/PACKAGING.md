# Packages

ChemLarp is packaged eight ways. Every one of them is built and installed by the `packages` workflow
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
| macOS app `ChemLarp-<v>-macos.dmg` | ChemLarp.app (universal arm64 + x86_64, macOS 12+): the web app in a WKWebView, served by a server inside the app on 127.0.0.1 | `scripts/build_macos_app.sh` (swiftc, `release/macos/main.swift`) | `--selftest` (server, MIME types, path traversal refused) and `--smoke` (an engine call through the app's worker, in its window); the disk image is mounted and checked |
| Windows installer `ChemLarp-<v>-windows-x64-setup.exe` | ChemLarp.exe (WebView2 window over the bundled web app, mapped to https://chemlarp.local) and `bin\chemlarp.exe`; per-user install, Start menu entry, uninstaller | `release/windows/main.cpp` with cl.exe, Inno Setup (`release/windows/chemlarp.iss`) | `--smoke` on the built program, then a silent install and `--smoke` on the installed one |
| Android app `ChemLarp-<v>-android.apk` | the web app in a WebView, served by an in-app server on 127.0.0.1 (Android 8+) | `scripts/build_android.sh` (Gradle, `release/android`) | installed on an emulator (API 34) and run with the `smoke` extra; the result is read from logcat |

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

# Native apps, from a built web app (scripts/build_web.sh output, or the `web-app` artifact of the workflow)
scripts/build_macos_app.sh build/web            # ChemLarp.app and the .dmg (Xcode command-line tools)
scripts/build_android.sh build/web              # the .apk (JDK 17, Gradle 8.9+, Android SDK 35)
```

The Windows app needs Visual Studio's C++ tools and the WebView2 SDK (NuGet package Microsoft.Web.WebView2);
the commands are in the `installer-windows` job of the workflow.

All three apps show the same web app as https://ayuroy01.github.io/ChemLarp/, built from the same commit, and
work offline. They are thin shells: the engine, the interface and every calculation are the WebAssembly
build. They do not update themselves; a newer version is a new download.

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

No code signing: the macOS app and command-line archive are signed ad hoc only and not notarized (Gatekeeper
asks before the first launch of a downloaded copy: right-click > Open; the wheel and Homebrew builds are not
affected), the Windows installer and executables are not Authenticode-signed (SmartScreen may warn), and the
Android APK is signed with a throwaway key made in CI, so each build is a different signer: an APK from a
later build will not install over an earlier one without uninstalling first, and Google Play needs the owner's
key and an app bundle (.aab) instead. On Android, exporting a drawing or file (SVG, PNG, MOL) does not save
anything yet (the WebView does not handle the page's blob downloads). The Windows app needs the WebView2
Runtime (part of Windows 11 and current Windows 10; the installer does not bundle it). No iOS app (needs an
Apple developer account and the App Store), no Linux desktop package (the installable web app covers Linux),
no Windows on Arm build. No Windows on Arm, musl Linux or 32-bit wheels, no conda-forge
recipe, and the npm package exposes 16 of the C API's functions (those of the web app's worker plus MBAR).
