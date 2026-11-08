# Build OmniChem Core on Windows and run the test suite.
# Requires CMake 3.20+ and a C++20 compiler (Visual Studio 2022 "Desktop development with C++",
# or clang-cl / MinGW-w64). Run from any directory:
#   powershell -ExecutionPolicy Bypass -File scripts\build.ps1            # Release
#   powershell -ExecutionPolicy Bypass -File scripts\build.ps1 -Config Debug
param(
  [ValidateSet("Release", "Debug", "RelWithDebInfo")]
  [string]$Config = "Release",
  [string]$Generator = ""
)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)

$configureArgs = @("-S", ".", "-B", "build")
if ($Generator -ne "") { $configureArgs += @("-G", $Generator) }
cmake @configureArgs
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }

cmake --build build --config $Config --parallel
if ($LASTEXITCODE -ne 0) { throw "Build failed" }

ctest --test-dir build -C $Config --output-on-failure
if ($LASTEXITCODE -ne 0) { throw "Tests failed" }

Write-Host "`nBuilt: build\omnichem.exe, build\omnichem.dll (Python: set OMNICHEM_LIB if you move it)"
