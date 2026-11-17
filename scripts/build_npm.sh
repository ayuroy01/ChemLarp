#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
# Node.js package: the C++ engine compiled to WebAssembly for Node, with the wrapper in release/npm.
# Output: build/npm and build/npm/chemlarp-<version>.tgz (npm pack); the package's smoke test is run first.
# Requires the Emscripten SDK (em++ on PATH) and Node 18+.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=build/npm
rm -rf "$OUT"
mkdir -p "$OUT"
EXPORTS='["_cl_version","_cl_parse_smiles_json","_cl_canonical_smiles32_json","_cl_inchi_json","_cl_iupac_name_json","_cl_depict_json","_cl_embed32_json","_cl_read_molecules_json","_cl_smarts_match_json","_cl_similarity_json","_cl_check_reaction_json","_cl_equilibrate_text_json","_cl_qc_json","_cl_mmff_json","_cl_pure_component_json","_cl_mbar_json","_cl_free_string","_malloc","_free"]'
em++ -std=c++20 -O3 -fexceptions -Iinclude -DCHEMLARP_STATIC src/*.cpp -o "$OUT/chemlarp.js" \
  -sMODULARIZE=1 -sEXPORT_NAME=createChemLarp -sENVIRONMENT=node \
  -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=4194304 \
  -sEXPORTED_FUNCTIONS="$EXPORTS" -sEXPORTED_RUNTIME_METHODS='["ccall","UTF8ToString"]'
cp release/npm/package.json release/npm/index.js release/npm/index.d.ts release/npm/test.js \
   release/npm/README.md LICENSE NOTICE "$OUT/"
(cd "$OUT" && node test.js && npm pack)
echo "npm package built in $OUT"
