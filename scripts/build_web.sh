#!/usr/bin/env bash
# Build the ChemLarp web app: the C++ engine compiled to WebAssembly plus the
# static UI in web/. Output: build/web (serve it with any static file server).
# Requires the Emscripten SDK (em++ on PATH).
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=build/web
rm -rf "$OUT"
mkdir -p "$OUT/data"
SOURCES=$(ls src/*.cpp)
EXPORTS='["_cl_flash_text_json","_cl_saturation_pressure_json","_cl_vle_text_json","_cl_flowsheet_json","_cl_aqueous_json","_cl_md_json","_cl_lj_json","_cl_mmff_json","_cl_pure_component_json","_cl_version","_cl_parse_smiles_json","_cl_canonical_smiles32_json","_cl_check_reaction_json","_cl_equilibrate_text_json","_cl_kinetics_text_json","_cl_embed32_json","_cl_depict_json","_cl_inchi_json","_cl_iupac_name_json","_cl_read_molecules_json","_cl_smarts_match_json","_cl_similarity_json","_cl_search_json","_cl_qc_json","_cl_free_string","_malloc","_free"]'
em++ -std=c++20 -O3 -fexceptions -Iinclude -DCHEMLARP_STATIC $SOURCES -o "$OUT/chemlarp.js" \
  -fexceptions -sMODULARIZE=1 -sEXPORT_NAME=createChemLarp -sENVIRONMENT=web,worker \
  -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=4194304 \
  -sEXPORTED_FUNCTIONS="$EXPORTS" -sEXPORTED_RUNTIME_METHODS='["ccall","UTF8ToString"]'
cp web/*.html web/*.js web/*.css web/manifest.webmanifest "$OUT/"
cp -r web/icons "$OUT/icons"
cp -r data/mechanisms "$OUT/data/mechanisms"
cp data/roadmap.json "$OUT/data/roadmap.json"
BUILD_ID="${GITHUB_SHA:-local}-$(date -u +%Y%m%d%H%M%S)"
sed -i.bak "s/__BUILD_ID__/${BUILD_ID}/g" "$OUT/sw.js" "$OUT/app.js" && rm -f "$OUT"/*.bak
echo "web app built in $OUT (build $BUILD_ID)"
