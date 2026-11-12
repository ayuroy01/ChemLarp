#!/usr/bin/env bash
# Build the OmniChem web app: the C++ engine compiled to WebAssembly plus the
# static UI in web/. Output: build/web (serve it with any static file server).
# Requires the Emscripten SDK (em++ on PATH).
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=build/web
rm -rf "$OUT"
mkdir -p "$OUT/data"
SOURCES=$(ls src/*.cpp)
EXPORTS='["_oc_flash_text_json","_oc_saturation_pressure_json","_oc_version","_oc_parse_smiles_json","_oc_canonical_smiles32_json","_oc_check_reaction_json","_oc_equilibrate_text_json","_oc_kinetics_text_json","_oc_embed32_json","_oc_depict_json","_oc_inchi_json","_oc_read_molecules_json","_oc_smarts_match_json","_oc_similarity_json","_oc_search_json","_oc_free_string","_malloc","_free"]'
em++ -std=c++20 -O3 -fexceptions -Iinclude -DOMNICHEM_STATIC $SOURCES -o "$OUT/omnichem.js" \
  -fexceptions -sMODULARIZE=1 -sEXPORT_NAME=createOmniChem -sENVIRONMENT=web,worker \
  -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=4194304 \
  -sEXPORTED_FUNCTIONS="$EXPORTS" -sEXPORTED_RUNTIME_METHODS='["ccall","UTF8ToString"]'
cp web/*.html web/*.js web/*.css web/manifest.webmanifest "$OUT/"
cp -r web/icons "$OUT/icons"
cp -r data/mechanisms "$OUT/data/mechanisms"
cp data/roadmap.json "$OUT/data/roadmap.json"
BUILD_ID="${GITHUB_SHA:-local}-$(date -u +%Y%m%d%H%M%S)"
sed -i.bak "s/__BUILD_ID__/${BUILD_ID}/g" "$OUT/sw.js" "$OUT/app.js" && rm -f "$OUT"/*.bak
echo "web app built in $OUT (build $BUILD_ID)"
