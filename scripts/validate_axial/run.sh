#!/bin/sh
# Axial chirality against CDK and the centres CIP labeller (Hanson, Musgrave, Mayfield et al.).
#
#   scripts/validate_axial/run.sh <cdk.jar> <centres-core.jar> <centres-cdk.jar> [work dir]
#
# Jars: CDK 2.13 (https://github.com/cdk/cdk/releases/tag/cdk-2.13), centres 1.5 from Maven Central
# (com.simolecule.centres:centres-core and centres-cdk). Build ChemLarp first (make).
#
# 1  367 allenes (random substituents, both configurations, cyclic and longer cumulenes): Ra/Sa from
#    ChemLarp's SMILES reading against centres (M = Ra, P = Sa) on CDK's reading of the same SMILES;
#    ChemLarp's canonical and three random spellings read back by CDK as the same stereoisomer.
# 2  The SMILES CDK writes for those allenes, read by ChemLarp: same descriptor.
# 3  124 biaryls (BINOL and random tri- and tetra-ortho-substituted ones) embedded in 3D by ChemLarp, written
#    as SDF; ChemLarp's Ra/Sa from that file against CDK's 3D atropisomer perception labelled by centres.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); ROOT="$HERE/../.."
CP="$1:$2:$3"; W=${4:-$(mktemp -d)}; mkdir -p "$W"
cc() { c++ -std=c++20 -O1 -I"$ROOT/include" "$HERE/$1" "$ROOT/build/libchemlarp.a" -o "$W/$2"; }
cc gen_allenes.cpp gen_allenes; cc gen_atropisomers.cpp gen_atropisomers; cc read_cdk_smiles.cpp read_cdk_smiles
"$W/gen_allenes" > "$W/allenes.tsv"
java -cp "$CP" "$HERE/Allenes.java" "$W/allenes.tsv" 2>/dev/null | grep -v "^RANDOM\|^CANON" || true
java -cp "$CP" "$HERE/WriteCdk.java" "$W/allenes.tsv" 2>/dev/null > "$W/cdk_written.tsv"
"$W/read_cdk_smiles" "$W/cdk_written.tsv" | tail -1
"$W/gen_atropisomers" "$W/atropisomers.sdf"
java -cp "$CP" "$HERE/Atropisomers.java" "$W/atropisomers.sdf" 2>/dev/null | tail -1
