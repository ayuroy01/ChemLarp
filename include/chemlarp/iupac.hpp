#pragma once
// Systematic IUPAC names (structure -> name), after the 2013 Recommendations
// (Favre & Powell, Nomenclature of Organic Chemistry, RSC 2014; rule numbers P-*).
//
// Substitutive nomenclature for organic compounds of C, H, N, O, S, P and the
// halogens:
//  - parents: acyclic chains, monocyclic carbocycles, benzene, retained and
//    Hantzsch-Widman heteromonocycles and fused systems with fixed numbering
//    (pyridine ... purine, naphthalene, indole, quinoline, carbazole,
//    cyclopenta[a]phenanthrene, ...) whose saturation is expressed by indicated
//    hydrogen and hydro prefixes; ortho- and peri-fused systems of any number
//    of rings by fusion nomenclature (FR-3 to FR-5); von Baeyer names for
//    bridged polycycles; Hantzsch-Widman heteromonocycles; bicyclic
//    von Baeyer and monospiro systems (bicyclo[2.2.1]heptane, spiro[4.5]decane);
//  - principal characteristic groups as suffixes in the order of P-41
//    (carboxylic and sulfonic acids, esters, amides and sulfonamides, nitriles,
//    aldehydes, ketones, alcohols, thiols, amines, imines), the rest as
//    prefixes; thioamides, hydrazides, thiones; functional class names for
//    esters, carbonates, carbamates, phosphate monoesters and phosphonic acids;
//    ureas, thioureas and guanidines; amidines; anions, cations, zwitterions, and salts, hydrogen
//    halide adducts and hydrates named component by component;
//  - parent selection by P-44 (rings before chains, then longest chain, then
//    unsaturation), numbering by the lowest-locant rules of P-31.1.4,
//    alphanumerical order and multiplying prefixes (P-14.5, P-16.3), enclosing
//    marks nested ( [ { (P-16.5);
//  - retained names that are preferred IUPAC names where they apply (acetic
//    acid, benzoic acid, phenol, aniline, benzaldehyde, benzonitrile, ...);
//  - CIP stereodescriptors (R/S, E/Z) on parent and substituent atoms.
//
// A structure outside this scope (spiro systems with fused components,
// uncommon ions, metal complexes, hydrazines, polyphosphates, ...) gets
// no name and a NamingUnsupported diagnostic saying why: the namer never guesses.

#include <cstdint>
#include <string>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

struct NameResult {
  bool ok = false;
  std::string name;
};

NameResult iupac_name(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);
NameResult iupac_name_from_smiles(const std::string& smiles, Diagnostics& diag);

}  // namespace chemlarp
