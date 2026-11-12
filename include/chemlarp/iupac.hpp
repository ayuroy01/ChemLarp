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
//    hydrogen and hydro prefixes, and bicyclic von Baeyer systems
//    (bicyclo[2.2.1]heptane, 8-azabicyclo[3.2.1]octane, ...);
//  - principal characteristic groups as suffixes in the order of P-41
//    (carboxylic and sulfonic acids, esters, amides and sulfonamides, nitriles,
//    aldehydes, ketones, alcohols, thiols, amines, imines), the rest as
//    prefixes; functional class names for esters, carbonates, carbamates,
//    phosphate monoesters and phosphonic acids; ureas;
//  - parent selection by P-44 (rings before chains, then longest chain, then
//    unsaturation), numbering by the lowest-locant rules of P-31.1.4,
//    alphanumerical order and multiplying prefixes (P-14.5, P-16.3), enclosing
//    marks nested ( [ { (P-16.5);
//  - retained names that are preferred IUPAC names where they apply (acetic
//    acid, benzoic acid, phenol, aniline, benzaldehyde, benzonitrile, ...);
//  - CIP stereodescriptors (R/S, E/Z) on parent and substituent atoms.
//
// A structure outside this scope (ions and salts, fused ring systems without a
// built-in parent, spiro and tricyclic systems, amidines, hydrazines, ...) gets
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
