// IUPAC names (2013 Recommendations). Every expected name was checked by parsing it back with
// OPSIN 2.8 (py2opsin) and comparing the InChIKey with the input structure (RDKit); the corpus
// harness does the same on 7,000 PubChem structures (see ARCHITECTURE.md).
#include <string>
#include <utility>
#include <vector>

#include "chemlarp/iupac.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {
std::string name_of(const char* smiles) {
  Diagnostics d;
  const NameResult r = iupac_name_from_smiles(smiles, d);
  return r.ok ? r.name : "!" + (d.items().empty() ? std::string() : d.items().back().message);
}
}  // namespace

TEST(iupac_names_round_trip_through_opsin) {
  const std::vector<std::pair<const char*, const char*>> cases = {
      {R"(CCO)", R"(ethanol)"},
      {R"(CCCO)", R"(propan-1-ol)"},
      {R"(CC(C)C)", R"(2-methylpropane)"},
      {R"(CC(=O)O)", R"(acetic acid)"},
      {R"(C(=O)O)", R"(formic acid)"},
      {R"(CCOC(C)=O)", R"(ethyl acetate)"},
      {R"(c1ccccc1O)", R"(phenol)"},
      {R"(Nc1ccccc1)", R"(aniline)"},
      {R"(CC(=O)c1ccccc1)", R"(1-phenylethan-1-one)"},
      {R"(CN(C)C=O)", R"(N,N-dimethylformamide)"},
      {R"(C[C@H](N)C(=O)O)", R"((2S)-2-aminopropanoic acid)"},
      {R"(OC(=O)CC(O)(CC(=O)O)C(=O)O)", R"(2-hydroxypropane-1,2,3-tricarboxylic acid)"},
      {R"(CC(C)Cc1ccc(cc1)C(C)C(=O)O)", R"(2-[4-(2-methylpropyl)phenyl]propanoic acid)"},
      {R"(CC(=O)Nc1ccc(O)cc1)", R"(N-(4-hydroxyphenyl)acetamide)"},
      {R"(C/C=C/C)", R"((2E)-but-2-ene)"},
      {R"(ClC(Br)F)", R"(bromo(chloro)(fluoro)methane)"},
      {R"(CCN(CC)CC)", R"(N,N-diethylethanamine)"},
      {R"(O=C1CCCN1C)", R"(1-methylpyrrolidin-2-one)"},
      {R"(CC(=O)OCCOC(C)=O)", R"(2-(acetyloxy)ethyl acetate)"},
      {R"(c1ccc2[nH]ccc2c1)", R"(1H-indole)"},
      {R"(O=c1cc[nH]c(=O)[nH]1)", R"(1,2,3,4-tetrahydropyrimidine-2,4-dione)"},
      {R"(C1CCc2ccccc2C1)", R"(1,2,3,4-tetrahydronaphthalene)"},
      {R"(O=c1ccc2ccccc2o1)", R"(2H-1-benzopyran-2-one)"},
      {R"(CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O)", R"(17-hydroxy-10,13-dimethyl-2,3,6,7,8,9,10,11,12,13,14,15,16,17-tetradecahydro-1H-cyclopenta[a]phenanthren-3-one)"},
      {R"(Clc1ccc(Oc2ccccc2)cc1)", R"(1-chloro-4-phenoxybenzene)"},
      {R"(CS(C)=O)", R"((methanesulfinyl)methane)"},
      {R"(Nc1ccc(S(N)(=O)=O)cc1)", R"(4-aminobenzene-1-sulfonamide)"},
      {R"(OCC(O)COP(=O)(O)O)", R"(2,3-dihydroxypropyl dihydrogen phosphate)"},
      {R"(CP(=O)(O)O)", R"(methylphosphonic acid)"},
      {R"(CNC(=O)Nc1ccccc1)", R"(N-methyl-N'-phenylurea)"},
      {R"(CCOC(=O)NC)", R"(ethyl N-methylcarbamate)"},
      {R"(COC(=O)OC)", R"(dimethyl carbonate)"},
      {R"(CC(C)=NO)", R"(N-hydroxypropan-2-imine)"},
      {R"(C1CC2CCC1C2)", R"(bicyclo[2.2.1]heptane)"},
      {R"(CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O)", R"((2S,5R,6R)-3,3-dimethyl-7-oxo-6-(2-phenylacetamido)-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid)"},
      {R"(CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2)", R"(8-methyl-8-azabicyclo[3.2.1]octan-3-yl 3-hydroxy-2-phenylpropanoate)"},
      {R"(C=CC#N)", R"(prop-2-enenitrile)"},
      {R"(CC(C)(C)c1ccc(O)cc1)", R"(4-(2-methylpropan-2-yl)phenol)"},
      {R"(O=C(O)c1ccncc1)", R"(pyridine-4-carboxylic acid)"},
      {R"(OC1CCCCC1)", R"(cyclohexanol)"},
      {R"(CC1=CC(=O)CC(C)(C)C1)", R"(3,5,5-trimethylcyclohex-2-en-1-one)"},
      {R"(C1CCOC1)", R"(oxolane)"},
      {R"(CCOC(=O)CCC(=O)OC)", R"(1-ethyl 4-methyl butanedioate)"},
      {R"(C/C(=N/O)c1ccccc1)", R"((1Z)-N-hydroxy-1-phenylethan-1-imine)"},
      {R"(C1CN2CCC1CC2)", R"(1-azabicyclo[2.2.2]octane)"},
      {R"(CC(C)NCC(O)COc1cccc2ccccc12)", R"(1-[(naphthalen-1-yl)oxy]-3-[(propan-2-yl)amino]propan-2-ol)"},
      {R"(CN1CCC[C@H]1c1cccnc1)", R"(3-[(2S)-1-methylpyrrolidin-2-yl]pyridine)"},
      {R"(NC(=O)c1cccnc1)", R"(pyridine-3-carboxamide)"},
      {R"(O=C1NC(=O)C(N1)c1ccccc1)", R"(5-phenylimidazolidine-2,4-dione)"},
      {R"(CC(=O)OC1CC2CCC1(C)C2(C)C)", R"(1,7,7-trimethylbicyclo[2.2.1]heptan-2-yl acetate)"},
      {R"(CC(C)C[C@H](N)C(=O)O)", R"((2S)-2-amino-4-methylpentanoic acid)"},
      {R"(OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O)", R"((2R,3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetraol)"},
      {R"(Cn1cnc2c1c(=O)n(C)c(=O)n2C)", R"(1,3,7-trimethyl-2,3,6,7-tetrahydro-1H-purine-2,6-dione)"},
  };
  for (const auto& [smiles, expected] : cases) CHECK_EQ(name_of(smiles), std::string(expected));
}

TEST(iupac_name_is_independent_of_input_order) {
  // the same structure written three ways
  CHECK_EQ(name_of("OC(=O)c1ccccc1OC(C)=O"), name_of("CC(=O)Oc1ccccc1C(=O)O"));
  CHECK_EQ(name_of("C(C(=O)O)(C)N"), name_of("CC(N)C(=O)O"));
  CHECK_EQ(name_of("c1ccc2c(c1)CCCC2"), name_of("C1CCc2ccccc2C1"));
}

TEST(iupac_unsupported_structures_are_reported_not_guessed) {
  for (const char* smiles : {"CC[Si](C)C", "[Na+].CC(=O)[O-]", "CC(=O)Cl", "c1ccc2cc3ccccc3cc2c1C1CC12CC2", "NN", "CC(=S)N",
                              "c1ccc2nc3ccccc3cc2c1.O", "C1CC2(C1)CCC2"}) {
    Diagnostics d;
    const NameResult r = iupac_name_from_smiles(smiles, d);
    CHECK(!r.ok);
    CHECK(d.contains(DiagCode::NamingUnsupported));
  }
}
