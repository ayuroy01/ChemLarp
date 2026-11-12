// IUPAC names (2013 Recommendations). Every expected name was checked by parsing it back with
// OPSIN 2.8 (py2opsin) and comparing the InChIKey with the input structure (RDKit); the corpus
// harness does the same on 15,000 PubChem structures, including salts and zwitterions (see ARCHITECTURE.md).
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
      {R"(CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21)", R"(7-chloro-1-methyl-5-phenyl-2,3-dihydro-1H-1,4-benzodiazepin-2-one)"},
      {R"(NC(=O)N1c2ccccc2C=Cc2ccccc21)", R"(5H-dibenzo[b,f]azepine-5-carboxamide)"},
      {R"(NC12CC3CC(CC(C3)C1)C2)", R"(adamantan-1-amine)"},
      {R"(CN(C)C(=O)Cc1c(-c2ccc(C)cc2)nc2ccc(C)cn12)", R"(N,N-dimethyl-2-[6-methyl-2-(4-methylphenyl)imidazo[1,2-a]pyridin-3-yl]acetamide)"},
      {R"(O=C1CCCCCCCCCCCCCCO1)", R"(1-oxacyclohexadecan-2-one)"},
      {R"(NC(=S)N)", R"(thiourea)"},
      {R"(CC(N)=S)", R"(ethanethioamide)"},
      {R"(O=C1CSC(=S)N1)", R"(2-sulfanylidene-1,3-thiazolidin-4-one)"},
      {R"(CC(=O)[O-].[Na+])", R"(sodium acetate)"},
      {R"(C[N+](C)(C)CC(=O)[O-])", R"(2-(trimethylazaniumyl)acetate)"},
      {R"(C[n+]1ccccc1.[I-])", R"(1-methylpyridin-1-ium iodide)"},
      {R"(CN(C)CCOC(c1ccccc1)c1ccccc1.Cl)", R"(2-(diphenylmethoxy)-N,N-dimethylethan-1-amine hydrochloride)"},
      {R"(NNC(=O)c1ccncc1)", R"(pyridine-4-carbohydrazide)"},
      {R"(CC(=O)NO)", R"(N-hydroxyacetamide)"},
      {R"(O=C(N/N=C/c1ccccc1)c1ccccc1)", R"(N'-[(E)-phenylmethylidene]benzohydrazide)"},
      {R"(CC(=O)[O-].CC(=O)[O-].[Ca+2])", R"(calcium diacetate)"},
      {R"(CC(=N)N)", R"(ethanimidamide)"},
      {R"(NC(=N)c1ccccc1)", R"(benzenecarboximidamide)"},
      {R"(CC(=NCc1ccccc1)N)", R"(N'-benzylethanimidamide)"},
      {R"(NC(=N)NCCCC(N)C(=O)O)", R"(2-amino-5-carbamimidamidopentanoic acid)"},
      {R"(NC(N)=NCCCC(N)C(=O)O)", R"(2-amino-5-carbamimidamidopentanoic acid)"},
      {R"(CNC(=N)N)", R"(N-methylguanidine)"},
      {R"(CN(CC(=O)O)C(=N)N)", R"(2-[carbamimidoyl(methyl)amino]acetic acid)"},
      {R"(C1SCCCS1)", R"(1,3-dithiane)"},
      {R"(C1CC2(CCN1)CCOC2)", R"(2-oxa-8-azaspiro[4.5]decane)"},
      {R"(O=C1CCC2(CCCCC2)O1)", R"(1-oxaspiro[4.5]decan-2-one)"},
      {R"(CN1CN(C(=S)SC1)C)", R"(3,5-dimethyl-1,3,5-thiadiazinane-2-thione)"},
      {R"(c1cnc2sccc2c1)", R"(thieno[2,3-b]pyridine)"},
      {R"(c1cnc2ccsc2c1)", R"(thieno[3,2-b]pyridine)"},
      {R"(CC1=CN2C=CSC2=N1)", R"(6-methylimidazo[2,1-b][1,3]thiazole)"},
      {R"(c1ccc2c(c1)CCCCC2)", R"(6,7,8,9-tetrahydro-5H-benzo[7]annulene)"},
      {R"(c1ccc2c(c1)COCN2)", R"(1,4-dihydro-2H-3,1-benzoxazine)"},
      {R"(c1ncc2[nH]ccc2n1)", R"(5H-pyrrolo[3,2-d]pyrimidine)"},
      {R"(c1ccc2nonc2c1)", R"(2,1,3-benzoxadiazole)"},
      {R"(CCOC(=O)N1CCN(CC1)C(=S)NC2=CC3=NSN=C3C=C2)", R"(ethyl 4-[(2,1,3-benzothiadiazol-5-yl)carbamothioyl]piperazine-1-carboxylate)"},
      {R"(c1ccc2c(c1)[nH]c1ncccc12)", R"(9H-pyrido[2,3-b]indole)"},
      {R"(c1ccc2c(c1)nc1ccccn12)", R"(pyrido[1,2-a]benzimidazole)"},
      {R"(O=c1ccc2cc3ccoc3cc2o1)", R"(7H-furo[3,2-g][1]benzopyran-7-one)"},
      {R"(c1ccc2c(c1)sc1ncncc12)", R"([1]benzothieno[2,3-d]pyrimidine)"},
      {R"(c1ccc2cc3cc4ccccc4cc3cc2c1)", R"(tetracene)"},
      {R"(c1ccc2c(c1)ccc1c3ccccc3ccc21)", R"(chrysene)"},
      {R"(C1=CC=CC2=C1C=1C(=CC=N2)C2=CC=CC=C2N1)", R"(indolo[3,2-d][1]benzazepine)"},
      {R"(C1=CC=C2C(=C1)OC3=CC=CC=C3C4=C(O2)C=CC=C4)", R"(dibenzo[b,e][1,4]benzodioxocine)"},
      {R"(C1=NN2C=NN=C2SC1)", R"(7H-[1,2,4]triazolo[3,4-b][1,3,4]thiadiazine)"},
      {R"(C1Cc2c(sc3ncnc(N)c23)C1)", R"(6,7-dihydro-5H-cyclopenta[4,5]thieno[2,3-d]pyrimidin-4-amine)"},
      {R"(Clc1ccc2c(c1)C(c1ccccc1)=NCc1nncn1-2)", R"(8-chloro-6-phenyl-4H-[1,2,4]triazolo[4,3-a][1,4]benzodiazepine)"},
      {R"(CN1CCN(CC1)C1=Nc2ccccc2Sc2ccc(C)cc12)", R"(2-methyl-11-(4-methylpiperazin-1-yl)benzo[b][1,4]benzothiazepine)"},
      {R"(C[N+]1(C)CCCCC1)", R"(1,1-dimethylpiperidin-1-ium)"},
      {R"(C[NH+]1CCCC1C(=O)[O-])", R"(1-methylpyrrolidin-1-ium-2-carboxylate)"},
      {R"(C[n+]1ccc(cc1)-c1cc[n+](C)cc1)", R"(1-methyl-4-(1-methylpyridin-1-ium-4-yl)pyridin-1-ium)"},
      {R"(NNc1ccccc1)", R"(phenylhydrazine)"},
      {R"(CC(C)=NNc1ccccc1)", R"(1-phenyl-2-(propan-2-ylidene)hydrazine)"},
      {R"(OC(=O)CC=NNc1ccccc1)", R"(3-(2-phenylhydrazinylidene)propanoic acid)"},
      {R"(CON(C)C(=O)c1ccccc1)", R"(N-methoxy-N-methylbenzamide)"},
      {R"(OC(=O)CCN(C)O)", R"(3-[hydroxy(methyl)amino]propanoic acid)"},
      {R"(O=S1(=O)CCCC1)", R"(1,1-dioxo-1λ6-thiolane)"},
      {R"(O=C1NS(=O)(=O)c2ccccc12)", R"(1,1-dioxo-2,3-dihydro-1λ6,2-benzothiazol-3-one)"},
      {R"(NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O)", R"(6-chloro-1,1-dioxo-3,4-dihydro-2H-1λ6,2,4-benzothiadiazine-7-sulfonamide)"},
      {R"(ClC(=O)c1ccccc1)", R"(benzoyl chloride)"},
      {R"(CS(=O)(=O)Cl)", R"(methanesulfonyl chloride)"},
      {R"(Cc1ccc(S(=O)(=O)OCC)cc1)", R"(ethyl 4-methylbenzene-1-sulfonate)"},
      {R"(CC(=O)SCCN)", R"(S-(2-aminoethyl) ethanethioate)"},
      {R"(COP(C)(=O)OC)", R"(dimethyl methylphosphonate)"},
      {R"(CCOP(=O)(O)OCC)", R"(diethyl hydrogen phosphate)"},
      {R"(CC(=O)NC(C)=O)", R"(N-acetylacetamide)"},
      {R"(NC(=NO)c1ccccc1)", R"(N'-hydroxybenzenecarboximidamide)"},
      {R"(OC(=O)c1ccc(OS(C)(=O)=O)cc1)", R"(4-[(methanesulfonyl)oxy]benzoic acid)"},
      {R"(OC(=O)CCOP(=O)(OC)OC)", R"(3-[(dimethoxyphosphoryl)oxy]propanoic acid)"},
      {R"(CC(=O)OC(C)=O)", R"(acetic anhydride)"},
      {R"(CCC(=O)OC(C)=O)", R"(acetic propanoic anhydride)"},
      {R"(COP(=O)(O)OP(=O)(O)O)", R"(methyl trihydrogen diphosphate)"},
      {R"(Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O)", R"([5-(6-amino-9H-purin-9-yl)-3,4-dihydroxyoxolan-2-yl]methyl tetrahydrogen triphosphate)"},
      {R"(OC(=O)CCOP(=O)(O)OP(=O)(O)O)", R"(3-{[hydroxy(phosphonooxy)phosphoryl]oxy}propanoic acid)"},
      {R"(COC1=C(C=C(C=C1)C(=O)NNC(=O)C2=CC(=CC=C2)O)OC)", R"(N'-(3-hydroxybenzoyl)-3,4-dimethoxybenzohydrazide)"},
      {R"(CCOC1=CC=C(C=C1)NC(=S)NC(=O)C(C)(C)C)", R"(N-[(4-ethoxyphenyl)carbamothioyl]-2,2-dimethylpropanamide)"},
      {R"(C[Si](C)(C)c1ccccc1)", R"(trimethyl(phenyl)silane)"},
      {R"(CC(C)(C)[Si](C)(C)OCC(=O)O)", R"(2-{[dimethyl(2-methylpropan-2-yl)silyl]oxy}acetic acid)"},
      {R"(OCCOS(=O)(=O)O)", R"(2-hydroxyethyl hydrogen sulfate)"},
      {R"(C(C(C(=O)O)N)S(=O)O)", R"(2-amino-3-sulfinopropanoic acid)"},
      {R"(c1ccc(/N=N/c2ccccc2)cc1)", R"((E)-diphenyldiazene)"},
      {R"(OC(=O)c1ccc(N=Nc2ccc(N(C)C)cc2)cc1)", R"(4-{[4-(dimethylamino)phenyl]diazenyl}benzoic acid)"},
      {R"([O-][n+]1ccccc1)", R"(1-oxidopyridin-1-ium)"},
      {R"(C[N+](C)(C)[O-])", R"(N,N-dimethyl-N-oxidomethanaminium)"},
      {R"(Nn1cnnc1)", R"(4H-1,2,4-triazol-4-amine)"},
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
  for (const char* smiles : {"CC[B](C)C", "CCO.CCCO", "CC(=O)OP(=O)(O)O", "c1cc2ccc3cccc4ccc(c1)c2c34", "CN=[N+]=[N-]", "COP(=S)(OC)OC",
                              "[Na+].[Cl-].c1ccccc1", "C1CC2(C1)c1ccccc12"}) {
    Diagnostics d;
    const NameResult r = iupac_name_from_smiles(smiles, d);
    CHECK(!r.ok);
    CHECK(d.contains(DiagCode::NamingUnsupported));
  }
}
