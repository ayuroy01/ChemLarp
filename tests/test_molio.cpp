// Structure-file readers (MOL V2000/V3000, SDF, PDB, CML).
// Expected SMILES for the PDB fixtures were checked against RDKit 2026.09 reading
// the same text (same canonical SMILES and InChIKey).
#include <cstdio>
#include <string>
#include <vector>

#include "omnichem/canonical.hpp"
#include "omnichem/conformer.hpp"
#include "omnichem/depict.hpp"
#include "omnichem/molio.hpp"
#include "omnichem/smiles.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

#include "data/molio_pdb.inc"

struct A {
  const char* sym;
  double x, y, z;
  int ccode = 0;  // V2000 charge code
  int vvv = 0;    // V2000 valence field
};
struct B {
  int a, b, type, stereo = 0;
};

std::string v2000(const std::vector<A>& atoms, const std::vector<B>& bonds, const std::string& props = "",
                  const std::string& title = "test") {
  std::string out = title + "\n  OmniChem test\n\n";
  char buf[128];
  std::snprintf(buf, sizeof buf, "%3zu%3zu  0  0  0  0  0  0  0  0999 V2000\n", atoms.size(), bonds.size());
  out += buf;
  for (const A& a : atoms) {
    std::snprintf(buf, sizeof buf, "%10.4f%10.4f%10.4f %-3s 0%3d  0  0  0%3d  0  0  0  0  0  0\n", a.x, a.y, a.z, a.sym,
                  a.ccode, a.vvv);
    out += buf;
  }
  for (const B& b : bonds) {
    std::snprintf(buf, sizeof buf, "%3d%3d%3d%3d\n", b.a, b.b, b.type, b.stereo);
    out += buf;
  }
  return out + props + "M  END\n";
}

// Reads one record and returns its canonical SMILES ("" on failure).
std::string read1(const std::string& text, Diagnostics* out_diag = nullptr, ReadOptions opt = {}) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto recs = read_molecules(st, text, d, opt);
  if (out_diag) *out_diag = d;
  if (recs.empty() || recs[0].species < 0) return "";
  return st.species[static_cast<std::size_t>(recs[0].species)].canonical_smiles;
}

std::string canon(const char* smiles) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = parse_smiles(st, smiles, d);
  return sp < 0 ? "" : st.species[static_cast<std::size_t>(sp)].canonical_smiles;
}

// Alanine drawn as in ChemDraw: C1 methyl, C2 centre, N3, C4 carboxyl, O5 (=), O6.
std::string alanine(int stereo) {
  return v2000({{"C", 0, 0.375, 0}, {"C", 1.299, -0.375, 0}, {"N", 1.299, -1.875, 0}, {"C", 2.598, 0.375, 0},
                {"O", 2.598, 1.875, 0}, {"O", 3.897, -0.375, 0}},
               {{2, 1, 1, stereo}, {2, 3, 1}, {2, 4, 1}, {4, 5, 2}, {4, 6, 1}});
}

}  // namespace

TEST(molio_v2000_wedges) {
  CHECK_EQ(read1(alanine(1)), std::string("C[C@H](N)C(=O)O"));
  CHECK_EQ(read1(alanine(6)), std::string("C[C@@H](N)C(=O)O"));
  CHECK_EQ(read1(alanine(4)), std::string("CC(N)C(=O)O"));  // wavy: unknown
  CHECK_EQ(read1(alanine(0)), std::string("CC(N)C(=O)O"));  // no wedge in 2D: no stereo
  // A wedge to an explicit hydrogen: H up means the methyl side is down.
  const std::string h_wedge =
      v2000({{"C", 0, 0.375, 0}, {"C", 1.299, -0.375, 0}, {"N", 1.299, -1.875, 0}, {"C", 2.598, 0.375, 0},
             {"O", 2.598, 1.875, 0}, {"O", 3.897, -0.375, 0}, {"H", 1.299, 1.125, 0}},
            {{2, 1, 1}, {2, 3, 1}, {2, 4, 1}, {4, 5, 2}, {4, 6, 1}, {2, 7, 1, 1}});
  CHECK_EQ(read1(h_wedge), canon("C[C@@H](N)C(=O)O"));
}

TEST(molio_v2000_double_bonds) {
  const std::vector<B> bonds{{1, 2, 1}, {2, 3, 2}, {3, 4, 1}};
  CHECK_EQ(read1(v2000({{"F", 0, 1, 0}, {"C", 0, 0, 0}, {"C", 1.3, 0, 0}, {"F", 1.3, 1, 0}}, bonds)), canon("F/C=C\\F"));
  CHECK_EQ(read1(v2000({{"F", 0, 1, 0}, {"C", 0, 0, 0}, {"C", 1.3, 0, 0}, {"F", 1.3, -1, 0}}, bonds)), canon("F/C=C/F"));
  CHECK_EQ(read1(v2000({{"F", 0, 1, 0}, {"C", 0, 0, 0}, {"C", 1.3, 0, 0}, {"F", 1.3, -1, 0}}, {{1, 2, 1}, {2, 3, 2, 3}, {3, 4, 1}})),
           std::string("FC=CF"));  // crossed: either
  // Cumulated: an allene carries no E/Z even when drawn bent.
  CHECK_EQ(read1(v2000({{"C", 0, 1, 0}, {"C", 0, 0, 0}, {"C", 1.3, 0, 0}, {"C", 2.6, 0, 0}, {"C", 2.6, 1, 0}},
                       {{1, 2, 1}, {2, 3, 2}, {3, 4, 2}, {4, 5, 1}})),
           std::string("CC=C=CC"));
}

TEST(molio_v2000_properties) {
  // Atom-block charge code 3 = +1; M  CHG replaces atom-block charges.
  CHECK_EQ(read1(v2000({{"N", 0, 0, 0, 3}}, {})), std::string("[NH4+]"));
  CHECK_EQ(read1(v2000({{"N", 0, 0, 0, 3}, {"O", 1, 0, 0}}, {}, "M  CHG  1   2  -1\n")), std::string("N.[OH-]"));
  CHECK_EQ(read1(v2000({{"C", 0, 0, 0}}, {}, "M  ISO  1   1  13\n")), std::string("[13CH4]"));
  CHECK_EQ(read1(v2000({{"C", 0, 0, 0}}, {}, "M  RAD  1   1   2\n")), std::string("[CH3]"));
  CHECK_EQ(read1(v2000({{"D", 0, 0, 0}, {"O", 1, 0, 0}}, {{1, 2, 1}})), std::string("[2H]O"));
  // Valence field: Si is outside the implicit-H set, so the field gives SiH2.
  CHECK_EQ(read1(v2000({{"C", 0, 0, 0}, {"Si", 1, 0, 0, 0, 4}, {"C", 2, 0, 0}}, {{1, 2, 1}, {2, 3, 1}})),
           std::string("C[SiH2]C"));
  CHECK_EQ(read1(v2000({{"C", 0, 0, 0}, {"Si", 1, 0, 0}, {"C", 2, 0, 0}}, {{1, 2, 1}, {2, 3, 1}})), std::string("C[Si]C"));
  CHECK_EQ(read1(v2000({{"Na", 0, 0, 0, 3}, {"Cl", 3, 0, 0, 5}}, {})), std::string("[Na+].[Cl-]"));
  CHECK_EQ(mdl_valence_field(6, 0, 4, 0), 0);
  CHECK_EQ(mdl_valence_field(6, 0, 1, 0), 1);   // [CH0] radical-like carbon
  CHECK_EQ(mdl_valence_field(14, 0, 2, 2), 4);  // SiH2: always explicit
  CHECK_EQ(mdl_valence_field(11, 1, 0, 0), 15);
}

TEST(molio_v3000) {
  const std::string mol =
      "bromochlorofluoromethane\n  test\n\n  0  0  0     0  0            999 V3000\n"
      "M  V30 BEGIN CTAB\n"
      "M  V30 COUNTS 5 4 0 0 0\n"
      "M  V30 BEGIN ATOM\n"
      "M  V30 1 C 0 0 0 0\n"
      "M  V30 2 F -1 0.5 0 0\n"
      "M  V30 3 Cl 1 0.5 0 0\n"
      "M  V30 4 Br 0 -1.2 0 0\n"
      "M  V30 5 N 0 3 0 0 CHG=1\n"
      "M  V30 END ATOM\n"
      "M  V30 BEGIN BOND\n"
      "M  V30 1 1 1 2 CFG=1\n"
      "M  V30 2 1 1 3\n"
      "M  V30 3 1 1 -\n"
      "M  V30 4\n"
      "M  V30 END BOND\n"
      "M  V30 END CTAB\n"
      "M  END\n";
  // Bond 3 continues on the next line ('-'), giving "3 1 1 4". F is wedged up and the
  // implicit H points down; Br -> Cl -> F turns counterclockwise from above: (S),
  // which RDKit's CIP labeller gives for F[C@@H](Cl)Br. N+ is a separate ammonium ion.
  Diagnostics d;
  const std::string got = read1(mol, &d);
  CHECK(!d.has_errors());
  CHECK_EQ(got, canon("F[C@@H](Cl)Br.[NH4+]"));
}

TEST(molio_sdf_records) {
  std::string sdf = alanine(1) + "> <NAME>\nL-alanine\n\n> <NOTE>\nline one\nline two\n\n$$$$\n";
  sdf += v2000({{"A", 0, 0, 0}}, {}, "", "query") + "$$$$\n";  // query atom: error, record skipped
  sdf += v2000({{"O", 0, 0, 0}}, {}, "", "water") + "$$$$\n";
  ChemicalSystemState st;
  Diagnostics d;
  const auto recs = read_molecules(st, sdf, d);
  CHECK_EQ(recs.size(), std::size_t{3});
  CHECK_EQ(recs[0].title, std::string("test"));
  CHECK_EQ(recs[0].properties.size(), std::size_t{2});
  CHECK_EQ(recs[0].properties[0].second, std::string("L-alanine"));
  CHECK_EQ(recs[0].properties[1].second, std::string("line one\nline two"));
  CHECK_EQ(recs[0].dimensions, 2);
  CHECK(recs[1].species < 0);
  CHECK(d.contains(DiagCode::UnsupportedFileFeature));
  CHECK_EQ(recs[2].title, std::string("water"));
  CHECK_EQ(st.species[static_cast<std::size_t>(recs[2].species)].canonical_smiles, std::string("O"));
  CHECK_EQ(st.species.size(), std::size_t{2});  // the failed record left nothing behind
  ReadOptions one;
  one.max_records = 1;
  ChemicalSystemState st2;
  Diagnostics d2;
  CHECK_EQ(read_molecules(st2, sdf, d2, one).size(), std::size_t{1});
}

TEST(molio_explicit_hydrogens) {
  const std::string methane = v2000({{"C", 0, 0, 0}, {"H", 1, 0, 0}, {"H", -1, 0, 0}, {"H", 0, 1, 0}, {"H", 0, -1, 0}},
                                    {{1, 2, 1}, {1, 3, 1}, {1, 4, 1}, {1, 5, 1}});
  ChemicalSystemState st;
  Diagnostics d;
  auto recs = read_molecules(st, methane, d);
  CHECK_EQ(st.species[static_cast<std::size_t>(recs[0].species)].atom_count, 1u);
  ReadOptions keep;
  keep.keep_hydrogens = true;
  ChemicalSystemState st2;
  recs = read_molecules(st2, methane, d, keep);
  CHECK_EQ(st2.species[static_cast<std::size_t>(recs[0].species)].atom_count, 5u);
  CHECK_EQ(st2.species[static_cast<std::size_t>(recs[0].species)].formula, std::string("CH4"));
}

// Our own writers read back to the same structure (2D wedges / 3D coordinates).
TEST(molio_round_trip_writers) {
  const char* cases[] = {
      "C[C@H](N)C(=O)O",         "F/C=C/[C@@H](Cl)Br",      "C/C=C\\C(=O)O",
      "C[S@](=O)c1ccccc1",       "O[C@H]1CCCC[C@@H]1O",     "[Na+].[O-]C(=O)c1ccccc1",
      "[13CH3][C@H](O)CC",       "C[SiH2]C",                "[B]c1ccccc1",
      "N[C@@H](Cc1ccccc1)C(=O)O", "CC(C)=CC/C=C(\\C)CO",    "C[N+](C)(C)CC([O-])=O",
  };
  for (const char* smi : cases) {
    Diagnostics d;
    const Depiction dp = depict_smiles(smi, d);
    CHECK(dp.ok);
    CHECK_EQ(read1(dp.molfile), canon(smi));
  }
  for (const char* smi : {"C[C@H](N)C(=O)O", "F/C=C/[C@@H](Cl)Br", "O[C@H]1CCCC[C@@H]1O", "C[S@](=O)CC"}) {
    ChemicalSystemState st;
    Diagnostics d;
    const auto sp = embed_smiles(st, smi, d);
    CHECK(sp >= 0);
    if (sp < 0) continue;
    CHECK_EQ(read1(write_sdf(st, static_cast<std::uint32_t>(sp))), canon(smi));
  }
  // 3D: an azide's N=N+ is cumulated and gets no E/Z.
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = embed_smiles(st, "CCN=[N+]=[N-]", d);
  CHECK(sp >= 0);
  if (sp >= 0) {
    ChemicalSystemState back;
    const auto recs = read_molecules(back, write_sdf(st, static_cast<std::uint32_t>(sp)), d);
    CHECK(back.double_bonds.empty());
    CHECK_EQ(recs[0].dimensions, 3);
  }
}

TEST(molio_pdb) {
  Diagnostics d;
  // Same canonical SMILES and InChIKey as RDKit for this text.
  CHECK_EQ(read1(kCrambin5, &d),
           std::string("C[C@@H](O)[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CS)C(=O)N[C@@H](CS)C(=O)N1CCC[C@H]1C=O)[C@@H](C)O"));
  CHECK(!d.has_errors());
  // Hydrogens listed: the N-terminus is NH3+, and the cut C-terminus keeps no H.
  Diagnostics d2;
  CHECK_EQ(read1(kTrpCage3, &d2), std::string("CC(C)C[C@H](NC(=O)[C@@H]([NH3+])CC(N)=O)C(=O)N[C@H]([C]=O)Cc1ccc(O)cc1"));
  CHECK(d2.contains(DiagCode::FileContentIgnored));  // second MODEL
  CHECK_EQ(detect_format(kCrambin5), MolFormat::Pdb);
}

TEST(molio_cml) {
  const std::string cml = R"(<?xml version="1.0"?>
<cml xmlns="http://www.xml-cml.org/schema">
 <molecule id="m1" title="L-alanine">
  <atomArray>
   <atom id="a1" elementType="C" x2="0" y2="0.375"/>
   <atom id="a2" elementType="C" x2="1.299" y2="-0.375"/>
   <atom id="a3" elementType="N" x2="1.299" y2="-1.875"/>
   <atom id="a4" elementType="C" x2="2.598" y2="0.375"/>
   <atom id="a5" elementType="O" x2="2.598" y2="1.875"/>
   <atom id="a6" elementType="O" x2="3.897" y2="-0.375"/>
  </atomArray>
  <bondArray>
   <bond atomRefs2="a2 a1" order="1"><bondStereo>W</bondStereo></bond>
   <bond atomRefs2="a2 a3" order="1"/>
   <bond atomRefs2="a2 a4" order="S"/>
   <bond atomRefs2="a4 a5" order="D"/>
   <bond atomRefs2="a4 a6" order="1"/>
  </bondArray>
 </molecule>
 <molecule id="m2">
  <name>ammonium &amp; friend</name>
  <atomArray atomID="n1 o1" elementType="N O" formalCharge="1 -1" hydrogenCount="4 1"/>
 </molecule>
</cml>)";
  ChemicalSystemState st;
  Diagnostics d;
  const auto recs = read_molecules(st, cml, d);
  CHECK(!d.has_errors());
  CHECK_EQ(recs.size(), std::size_t{2});
  CHECK_EQ(recs[0].title, std::string("L-alanine"));
  CHECK_EQ(st.species[static_cast<std::size_t>(recs[0].species)].canonical_smiles, std::string("C[C@H](N)C(=O)O"));
  CHECK_EQ(recs[1].title, std::string("ammonium & friend"));
  CHECK_EQ(st.species[static_cast<std::size_t>(recs[1].species)].canonical_smiles, std::string("[NH4+].[OH-]"));
  CHECK_EQ(detect_format(cml), MolFormat::Cml);
}

TEST(molio_detect_and_errors) {
  CHECK_EQ(detect_format("", "x.SDF"), MolFormat::Molfile);
  CHECK_EQ(detect_format("", "1abc.pdb"), MolFormat::Pdb);
  CHECK_EQ(detect_format(alanine(0)), MolFormat::Molfile);
  CHECK_EQ(detect_format("CCO"), MolFormat::Auto);
  ChemicalSystemState st;
  Diagnostics d;
  CHECK(read_molecules(st, "CCO", d).empty());
  CHECK(d.contains(DiagCode::UnknownFileFormat));
  // Truncated atom block.
  Diagnostics d2;
  const std::string cut = alanine(0).substr(0, 120);
  CHECK_EQ(read1(cut, &d2), std::string(""));
  CHECK(d2.contains(DiagCode::MalformedRecord));
  // Valence error is reported with the file line and atom number.
  Diagnostics d3;
  CHECK_EQ(read1(v2000({{"C", 0, 0, 0}, {"O", 1, 0, 0}}, {{1, 2, 3}, {1, 2, 2}}), &d3), std::string(""));
  CHECK(d3.has_errors());
}
