// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Hostile and extreme input through the public C API (docs/AUDITS.md, audits 1.7 and 2.6).
//
// Every entry point must come back with a JSON object and an error vector, never crash, hang, read out of
// bounds or loop forever. Null pointers return null. Nothing here checks chemistry: only that malformed,
// empty, enormous, deeply nested, non-UTF-8 and out-of-range input is rejected in an orderly way. Run it
// under the sanitizers as well as normally.
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "chemlarp/chemlarp_c.h"
#include "test_framework.hpp"

namespace {

// Calls f, requires a non-null JSON object back, and frees it.
std::string take(char* s) {
  CHECK(s != nullptr);
  std::string out = s ? s : "";
  cl_free_string(s);
  CHECK(!out.empty());
  CHECK(out.front() == '{' || out.front() == '[');
  return out;
}

// Inputs that have broken parsers before: empty, blank, huge, deeply nested, unterminated, binary.
std::vector<std::string> nasty_strings() {
  std::vector<std::string> v = {
      "", " ", "\n", "\t\t\t", "C" + std::string(20000, 'C'), std::string(5000, '('),
      std::string(5000, ')'), std::string(2000, '['), "[", "]", "(", ")", "%", "%%", "%99", "@", "@@",
      "/\\", "$", "&", "*", "~", "#", "=", ".", "..", "....", "-]", "[*]", "[]", "[C@@H", "C1CC", "C%",
      "C1" + std::string(500, 'C'), std::string(1000, '1'), "\xff\xfe\xfd", "\xc3\x28", "\xe2\x82",
      "C\xff", "<<<>>>", "{\"a\":", "\"", "'", "\\", "\\\\", "N[C@@H](C)C(=O)O\n\n\n", "0", "-1", "1e400",
      "NaN", "inf", "null", "undefined", std::string(20000, 'x'), std::string(200, '\n'),
      "C(((((((((((((((((((((C)))))))))))))))))))))", "c1ccccc1" + std::string(100, '1'),
  };
  return v;
}

}  // namespace

TEST(c_api_rejects_null_pointers) {
  CHECK(cl_parse_smiles_json(nullptr) == nullptr);
  CHECK(cl_canonical_smiles_json(nullptr, 0, 0) == nullptr);
  CHECK(cl_inchi_json(nullptr) == nullptr);
  CHECK(cl_iupac_name_json(nullptr) == nullptr);
  CHECK(cl_check_reaction_json(nullptr, 298.15) == nullptr);
  CHECK(cl_aqueous_json(nullptr) == nullptr);
  CHECK(cl_flowsheet_json(nullptr) == nullptr);
  CHECK(cl_lab_mix_json(nullptr) == nullptr);
  CHECK(cl_lab_titration_json(nullptr) == nullptr);
  CHECK(cl_mbar_json(nullptr) == nullptr);
  CHECK(cl_ghs_json(nullptr) == nullptr);
  CHECK(cl_hydration_json(nullptr, 0, 0, 0, 0, 0) == nullptr);
  CHECK(cl_dock_json(nullptr, "C", nullptr, 0, 0, 0, 10, 10, 10, 1, 1, 1, 7.4, 1, 0) == nullptr);
  CHECK(cl_dock_json("ATOM\n", nullptr, nullptr, 0, 0, 0, 10, 10, 10, 1, 1, 1, 7.4, 1, 0) == nullptr);
  CHECK(cl_smarts_match_json(nullptr, "C", 0, 10) == nullptr);
  CHECK(cl_smarts_match_json("C", nullptr, 0, 10) == nullptr);
  CHECK(cl_similarity_json(nullptr, "C", 2, 2048) == nullptr);
  CHECK(cl_read_molecules_json(nullptr, nullptr, 0, 0) == nullptr);
  cl_free_string(nullptr);  // must be a no-op
}

TEST(structure_parsers_survive_hostile_input) {
  for (const auto& s : nasty_strings()) {
    take(cl_parse_smiles_json(s.c_str()));  // the parser sees everything, however large
    if (s.size() > 2000) continue;          // canonicalisation and identity are superlinear
    take(cl_canonical_smiles_json(s.c_str(), 2, 1));
    take(cl_inchi_json(s.c_str()));
    take(cl_ghs_json(s.c_str()));
    take(cl_smarts_match_json("CCO", s.c_str(), 0, 100));
    take(cl_smarts_match_json(s.c_str(), "[#6]", 0, 100));
    take(cl_similarity_json(s.c_str(), "CCO", 2, 2048));
    take(cl_read_molecules_json(s.c_str(), nullptr, 0, 0));
  }
}

TEST(naming_and_depiction_survive_hostile_input) {
  for (const auto& s : nasty_strings()) {
    if (s.size() > 4096) continue;  // naming is superlinear; the parser is covered above
    take(cl_iupac_name_json(s.c_str()));
    take(cl_depict_json(s.c_str(), 30.0, 0));
  }
}

TEST(file_readers_survive_malformed_files) {
  const std::vector<std::string> files = {
      "", "\n\n\n", "$$$$", "M  END", "garbage\nnot a molfile\n$$$$\n",
      // A V2000 header promising atoms and bonds the file does not have.
      "x\n\n\n 99 99  0  0  0  0  0  0  0  0999 V2000\nM  END\n",
      "x\n\n\n-1 -1  0  0  0  0  0  0  0  0999 V2000\nM  END\n",
      "x\n\n\n  1  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n",
      "ATOM  xxxx\n", "ATOM      1  CA  ALA A   1       x       y       z\n",
      "HETATM\nATOM\nTER\nEND\n", "data_x\nloop_\n_atom_site.id\n",
      std::string(100000, 'A') + "\n",
  };
  for (const auto& f : files) {
    take(cl_read_molecules_json(f.c_str(), nullptr, 0, 0));
    take(cl_read_molecules_json(f.c_str(), "sdf", 1, 5));
    take(cl_protein_json(f.c_str(), 7.4));
  }
}

TEST(numeric_entry_points_reject_out_of_range_conditions) {
  const double bad_T[] = {0.0, -1.0, -1e30, 1e30, 1e-30};
  for (double T : bad_T) {
    take(cl_check_reaction_json("1 C + 1 O=O -> 1 O=C=O", T));
    take(cl_equilibrate_text_json("O 1\n[H][H] 1", T, 1e5));
    take(cl_saturation_pressure_json("CCO", T, "pr"));
    take(cl_pure_component_json("CCO", T, 1e5));
    take(cl_flash_text_json("CCO 1", T, 1e5, "pr"));
  }
  const double bad_P[] = {0.0, -1.0, -1e20, 1e30};
  for (double P : bad_P) {
    take(cl_equilibrate_text_json("O 1\n[H][H] 1", 1000.0, P));
    take(cl_flash_text_json("CCO 0.5\nCCCC 0.5", 300.0, P, "srk"));
  }
  // Zero and negative amounts, and a mixture that sums to nothing.
  take(cl_equilibrate_text_json("O 0\n[H][H] 0", 1000.0, 1e5));
  take(cl_equilibrate_text_json("O -1", 1000.0, 1e5));
  take(cl_equilibrate_text_json("", 1000.0, 1e5));
  take(cl_flash_text_json("CCO 0\nCCCC 0", 300.0, 1e5, "pr"));
  take(cl_vle_text_json("CCO 0.5\nO 0.5", "txy", -5.0));
  take(cl_vle_text_json("", "txy", 1.0));
}

TEST(simulation_entry_points_reject_absurd_settings) {
  // Zero, negative and enormous step counts, step sizes and system sizes.
  take(cl_md_json("CCO", 300.0, 0, 0.0, nullptr, 0.0, 0, 0, 0, nullptr));
  take(cl_md_json("CCO", -300.0, 10, -1.0, "bussi", -1.0, 0, 1, 1, "uff"));
  take(cl_md_json("", 300.0, 10, 1.0, nullptr, 100.0, 1, 1, 0, nullptr));
  take(cl_md_json("CCO", 300.0, 5, 1.0, "no-such-thermostat", 100.0, 1, 1, 0, "no-such-field"));
  take(cl_lj_json(0, 0.0, 0.0, 0.0, 0.0, 0, 0, nullptr, 0));
  take(cl_lj_json(-5, -1.0, -1.0, -1.0, -1.0, 0, 0, "none", 1));
  take(cl_mmff_json("", "mmff94", 1));
  take(cl_mmff_json("CCO", "not-a-variant", 1));
  take(cl_embed_json("", 1, 1, 1));
  take(cl_embed_json("CCO", 1, 0, 0));
  take(cl_embed_json("CCO", 1, -3, 1));
  take(cl_qc_json("", "hf", "STO-3G", 0, 2, 0));
  take(cl_qc_json("CCO", "no-such-method", "no-such-basis", -2, -1, 0));
  // Boxes too small for the cutoff, and no molecules at all.
  take(cl_water_box_json(1, 298.15, 1.0, 0.01, 0.01, 1));
}

// An argument out of range must be reported at once. These used to fall back to the defaults, so a caller
// passing -1 started a 400-water, 120 ps simulation instead of getting an error (docs/AUDITS.md, audit 2.6).
TEST(expensive_entry_points_reject_out_of_range_arguments_at_once) {
  const auto rejected = [](char* s) {
    const std::string j = take(s);
    CHECK(j.rfind("{\"ok\":false", 0) == 0);
  };
  rejected(cl_water_box_json(-10, 298.15, 1.0, 1.0, 1.0, 1));
  rejected(cl_water_box_json(100, -5.0, 1.0, 1.0, 1.0, 1));
  rejected(cl_water_box_json(100, 298.15, 1.0, -1.0, 1.0, 1));
  rejected(cl_water_box_json(100, 298.15, 1.0, 1.0, -1.0, 1));
  rejected(cl_water_box_json(1 << 20, 298.15, 1.0, 1.0, 1.0, 1));
  rejected(cl_water_box_json(100, 1e9, 1.0, 1.0, 1.0, 1));
  rejected(cl_water_box_json(100, 298.15, 1.0, 1e9, 1.0, 1));
  rejected(cl_hydration_json("CCO", -1, 298.15, 1.0, 1.0, 1));
  rejected(cl_hydration_json("CCO", 100, -1.0, 1.0, 1.0, 1));
  rejected(cl_hydration_json("CCO", 100, 298.15, -1.0, 1.0, 1));
  rejected(cl_hydration_json("CCO", 100, 298.15, 1.0, -1.0, 1));
  // Docking: a negative or absurd box, run count, step count, pose count or pH is a mistake by the caller,
  // not a request for the default, which could start a search of hours.
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, -10, 10, 10, 1, 1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 1e9, 10, 10, 1, 1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 0, 10, 1, 1, 1, 7.4, 1, 0));  // two edges of three
  rejected(cl_dock_json("ATOM\n", "C", nullptr, std::nan(""), 0, 0, 10, 10, 10, 1, 1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, HUGE_VAL, 0, 0, 10, 10, 10, 1, 1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, -1, 1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1 << 20, 1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1, -1, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1, 1 << 24, 1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1, 1, -1, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1, 1, 10000, 7.4, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1, 1, 1, -5.0, 1, 0));
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 10, 10, 10, 1, 1, 1, 100.0, 1, 0));
  // No box and no reference ligand: there is nothing to say where to search.
  rejected(cl_dock_json("ATOM\n", "C", nullptr, 0, 0, 0, 0, 0, 0, 1, 1, 1, 7.4, 1, 0));
}

TEST(docking_survives_hostile_receptors_and_ligands) {
  // A receptor of one atom, so a search is cheap, and every nasty string in turn as the ligand.
  const char* receptor =
      "HETATM    1  C1  PKT A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n";
  for (const auto& s : nasty_strings()) {
    // Edge cases, not long inputs: a docking call parses its ligand, builds a grid and searches, and the
    // long strings are already put through every parser by the tests above. One oversized ligand is checked
    // on its own below.
    if (s.size() > 500) continue;
    take(cl_dock_json(receptor, s.c_str(), nullptr, 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
    take(cl_dock_json(s.c_str(), "C", nullptr, 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
    take(cl_dock_json(receptor, "C", s.c_str(), 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
  }
  // A "ligand" of hundreds of atoms is refused on sight: generating a conformation for it would take hours.
  {
    const std::string huge = "C" + std::string(400, 'C');
    const std::string j = take(cl_dock_json(receptor, huge.c_str(), nullptr, 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
    CHECK(j.rfind("{\"ok\":false", 0) == 0);
    CHECK(j.find("heavy atoms") != std::string::npos);
  }
  // A receptor whose coordinates are absurd: the grid and the cell lists must clamp, not overflow.
  const char* far =
      "HETATM    1  C1  PKT A   1    99999.9 99999.9 99999.9  1.00  0.00           C\n"
      "HETATM    2  C2  PKT A   1   -99999.9-99999.9-99999.9  1.00  0.00           C\nEND\n";
  take(cl_dock_json(far, "C", nullptr, 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
  take(cl_dock_json(far, "C", nullptr, 99999.9, 99999.9, 99999.9, 200, 200, 200, 1, 1, 1, 7.4, 1, 0));
}

TEST(text_formats_survive_hostile_input) {
  for (const auto& s : nasty_strings()) {
    if (s.size() > 70000) continue;
    take(cl_aqueous_json(s.c_str()));
    take(cl_flowsheet_json(s.c_str()));
    take(cl_lab_mix_json(s.c_str()));
    take(cl_lab_titration_json(s.c_str()));
    take(cl_mbar_json(s.c_str()));
    take(cl_check_reaction_json(s.c_str(), 298.15));
    take(cl_equilibrate_text_json(s.c_str(), 1000.0, 1e5));
    take(cl_vle_text_json(s.c_str(), "txy", 1.0));
    take(cl_flash_text_json(s.c_str(), 300.0, 1e5, "pr"));
  }
}

TEST(mechanism_readers_survive_malformed_mechanisms) {
  const std::vector<std::string> mechs = {
      "", "ELEMENTS\n", "REACTIONS\nEND\n", "species:\n", "phases:\n  - name: x\n",
      "reaction H2 + O2 => H2O  1 0 0\n", "reaction  => \n", "reaction A => B 1e300 1e300 1e300\n",
      "species A\nreaction A => A 1 0 0\n", "{\"not\": \"a mechanism\"}",
      std::string(50000, 'r') + "\n",
  };
  for (const auto& m : mechs) {
    take(cl_kinetics_text_json(m.c_str(), nullptr, nullptr, "adiabatic-v", 1000.0, 1e5, "", 1e-6, 10, 1e-8));
    take(cl_kinetics_text_json(m.c_str(), m.c_str(), m.c_str(), "no-such-model", -1.0, -1.0, "H2=1", -1.0, 0, 0.0));
    take(cl_kinetics_text_json(m.c_str(), nullptr, nullptr, "iso-p", 300.0, 1e5, "not a feed line", 1e-6, 1, 1e-8));
  }
  // `mechanism_text` and `model` are required: NULL in, NULL out (the header says so).
  CHECK(cl_kinetics_text_json(nullptr, nullptr, nullptr, "iso-v", 300.0, 1e5, "", 1e-6, 1, 1e-8) == nullptr);
  CHECK(cl_kinetics_text_json("", nullptr, nullptr, nullptr, 300.0, 1e5, "", 1e-6, 1, 1e-8) == nullptr);
}

TEST(search_entry_point_survives_hostile_input) {
  take(cl_search_json("", "[#6]", "CCO", 2, 2048));
  take(cl_search_json("CCO\nCCC\n", "", "", 2, 2048));
  take(cl_search_json("CCO\n\n\nnot-a-smiles\n", "[", "](", 0, 0));
  take(cl_search_json("CCO", "[#6]", "CCO", -1, -1));
  take(cl_search_json(std::string(20000, 'C').c_str(), "[#6]", "C", 2, 2048));
}

// A long chain used to recurse once per atom in the canonical SMILES writer, so it overflowed the stack and
// crashed the process; AddressSanitizer caught it in CI on a 20,000-atom chain (docs/AUDITS.md, audit 1.3).
// The writer now follows the last child in a loop. Deeply nested *branches* still recurse, so they are
// refused with a diagnostic instead of crashing. 4,000 atoms keeps this test to a couple of seconds while
// staying far past the depth that overflowed under a sanitizer's larger stack frames.
TEST(long_chains_and_deep_branches_do_not_overflow_the_stack) {
  const std::string chain(4000, 'C');
  const std::string out = take(cl_parse_smiles_json(chain.c_str()));
  CHECK(out.rfind("{\"ok\":true", 0) == 0);
  CHECK(out.find(std::string(64, 'C')) != std::string::npos);  // the chain came back as a chain
  take(cl_canonical_smiles_json(chain.c_str(), 0, 1));

  // Branches nested past the writer's limit: an error, not a crash. "C(C(C(...)C)C)C" nests one per atom.
  const std::size_t deep = 6000;
  std::string nested;
  nested.reserve(4 * deep);
  for (std::size_t i = 0; i < deep; ++i) nested += "C(";
  nested += "C";
  for (std::size_t i = 0; i < deep; ++i) nested += ")C";
  take(cl_parse_smiles_json(nested.c_str()));  // parses or is rejected; either way it must not crash
}

// A PDB coordinate far outside the usual range used to be converted to int when hashing atoms into cells,
// which is undefined behaviour for a value that does not fit. The fuzzer found it with 2.16e8 in a mangled
// ATOM record (docs/AUDITS.md, audit 1.5); the index is now clamped, which also covers NaN and infinity.
TEST(protein_reader_survives_absurd_coordinates) {
  const char* fmt =
      "ATOM      1  N   ASN A   1    %8s%8s%8s  1.00  0.00           N\n"
      "ATOM      2  CA  ASN A   1       1.000   2.000   3.000  1.00  0.00           C\n"
      "END\n";
  for (const char* v : {"2.16e+08", "-2.16e+08", "1e+300", "-1e+300", "     nan", "     inf", "    -inf",
                        "99999999", "      ##"}) {
    char pdb[512];
    std::snprintf(pdb, sizeof pdb, fmt, v, v, v);
    take(cl_protein_json(pdb, 7.4));
  }
}

// Numbers inside a bracket atom were accumulated digit by digit with no cap, so a long run of digits
// overflowed before the range check could reject it — undefined behaviour for the signed ones. The fuzzer
// found it with a hydrogen count (docs/AUDITS.md, audit 1.5); the same shape was in the isotope, the charge
// magnitude, the atom class and SMARTS's number(). All are capped now, and still rejected as out of range.
TEST(bracket_atom_numbers_do_not_overflow) {
  for (const char* s : {"[CH999999999999]", "[999999999999CH4]", "[C+999999999999]", "[C-999999999999]",
                        "[12345678901234567890CH4]", "[C+++++++++++++++++++++++++++++++++]"}) {
    const std::string j = take(cl_parse_smiles_json(s));
    CHECK(j.rfind("{\"ok\":false", 0) == 0);  // out of range, reported rather than wrapped
  }
  // The atom class is clamped to 65535 instead of rejected, which is what it did before; what matters here
  // is that the digits no longer overflow on the way.
  take(cl_parse_smiles_json("[CH4:999999999999]"));
  // The legal forms still parse, so the cap did not narrow what SMILES accepts.
  for (const char* s : {"[13CH3]", "[235U]", "[NH4+]", "[O-2]", "[Fe+3]", "[CH4:99]", "[CH3:65535]"})
    CHECK(take(cl_parse_smiles_json(s)).rfind("{\"ok\":true", 0) == 0);
  take(cl_smarts_match_json("CCO", "[R999999999999]", 0, 10));
  take(cl_smarts_match_json("CCO", "[v999999999999]", 0, 10));
}
