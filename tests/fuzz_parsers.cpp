// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Coverage-guided fuzzing of every text input (docs/AUDITS.md, audit 1.5).
//
// One libFuzzer target for all of them: the first byte of the input chooses which parser sees the rest, so a
// single corpus and a single run cover the lot. Everything goes through the C API, which is the boundary an
// attacker reaches through the web app, the command line or a binding, and which must never crash, hang or
// read out of bounds, whatever it is given.
//
// Build and run (Linux, or any clang with libFuzzer; Apple's clang does not ship it):
//
//   clang++ -std=c++20 -O1 -g -fsanitize=fuzzer,address,undefined -fno-sanitize-recover=undefined \
//     -Iinclude -DCHEMLARP_STATIC src/*.cpp tests/fuzz_parsers.cpp -o fuzz_parsers
//   ./fuzz_parsers corpus/ -max_total_time=300 -max_len=4096
//
// The `fuzz` job in .github/workflows/ci.yml runs it on every push for a short time. A crash is written to
// `crash-*`; replay it with `./fuzz_parsers crash-<hash>`.
#include <cstddef>
#include <cstdint>
#include <string>

#include "chemlarp/chemlarp_c.h"

namespace {

void consume(char* s) {
  if (s) cl_free_string(s);
}

// A receptor of one atom, so that fuzzing the docking entry point costs a parse and a one-step search
// rather than a real docking run.
const char* kTinyReceptor =
    "HETATM    1  C1  PKT A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n";

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  if (size < 2) return 0;
  const uint8_t which = data[0];
  // Keep inputs small: the point is parser edge cases, not how long a 100,000-atom molecule takes.
  const std::string in(reinterpret_cast<const char*>(data) + 1, size > 8193 ? 8192 : size - 1);
  switch (which % 16) {
    case 0: consume(cl_parse_smiles_json(in.c_str())); break;
    case 1: consume(cl_canonical_smiles_json(in.c_str(), 1, 1)); break;
    case 2: consume(cl_read_molecules_json(in.c_str(), nullptr, 0, 4)); break;
    case 3: consume(cl_protein_json(in.c_str(), 7.4)); break;
    case 4: consume(cl_smarts_match_json("CC(=O)Oc1ccccc1C(=O)O", in.c_str(), 0, 50)); break;
    case 5: consume(cl_check_reaction_json(in.c_str(), 298.15)); break;
    case 6: consume(cl_equilibrate_text_json(in.c_str(), 1000.0, 1e5)); break;
    case 7: consume(cl_aqueous_json(in.c_str())); break;
    case 8: consume(cl_flowsheet_json(in.c_str())); break;
    case 9: consume(cl_lab_mix_json(in.c_str())); break;
    case 10: consume(cl_lab_titration_json(in.c_str())); break;
    case 11: consume(cl_mbar_json(in.c_str())); break;
    case 12:
      // Mechanisms are the largest grammar: .ocm, CHEMKIN-II and Cantera YAML, detected from content.
      consume(cl_kinetics_text_json(in.c_str(), nullptr, nullptr, "iso-v", 1000.0, 1e5, "", 1e-9, 2, 1e-6));
      break;
    case 13: consume(cl_vle_text_json(in.c_str(), "txy", 1.0)); break;
    case 14:
      // Docking, with the input as the receptor: the PDB/mmCIF reader again, then the receptor typing, the
      // grid and the search, all of which convert coordinates to grid indices.
      consume(cl_dock_json(in.c_str(), "C", nullptr, 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
      break;
    case 15:
      // The input as the ligand (SMILES or MOL/SDF, detected) and as the reference, against a receptor of
      // one atom: the torsion tree, the isomorphism search and the pose code.
      consume(cl_dock_json(kTinyReceptor, in.c_str(), in.c_str(), 0, 0, 0, 8, 8, 8, 1, 1, 1, 7.4, 1, 0));
      break;
    default: break;
  }
  return 0;
}
