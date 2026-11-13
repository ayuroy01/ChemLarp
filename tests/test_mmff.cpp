// MMFF94 / MMFF94s: atom types, partial charges, energy and gradient.
//
// data/mmff_reference.sdf holds twelve molecules (ChemLarp-generated coordinates) with the atom types,
// partial charges and MMFF94 / MMFF94s energies that RDKit's independent implementation computes for
// exactly those coordinates. The official MMFF94 validation suite (759 molecules against Merck's OPTIMOL)
// is run by scripts/validate_mmff.py; see "MMFF94: validation" in docs/ARCHITECTURE.md.
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "chemlarp/chemlarp_c.h"
#include "chemlarp/conformer.hpp"
#include "chemlarp/md.hpp"
#include "chemlarp/mmff.hpp"
#include "chemlarp/molio.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

std::vector<double> coords(const ChemicalSystemState& st, std::uint32_t sp) {
  const auto& s = st.species[sp];
  std::vector<double> x;
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    x.push_back(st.atoms.x[s.atom_begin + i]);
    x.push_back(st.atoms.y[s.atom_begin + i]);
    x.push_back(st.atoms.z[s.atom_begin + i]);
  }
  return x;
}

std::string prop(const MolRecord& r, const std::string& key) {
  for (const auto& [k, v] : r.properties)
    if (k == key) return v;
  return "";
}

template <class T>
std::vector<T> numbers(const std::string& s) {
  std::istringstream in(s);
  std::vector<T> v;
  for (T x; in >> x;) v.push_back(x);
  return v;
}

}  // namespace

TEST(mmff_matches_rdkit_reference_molecules) {
  std::ifstream f(cltest::data_path("mmff_reference.sdf"));
  CHECK(f.good());
  std::stringstream ss;
  ss << f.rdbuf();
  ChemicalSystemState st;
  Diagnostics d;
  ReadOptions ro;
  ro.keep_hydrogens = true;
  const auto recs = read_molecules(st, ss.str(), d, ro);
  CHECK_EQ(recs.size(), std::size_t{12});
  for (const auto& r : recs) {
    CHECK(r.species >= 0);
    const auto sp = std::uint32_t(r.species);
    const auto x = coords(st, sp);
    auto ff = MmffForceField::build(st, sp, d);
    CHECK(ff.has_value());
    const auto types = numbers<int>(prop(r, "RDKIT_MMFF94_TYPES"));
    const auto charges = numbers<double>(prop(r, "RDKIT_MMFF94_CHARGES"));
    CHECK_EQ(ff->types().size(), types.size());
    for (std::size_t i = 0; i < types.size(); ++i) {
      CHECK_EQ(ff->types()[i], types[i]);
      CHECK_NEAR(ff->charges()[i], charges[i], 6e-5);
    }
    CHECK_NEAR(ff->energy(x), std::stod(prop(r, "RDKIT_MMFF94_ENERGY")), 2e-4);
    auto ffs = MmffForceField::build(st, sp, d, MmffVariant::Mmff94s);
    CHECK(ffs.has_value());
    CHECK_NEAR(ffs->energy(x), std::stod(prop(r, "RDKIT_MMFF94S_ENERGY")), 2e-4);
    // Charges sum to the net formal charge.
    double q = 0.0, fq = 0.0;
    for (double c : ff->charges()) q += c;
    for (std::uint32_t i = 0; i < st.species[sp].atom_count; ++i) fq += st.atoms.formal_charge[st.species[sp].atom_begin + i];
    CHECK_NEAR(q, fq, 1e-9);
  }
}

TEST(mmff_gradient_matches_finite_differences) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = embed_smiles(st, "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", d);  // sulfamethoxazole
  CHECK(sp >= 0);
  auto ff = MmffForceField::build(st, std::uint32_t(sp), d);
  CHECK(ff.has_value());
  auto x = coords(st, std::uint32_t(sp));
  std::vector<double> g;
  ff->energy_gradient(x, g);
  const double h = 1e-5;
  for (std::size_t k = 0; k < x.size(); k += 7) {
    const double x0 = x[k];
    x[k] = x0 + h;
    const double ep = ff->energy(x);
    x[k] = x0 - h;
    const double em = ff->energy(x);
    x[k] = x0;
    CHECK_NEAR(g[k], (ep - em) / (2 * h), 1e-4 * std::max(1.0, std::fabs(g[k])));
  }
}

TEST(mmff_minimization_and_dynamics) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = embed_smiles(st, "CC(=O)Nc1ccc(O)cc1", d);  // paracetamol
  CHECK(sp >= 0);
  auto r = mmff_minimize(st, std::uint32_t(sp), d, MmffVariant::Mmff94, LbfgsOptions{5000, 10, 1e-5, 1e-15, 0.3});
  CHECK(r.ok && r.converged);
  CHECK(r.energy_final <= r.energy_initial);
  CHECK(r.grad_max < 1e-5);
  // Constant-energy MMFF94 dynamics conserves energy.
  MdOptions o;
  o.thermostat = Thermostat::None;
  o.steps = 2000;
  o.dt_fs = 0.5;
  o.sample_every = 10;
  auto md = molecule_md("CC(=O)Nc1ccc(O)cc1", o, false, d, ForceFieldChoice::Mmff94);
  CHECK(md.has_value());
  CHECK(md->md.fluctuation < 0.05);
  CHECK(std::fabs(md->md.drift_kcal_mol_ps) < 0.05);
}

TEST(mmff_errors_and_c_api) {
  ChemicalSystemState st;
  Diagnostics d;
  // An element MMFF94 does not cover.
  const std::int32_t sp = embed_smiles(st, "C[Se]C", d);
  if (sp >= 0) {
    Diagnostics e;
    CHECK(!MmffForceField::build(st, std::uint32_t(sp), e).has_value());
    CHECK(e.has_errors());
  }
  char* ok = cl_mmff_json("c1ccccc1O", "MMFF94", 1);
  const std::string s(ok);
  cl_free_string(ok);
  CHECK(s.rfind("{\"ok\":true", 0) == 0);
  CHECK(s.find("\"types\":[37,37,37,37,37,37,6") != std::string::npos);
  CHECK(s.find("\"electrostatic\":") != std::string::npos);
}
