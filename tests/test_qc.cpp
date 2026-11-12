// Hartree-Fock and DFT single points against PySCF 2.15 on the same geometries
// (cart=True, so 6-31G* has six Cartesian d functions; DFT reference grid level 9).
#include <cmath>
#include <string>

#include "chemlarp/molio.hpp"
#include "chemlarp/qc.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

const char* kWater =
    "water\n  test 3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    "   -0.0023    0.3888   -0.0160 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7741   -0.2259   -0.0200 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.7914   -0.2084    0.0209 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\n  1  3  1  0\nM  END\n";
const char* kMethyl =
    "methyl\n  test 3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0067    0.0128   -0.0015 C   0  0  0  0  0  3  0  0  0  0  0  0\n"
    "    0.7426    0.8090    0.0181 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.0510    0.2521   -0.0071 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.3284   -1.0227   -0.0155 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\n  1  3  1  0\n  1  4  1  0\nM  END\n";
const char* kO2 =
    "dioxygen\n  test 3D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
    "   -0.5734   -0.0203   -0.0017 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.5775   -0.0173   -0.0071 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  2  0\nM  END\n";

QcResult run(const char* mol, QcMethod m, const char* basis, int mult = 0, Diagnostics* out = nullptr) {
  ChemicalSystemState st;
  Diagnostics d;
  ReadOptions ro;
  ro.keep_hydrogens = true;
  const auto recs = read_molecules(st, mol, d, ro);
  QcOptions o;
  o.method = m;
  o.basis = basis;
  o.multiplicity = mult;
  o.grid_level = 3;
  QcResult r = run_qc(st, static_cast<std::uint32_t>(recs.at(0).species), d, o);
  if (out) *out = d;
  return r;
}

}  // namespace

TEST(qc_hartree_fock_matches_pyscf) {
  QcResult r = run(kWater, QcMethod::HF, "STO-3G");
  CHECK(r.ok && r.converged);
  CHECK_EQ(r.n_basis, 7);
  CHECK_NEAR(r.energy, -74.9650023760, 1e-6);
  CHECK_NEAR(r.dipole[1], -1.68128, 1e-4);
  r = run(kWater, QcMethod::HF, "6-31G*");
  CHECK_EQ(r.n_basis, 19);  // six Cartesian d functions on O
  CHECK_NEAR(r.energy, -76.0071681369, 1e-6);
  CHECK_NEAR(r.dipole[1], -2.25082, 1e-4);
  // UHF doublet and triplet
  r = run(kMethyl, QcMethod::HF, "6-31G");
  CHECK(r.unrestricted && r.converged);
  CHECK_EQ(r.multiplicity, 2);
  CHECK_NEAR(r.energy, -39.5463173334, 1e-6);
  CHECK_NEAR(r.s_squared, 0.762257, 1e-4);
  r = run(kO2, QcMethod::HF, "STO-3G", 3);
  CHECK_NEAR(r.energy, -147.6223422010, 1e-6);
  CHECK_NEAR(r.s_squared, 2.003366, 1e-4);
}

TEST(qc_dft_matches_pyscf) {
  QcResult r = run(kWater, QcMethod::SVWN5, "6-31G");
  CHECK(r.converged);
  CHECK_NEAR(r.energy, -75.8195276647, 5e-6);
  CHECK_NEAR(r.electrons_on_grid, 10.0, 1e-4);
  r = run(kMethyl, QcMethod::PBE, "6-31G");
  CHECK_NEAR(r.energy, -39.7621291674, 5e-6);
  CHECK_NEAR(r.s_squared, 0.753161, 1e-4);
  r = run(kWater, QcMethod::B3LYP, "6-31G*");
  CHECK_NEAR(r.energy, -76.4080973193, 5e-6);
  CHECK_NEAR(r.dipole[1], -2.08689, 5e-4);
  r = run(kMethyl, QcMethod::B3LYP, "STO-3G");  // open shell: VWN-RPA with the f(zeta) interpolation
  CHECK_NEAR(r.energy, -39.3440930533, 5e-6);
}

TEST(qc_errors) {
  Diagnostics d;
  QcResult r = run(kWater, QcMethod::HF, "cc-pVTZ", 0, &d);
  CHECK(!r.ok);
  CHECK(d.has_errors());
  Diagnostics d2;
  r = run(kWater, QcMethod::HF, "STO-3G", 2, &d2);  // 10 electrons cannot be a doublet
  CHECK(!r.ok);
  CHECK(d2.has_errors());
  QcMethod m;
  CHECK(parse_qc_method("B3LYP", m) && m == QcMethod::B3LYP);
  CHECK(!parse_qc_method("mp2", m));
}
