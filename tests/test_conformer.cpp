#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "omnichem/cip.hpp"
#include "omnichem/conformer.hpp"
#include "omnichem/optimize.hpp"
#include "omnichem/smiles.hpp"
#include "omnichem/uff.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

struct Embedded {
  ChemicalSystemState st;
  Diagnostics d;
  EmbedResult r;
  std::int32_t sp = -1;
  const Species& s() const { return st.species.at(static_cast<std::size_t>(sp)); }
};

Embedded embed(const std::string& smi, std::uint64_t seed = 42) {
  Embedded e;
  EmbedOptions o;
  o.seed = seed;
  e.sp = embed_smiles(e.st, smi, e.d, &e.r, o);
  if (e.sp < 0) throw octest::Failure("embedding failed for " + smi + ": " + (e.d.empty() ? "" : e.d.items().back().message));
  return e;
}

double dist(const ChemicalSystemState& st, std::uint32_t a, std::uint32_t b) {
  const double dx = st.atoms.x[a] - st.atoms.x[b], dy = st.atoms.y[a] - st.atoms.y[b], dz = st.atoms.z[a] - st.atoms.z[b];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double angle_deg(const ChemicalSystemState& st, std::uint32_t a, std::uint32_t b, std::uint32_t c) {
  const double ux = st.atoms.x[a] - st.atoms.x[b], uy = st.atoms.y[a] - st.atoms.y[b], uz = st.atoms.z[a] - st.atoms.z[b];
  const double vx = st.atoms.x[c] - st.atoms.x[b], vy = st.atoms.y[c] - st.atoms.y[b], vz = st.atoms.z[c] - st.atoms.z[b];
  const double cs = (ux * vx + uy * vy + uz * vz) / std::sqrt((ux * ux + uy * uy + uz * uz) * (vx * vx + vy * vy + vz * vz));
  return std::acos(std::clamp(cs, -1.0, 1.0)) * 180.0 / 3.14159265358979323846;
}

// Mean length of bonds between elements za and zb (and orders) in the species.
double mean_bond(const Embedded& e, std::uint8_t za, std::uint8_t zb, int order) {
  double sum = 0;
  int n = 0;
  for (std::uint32_t b = e.s().bond_begin; b < e.s().bond_begin + e.s().bond_count; ++b) {
    const std::uint32_t i = e.st.bonds.begin[b], j = e.st.bonds.end[b];
    const std::uint8_t a = e.st.atoms.element[i], c = e.st.atoms.element[j];
    if (!((a == za && c == zb) || (a == zb && c == za))) continue;
    if (order > 0 && e.st.bonds.order[b] != order && !(order == 4 && (e.st.bonds.flags[b] & kBondAromatic))) continue;
    sum += dist(e.st, i, j);
    ++n;
  }
  if (n == 0) throw octest::Failure("no such bond");
  return sum / n;
}

// R/S straight from the coordinates: R iff V(lowest, p1, p2, p3) < 0, where V(a,b,c,d) = (b-a)·[(c-a)×(d-a)]
// and p1 > p2 > p3 are CIP priorities. Independent of the parity bookkeeping used by the embedder.
CipLabel geometric_cip(const ChemicalSystemState& st, std::uint32_t centre) {
  std::vector<std::uint32_t> nb(st.adjacency.neighbors(centre).begin(), st.adjacency.neighbors(centre).end());
  if (nb.size() != 4) return CipLabel::None;
  std::sort(nb.begin(), nb.end(), [&](std::uint32_t a, std::uint32_t b) { return cip_compare(st, centre, a, b) > 0; });
  auto P = [&](std::uint32_t i) { return std::array<double, 3>{st.atoms.x[i], st.atoms.y[i], st.atoms.z[i]}; };
  const auto L = P(nb[3]), a = P(nb[0]), b = P(nb[1]), c = P(nb[2]);
  const double v1[3] = {a[0] - L[0], a[1] - L[1], a[2] - L[2]}, v2[3] = {b[0] - L[0], b[1] - L[1], b[2] - L[2]},
               v3[3] = {c[0] - L[0], c[1] - L[1], c[2] - L[2]};
  const double V = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) + v1[1] * (v2[2] * v3[0] - v2[0] * v3[2]) +
                   v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
  return V < 0 ? CipLabel::R : CipLabel::S;
}

}  // namespace

TEST(lbfgs_minimises_rosenbrock) {
  std::vector<double> x = {-1.2, 1.0};
  LbfgsOptions o;
  o.grad_tol = 1e-8;
  o.max_step = 1.0;
  auto r = minimize_lbfgs([](const std::vector<double>& v, std::vector<double>& g) {
    const double a = 1 - v[0], b = v[1] - v[0] * v[0];
    g = {-2 * a - 400 * v[0] * b, 200 * b};
    return a * a + 100 * b * b;
  }, x, o);
  CHECK(r.converged);
  CHECK_NEAR(x[0], 1.0, 1e-6);
  CHECK_NEAR(x[1], 1.0, 1e-6);
}

TEST(uff_atom_typing) {
  auto types_of = [](const std::string& smi) {
    ChemicalSystemState st;
    Diagnostics d;
    const auto sp = parse_smiles(st, smi, d);
    std::vector<std::string> out;
    for (int t : assign_uff_types(st, static_cast<std::uint32_t>(sp), d)) out.push_back(t < 0 ? "?" : uff_parameters(t).label);
    return out;
  };
  CHECK(types_of("c1ccccc1")[0] == "C_R");
  CHECK((types_of("CC(C)=O") == std::vector<std::string>{"C_3", "C_2", "C_3", "O_2"}));
  CHECK((types_of("CC#N") == std::vector<std::string>{"C_3", "C_1", "N_1"}));
  CHECK((types_of("CC(N)=O") == std::vector<std::string>{"C_3", "C_2", "N_R", "O_2"}));  // amide N
  CHECK(types_of("N")[0] == "N_3");
  CHECK(types_of("c1ccncc1")[3] == "N_R");
  CHECK((types_of("CS(C)=O")[1] == "S_3+4"));
  CHECK((types_of("OS(=O)(=O)O")[1] == "S_3+6"));
  CHECK(types_of("C[Si](C)(C)C")[1] == "Si3");
  Diagnostics d;
  ChemicalSystemState st;
  const auto sp = parse_smiles(st, "[Fe]", d);
  assign_uff_types(st, static_cast<std::uint32_t>(sp), d);
  CHECK(d.contains(DiagCode::NoForceFieldType));
}

TEST(uff_rest_lengths) {
  const int C3 = uff_type_by_label("C_3"), CR = uff_type_by_label("C_R"), C2 = uff_type_by_label("C_2"),
            C1 = uff_type_by_label("C_1");
  CHECK_NEAR(uff_rest_length(C3, C3, 1), 1.514, 1e-12);
  CHECK_NEAR(uff_rest_length(CR, CR, 1.5), 1.458 - 0.1332 * 1.458 * std::log(1.5), 1e-12);
  CHECK_NEAR(uff_rest_length(C2, C2, 2), 1.464 - 0.1332 * 1.464 * std::log(2.0), 1e-12);
  CHECK_NEAR(uff_rest_length(C1, C1, 3), 1.412 - 0.1332 * 1.412 * std::log(3.0), 1e-12);
}

TEST(uff_gradient_matches_finite_differences) {
  for (const char* smi : {"CC", "C=C", "CC=C", "OO", "c1ccccc1", "CC(=O)N", "CCCC", "C[C@H](N)C(=O)O", "CS(C)=O"}) {
    ChemicalSystemState st;
    Diagnostics d;
    EmbedOptions o;
    o.minimize = false;
    const auto sp = embed_smiles(st, smi, d, nullptr, o);
    CHECK(sp >= 0);
    auto ff = UffForceField::build(st, static_cast<std::uint32_t>(sp), d);
    CHECK(ff.has_value());
    const Species& s = st.species[static_cast<std::size_t>(sp)];
    std::vector<double> x(3 * s.atom_count), g;
    for (std::uint32_t i = 0; i < s.atom_count; ++i) {
      x[3 * i] = st.atoms.x[s.atom_begin + i] + 0.05 * std::sin(1.7 * i);
      x[3 * i + 1] = st.atoms.y[s.atom_begin + i] + 0.05 * std::cos(2.3 * i);
      x[3 * i + 2] = st.atoms.z[s.atom_begin + i];
    }
    ff->energy_gradient(x, g);
    for (std::size_t i = 0; i < x.size(); ++i) {
      auto xp = x, xm = x;
      xp[i] += 1e-6;
      xm[i] -= 1e-6;
      const double fd = (ff->energy(xp) - ff->energy(xm)) / 2e-6;
      if (std::fabs(fd - g[i]) > 1e-5 * (std::fabs(fd) + 1.0))
        throw octest::Failure(std::string(smi) + ": coordinate " + std::to_string(i) + " analytic " +
                              std::to_string(g[i]) + " vs finite difference " + std::to_string(fd));
    }
  }
}

// Experimental reference geometries (gas phase, CRC Handbook). UFF is a generic force
// field: agreement within a few hundredths of an Å and a few degrees is expected.
TEST(uff_geometries_against_experiment) {
  CHECK_NEAR(mean_bond(embed("C"), 6, 1, 1), 1.087, 0.03);
  CHECK_NEAR(mean_bond(embed("CC"), 6, 6, 1), 1.535, 0.03);
  CHECK_NEAR(mean_bond(embed("C=C"), 6, 6, 2), 1.339, 0.03);
  CHECK_NEAR(mean_bond(embed("C#C"), 6, 6, 3), 1.203, 0.03);
  CHECK_NEAR(mean_bond(embed("c1ccccc1"), 6, 6, 4), 1.397, 0.03);
  CHECK_NEAR(mean_bond(embed("C=O"), 6, 8, 2), 1.205, 0.04);
  auto water = embed("O");
  CHECK_NEAR(mean_bond(water, 8, 1, 1), 0.958, 0.04);
  CHECK_NEAR(angle_deg(water.st, water.s().atom_begin + 1, water.s().atom_begin, water.s().atom_begin + 2), 104.5, 3.0);
  auto methane = embed("C");
  CHECK_NEAR(angle_deg(methane.st, methane.s().atom_begin + 1, methane.s().atom_begin, methane.s().atom_begin + 2), 109.47, 1.0);
}

TEST(embedding_corpus_is_clean) {
  for (const char* smi : {"CCO", "c1ccccc1", "C1CCCCC1", "CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                          "C1CC2CCC1CC2", "CC(C)(C)c1ccc(O)cc1", "C[S@](=O)CC", "[O-][N+](=O)c1ccccc1", "[NH4+].[Cl-]",
                          "CC(C)C[C@H](NC(=O)[C@@H](N)Cc1ccccc1)C(=O)O"}) {
    auto e = embed(smi);
    CHECK(e.r.ok);
    CHECK(e.r.uff.converged);
    CHECK_EQ(e.r.stereo_centres_ok, e.r.stereo_centres);
    const Species& s = e.s();
    // No clashes: atoms not within 1-3 of each other stay > 1.6 Å apart (H...H contacts included).
    std::vector<std::vector<int>> topo(s.atom_count, std::vector<int>(s.atom_count, 9));
    for (std::uint32_t i = 0; i < s.atom_count; ++i) topo[i][i] = 0;
    for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
      const std::uint32_t i = e.st.bonds.begin[b] - s.atom_begin, j = e.st.bonds.end[b] - s.atom_begin;
      topo[i][j] = topo[j][i] = 1;
    }
    for (std::uint32_t k = 0; k < s.atom_count; ++k)
      for (std::uint32_t i = 0; i < s.atom_count; ++i)
        for (std::uint32_t j = 0; j < s.atom_count; ++j) topo[i][j] = std::min(topo[i][j], topo[i][k] + topo[k][j]);
    for (std::uint32_t i = 0; i < s.atom_count; ++i)
      for (std::uint32_t j = i + 1; j < s.atom_count; ++j)
        if (topo[i][j] > 2 && dist(e.st, s.atom_begin + i, s.atom_begin + j) < 1.6)
          throw octest::Failure(std::string(smi) + ": atoms " + std::to_string(i) + " and " + std::to_string(j) +
                                " are only " + std::to_string(dist(e.st, s.atom_begin + i, s.atom_begin + j)) + " Å apart");
  }
}

TEST(embedding_reproduces_cip_labels_geometrically) {
  struct Case { const char* smiles; std::vector<CipLabel> expected; };
  const std::vector<Case> cases = {
      {"C[C@@H](C(=O)O)N", {CipLabel::S}},                  // L-alanine
      {"N[C@H](C)C(=O)O", {CipLabel::R}},                   // D-alanine
      {"CC[C@@H](C)O", {CipLabel::R}},                      // (R)-butan-2-ol
      {"F[C@H](Cl)Br", {CipLabel::R}},
      {"C[C@@H](O)[C@H](C)O", {CipLabel::R, CipLabel::S}},  // meso-butane-2,3-diol
      {"C[C@@H](O)[C@@H](C)O", {CipLabel::R, CipLabel::R}}, // (2R,3R)-butane-2,3-diol
  };
  for (const auto& c : cases) {
    for (std::uint64_t seed : {1u, 7u, 42u}) {
      auto e = embed(c.smiles, seed);
      std::vector<CipLabel> got;
      for (const auto& t : e.st.tetrahedral)
        if (e.st.atoms.species[t.atom] == static_cast<std::uint32_t>(e.sp)) got.push_back(geometric_cip(e.st, t.atom));
      if (got != c.expected) throw octest::Failure(std::string(c.smiles) + ": geometry does not give the expected R/S");
      for (const auto& t : e.st.tetrahedral)  // and the label from the parsed stereo agrees
        if (e.st.atoms.species[t.atom] == static_cast<std::uint32_t>(e.sp)) CHECK(cip_label(e.st, t) == geometric_cip(e.st, t.atom));
    }
  }
}

TEST(embedding_reproduces_double_bond_geometry) {
  auto trans = embed("C/C=C/C");
  auto cis = embed("C/C=C\\C");
  const auto& d1 = trans.st.double_bonds.back();
  const auto& d2 = cis.st.double_bonds.back();
  CHECK_NEAR(std::fabs(dihedral_degrees(trans.st, d1.ref_begin, trans.st.bonds.begin[d1.bond], trans.st.bonds.end[d1.bond], d1.ref_end)), 180.0, 5.0);
  CHECK_NEAR(std::fabs(dihedral_degrees(cis.st, d2.ref_begin, cis.st.bonds.begin[d2.bond], cis.st.bonds.end[d2.bond], d2.ref_end)), 0.0, 5.0);
  auto diene = embed("C/C=C/C=C\\C");
  CHECK_EQ(diene.r.stereo_bonds_ok, 2u);
}

TEST(cyclohexane_chair_and_planar_benzene) {
  auto chx = embed("C1CCCCC1");
  const std::uint32_t b = chx.s().atom_begin;
  double prev = 0;
  for (std::uint32_t k = 0; k < 6; ++k) {
    const double phi = dihedral_degrees(chx.st, b + k, b + (k + 1) % 6, b + (k + 2) % 6, b + (k + 3) % 6);
    CHECK(std::fabs(phi) > 45.0 && std::fabs(phi) < 65.0);  // chair: ±(50-60)°
    if (k > 0) CHECK(phi * prev < 0);                       // alternating signs
    prev = phi;
  }
  auto bz = embed("c1ccccc1");
  const std::uint32_t c = bz.s().atom_begin;
  for (std::uint32_t k = 0; k < 6; ++k)
    CHECK(std::fabs(dihedral_degrees(bz.st, c + k, c + (k + 1) % 6, c + (k + 2) % 6, c + (k + 3) % 6)) < 1.0);
}

TEST(explicit_hydrogens_preserve_identity) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = parse_smiles(st, "N[C@@H](C)C(=O)O", d);
  const auto sh = add_explicit_hydrogens(st, static_cast<std::uint32_t>(sp), d);
  CHECK(sh >= 0);
  const Species& a = st.species[static_cast<std::size_t>(sp)];
  const Species& b = st.species[static_cast<std::size_t>(sh)];
  CHECK_EQ(b.canonical_smiles, a.canonical_smiles);
  CHECK(b.composition == a.composition);
  CHECK_EQ(b.atom_count, a.atom_count + 7u);
  for (std::uint32_t g = b.atom_begin; g < b.atom_begin + b.atom_count; ++g) CHECK_EQ(st.atoms.implicit_h[g], 0);
  for (const auto& t : st.tetrahedral)
    if (st.atoms.species[t.atom] == static_cast<std::uint32_t>(sh))
      for (auto r : t.refs) CHECK(r != kImplicitHRef);
}

TEST(embedding_is_deterministic) {
  auto a = embed("CC(=O)Oc1ccccc1C(=O)O", 123), b = embed("CC(=O)Oc1ccccc1C(=O)O", 123);
  for (std::uint32_t i = 0; i < a.s().atom_count; ++i) {
    CHECK_EQ(a.st.atoms.x[a.s().atom_begin + i], b.st.atoms.x[b.s().atom_begin + i]);
    CHECK_EQ(a.st.atoms.z[a.s().atom_begin + i], b.st.atoms.z[b.s().atom_begin + i]);
  }
}

TEST(structure_output_formats) {
  auto e = embed("[O-][N+](=O)c1ccccc1");
  const std::string sdf = write_sdf(e.st, static_cast<std::uint32_t>(e.sp));
  CHECK(sdf.find(" 14 14  0  0  0  0  0  0  0  0999 V2000") != std::string::npos);  // C6H5NO2
  CHECK(sdf.find("M  CHG  2") != std::string::npos);
  CHECK(sdf.size() >= 12 && sdf.compare(sdf.size() - 12, 12, "M  END\n$$$$\n") == 0);
  const std::string xyz = write_xyz(e.st, static_cast<std::uint32_t>(e.sp));
  CHECK(xyz.rfind("14\n", 0) == 0);
}

TEST(conformer_error_vectors) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = parse_smiles(st, "CCO", d);
  embed_species(st, static_cast<std::uint32_t>(sp), d);  // still has implicit hydrogens
  CHECK(d.contains(DiagCode::ImplicitHydrogens));
  Diagnostics d2;
  CHECK(embed_smiles(st, "[Fe+2].[Cl-].[Cl-]", d2) < 0);
  CHECK(d2.contains(DiagCode::NoForceFieldType));
}
