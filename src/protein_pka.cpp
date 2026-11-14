// Structure-based pKa prediction: geometric features with weights fitted to measured values.
#include <cstring>
#include <set>

#include "protein_internal.hpp"

namespace chemlarp {

using namespace detail;

namespace {

// Group order of kPkaOffset; weight classes: 0 carboxyl (ASP GLU C-), 1 HIS, 2 amine (LYS N+), 3 TYR CYS.
const char* const kGroups[9] = {"ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG", "N+", "C-"};
const int kClass[9] = {0, 0, 1, 3, 3, 2, -1, 2, 0};
const bool kAcid[9] = {true, true, false, true, true, false, false, false, true};
#include "protein_pka.inc"

struct Site {
  std::size_t res;
  int group;
  std::vector<Vec> atoms;  // titratable atoms
  Vec centre;
  double burial = 0.0;
};

double ramp(double d, double d1, double d2) { return d <= d1 ? 1.0 : d >= d2 ? 0.0 : (d2 - d) / (d2 - d1); }

int metal_charge(const std::string& e) {
  if (e == "Zn" || e == "Ca" || e == "Mg" || e == "Mn" || e == "Fe" || e == "Cu" || e == "Co" || e == "Ni" || e == "Cd") return 2;
  if (e == "Na" || e == "K") return 1;
  return 0;
}

}  // namespace

std::vector<PkaSite> predict_pka(const Protein& p) {
  const Pka K;
  const double ref[9] = {K.asp, K.glu, K.his, K.cys, K.tyr, K.lys, K.arg, K.nterm, K.cterm};
  const auto& defs = residue_defs();
  std::set<std::size_t> in_ss;
  for (const auto& [a, b] : p.disulfides) in_ss.insert(a), in_ss.insert(b);
  // Sites.
  std::vector<Site> sites;
  auto add_site = [&](std::size_t ri, int g, std::initializer_list<const char*> names, std::initializer_list<const char*> centre_names) {
    Site s;
    s.res = ri;
    s.group = g;
    for (const char* n : names)
      if (const ProteinAtom* a = p.residues[ri].atom(n)) s.atoms.push_back(pos(*a));
    if (s.atoms.empty()) return;
    Vec c{0, 0, 0};
    int k = 0;
    for (const char* n : centre_names)
      if (const ProteinAtom* a = p.residues[ri].atom(n)) c = add(c, pos(*a)), ++k;
    if (!k) {
      for (const auto& x : s.atoms) c = add(c, x), ++k;
    }
    s.centre = mul(c, 1.0 / k);
    sites.push_back(std::move(s));
  };
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    const auto& r = p.residues[i];
    const std::string base = canonical_name(r.name);
    if (r.hetero || !defs.count(base)) continue;
    if (r.n_terminal) add_site(i, 7, {"N"}, {});
    if (base == "ASP") add_site(i, 0, {"OD1", "OD2"}, {});
    else if (base == "GLU") add_site(i, 1, {"OE1", "OE2"}, {});
    else if (base == "HIS") add_site(i, 2, {"ND1", "NE2"}, {"CG", "ND1", "CD2", "CE1", "NE2"});
    else if (base == "CYS" && !in_ss.count(i)) add_site(i, 3, {"SG"}, {});
    else if (base == "TYR") add_site(i, 4, {"OH"}, {});
    else if (base == "LYS") add_site(i, 5, {"NZ"}, {});
    else if (base == "ARG") add_site(i, 6, {"NE", "NH1", "NH2"}, {"CZ"});
    if (r.c_terminal) add_site(i, 8, {"O", "OXT"}, {"C"});
  }
  // Atoms: all heavy atoms except water (burial); polar partners by type.
  struct Polar {
    std::uint32_t res;
    int type;  // 0 backbone N, 1 backbone O, 2 side-chain NH, 3 hydroxyl, 4 amide O
  };
  std::vector<std::uint32_t> heavy_res;
  std::vector<Polar> polar;
  struct Ion {
    Vec x;
    int q;
  };
  std::vector<Ion> ions;
  Grid heavy(6.0), pgrid(4.0);
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    const auto& r = p.residues[i];
    if (r.name == "HOH" || r.name == "WAT" || r.name == "DOD") continue;
    const std::string base = canonical_name(r.name);
    const bool std_res = !r.hetero && defs.count(base);
    for (const auto& a : r.atoms) {
      if (a.element == "H" || a.element == "D") continue;
      heavy.add(pos(a), std::uint32_t(heavy_res.size()));
      heavy_res.push_back(std::uint32_t(i));
      if (r.hetero) {
        if (const int q = metal_charge(a.element)) ions.push_back({pos(a), q});
        continue;
      }
      if (!std_res) continue;
      int type = -1;
      if (a.name == "N" && base != "PRO" && !r.n_terminal) type = 0;
      else if (a.name == "O" && !r.c_terminal) type = 1;
      else if ((base == "ASN" && a.name == "ND2") || (base == "GLN" && a.name == "NE2") || (base == "TRP" && a.name == "NE1")) type = 2;
      else if ((base == "SER" && a.name == "OG") || (base == "THR" && a.name == "OG1")) type = 3;
      else if ((base == "ASN" && a.name == "OD1") || (base == "GLN" && a.name == "OE1")) type = 4;
      if (type >= 0) {
        pgrid.add(pos(a), std::uint32_t(polar.size()));
        polar.push_back({std::uint32_t(i), type});
      }
    }
  }
  Grid sgrid(10.0);
  for (std::size_t k = 0; k < sites.size(); ++k) sgrid.add(sites[k].centre, std::uint32_t(k));
  std::vector<PkaSite> out(sites.size());
  std::vector<double> local(sites.size());
  for (std::size_t k = 0; k < sites.size(); ++k) {
    Site& s = sites[k];
    int n15 = 0, n6 = 0;
    heavy.query(s.centre, 15.5, [&](std::uint32_t id, Vec, double d2) {
      if (heavy_res[id] == s.res) return;
      ++n15;
      if (d2 < 36.0) ++n6;
    });
    s.burial = std::clamp((n15 - 150.0) / 400.0, 0.0, 1.0);
    local[k] = n6 / 20.0;
  }
  for (std::size_t k = 0; k < sites.size(); ++k) {
    const Site& s = sites[k];
    PkaSite& o = out[k];
    o.residue = s.res;
    o.group = kGroups[s.group];
    o.reference = ref[s.group];
    o.buried = s.burial;
    double* f = o.features;
    f[0] = s.burial;
    f[1] = local[k];
    // Hydrogen-bond partners: nearest titratable atom to each polar atom.
    std::map<std::uint32_t, double> nearest;
    for (const Vec& t : s.atoms)
      pgrid.query(t, 4.0, [&](std::uint32_t id, Vec, double d2) {
        if (polar[id].res == s.res) return;
        auto [it, fresh] = nearest.emplace(id, d2);
        if (!fresh) it->second = std::min(it->second, d2);
      });
    for (const auto& [id, d2] : nearest) f[2 + polar[id].type] += ramp(std::sqrt(d2), 3.0, 4.0);
    // Charges of the other groups at the pH where this one titrates, and metal ions.
    auto coulomb = [&](double q, double r, double bij) {
      if (r >= 10.0) return;
      const double w = q * (1.0 / std::max(r, 3.0) - 0.1);
      f[7] += w * (1.0 - bij);
      f[8] += w * bij;
    };
    sgrid.query(s.centre, 16.0, [&](std::uint32_t j, Vec, double) {
      if (j == k) return;
      const Site& t = sites[j];
      if (t.res == s.res && !(s.group >= 7 || t.group >= 7)) return;
      double r = 1e9;
      for (const Vec& a : s.atoms)
        for (const Vec& b : t.atoms) r = std::min(r, norm(sub(a, b)));
      const double dpk = ref[t.group] - ref[s.group];
      const double q = kAcid[t.group] ? -1.0 / (1.0 + std::pow(10.0, dpk)) : 1.0 / (1.0 + std::pow(10.0, -dpk));
      coulomb(q, r, 0.5 * (s.burial + t.burial));
      f[9] += q * ramp(r, 3.0, 4.0);  // salt bridge or hydrogen bond between titratable groups
    });
    for (const auto& ion : ions) {
      double r = 1e9;
      for (const Vec& a : s.atoms) r = std::min(r, norm(sub(a, ion.x)));
      coulomb(ion.q, r, s.burial);
      f[9] += ion.q * ramp(r, 2.2, 3.2);
    }
    const int c = kClass[s.group];
    if (c >= 0) {
      const double* w = kPkaWeights[c];
      o.desolvation = w[0] * f[0] + w[1] * f[1];
      for (int i = 2; i < 7; ++i) o.hydrogen_bonds += w[i] * f[i];
      o.charges = w[7] * f[7] + w[8] * f[8] + w[9] * f[9];
      o.reference += kPkaOffset[s.group];
    }
    o.pka = o.reference + o.desolvation + o.hydrogen_bonds + o.charges;
  }
  return out;
}

}  // namespace chemlarp
