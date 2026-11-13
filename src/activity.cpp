#include "chemlarp/activity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "chemlarp/smarts.hpp"
#include "chemlarp/smiles.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {
namespace {

struct UnifacSubgroup {
  int id, main;
  double R, Q;
  const char* name;
};
struct UnifacInteraction {
  int m, n;
  double a;
};
#include "unifac_data.inc"

const std::map<int, UnifacSubgroupInfo>& subgroup_table() {
  static const std::map<int, UnifacSubgroupInfo> t = [] {
    std::map<int, UnifacSubgroupInfo> m;
    for (const auto& g : kUnifacSubgroups) m[g.id] = {g.id, g.main, g.R, g.Q, g.name};
    return m;
  }();
  return t;
}

const std::map<std::pair<int, int>, double>& interaction_table() {
  static const std::map<std::pair<int, int>, double> t = [] {
    std::map<std::pair<int, int>, double> m;
    for (const auto& p : kUnifacInteractions) m[{p.m, p.n}] = p.a;
    return m;
  }();
  return t;
}

// Mole fractions: finite, non-negative, positive sum; returned normalised.
std::optional<std::vector<double>> normalise(const std::vector<double>& x, std::size_t n, Diagnostics& diag) {
  if (x.size() != n || n == 0) {
    diag.error(DiagCode::InvalidMixtureInput, "expected " + std::to_string(n) + " mole fractions, got " + std::to_string(x.size()));
    return std::nullopt;
  }
  double s = 0.0;
  for (double v : x) {
    if (!std::isfinite(v) || v < 0.0) {
      diag.error(DiagCode::InvalidMixtureInput, "mole fractions must be finite and non-negative");
      return std::nullopt;
    }
    s += v;
  }
  if (!(s > 0.0)) {
    diag.error(DiagCode::InvalidMixtureInput, "mole fractions sum to zero");
    return std::nullopt;
  }
  std::vector<double> out(x);
  for (double& v : out) v /= s;
  return out;
}

bool valid_temperature(double T, Diagnostics& diag) {
  if (std::isfinite(T) && T > 0.0) return true;
  diag.error(DiagCode::InvalidMixtureInput, "temperature must be a positive number of kelvin");
  return false;
}

double at(const std::vector<double>& m, std::size_t n, std::size_t i, std::size_t j) { return m.empty() ? 0.0 : m[i * n + j]; }

}  // namespace

std::optional<std::vector<double>> nrtl_ln_gamma(const NrtlParameters& p, double T, const std::vector<double>& x_in,
                                                 Diagnostics& diag) {
  const std::size_t n = p.n;
  if (!valid_temperature(T, diag)) return std::nullopt;
  for (const auto* m : {&p.a, &p.b, &p.e, &p.f, &p.alpha_c, &p.alpha_d})
    if (!m->empty() && m->size() != n * n) {
      diag.error(DiagCode::InvalidMixtureInput, "NRTL parameter matrices must be n*n");
      return std::nullopt;
    }
  if (p.alpha_c.empty() && n > 1) {
    diag.error(DiagCode::InvalidMixtureInput, "NRTL non-randomness parameters alpha are required");
    return std::nullopt;
  }
  auto x = normalise(x_in, n, diag);
  if (!x) return std::nullopt;
  std::vector<double> tau(n * n, 0.0), G(n * n, 1.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j) continue;
      const double t = at(p.a, n, i, j) + at(p.b, n, i, j) / T + at(p.e, n, i, j) * std::log(T) + at(p.f, n, i, j) * T;
      const double alpha = at(p.alpha_c, n, i, j) + at(p.alpha_d, n, i, j) * (T - 273.15);
      if (!std::isfinite(t) || !std::isfinite(alpha)) {
        diag.error(DiagCode::InvalidMixtureInput, "NRTL parameters give a non-finite tau or alpha");
        return std::nullopt;
      }
      tau[i * n + j] = t;
      G[i * n + j] = std::exp(-alpha * t);
    }
  // S_j = Σ_k x_k G_kj, C_j = Σ_m x_m τ_mj G_mj
  std::vector<double> S(n, 0.0), C(n, 0.0);
  for (std::size_t j = 0; j < n; ++j)
    for (std::size_t k = 0; k < n; ++k) {
      S[j] += (*x)[k] * G[k * n + j];
      C[j] += (*x)[k] * tau[k * n + j] * G[k * n + j];
    }
  std::vector<double> lng(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double v = C[i] / S[i];
    for (std::size_t j = 0; j < n; ++j) v += (*x)[j] * G[i * n + j] / S[j] * (tau[i * n + j] - C[j] / S[j]);
    lng[i] = v;
  }
  return lng;
}

std::optional<UnifacSubgroupInfo> unifac_subgroup(int id) {
  auto it = subgroup_table().find(id);
  if (it == subgroup_table().end()) return std::nullopt;
  return it->second;
}

std::optional<double> unifac_interaction(int m, int n) {
  if (m == n) return 0.0;
  auto it = interaction_table().find({m, n});
  if (it == interaction_table().end()) return std::nullopt;
  return it->second;
}

std::optional<std::vector<double>> unifac_ln_gamma(const std::vector<UnifacGroups>& comps, double T,
                                                   const std::vector<double>& x_in, Diagnostics& diag) {
  const std::size_t nc = comps.size();
  if (!valid_temperature(T, diag)) return std::nullopt;
  auto xo = normalise(x_in, nc, diag);
  if (!xo) return std::nullopt;
  const auto& x = *xo;
  // the subgroups present, in a fixed order
  std::vector<int> ids;
  for (const auto& c : comps) {
    if (c.empty()) {
      diag.error(DiagCode::InvalidMixtureInput, "a UNIFAC component has no subgroups");
      return std::nullopt;
    }
    for (const auto& [id, cnt] : c) {
      if (!unifac_subgroup(id)) {
        diag.error(DiagCode::InvalidMixtureInput, "unknown UNIFAC subgroup " + std::to_string(id));
        return std::nullopt;
      }
      if (cnt <= 0) {
        diag.error(DiagCode::InvalidMixtureInput, "UNIFAC subgroup counts must be positive");
        return std::nullopt;
      }
      if (std::find(ids.begin(), ids.end(), id) == ids.end()) ids.push_back(id);
    }
  }
  std::sort(ids.begin(), ids.end());
  const std::size_t ng = ids.size();
  std::vector<double> R(ng), Q(ng);
  std::vector<int> mg(ng);
  for (std::size_t k = 0; k < ng; ++k) {
    const auto g = *unifac_subgroup(ids[k]);
    R[k] = g.R, Q[k] = g.Q, mg[k] = g.main_group;
  }
  // Ψ_mn between the subgroups' main groups; a missing published pair is an error
  std::vector<double> Psi(ng * ng, 1.0);
  for (std::size_t a = 0; a < ng; ++a)
    for (std::size_t b = 0; b < ng; ++b) {
      auto amn = unifac_interaction(mg[a], mg[b]);
      if (!amn) {
        diag.error(DiagCode::NoInteractionParameter, "no published UNIFAC interaction parameter between main groups " +
                                                         std::to_string(mg[a]) + " and " + std::to_string(mg[b]));
        return std::nullopt;
      }
      Psi[a * ng + b] = std::exp(-*amn / T);
    }
  // ν_k^(i)
  std::vector<std::vector<double>> nu(nc, std::vector<double>(ng, 0.0));
  for (std::size_t i = 0; i < nc; ++i)
    for (const auto& [id, cnt] : comps[i]) {
      const std::size_t k = static_cast<std::size_t>(std::find(ids.begin(), ids.end(), id) - ids.begin());
      nu[i][k] += cnt;
    }
  // ln Γ_k for a group mixture with group mole numbers v
  auto ln_Gamma = [&](const std::vector<double>& v) {
    double tot = 0.0;
    for (double w : v) tot += w;
    std::vector<double> X(ng), Theta(ng);
    double sQX = 0.0;
    for (std::size_t k = 0; k < ng; ++k) X[k] = v[k] / tot, sQX += Q[k] * X[k];
    for (std::size_t k = 0; k < ng; ++k) Theta[k] = Q[k] * X[k] / sQX;
    std::vector<double> sum_m(ng, 0.0);  // Σ_m Θ_m Ψ_mk
    for (std::size_t k = 0; k < ng; ++k)
      for (std::size_t m = 0; m < ng; ++m) sum_m[k] += Theta[m] * Psi[m * ng + k];
    std::vector<double> out(ng, 0.0);
    for (std::size_t k = 0; k < ng; ++k) {
      double s2 = 0.0;
      for (std::size_t m = 0; m < ng; ++m) s2 += Theta[m] * Psi[k * ng + m] / sum_m[m];
      out[k] = Q[k] * (1.0 - std::log(sum_m[k]) - s2);
    }
    return out;
  };
  // combinatorial part
  std::vector<double> r(nc, 0.0), q(nc, 0.0);
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t k = 0; k < ng; ++k) r[i] += nu[i][k] * R[k], q[i] += nu[i][k] * Q[k];
  double sr = 0.0, sq = 0.0, sxl = 0.0;
  std::vector<double> l(nc);
  for (std::size_t i = 0; i < nc; ++i) {
    sr += x[i] * r[i], sq += x[i] * q[i];
    l[i] = 5.0 * (r[i] - q[i]) - (r[i] - 1.0);
    sxl += x[i] * l[i];
  }
  // residual part
  std::vector<double> vmix(ng, 0.0);
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t k = 0; k < ng; ++k) vmix[k] += x[i] * nu[i][k];
  const auto lnG_mix = ln_Gamma(vmix);
  std::vector<double> lng(nc, 0.0);
  for (std::size_t i = 0; i < nc; ++i) {
    // Φ_i/x_i = r_i / Σ x r and θ_i/Φ_i = (q_i/Σ x q) / (r_i/Σ x r): finite also at x_i = 0 (infinite dilution)
    const double phi_x = r[i] / sr, theta_phi = (q[i] / sq) / (r[i] / sr);
    double lc = std::log(phi_x) + 5.0 * q[i] * std::log(theta_phi) + l[i] - phi_x * sxl;
    const auto lnG_pure = ln_Gamma(nu[i]);
    double lr = 0.0;
    for (std::size_t k = 0; k < ng; ++k)
      if (nu[i][k] > 0.0) lr += nu[i][k] * (lnG_mix[k] - lnG_pure[k]);
    lng[i] = lc + lr;
  }
  return lng;
}

namespace {

// Subgroup patterns in priority order (larger and more specific groups first, as in the DDBST assignments).
// The atoms a pattern matches are the group's heavy atoms; $(...) only constrains the surroundings.
// `whole`: the group is the entire molecule (water, methanol, CS2, ...).
struct GroupRule {
  int id;
  const char* smarts;
  bool whole;
  int pass = 0;  // 0: always; 1: only in the first ordering; 2: only in the fallback ordering
};
const GroupRule kGroupRules[] = {
    {16, "[OH2]", true}, {15, "[CH3][OH]", true}, {58, "S=C=S", true}, {43, "[CH1](=O)[OH]", true},
    {52, "ClC(Cl)(Cl)Cl", true}, {50, "Cl[CH1](Cl)Cl", true}, {47, "Cl[CH2]Cl", true},
    {61, "O=[CH1]c1ccco1", true}, {62, "[OH][CH2][CH2][OH]", true}, {67, "[CH3][SX3](=O)[CH3]", true},
    {68, "[CH2]=[CH1]C#N", true}, {72, "[CH3]N([CH3])[CH1]=O", true}, {85, "[CH3]N1[CH2][CH2][CH2]C1=O", true},
    {105, "[CH2]1[CH2]O[CH2][CH2][NH1]1", true},
    // nitro, nitriles
    {57, "c[N+](=O)[O-]", false}, {54, "[CH3][N+](=O)[O-]", false}, {55, "[CH2][N+](=O)[O-]", false}, {56, "[CH1][N+](=O)[O-]", false},
    {40, "[CH3]C#N", false}, {41, "[CH2]C#N", false},
    // esters, acids, ketones, aldehydes
    {42, "[CX3](=O)[OH]", false},
    {21, "[CH3][CX3](=O)[OX2H0;$(O[#6])]", false}, {22, "[CH2][CX3](=O)[OX2H0;$(O[#6])]", false},
    {23, "[CH1](=O)[OX2H0;$(O[#6])]", false}, {77, "[CX3](=O)[OX2H0;$(O[#6])]", false},
    {18, "[CH3][CX3H0;$(C(=O)([#6])[#6])]=O", false}, {19, "[CH2][CX3H0;$(C(=O)([#6])[#6])]=O", false},
    {20, "[CX3H1;$(C(=O)[#6])]=O", false},
    // hydroxy ethers: HO-CH2-CH2-O- (C2H5O2), HO-CH-CH2-O- (C2H4O2)
    {100, "[OH][CH2;X4;!R][CH2;X4;!R][OX2H0;!R;$(O([#6])[#6])]", false},
    {101, "[OH][CH1;X4;!R][CH2;X4;!R][OX2H0;!R;$(O([#6])[#6])]", false},
    {101, "[OH][CH2;X4;!R][CH1;X4;!R][OX2H0;!R;$(O([#6])[#6])]", false},
    // sulfides take benzylic carbons before ACCH2 does (DDBST: Ph-CH2-S- is AC + CH2S)
    {102, "[CH3][SX2;$(S([#6])[#6])]", false}, {103, "[CH2;X4][SX2;$(S([#6])[#6])]", false}, {104, "[CH1;X4][SX2;$(S([#6])[#6])]", false},
    // aromatic carbon with an alkyl carbon: first ordering before ethers and amines (benzyl ethers: ACCH2),
    // fallback after them when the first ordering leaves an atom without a group (Ph-O-CH2-Ph)
    {11, "c[CH3]", false, 1}, {12, "c[CH2;X4]", false, 1}, {13, "c[CH1;X4]", false, 1},
    // ethers (cyclic CH2-O as THF)
    {27, "[CH2;R][OX2;R;$(O([#6])[#6])]", false},
    {24, "[CH3][OX2;$(O([#6])[#6])]", false}, {25, "[CH2][OX2;$(O([#6])[#6])]", false}, {26, "[CH1][OX2;$(O([#6])[#6])]", false},
    // amines
    {28, "[CH3][NH2]", false}, {29, "[CH2;X4][NH2]", false}, {30, "[CH1;X4][NH2]", false},
    {31, "[CH3][NH1;$(N([CX4])[CX4])]", false}, {32, "[CH2;X4][NH1;$(N([CX4])[CX4])]", false}, {33, "[CH1;X4][NH1;$(N([CX4])[CX4])]", false},
    {34, "[CH3][NX3H0;$(N([CX4])([CX4])[CX4])]", false}, {35, "[CH2;X4][NX3H0;$(N([CX4])([CX4])[CX4])]", false},
    {36, "c[NH2]", false}, {17, "c[OH]", false}, {14, "[OH;$(O[CX4])]", false},
    // halogens
    {51, "[CX4H0](Cl)(Cl)Cl", false}, {48, "[CH1;X4](Cl)Cl", false}, {49, "[CX4H0](Cl)Cl", false},
    {44, "[CH2;X4]Cl", false}, {45, "[CH1;X4]Cl", false}, {46, "[CX4H0]Cl", false}, {53, "cCl", false},
    {69, "[Cl;$(Cl[CX3]=[CX3])]", false}, {63, "[I;$(I[#6])]", false}, {64, "[Br;$(Br[#6])]", false},
    {71, "cF", false}, {74, "[CX4](F)(F)F", false}, {75, "[CX4](F)F", false}, {76, "[CX4]F", false},
    // sulfur
    {59, "[CH3][SH1]", false}, {60, "[CH2;X4][SH1]", false},
    // alkynes, alkenes
    {65, "[CH1]#C", false}, {66, "C#C", false},
    {5, "[CH2]=[CH1]", false}, {6, "[CH1]=[CH1]", false}, {7, "[CH2]=[CX3H0]", false}, {8, "[CH1]=[CX3H0]", false},
    {70, "[CX3H0]=[CX3H0]", false},
    // aromatic carbon with alkyl, aromatic carbon, alkyl
    {11, "c[CH3]", false, 2}, {12, "c[CH2;X4]", false, 2}, {13, "c[CH1;X4]", false, 2},
    {9, "[cH]", false}, {10, "[cH0]", false},
    {1, "[CH3;X4]", false}, {2, "[CH2;X4]", false}, {3, "[CH1;X4]", false}, {4, "[CX4H0]", false},
};

const std::vector<std::pair<GroupRule, SmartsQuery>>& compiled_rules() {
  static const std::vector<std::pair<GroupRule, SmartsQuery>> v = [] {
    std::vector<std::pair<GroupRule, SmartsQuery>> out;
    for (const auto& r : kGroupRules) {
      Diagnostics d;
      out.push_back({r, parse_smarts(r.smarts, d)});
    }
    return out;
  }();
  return v;
}

}  // namespace

std::optional<UnifacGroups> unifac_groups_from_smiles(const std::string& smiles, Diagnostics& diag) {
  ChemicalSystemState st;
  const std::int32_t sp = parse_smiles(st, smiles, diag);
  if (sp < 0) return std::nullopt;
  const Species& s = st.species[static_cast<std::size_t>(sp)];
  std::vector<std::uint32_t> heavy;
  for (std::uint32_t i = 0; i < s.atom_count; ++i)
    if (st.atoms.element[s.atom_begin + i] != 1) heavy.push_back(s.atom_begin + i);
  // heavy-atom count of each connected component (a whole-molecule group must cover its component)
  std::map<std::uint32_t, int> comp_of;
  std::map<int, int> comp_size;
  {
    int c = 0;
    for (auto g0 : heavy) {
      if (comp_of.count(g0)) continue;
      std::vector<std::uint32_t> stack = {g0};
      comp_of[g0] = c;
      while (!stack.empty()) {
        const auto g = stack.back();
        stack.pop_back();
        ++comp_size[c];
        for (auto b : st.adjacency.neighbors(g))
          if (st.atoms.element[b] != 1 && !comp_of.count(b)) comp_of[b] = c, stack.push_back(b);
      }
      ++c;
    }
  }
  // one ordering of the rules; false when an atom is left without a group
  auto attempt = [&](int pass, std::map<int, int>& count, std::uint32_t& uncovered) {
    std::map<std::uint32_t, bool> used;
    for (auto g : heavy) used[g] = false;
    // pyridine and thiophene rings: one subgroup per ring, by the hydrogens on its carbons
    {
      Diagnostics d;
      static const SmartsQuery pyr = parse_smarts("n1ccccc1", d), thio = parse_smarts("s1cccc1", d);
      for (const auto& [q, ids, n_atoms] : {std::tuple{&pyr, std::array<int, 3>{37, 38, 39}, 6}, std::tuple{&thio, std::array<int, 3>{106, 107, 108}, 5}}) {
        for (const auto& m : match_smarts(st, static_cast<std::uint32_t>(sp), *q)) {
          if (std::any_of(m.begin(), m.end(), [&](std::uint32_t at) { return used[at]; })) continue;
          int h = 0;
          for (auto at : m)
            if (st.atoms.element[at] == 6) h += st.atoms.implicit_h[at];
          const int nc = n_atoms - 1, k = nc - h;  // substituted ring carbons
          if (k < 0 || k > 2) continue;
          for (auto at : m) used[at] = true;
          ++count[ids[static_cast<std::size_t>(k)]];
        }
      }
    }
    for (const auto& [rule, q] : compiled_rules()) {
      if (!q.valid() || (rule.pass != 0 && rule.pass != pass)) continue;
      for (const auto& m : match_smarts(st, static_cast<std::uint32_t>(sp), q)) {
        if (rule.whole && static_cast<int>(m.size()) != comp_size[comp_of[m[0]]]) continue;
        if (std::any_of(m.begin(), m.end(), [&](std::uint32_t at) { return used[at]; })) continue;
        for (auto at : m) used[at] = true;
        ++count[rule.id];
      }
    }
    for (auto g : heavy)
      if (!used[g]) {
        uncovered = g;
        return false;
      }
    return true;
  };
  std::map<int, int> count;
  std::uint32_t uncovered = 0;
  if (!attempt(1, count, uncovered)) {
    count.clear();
    if (!attempt(2, count, uncovered)) {
      diag.error(DiagCode::NoUnifacGroups,
                 "no supported UNIFAC subgroup covers atom " + std::to_string(uncovered - s.atom_begin + 1) + " of " + smiles);
      return std::nullopt;
    }
  }
  UnifacGroups out(count.begin(), count.end());
  return out;
}


// ---- vapour pressure and vapour-liquid equilibrium ------------------------------------------------------

namespace {
struct AntoineRow {
  const char* smiles;
  const char* name;
  double A, B, C, Tmin, Tmax;
};
#include "antoine_data.inc"

std::string canonical(const std::string& smiles) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  if (sp < 0) return {};
  return st.species[static_cast<std::size_t>(sp)].canonical_smiles;
}

const std::map<std::string, AntoineCoefficients>& antoine_table() {
  static const std::map<std::string, AntoineCoefficients> t = [] {
    std::map<std::string, AntoineCoefficients> m;
    for (const auto& r : kAntoine) {
      const std::string key = canonical(r.smiles);
      if (key.empty() || m.count(key)) continue;  // first entry wins for structures listed twice
      m[key] = AntoineCoefficients{r.name, r.A, r.B, r.C, r.Tmin, r.Tmax};
    }
    return m;
  }();
  return t;
}

// ln P^sat_i(T) for every component; flags use outside the fitted range
std::vector<double> ln_psat(const ActivitySystem& sys, double T, bool& extrapolated) {
  std::vector<double> v(sys.size());
  for (std::size_t i = 0; i < sys.size(); ++i) {
    const auto& c = sys.antoine[i];
    v[i] = (c.A - c.B / (T + c.C)) * std::log(10.0);
    if (T < c.Tmin_K || T > c.Tmax_K) extrapolated = true;
  }
  return v;
}

bool check_system(const ActivitySystem& sys, Diagnostics& diag) {
  const std::size_t n = sys.size();
  if (n == 0 || (sys.model == ActivityModel::Unifac && sys.groups.size() != n) || (sys.model == ActivityModel::Nrtl && sys.nrtl.n != n)) {
    diag.error(DiagCode::InvalidMixtureInput, "activity system: component counts of the model and vapour pressures differ");
    return false;
  }
  return true;
}

void warn_extrapolated(const VleResult& r, Diagnostics& diag) {
  if (r.extrapolated)
    diag.warn(DiagCode::NoVaporPressureData, "a vapour pressure was extrapolated outside its Antoine range at " + std::to_string(r.T_K) + " K");
}

// Root of f on [lo, hi] with f(lo) < 0 < f(hi) (Illinois false position).
template <class F>
std::optional<double> solve_bracketed(F f, double lo, double hi, double flo, double fhi, int& iters) {
  int side = 0;
  for (iters = 0; iters < 200; ++iters) {
    const double x = (lo * fhi - hi * flo) / (fhi - flo);
    const auto fx = f(x);
    if (!fx) return std::nullopt;
    if (std::abs(*fx) < 1e-12 || hi - lo < 1e-10) return x;
    if (*fx < 0) {
      lo = x, flo = *fx;
      if (side == -1) fhi *= 0.5;
      side = -1;
    } else {
      hi = x, fhi = *fx;
      if (side == 1) flo *= 0.5;
      side = 1;
    }
  }
  return std::nullopt;
}

// Temperature where g(T) = 0 for g increasing in T, starting from T0.
template <class G>
std::optional<double> solve_temperature(G g, double T0, int& iters, Diagnostics& diag) {
  double lo = T0, hi = T0;
  auto glo = g(lo);
  if (!glo) return std::nullopt;
  auto ghi = glo;
  double step = 5.0;
  for (int k = 0; k < 60 && *glo > 0; ++k) {
    hi = lo, ghi = glo;
    lo = std::max(1.0, lo - step), step *= 1.5;
    glo = g(lo);
    if (!glo) return std::nullopt;
  }
  step = 5.0;
  for (int k = 0; k < 60 && *ghi < 0; ++k) {
    lo = hi, glo = ghi;
    hi += step, step *= 1.5;
    ghi = g(hi);
    if (!ghi) return std::nullopt;
  }
  if (!(*glo <= 0 && *ghi >= 0)) {
    diag.error(DiagCode::VleNotConverged, "no temperature bracket found for the bubble/dew point");
    return std::nullopt;
  }
  if (*glo == 0) return lo;
  if (*ghi == 0) return hi;
  auto T = solve_bracketed(g, lo, hi, *glo, *ghi, iters);
  if (!T) diag.error(DiagCode::VleNotConverged, "bubble/dew temperature iteration did not converge");
  return T;
}

// Dew pressure at T for vapour y: successive substitution on the liquid composition.
std::optional<VleResult> dew_pressure_quiet(const ActivitySystem& sys, double T, const std::vector<double>& y, Diagnostics& diag) {
  VleResult r;
  r.T_K = T;
  r.y = y;
  const auto lps = ln_psat(sys, T, r.extrapolated);
  std::vector<double> x = y;
  for (int it = 0; it < 1000; ++it) {
    auto lg = activity_ln_gamma(sys, T, x, diag);
    if (!lg) return std::nullopt;
    double inv = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) inv += y[i] / std::exp((*lg)[i] + lps[i]);
    const double P = 1.0 / inv;
    double change = 0.0;
    std::vector<double> xn(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
      xn[i] = y[i] * P / std::exp((*lg)[i] + lps[i]);
      change = std::max(change, std::abs(xn[i] - x[i]));
    }
    x = xn;
    r.iterations = it + 1;
    if (change < 1e-12) {
      r.P_Pa = P;
      r.x = x;
      r.gamma.resize(x.size());
      auto lg2 = activity_ln_gamma(sys, T, x, diag);
      if (!lg2) return std::nullopt;
      for (std::size_t i = 0; i < x.size(); ++i) r.gamma[i] = std::exp((*lg2)[i]);
      return r;
    }
  }
  diag.error(DiagCode::VleNotConverged, "dew-point liquid composition did not converge (liquid-liquid split?)");
  return std::nullopt;
}

std::optional<std::vector<double>> fractions(const std::vector<double>& v, std::size_t n, Diagnostics& diag) {
  return normalise(v, n, diag);
}

// Saturation temperature of a pure component at P (Antoine inverted), used as a starting point.
double t_sat(const AntoineCoefficients& c, double P) { return c.B / (c.A - std::log10(P)) - c.C; }

}  // namespace

double AntoineCoefficients::pressure(double T) const { return std::pow(10.0, A - B / (T + C)); }

std::optional<AntoineCoefficients> antoine_coefficients(const std::string& smiles) {
  const auto key = canonical(smiles);
  auto it = antoine_table().find(key);
  if (key.empty() || it == antoine_table().end()) return std::nullopt;
  return it->second;
}

std::optional<ActivitySystem> unifac_system_from_smiles(const std::vector<std::string>& smiles, Diagnostics& diag) {
  ActivitySystem sys;
  sys.model = ActivityModel::Unifac;
  bool ok = true;
  for (const auto& s : smiles) {
    auto g = unifac_groups_from_smiles(s, diag);
    auto a = antoine_coefficients(s);
    if (!a) diag.error(DiagCode::NoVaporPressureData, "no Antoine vapour-pressure coefficients for " + s);
    if (!g || !a) {
      ok = false;
      continue;
    }
    sys.groups.push_back(*g);
    sys.antoine.push_back(*a);
  }
  if (!ok) return std::nullopt;
  return sys;
}

std::optional<std::vector<double>> activity_ln_gamma(const ActivitySystem& sys, double T, const std::vector<double>& x,
                                                     Diagnostics& diag) {
  if (sys.model == ActivityModel::Unifac) return unifac_ln_gamma(sys.groups, T, x, diag);
  return nrtl_ln_gamma(sys.nrtl, T, x, diag);
}

std::optional<VleResult> bubble_pressure(const ActivitySystem& sys, double T, const std::vector<double>& x_in, Diagnostics& diag) {
  if (!check_system(sys, diag) || !valid_temperature(T, diag)) return std::nullopt;
  auto x = fractions(x_in, sys.size(), diag);
  if (!x) return std::nullopt;
  auto lg = activity_ln_gamma(sys, T, *x, diag);
  if (!lg) return std::nullopt;
  VleResult r;
  r.T_K = T;
  r.x = *x;
  const auto lps = ln_psat(sys, T, r.extrapolated);
  std::vector<double> partial(x->size());
  for (std::size_t i = 0; i < x->size(); ++i) partial[i] = (*x)[i] * std::exp((*lg)[i] + lps[i]), r.P_Pa += partial[i];
  r.y.resize(x->size());
  r.gamma.resize(x->size());
  for (std::size_t i = 0; i < x->size(); ++i) r.y[i] = partial[i] / r.P_Pa, r.gamma[i] = std::exp((*lg)[i]);
  r.iterations = 1;
  warn_extrapolated(r, diag);
  return r;
}

std::optional<VleResult> dew_pressure(const ActivitySystem& sys, double T, const std::vector<double>& y_in, Diagnostics& diag) {
  if (!check_system(sys, diag) || !valid_temperature(T, diag)) return std::nullopt;
  auto y = fractions(y_in, sys.size(), diag);
  if (!y) return std::nullopt;
  auto r = dew_pressure_quiet(sys, T, *y, diag);
  if (r) warn_extrapolated(*r, diag);
  return r;
}

std::optional<VleResult> bubble_temperature(const ActivitySystem& sys, double P, const std::vector<double>& x_in, Diagnostics& diag) {
  if (!check_system(sys, diag)) return std::nullopt;
  if (!(std::isfinite(P) && P > 0)) {
    diag.error(DiagCode::InvalidMixtureInput, "pressure must be a positive number of pascals");
    return std::nullopt;
  }
  auto x = fractions(x_in, sys.size(), diag);
  if (!x) return std::nullopt;
  double T0 = 0.0;
  for (std::size_t i = 0; i < x->size(); ++i) T0 += (*x)[i] * t_sat(sys.antoine[i], P);
  Diagnostics inner;
  auto g = [&](double T) -> std::optional<double> {
    auto r = bubble_pressure(sys, T, *x, inner);
    if (!r) return std::nullopt;
    return std::log(r->P_Pa / P);
  };
  int iters = 0;
  auto T = solve_temperature(g, T0, iters, diag);
  if (!T) {
    for (const auto& it : inner.items())
      if (it.severity == Severity::Error) diag.error(it.code, it.message);
    return std::nullopt;
  }
  Diagnostics last;
  auto r = bubble_pressure(sys, *T, *x, last);
  if (!r) return std::nullopt;
  r->iterations = iters;
  warn_extrapolated(*r, diag);
  return r;
}

std::optional<VleResult> dew_temperature(const ActivitySystem& sys, double P, const std::vector<double>& y_in, Diagnostics& diag) {
  if (!check_system(sys, diag)) return std::nullopt;
  if (!(std::isfinite(P) && P > 0)) {
    diag.error(DiagCode::InvalidMixtureInput, "pressure must be a positive number of pascals");
    return std::nullopt;
  }
  auto y = fractions(y_in, sys.size(), diag);
  if (!y) return std::nullopt;
  double T0 = 0.0;
  for (std::size_t i = 0; i < y->size(); ++i) T0 += (*y)[i] * t_sat(sys.antoine[i], P);
  Diagnostics inner;
  auto g = [&](double T) -> std::optional<double> {
    auto r = dew_pressure_quiet(sys, T, *y, inner);
    if (!r) return std::nullopt;
    return std::log(r->P_Pa / P);
  };
  int iters = 0;
  auto T = solve_temperature(g, T0, iters, diag);
  if (!T) {
    for (const auto& it : inner.items())
      if (it.severity == Severity::Error) diag.error(it.code, it.message);
    return std::nullopt;
  }
  Diagnostics last;
  auto r = dew_pressure_quiet(sys, *T, *y, last);
  if (!r) return std::nullopt;
  r->iterations = iters;
  warn_extrapolated(*r, diag);
  return r;
}

}  // namespace chemlarp
