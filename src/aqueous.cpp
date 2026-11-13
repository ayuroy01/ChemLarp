#include "chemlarp/aqueous.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>

#include "chemlarp/json.hpp"

namespace chemlarp {
namespace {

struct AqTerm {
  int index = 0;  // >= 0 component, -1 H+, -2 H2O
  double coef = 0.0;
};
struct AqMaster {
  const char* name;
  const char* species;
  int charge;
  double alk;
  const char* element;
  double gfw;
};
struct AqDefault {
  const char* element;
  const char* component;
};
struct AqSpeciesData {
  const char* name;
  int charge;
  AqTerm st[6];  // unused slots have coef 0
  double an[6];
  double lk, dh;
  bool has_gamma;
  double a, b;
};
struct AqGamma {
  bool has;
  double a, b;
};
struct AqPhaseData {
  const char* name;
  const char* formula;
  AqTerm st[6];
  double an[6];
  double lk, dh;
};

#include "aqueous_data.inc"

constexpr double kLn10 = 2.302585092994046;
constexpr double kRkJ = 8.314462618e-3;
constexpr std::size_t kNC = sizeof(kAqMasters) / sizeof(kAqMasters[0]);

double log_k(const double* a, double lk, double dh, double T) {
  return a[0] + a[1] * T + a[2] / T + a[3] * std::log10(T) + a[4] / (T * T) + a[5] * T * T + lk -
         dh / (kRkJ * kLn10) * (1.0 / T - 1.0 / 298.15);
}

// Debye–Hückel A and B at 1 atm from the dielectric constant (Malmberg & Maryott, J. Res. NBS 56 (1956) 1)
// and the density of water (Kell, J. Chem. Eng. Data 20 (1975) 97), t in °C.
void debye_huckel(double T, double& A, double& B) {
  const double t = T - 273.15;
  const double eps = 87.740 - 0.40008 * t + 9.398e-4 * t * t - 1.410e-6 * t * t * t;
  const double rho = (999.83952 + 16.945176 * t - 7.9870401e-3 * t * t - 46.170461e-6 * t * t * t +
                      105.56302e-9 * t * t * t * t - 280.54253e-12 * t * t * t * t * t) /
                     (1.0 + 16.879850e-3 * t) / 1000.0;  // g/cm³
  A = 1.82483e6 * std::sqrt(rho) / std::pow(eps * T, 1.5);
  B = 50.2916 * std::sqrt(rho) / std::sqrt(eps * T);
}

std::string num(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%g", v);
  return b;
}

int component_index(const std::string& name) {
  std::string n = name;
  // C(+4) -> C(4)
  const auto p = n.find("(+");
  if (p != std::string::npos) n.erase(p + 1, 1);
  for (std::size_t k = 0; k < kNC; ++k)
    if (n == kAqMasters[k].name) return int(k);
  for (const auto& d : kAqDefaults)
    if (n == d.element) return component_index(d.component);
  return -1;
}

int phase_index(const std::string& name) {
  for (std::size_t p = 0; p < sizeof(kAqPhases) / sizeof(kAqPhases[0]); ++p)
    if (name == kAqPhases[p].name) return int(p);
  // Case-insensitive fallback.
  auto lower = [](std::string s) {
    for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
    return s;
  };
  for (std::size_t p = 0; p < sizeof(kAqPhases) / sizeof(kAqPhases[0]); ++p)
    if (lower(name) == lower(kAqPhases[p].name)) return int(p);
  return -1;
}

// ---- "add <formula>": element counts and valence states ----

bool parse_formula(const std::string& f, std::size_t& i, std::map<std::string, double>& out, std::string& err) {
  while (i < f.size() && f[i] != ')') {
    if (f[i] == '(') {
      ++i;
      std::map<std::string, double> inner;
      if (!parse_formula(f, i, inner, err)) return false;
      if (i >= f.size() || f[i] != ')') return err = "unbalanced parentheses", false;
      ++i;
      double n = 1;
      if (i < f.size() && std::isdigit(static_cast<unsigned char>(f[i]))) n = std::strtod(f.c_str() + i, nullptr);
      while (i < f.size() && (std::isdigit(static_cast<unsigned char>(f[i])) || f[i] == '.')) ++i;
      for (auto& [e, c] : inner) out[e] += c * n;
    } else if (std::isupper(static_cast<unsigned char>(f[i]))) {
      std::string el(1, f[i++]);
      while (i < f.size() && std::islower(static_cast<unsigned char>(f[i]))) el += f[i++];
      double n = 1;
      if (i < f.size() && std::isdigit(static_cast<unsigned char>(f[i]))) n = std::strtod(f.c_str() + i, nullptr);
      while (i < f.size() && (std::isdigit(static_cast<unsigned char>(f[i])) || f[i] == '.')) ++i;
      out[el] += n;
    } else {
      return err = std::string("unexpected '") + f[i] + "'", false;
    }
  }
  return true;
}

// Component amounts per mole of a neutral formula unit. Hydrate water (after ':' or '·') is ignored.
bool formula_components(std::string f, std::vector<std::pair<int, double>>& out, std::string& err) {
  for (const char* sep : {":", "\xc2\xb7"}) {
    const auto p = f.find(sep);
    if (p != std::string::npos) f = f.substr(0, p);
  }
  // Common groups fix the valence of their central atom: NH4+ N(-3), SO4-2 S(6), CO3-2 C(4), NO3- N(5),
  // NO2- N(3). Each is counted (with a parenthesized multiplier) and removed before the rest is read.
  struct Group {
    const char* text;
    const char* component;
    int charge;
  };
  static const Group groups[] = {{"NH4", "N(-3)", 1}, {"SO4", "S(6)", -2}, {"CO3", "C(4)", -2}, {"NO3", "N(5)", -1}, {"NO2", "N(3)", -1}};
  double charge = 0.0;
  for (const Group& g : groups) {
    const std::string paren = std::string("(") + g.text + ")";
    double count = 0;
    for (std::size_t p; (p = f.find(paren)) != std::string::npos;) {
      std::size_t q = p + paren.size();
      double n = 1;
      if (q < f.size() && std::isdigit(static_cast<unsigned char>(f[q]))) n = std::strtod(f.c_str() + q, nullptr);
      while (q < f.size() && std::isdigit(static_cast<unsigned char>(f[q]))) ++q;
      count += n;
      f.erase(p, q - p);
    }
    for (std::size_t p; (p = f.find(g.text)) != std::string::npos;) {
      // A bare group takes no multiplier of its own ("CaCO3", "NH4NO3").
      count += 1;
      f.erase(p, std::string(g.text).size());
    }
    if (count > 0) {
      out.push_back({component_index(g.component), count});
      charge += g.charge * count;
    }
  }
  std::map<std::string, double> el;
  std::size_t i = 0;
  if (!parse_formula(f, i, el, err)) return false;
  if (i != f.size()) return err = "unbalanced parentheses", false;
  static const std::map<std::string, int> fixed = {
      {"H", 1},  {"O", -2}, {"F", -1}, {"Cl", -1}, {"Br", -1}, {"Li", 1}, {"Na", 1}, {"K", 1},  {"Mg", 2},
      {"Ca", 2}, {"Sr", 2}, {"Ba", 2}, {"Al", 3},  {"Zn", 2},  {"Cd", 2}, {"Pb", 2}, {"Si", 4}, {"B", 3}, {"P", 5}};
  std::string variable;
  for (const auto& [e, n] : el) {
    auto it = fixed.find(e);
    const bool oxyhalogen = (e == "Cl" || e == "Br" || e == "F") && el.count("O");  // ClO4-, BrO3-: not halides
    if (it != fixed.end() && !oxyhalogen) {
      charge += it->second * n;
    } else {
      if (!variable.empty()) return err = "more than one element with a variable valence (" + variable + ", " + e + "); give the components instead", false;
      variable = e;
    }
  }
  for (const auto& [e, n] : el) {
    if (e == "H" || e == "O") continue;
    std::string comp = e;
    if (e == variable) {
      const double v = -charge / n;
      if (std::fabs(v - std::round(v)) > 1e-9) return err = "no whole oxidation state for " + e, false;
      comp = e + "(" + std::to_string(static_cast<long>(std::lround(v))) + ")";
    }
    const int k = component_index(comp);
    if (k < 0) return err = "the database has no component " + comp, false;
    out.push_back({k, n});
  }
  return true;
}

struct AqSpeciesRow {
  std::string name;
  int charge = 0;
  std::vector<AqTerm> st;
  double logK = 0.0;
  bool has_gamma = false;
  double a = 0.0, b = 0.0;
  double alk = 0.0;
};

double log_gamma(const AqSpeciesRow& s, double I, double A, double B) {
  if (s.charge == 0) return 0.1 * I;
  const double z2 = double(s.charge * s.charge), sI = std::sqrt(I);
  if (s.has_gamma) return -A * z2 * sI / (1.0 + B * s.a * sI) + s.b * I;
  return -A * z2 * (sI / (1.0 + sI) - 0.3 * I);
}

// Dense Gaussian elimination with partial pivoting; returns false when singular.
bool aq_solve_dense(std::vector<double> a, std::vector<double>& b, std::size_t n) {
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t r = k + 1; r < n; ++r)
      if (std::fabs(a[r * n + k]) > std::fabs(a[p * n + k])) p = r;
    if (!(std::fabs(a[p * n + k]) > 1e-300)) return false;
    if (p != k) {
      for (std::size_t c = 0; c < n; ++c) std::swap(a[k * n + c], a[p * n + c]);
      std::swap(b[k], b[p]);
    }
    for (std::size_t r = k + 1; r < n; ++r) {
      const double f = a[r * n + k] / a[k * n + k];
      if (f == 0.0) continue;
      for (std::size_t c = k; c < n; ++c) a[r * n + c] -= f * a[k * n + c];
      b[r] -= f * b[k];
    }
  }
  for (std::size_t r = n; r-- > 0;) {
    double s = b[r];
    for (std::size_t c = r + 1; c < n; ++c) s -= a[r * n + c] * b[c];
    b[r] = s / a[r * n + r];
  }
  return true;
}

std::vector<AqTerm> terms(const AqTerm (&t)[6]) {
  std::vector<AqTerm> v;
  for (const auto& x : t)
    if (x.coef != 0.0) v.push_back(x);
  return v;
}

}  // namespace

std::vector<std::string> aqueous_components() {
  std::vector<std::string> v;
  for (const auto& m : kAqMasters) v.push_back(m.name);
  return v;
}

std::vector<std::string> aqueous_phases() {
  std::vector<std::string> v;
  for (const auto& p : kAqPhases) v.push_back(p.name);
  return v;
}

std::optional<AqueousInput> parse_aqueous(std::string_view text, Diagnostics& diag) {
  AqueousInput in;
  double unit = 1e-3;
  bool ok = true;
  int line = 0;
  std::map<int, double> totals;
  std::vector<int> order;
  auto fail = [&](const std::string& m) {
    diag.error(DiagCode::InvalidInput, "line " + std::to_string(line) + ": " + m);
    ok = false;
  };
  auto number = [](const std::string& s, double& v) {
    char* e = nullptr;
    v = std::strtod(s.c_str(), &e);
    return !s.empty() && e == s.c_str() + s.size() && std::isfinite(v);
  };
  auto add_total = [&](int k, double v) {
    if (!totals.count(k)) order.push_back(k);
    totals[k] += v;
  };
  std::istringstream ss{std::string(text)};
  for (std::string raw; std::getline(ss, raw);) {
    ++line;
    const auto h = raw.find('#');
    if (h != std::string::npos) raw = raw.substr(0, h);
    std::istringstream ls(raw);
    std::vector<std::string> t;
    for (std::string w; ls >> w;) t.push_back(w);
    if (t.empty()) continue;
    std::string kw = t[0];
    std::string low = kw;
    for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
    double v = 0;
    if (low == "temperature" || low == "temp") {
      if (t.size() != 2 || !number(t[1], in.T_C) || in.T_C < 0 || in.T_C > 100) fail("temperature must be 0–100 °C");
    } else if (low == "units") {
      if (t.size() == 2 && (t[1] == "mol/kgw" || t[1] == "mol/kg")) unit = 1.0;
      else if (t.size() == 2 && (t[1] == "mmol/kgw" || t[1] == "mmol/kg")) unit = 1e-3;
      else fail("units must be mol/kgw or mmol/kgw");
    } else if (low == "ph") {
      if (t.size() < 2 || t.size() > 3 || !number(t[1], in.pH) || in.pH < -2 || in.pH > 16) fail("pH <value> [charge]");
      in.charge_balance = t.size() == 3 && t[2] == "charge";
      if (t.size() == 3 && t[2] != "charge") fail("expected 'charge' after the pH");
    } else if (low == "add") {
      if (t.size() != 3 || !number(t[2], v) || v < 0) {
        fail("add <formula> <amount>");
        continue;
      }
      std::vector<std::pair<int, double>> comps;
      std::string err;
      if (!formula_components(t[1], comps, err)) {
        fail(t[1] + ": " + err);
        continue;
      }
      for (const auto& [k, n] : comps) add_total(k, n * v * unit);
    } else if (low == "equilibrium" || low == "equilibrium_phases") {
      if (t.size() < 2 || t.size() > 4) {
        fail("equilibrium <phase> [target SI] [moles]");
        continue;
      }
      AqueousPhaseTarget p;
      const int pi = phase_index(t[1]);
      if (pi < 0) {
        fail("no phase named " + t[1] + " in the database");
        continue;
      }
      p.name = kAqPhases[pi].name;
      if (t.size() >= 3 && !number(t[2], p.target_si)) fail("bad target saturation index " + t[2]);
      if (t.size() == 4 && (!number(t[3], p.amount_mol) || p.amount_mol < 0)) fail("bad amount " + t[3]);
      in.phases.push_back(p);
    } else {
      const int k = component_index(kw);
      if (k < 0) {
        fail("unknown component or statement '" + kw + "'");
        continue;
      }
      if (t.size() != 2 || !number(t[1], v) || v < 0) {
        fail(kw + " needs one non-negative total");
        continue;
      }
      add_total(k, v * unit);
    }
  }
  if (!ok) return std::nullopt;
  for (int k : order) in.totals.push_back({kAqMasters[k].name, totals[k]});
  return in;
}

std::optional<AqueousResult> aqueous_equilibrium(const AqueousInput& input, Diagnostics& diag) {
  const double T = input.T_C + 273.15;
  if (!(input.T_C >= 0 && input.T_C <= 100)) {
    diag.error(DiagCode::InvalidInput, "the temperature must be 0–100 °C");
    return std::nullopt;
  }
  AqueousResult res;
  res.T_K = T;
  debye_huckel(T, res.A_dh, res.B_dh);
  const double A = res.A_dh, B = res.B_dh;

  // Components in play: given totals, and those of the equilibrium phases.
  std::vector<double> total(kNC, 0.0);
  std::vector<bool> active(kNC, false);
  for (const auto& [name, v] : input.totals) {
    const int k = component_index(name);
    if (k < 0 || !(v >= 0) || !std::isfinite(v)) {
      diag.error(DiagCode::InvalidInput, "unknown component or bad total: " + name);
      return std::nullopt;
    }
    total[k] += v;
    if (v > 0) active[k] = true;
  }
  struct Phase {
    int index;
    std::vector<AqTerm> st;
    double logK, target, amount, n = 0.0;
    bool on = true;
  };
  std::vector<Phase> phases;
  for (const auto& p : input.phases) {
    const int pi = phase_index(p.name);
    if (pi < 0) {
      diag.error(DiagCode::InvalidInput, "no phase named " + p.name);
      return std::nullopt;
    }
    const auto& d = kAqPhases[pi];
    Phase ph{pi, terms(d.st), log_k(d.an, d.lk, d.dh, T), p.target_si, p.amount_mol};
    for (const auto& t : ph.st)
      if (t.index >= 0) active[t.index] = true;
    phases.push_back(ph);
  }
  // Index of each active component in the unknown vector.
  std::vector<int> col(kNC, -1);
  std::vector<int> comps;
  for (std::size_t k = 0; k < kNC; ++k)
    if (active[k]) col[k] = int(comps.size()), comps.push_back(int(k));
  const std::size_t nc = comps.size();

  // AqSpeciesRow: the masters, H+, and every species whose components are all active.
  std::vector<AqSpeciesRow> sp;
  {
    AqSpeciesRow h;
    h.name = "H+";
    h.charge = 1;
    h.st = {{-1, 1.0}};
    h.has_gamma = kAqHydrogenGamma.has;
    h.a = kAqHydrogenGamma.a;
    h.b = kAqHydrogenGamma.b;
    h.alk = -1.0;
    sp.push_back(h);
  }
  for (int k : comps) {
    AqSpeciesRow m;
    m.name = kAqMasters[k].species;
    m.charge = kAqMasters[k].charge;
    m.st = {{k, 1.0}};
    m.has_gamma = kAqMasterGamma[k].has;
    m.a = kAqMasterGamma[k].a;
    m.b = kAqMasterGamma[k].b;
    m.alk = kAqMasters[k].alk;
    sp.push_back(m);
  }
  for (const auto& d : kAqSpecies) {
    AqSpeciesRow s;
    s.name = d.name;
    s.charge = d.charge;
    s.st = terms(d.st);
    bool ok = true;
    for (const auto& t : s.st) {
      if (t.index >= 0 && !active[t.index]) ok = false;
      if (t.index >= 0) s.alk += t.coef * kAqMasters[t.index].alk;
      if (t.index == -1) s.alk -= t.coef;
    }
    if (!ok) continue;
    s.logK = log_k(d.an, d.lk, d.dh, T);
    s.has_gamma = d.has_gamma;
    s.a = d.a;
    s.b = d.b;
    sp.push_back(s);
  }
  const std::size_t ns = sp.size();

  // Unknowns: log10 a of each active component, log10 a(H+) when the pH is not fixed, phase amounts.
  const bool cb = input.charge_balance;
  std::vector<double> x(nc);
  for (std::size_t c = 0; c < nc; ++c) x[c] = std::log10(std::max(total[comps[c]], 1e-10) * 0.5);
  double xH = -input.pH, logaw = 0.0, I = 0.0;
  std::vector<double> m(ns), lg(ns, 0.0);
  auto molalities = [&]() {
    for (std::size_t j = 0; j < ns; ++j) {
      double e = sp[j].logK - lg[j];
      for (const auto& t : sp[j].st) e += t.coef * (t.index >= 0 ? x[col[t.index]] : t.index == -1 ? xH : logaw);
      m[j] = std::pow(10.0, std::clamp(e, -300.0, 300.0));
    }
  };
  auto update_activity = [&]() {
    I = 0.0;
    double sum = 0.0;
    for (std::size_t j = 0; j < ns; ++j) {
      I += 0.5 * sp[j].charge * sp[j].charge * m[j];
      sum += m[j];
    }
    I = std::min(I, 10.0);
    for (std::size_t j = 0; j < ns; ++j) lg[j] = log_gamma(sp[j], I, A, B);
    logaw = std::log10(std::max(1.0 - 0.017 * sum, 0.05));
  };

  const int max_it = 400;
  int it_total = 0;
  bool converged = false;
  // Start: a component that only a phase supplies takes the activity that puts that phase at its target
  // saturation index (with the other activities as they stand), so Newton begins near the solution.
  {
    std::vector<bool> fixed_by_phase(kNC, false);
    for (auto& ph : phases) {
      int pick = -1;
      double coef = 0.0;
      for (const auto& t : ph.st)
        if (t.index >= 0 && !(total[t.index] > 0) && !fixed_by_phase[t.index]) pick = t.index, coef = t.coef;
      if (pick < 0) continue;
      double rest = ph.logK + ph.target;
      for (const auto& t : ph.st)
        if (t.index != pick) rest -= t.coef * (t.index >= 0 ? x[col[t.index]] : t.index == -1 ? xH : 0.0);
      x[col[pick]] = std::clamp(rest / coef, -20.0, 1.0);
      fixed_by_phase[pick] = true;
    }
  }
  for (int pass = 0; pass < 20; ++pass) {
    std::vector<std::size_t> on;  // active phases
    for (std::size_t p = 0; p < phases.size(); ++p)
      if (phases[p].on) on.push_back(p);
    const std::size_t n = nc + (cb ? 1 : 0) + on.size();
    converged = false;
    bool restart = false;
    std::vector<int> at_cap(phases.size(), 0);
    molalities();
    update_activity();
    for (int it = 0; it < max_it; ++it, ++it_total) {
      molalities();
      // Residuals (unscaled) and their scales.
      std::vector<double> f(n, 0.0), scale(n, 1.0), J(n * n, 0.0);
      for (std::size_t c = 0; c < nc; ++c) {
        const int k = comps[c];
        double tot = total[k], s = 0.0;
        for (const auto& ph : phases)
          for (const auto& t : ph.st)
            if (t.index == k) tot += t.coef * ph.n;
        f[c] = -tot;
        scale[c] = std::max(std::fabs(tot), 1e-30);
        for (std::size_t j = 0; j < ns; ++j)
          for (const auto& t : sp[j].st)
            if (t.index == k) {
              f[c] += t.coef * m[j];
              s += std::fabs(t.coef) * m[j];
            }
        scale[c] = std::max(scale[c], s);
      }
      if (cb) {
        double s = 0.0;
        for (std::size_t j = 0; j < ns; ++j) {
          f[nc] += sp[j].charge * m[j];
          s += std::abs(sp[j].charge) * m[j];
        }
        scale[nc] = std::max(s, 1e-30);
      }
      const std::size_t p0 = nc + (cb ? 1 : 0);
      for (std::size_t q = 0; q < on.size(); ++q) {
        const Phase& ph = phases[on[q]];
        double e = -ph.logK - ph.target;
        for (const auto& t : ph.st) e += t.coef * (t.index >= 0 ? x[col[t.index]] : t.index == -1 ? xH : logaw);
        f[p0 + q] = e;
      }
      double err = 0.0;
      for (std::size_t r = 0; r < n; ++r) err = std::max(err, std::fabs(f[r]) / (r < p0 ? scale[r] : 1.0));
      // Jacobian: d m_j / d x_l = ln10 ν_jl m_j (activity coefficients lagged).
      auto var = [&](const AqTerm& t) -> int {
        if (t.index >= 0) return col[t.index];
        if (t.index == -1 && cb) return int(nc);
        return -1;
      };
      for (std::size_t j = 0; j < ns; ++j) {
        for (const auto& tr : sp[j].st) {
          const int r = tr.index >= 0 ? col[tr.index] : -1;
          for (const auto& tc : sp[j].st) {
            const int c2 = var(tc);
            if (c2 < 0) continue;
            if (r >= 0) J[std::size_t(r) * n + std::size_t(c2)] += kLn10 * tr.coef * tc.coef * m[j];
          }
        }
        if (cb)
          for (const auto& tc : sp[j].st) {
            const int c2 = var(tc);
            if (c2 >= 0) J[nc * n + std::size_t(c2)] += kLn10 * sp[j].charge * tc.coef * m[j];
          }
      }
      for (std::size_t q = 0; q < on.size(); ++q) {
        const Phase& ph = phases[on[q]];
        for (const auto& t : ph.st) {
          const int c2 = var(t);
          if (c2 >= 0) J[(p0 + q) * n + std::size_t(c2)] += t.coef;
          if (t.index >= 0) J[std::size_t(col[t.index]) * n + p0 + q] -= t.coef;  // d(mass balance)/d n
        }
      }
      const double Iold = I;
      if (err < 1e-13 && it > 0) {
        update_activity();
        if (std::fabs(I - Iold) <= 1e-12 * std::max(I, 1e-12)) {
          converged = true;
          break;
        }
        continue;
      }
      // Row scaling keeps mass balances of very different sizes comparable.
      for (std::size_t r = 0; r < p0; ++r) {
        for (std::size_t c2 = 0; c2 < n; ++c2) J[r * n + c2] /= scale[r];
        f[r] /= scale[r];
      }
      std::vector<double> dx(f);
      for (double& v : dx) v = -v;
      if (!aq_solve_dense(J, dx, n)) break;
      double big = 0.0;
      for (std::size_t c = 0; c < p0; ++c) big = std::max(big, std::fabs(dx[c]));
      const double step = big > 1.0 ? 1.0 / big : 1.0;  // at most one log unit per iteration
      for (std::size_t c = 0; c < nc; ++c) x[c] += step * dx[c];
      if (cb) xH += step * dx[nc];
      for (std::size_t q = 0; q < on.size(); ++q) {
        Phase& ph = phases[on[q]];
        // A phase cannot dissolve more than is present: held at its amount for a few steps, it is dissolved
        // completely and leaves the active set.
        const double want = ph.n + step * dx[p0 + q];
        ph.n = std::min(want, ph.amount);
        at_cap[on[q]] = want > ph.amount ? at_cap[on[q]] + 1 : 0;
        if (at_cap[on[q]] >= 3) {
          ph.n = ph.amount;
          ph.on = false;
          restart = true;
        }
      }
      update_activity();
      if (restart) break;
    }
    if (restart) continue;
    if (!converged) break;
    // Active set: a phase cannot dissolve more than is present; an absent phase only precipitates when
    // the solution is supersaturated.
    bool changed = false;
    for (auto& ph : phases) {
      if (ph.on && ph.n > ph.amount) {
        ph.n = ph.amount;
        ph.on = false;
        changed = true;
      }
    }
    for (auto& ph : phases) {
      if (ph.on) continue;
      double si = -ph.logK;
      for (const auto& t : ph.st) si += t.coef * (t.index >= 0 ? x[col[t.index]] : t.index == -1 ? xH : logaw);
      if (si > ph.target + 1e-9) {
        ph.on = true;
        changed = true;
      }
    }
    if (!changed) break;
  }
  res.iterations = it_total;
  if (!converged) {
    diag.error(DiagCode::NotConverged, "the speciation did not converge (" + std::to_string(it_total) + " iterations)");
    return std::nullopt;
  }

  // ---- Results ----
  res.pH = -xH;
  res.ionic_strength = I;
  res.water_activity = std::pow(10.0, logaw);
  for (std::size_t j = 0; j < ns; ++j) {
    res.charge_balance_eq += sp[j].charge * m[j];
    res.alkalinity_eq += sp[j].alk * m[j];
    AqueousSpecies out{sp[j].name, sp[j].charge, m[j], m[j] * std::pow(10.0, lg[j]), lg[j], {}};
    for (const auto& t : sp[j].st)
      if (t.index >= 0) out.components.push_back({kAqMasters[t.index].name, t.coef});
    res.species.push_back(std::move(out));
  }
  // OH- is formed from water: add it if not present as a species (it always is: H2O − H+).
  std::sort(res.species.begin(), res.species.end(), [](const auto& a, const auto& b) { return a.molality > b.molality; });
  for (std::size_t c = 0; c < nc; ++c) {
    double tot = total[comps[c]];
    for (const auto& ph : phases)
      for (const auto& t : ph.st)
        if (t.index == comps[c]) tot += t.coef * ph.n;
    res.totals.push_back({kAqMasters[comps[c]].name, tot});
  }
  for (const auto& d : kAqPhases) {
    const auto st = terms(d.st);
    bool ok = true;
    double iap = 0.0;
    for (const auto& t : st) {
      if (t.index >= 0 && !active[t.index]) ok = false;
      else iap += t.coef * (t.index >= 0 ? x[col[t.index]] : t.index == -1 ? xH : logaw);
    }
    if (!ok) continue;
    const double lk = log_k(d.an, d.lk, d.dh, T);
    res.saturation.push_back({d.name, d.formula, iap - lk, iap, lk});
  }
  for (const auto& ph : phases)
    res.phases.push_back({kAqPhases[ph.index].name, ph.target, ph.n, ph.amount - ph.n});
  if (I > 1.0)
    diag.warn(DiagCode::ThermoApproximation, "ionic strength " + num(I) +
                                                 " mol/kg: the Debye–Hückel and Davies equations are unreliable above about 0.5–1 mol/kg");
  return res;
}

std::string aqueous_json(std::string_view text) {
  Diagnostics d;
  auto in = parse_aqueous(text, d);
  std::optional<AqueousResult> r;
  if (in) r = aqueous_equilibrium(*in, d);
  auto jn = [](double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.12g", v);
    return std::isfinite(v) ? std::string(b) : std::string("null");
  };
  std::string o = std::string("{\"ok\":") + (r && !d.has_errors() ? "true" : "false");
  if (r) {
    o += ",\"T_K\":" + jn(r->T_K) + ",\"pH\":" + jn(r->pH) + ",\"ionic_strength\":" + jn(r->ionic_strength) +
         ",\"water_activity\":" + jn(r->water_activity) + ",\"A_dh\":" + jn(r->A_dh) + ",\"B_dh\":" + jn(r->B_dh) +
         ",\"charge_balance_eq\":" + jn(r->charge_balance_eq) + ",\"alkalinity_eq\":" + jn(r->alkalinity_eq) +
         ",\"iterations\":" + std::to_string(r->iterations) + ",\"totals\":[";
    for (std::size_t i = 0; i < r->totals.size(); ++i)
      o += (i ? "," : "") + std::string("{\"component\":") + json_escape(r->totals[i].first) + ",\"molality\":" + jn(r->totals[i].second) + "}";
    o += "],\"species\":[";
    for (std::size_t i = 0; i < r->species.size(); ++i) {
      const auto& s = r->species[i];
      o += (i ? "," : "") + std::string("{\"name\":") + json_escape(s.name) + ",\"charge\":" + std::to_string(s.charge) +
           ",\"molality\":" + jn(s.molality) + ",\"activity\":" + jn(s.activity) + ",\"log_gamma\":" + jn(s.log_gamma) +
           ",\"components\":{";
      for (std::size_t k = 0; k < s.components.size(); ++k)
        o += (k ? "," : "") + json_escape(s.components[k].first) + ":" + jn(s.components[k].second);
      o += "}}";
    }
    o += "],\"saturation\":[";
    for (std::size_t i = 0; i < r->saturation.size(); ++i) {
      const auto& s = r->saturation[i];
      o += (i ? "," : "") + std::string("{\"phase\":") + json_escape(s.phase) + ",\"formula\":" + json_escape(s.formula) +
           ",\"si\":" + jn(s.si) + ",\"log_iap\":" + jn(s.log_iap) + ",\"log_k\":" + jn(s.log_k) + "}";
    }
    o += "],\"phases\":[";
    for (std::size_t i = 0; i < r->phases.size(); ++i) {
      const auto& p = r->phases[i];
      o += (i ? "," : "") + std::string("{\"phase\":") + json_escape(p.phase) + ",\"target_si\":" + jn(p.target_si) +
           ",\"dissolved_mol\":" + jn(p.dissolved_mol) + ",\"remaining_mol\":" + jn(p.remaining_mol) + "}";
    }
    o += "]";
  }
  return o + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp
