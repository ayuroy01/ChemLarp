// Mechanism input: OmniChem .ocm, CHEMKIN-II (+ NASA-7 thermo cards), Cantera YAML.
//
// All three formats feed one ReactionBuilder, so unit conversion, reaction-kind
// validation, element balance and duplicate detection are identical for every
// source. Unit conversion: A is given per (volume/amount)^(order-1); with
// conc = m^3 per amount unit, A_SI = A * conc^(order-1).

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "omnichem/elements.hpp"
#include "omnichem/formula.hpp"
#include "omnichem/kinetics.hpp"
#include "omnichem/smiles.hpp"
#include "yaml_lite.hpp"

namespace omnichem {
namespace {

constexpr double kAtm = 101325.0;

// ---- text helpers ------------------------------------------------------------

struct Line {
  std::string text;
  std::size_t pos = 0;  // offset of the line in the source text
  std::size_t no = 0;   // 1-based line number
};

std::vector<Line> split_lines(std::string_view text) {
  std::vector<Line> out;
  std::size_t start = 0, no = 0;
  while (start <= text.size()) {
    std::size_t e = text.find('\n', start);
    if (e == std::string_view::npos) e = text.size();
    std::string s(text.substr(start, e - start));
    if (!s.empty() && s.back() == '\r') s.pop_back();
    out.push_back({std::move(s), start, ++no});
    if (e == text.size()) break;
    start = e + 1;
  }
  return out;
}

struct Tok {
  std::string text;
  std::size_t pos;
};

std::vector<Tok> split_ws(std::string_view s, std::size_t base) {
  std::vector<Tok> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    if (i >= s.size()) break;
    const std::size_t st = i;
    while (i < s.size() && !std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    out.push_back({std::string(s.substr(st, i - st)), base + st});
  }
  return out;
}

std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string trim(std::string_view s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return std::string(s.substr(a, b - a));
}

// Accepts Fortran exponents (1.0D+03).
bool parse_number(std::string s, double& v) {
  s = trim(s);
  if (s.empty()) return false;
  for (char& c : s)
    if (c == 'D' || c == 'd') c = 'E';
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  return end == s.c_str() + s.size() && std::isfinite(v);
}

// '!' starts a comment anywhere; '#' only at line start or after whitespace
// (it is the SMILES triple bond).
std::string strip_comment(std::string_view line, bool hash_comments) {
  std::size_t cut = line.find('!');
  if (hash_comments)
    for (std::size_t c = 0; c < line.size() && c < cut; ++c)
      if (line[c] == '#' && (c == 0 || std::isspace(static_cast<unsigned char>(line[c - 1])))) { cut = c; break; }
  return std::string(cut == std::string_view::npos ? line : line.substr(0, cut));
}

std::string at_line(std::size_t no) { return "line " + std::to_string(no) + ": "; }

// ---- units -------------------------------------------------------------------

struct Units {
  double ea = 4.184;   // J/mol per unit of Ea (CHEMKIN default cal/mol)
  double conc = 1e-6;  // m^3 per (amount unit) for A conversion (CHEMKIN default cm^3/mol)
  double pressure = kAtm;  // Pa per PLOG pressure unit (CHEMKIN: atm)
};

enum class UnitResult { Ok, Unknown, Unsupported };

UnitResult apply_unit(const std::string& tok, Units& u) {
  const std::string t = lower(tok);
  if (t == "cal/mol" || t == "cal/mole") u.ea = 4.184;
  else if (t == "kcal/mol" || t == "kcal/mole") u.ea = 4184.0;
  else if (t == "j/mol" || t == "joules/mole" || t == "j/mole") u.ea = 1.0;
  else if (t == "kj/mol" || t == "kjoules/mole" || t == "kj/mole") u.ea = 1000.0;
  else if (t == "k" || t == "kelvin" || t == "kelvins") u.ea = kGasConstant;
  else if (t == "ev" || t == "evolts") u.ea = 96485.33212;
  else if (t == "cm3" || t == "moles" || t == "mole" || t == "mol") u.conc = 1e-6;
  else if (t == "m3") u.conc = 1.0;
  else if (t == "molecules") return UnitResult::Unsupported;
  else return UnitResult::Unknown;
  return UnitResult::Ok;
}

// ---- NASA-7 thermo cards -------------------------------------------------------

struct Card {
  std::string name;
  Composition comp{};
  Nasa7 nasa;
  std::size_t line = 0;
};

std::string col(const std::string& s, std::size_t a, std::size_t n) {
  if (a >= s.size()) return "";
  return s.substr(a, n);
}

// CHEMKIN-II fixed format: card line 1 has '1' in column 80, then lines with 2, 3, 4.
// Line 1: name (1-18), 4 x (element 2 + count 3) at 25-44, phase 45, Tlow 46-55,
// Thigh 56-65, Tmid 66-73, optional 5th element at 74-78.
// Lines 2-4: 15-column coefficients: a1..a5 (high), a6 a7 (high) a1..a3 (low), a4..a7 (low).
void parse_cards(const std::vector<Line>& lines, std::map<std::string, Card>& out, Diagnostics& d) {
  double t_mid_default = 1000.0;
  for (std::size_t i = 0; i < lines.size();) {
    const std::string& L = lines[i].text;
    const std::string t = trim(L);
    if (t.empty() || t[0] == '!') { ++i; continue; }
    const std::string first = upper(split_ws(t, 0)[0].text);
    if (first.rfind("THER", 0) == 0 || first == "END") { ++i; continue; }
    if (L.size() < 80 || L[79] != '1') {
      // Global defaults line "Tlow Tmid Thigh".
      auto toks = split_ws(t, 0);
      double a, b, c;
      if (toks.size() == 3 && parse_number(toks[0].text, a) && parse_number(toks[1].text, b) &&
          parse_number(toks[2].text, c)) {
        t_mid_default = b;
        ++i;
        continue;
      }
      d.error(DiagCode::MechanismSyntax, at_line(lines[i].no) +
                                             "unrecognised thermo line (expected a NASA-7 card with '1' in column 80)")
          .at(lines[i].pos);
      ++i;
      continue;
    }
    if (i + 3 >= lines.size()) {
      d.error(DiagCode::MechanismSyntax, at_line(lines[i].no) + "incomplete NASA-7 card (needs 4 lines)").at(lines[i].pos);
      return;
    }
    Card card;
    card.line = lines[i].no;
    bool ok = true;
    for (int k = 1; k <= 3; ++k) {
      const std::string& Lk = lines[i + static_cast<std::size_t>(k)].text;
      if (Lk.size() >= 80 && Lk[79] != static_cast<char>('1' + k)) {
        d.error(DiagCode::MechanismSyntax, at_line(lines[i + static_cast<std::size_t>(k)].no) + "expected card line " +
                                               std::to_string(k + 1) + " (column 80)")
            .at(lines[i + static_cast<std::size_t>(k)].pos);
        ok = false;
      }
    }
    const auto name_toks = split_ws(col(L, 0, 18), 0);
    card.name = name_toks.empty() ? "" : name_toks[0].text;
    std::vector<std::pair<std::string, std::string>> elems;
    for (int k = 0; k < 4; ++k)
      elems.push_back({trim(col(L, 24 + 5 * static_cast<std::size_t>(k), 2)), trim(col(L, 26 + 5 * static_cast<std::size_t>(k), 3))});
    if (L.size() >= 78) elems.push_back({trim(col(L, 73, 2)), trim(col(L, 75, 3))});
    for (auto& [sym, cnt] : elems) {
      if (sym.empty() || cnt.empty() || sym == "0") continue;
      double n = 0;
      if (!parse_number(cnt, n)) {
        d.error(DiagCode::MechanismSyntax, at_line(card.line) + "bad element count '" + cnt + "' for " + card.name).at(lines[i].pos);
        ok = false;
        continue;
      }
      if (n == 0) continue;
      std::string s2 = sym;
      s2[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s2[0])));
      if (s2.size() > 1) s2[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(s2[1])));
      if (s2 == "E") {
        d.error(DiagCode::UnsupportedReactionType, at_line(card.line) + card.name +
                                                       ": electrons / ionic species are not supported")
            .at(lines[i].pos);
        ok = false;
        continue;
      }
      const ElementInfo* e = find_element_by_symbol(s2);
      if (!e) {
        d.error(DiagCode::UnknownElement, at_line(card.line) + "unknown element '" + sym + "' in " + card.name).at(lines[i].pos);
        ok = false;
        continue;
      }
      if (std::fabs(n - std::round(n)) > 1e-9) {
        d.error(DiagCode::MechanismSyntax, at_line(card.line) + card.name + ": non-integer element count").at(lines[i].pos);
        ok = false;
        continue;
      }
      card.comp[e->z] += static_cast<std::int32_t>(std::lround(n));
    }
    double tl = 0, th = 0, tm = t_mid_default;
    if (!parse_number(col(L, 45, 10), tl) || !parse_number(col(L, 55, 10), th)) {
      d.error(DiagCode::MechanismSyntax, at_line(card.line) + card.name + ": temperature range missing (columns 46-65)").at(lines[i].pos);
      ok = false;
    }
    if (!trim(col(L, 65, 8)).empty() && !parse_number(col(L, 65, 8), tm)) {
      d.error(DiagCode::MechanismSyntax, at_line(card.line) + card.name + ": bad common temperature (columns 66-73)").at(lines[i].pos);
      ok = false;
    }
    std::array<double, 15> c{};
    for (std::size_t k = 0; k < 3 && ok; ++k) {
      const std::string& Lk = lines[i + 1 + k].text;
      const std::size_t nfields = k < 2 ? 5 : 4;
      for (std::size_t f = 0; f < nfields; ++f) {
        if (!parse_number(col(Lk, 15 * f, 15), c[5 * k + f])) {
          d.error(DiagCode::MechanismSyntax, at_line(lines[i + 1 + k].no) + card.name + ": bad coefficient in columns " +
                                                 std::to_string(15 * f + 1) + "-" + std::to_string(15 * f + 15))
              .at(lines[i + 1 + k].pos);
          ok = false;
          break;
        }
      }
    }
    if (ok) {
      card.nasa.t_low = tl;
      card.nasa.t_high = th;
      card.nasa.t_mid = tm;
      card.nasa.high = {c[0], c[1], c[2], c[3], c[4], c[5], c[6]};
      card.nasa.low = {c[7], c[8], c[9], c[10], c[11], c[12], c[13]};
      const std::string key = upper(card.name);
      if (out.count(key)) {
        d.warn(DiagCode::DuplicateSpecies, at_line(card.line) + "second thermo card for " + card.name + " ignored").at(lines[i].pos);
      } else {
        out.emplace(key, std::move(card));
      }
    }
    i += 4;
  }
}

// ---- species -----------------------------------------------------------------

std::int32_t add_smiles_species(Mechanism& m, const std::string& name, const std::string& smiles, std::size_t pos,
                                const std::string& where, Diagnostics& d) {
  const std::int32_t sp = parse_smiles(m.state, smiles, d, pos);
  if (sp < 0) return -1;
  for (std::size_t i = 0; i + 1 < m.state.species.size(); ++i) {
    const auto& other = m.state.species[i];
    if (!other.canonical_smiles.empty() && other.canonical_smiles == m.state.species.back().canonical_smiles) {
      d.warn(DiagCode::DuplicateSpecies, where + "'" + name + "' has the same structure as '" + m.names[i] + "' (" +
                                             other.canonical_smiles + "); fine for electronic states such as CH2/CH2(S)")
          .at(pos);
    }
  }
  m.names.push_back(name);
  return sp;
}

std::int32_t add_formula_species(Mechanism& m, const std::string& name, const Composition& comp) {
  Species s;
  s.composition = comp;
  s.formula = hill_formula(comp, 0);
  s.molar_mass = molar_mass(comp);
  s.atom_begin = static_cast<std::uint32_t>(m.state.atoms.size());
  s.bond_begin = static_cast<std::uint32_t>(m.state.bonds.size());
  const std::uint32_t sp = m.state.append_species(std::move(s));
  m.names.push_back(name);
  return static_cast<std::int32_t>(sp);
}

using SmilesMap = std::map<std::string, std::pair<std::string, std::size_t>>;  // NAME -> (SMILES, offset)

SmilesMap parse_smiles_map(std::string_view text, Diagnostics& d) {
  SmilesMap out;
  for (const auto& ln : split_lines(text)) {
    auto toks = split_ws(strip_comment(ln.text, true), ln.pos);
    if (toks.empty()) continue;
    if (toks.size() != 2) {
      d.error(DiagCode::MechanismSyntax, "SMILES map " + at_line(ln.no) + "expected 'NAME SMILES'").at(toks[0].pos);
      continue;
    }
    out[upper(toks[0].text)] = {toks[1].text, toks[1].pos};
  }
  return out;
}

// ---- reaction builder --------------------------------------------------------

struct RawRate {
  double A = 0.0, b = 0.0, Ea_J = 0.0;
  bool set = false;
};

struct Draft {
  bool active = false, ok = true;
  std::string equation;
  std::size_t pos = Diagnostic::npos, line = 0;
  double conc = 1e-6;
  std::vector<std::pair<std::uint32_t, double>> reactants, products;
  bool reversible = false;
  bool plus_m = false;
  bool falloff = false;     // "(+M)" or "(+SPECIES)"
  std::int64_t collider = -1;  // specific collider species for "(+SPECIES)"
  RawRate main, low, high, rev;
  std::optional<TroeParameters> troe;
  std::vector<std::pair<double, RawRate>> plog;  // (P in Pa, rate)
  std::vector<std::pair<std::uint32_t, double>> eff;
  double default_eff = 1.0;
  bool have_eff = false;
  bool duplicate = false;
  bool allow_negative_A = false;
};

void merge(std::vector<std::pair<std::uint32_t, double>>& side) {
  std::sort(side.begin(), side.end());
  std::vector<std::pair<std::uint32_t, double>> out;
  for (auto& p : side) {
    if (!out.empty() && out.back().first == p.first) out.back().second += p.second;
    else out.push_back(p);
  }
  side = std::move(out);
}

class ReactionBuilder {
 public:
  ReactionBuilder(Mechanism& m, Diagnostics& d) : m_(m), d_(d) {}

  Draft& draft() { return r_; }

  void error(DiagCode c, const std::string& msg) {
    d_.error(c, at_line(r_.line) + msg).at(r_.pos);
    r_.ok = false;
  }

  // `main.Ea_J` already in J/mol; `main.A` in source units with `conc` m^3 per amount unit.
  void begin(const std::string& equation, RawRate main, double conc, std::size_t pos, std::size_t line) {
    finish();
    r_ = Draft{};
    r_.active = true;
    r_.pos = pos;
    r_.line = line;
    r_.conc = conc;
    r_.main = main;
    std::string eq;
    for (char c : equation)
      if (!std::isspace(static_cast<unsigned char>(c))) eq += c;
    r_.equation = equation;
    std::size_t arrow = eq.find("<=>"), len = 3;
    r_.reversible = true;
    if (arrow == std::string::npos) {
      arrow = eq.find("=>");
      len = 2;
      r_.reversible = false;
      if (arrow == std::string::npos) {
        arrow = eq.find('=');
        len = 1;
        r_.reversible = true;
      }
    }
    if (arrow == std::string::npos) {
      error(DiagCode::MechanismSyntax, "reaction has no '=', '=>' or '<=>'");
      return;
    }
    std::string lhs = eq.substr(0, arrow), rhs = eq.substr(arrow + len);
    std::string fl, fr;
    const bool hl = strip_falloff(lhs, fl), hr = strip_falloff(rhs, fr);
    if (!r_.ok) return;
    if (hl != hr || fl != fr) {
      error(DiagCode::MechanismSyntax, "'(+M)' / '(+species)' must appear identically on both sides");
      return;
    }
    if (hl) {
      r_.falloff = true;
      if (upper(fl) != "M") {
        const std::int32_t c = m_.find(fl);
        if (c < 0) {
          error(DiagCode::UnknownSpecies, "unknown collider species '" + fl + "' in (+" + fl + ")");
          return;
        }
        r_.collider = c;
      }
    }
    int ml = 0, mr = 0;
    parse_side(lhs, r_.reactants, ml);
    parse_side(rhs, r_.products, mr);
    if (ml != mr || ml > 1) {
      error(DiagCode::MechanismSyntax, "third body M must appear exactly once on each side");
      return;
    }
    r_.plus_m = ml == 1;
  }

  // One auxiliary line: KEYWORD / values / groups and NAME / efficiency / groups.
  void aux(const std::string& text, std::size_t pos, std::size_t line, const Units& u) {
    if (!r_.active) {
      d_.error(DiagCode::MechanismSyntax, at_line(line) + "auxiliary data without a preceding reaction").at(pos);
      return;
    }
    std::string s = text;
    const std::string t = trim(s);
    if (lower(t).rfind("efficiencies:", 0) == 0) {  // .ocm legacy: "efficiencies: H2:2.5 H2O:12"
      s = t.substr(13);
      for (char& c : s)
        if (c == ':') c = '/';
      std::string rebuilt;
      for (const auto& tk : split_ws(s, 0)) rebuilt += tk.text + "/ ";
      s = rebuilt;
    }
    std::size_t i = 0;
    while (i < s.size()) {
      while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
      if (i >= s.size()) break;
      const std::size_t st = i;
      while (i < s.size() && !std::isspace(static_cast<unsigned char>(s[i])) && s[i] != '/') ++i;
      const std::string name = s.substr(st, i - st);
      while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
      std::string vals;
      bool slash = false;
      if (i < s.size() && s[i] == '/') {
        const std::size_t j = s.find('/', i + 1);
        if (j == std::string::npos) {
          d_.error(DiagCode::MechanismSyntax, at_line(line) + "unterminated '/' after '" + name + "'").at(pos + st);
          r_.ok = false;
          return;
        }
        vals = s.substr(i + 1, j - i - 1);
        i = j + 1;
        slash = true;
      }
      keyword(name, vals, slash, pos + st, line, u);
    }
  }

  void finish() {
    if (!r_.active) return;
    r_.active = false;
    if (!r_.ok) return;
    Draft& r = r_;
    merge(r.reactants);
    merge(r.products);
    double nf = 0, nr = 0;
    for (auto& [s, nu] : r.reactants) nf += nu;
    for (auto& [s, nu] : r.products) nr += nu;

    KineticReaction k;
    k.reactants = r.reactants;
    k.products = r.products;
    k.reversible = r.reversible;
    k.equation = r.equation;
    k.line = r.line;
    k.duplicate = r.duplicate;

    if (r.falloff && r.plus_m) return error(DiagCode::MechanismSyntax, "use either '+ M' or '(+M)', not both");
    if (r.falloff) {
      if (r.low.set && r.high.set) return error(DiagCode::MechanismSyntax, "fall-off reaction has both LOW and HIGH");
      if (!r.low.set && !r.high.set)
        return error(DiagCode::MechanismSyntax, "'(+M)' reaction needs a LOW (fall-off) or HIGH (chemically activated) line");
      k.kind = r.low.set ? RateKind::FallOff : RateKind::ChemicallyActivated;
    } else if (!r.plog.empty()) {
      if (r.plus_m) return error(DiagCode::MechanismSyntax, "PLOG reactions cannot have a third body");
      k.kind = RateKind::PressureLog;
    } else {
      k.kind = r.plus_m ? RateKind::ThirdBody : RateKind::Elementary;
    }
    if ((r.low.set || r.high.set || r.troe) && !r.falloff)
      return error(DiagCode::MechanismSyntax, "LOW / HIGH / TROE require a '(+M)' fall-off reaction");
    if (r.rev.set && (r.falloff || !r.plog.empty()))
      return error(DiagCode::UnsupportedReactionType, "explicit REV on fall-off or PLOG reactions is not supported");
    if (r.rev.set && !r.reversible) return error(DiagCode::MechanismSyntax, "REV given for an irreversible '=>' reaction");
    if (r.have_eff && !(r.plus_m || (r.falloff && r.collider < 0)))
      return error(DiagCode::MechanismSyntax, "third-body efficiencies need '+ M' or '(+M)'");
    if (!r.main.set && k.kind != RateKind::PressureLog)
      return error(DiagCode::MechanismSyntax, "missing rate constant");
    const bool neg_ok = r.duplicate || r.allow_negative_A;
    if (k.kind != RateKind::PressureLog && r.main.A < 0 && !neg_ok)
      return error(DiagCode::InvalidRateParameter, "negative A is only allowed on DUPLICATE reactions");
    if ((r.low.set && r.low.A < 0) || (r.high.set && r.high.A < 0))
      return error(DiagCode::InvalidRateParameter, "LOW / HIGH A must be >= 0");

    auto conv = [&](const RawRate& x, double order) {
      return ArrheniusRate{x.A * std::pow(r.conc, order - 1.0), x.b, x.Ea_J};
    };
    switch (k.kind) {
      case RateKind::Elementary: k.forward = conv(r.main, nf); break;
      case RateKind::ThirdBody: k.forward = conv(r.main, nf + 1); break;
      case RateKind::FallOff:
        k.forward = conv(r.main, nf);
        k.low_or_high = conv(r.low, nf + 1);
        break;
      case RateKind::ChemicallyActivated:
        k.forward = conv(r.main, nf + 1);
        k.low_or_high = conv(r.high, nf);
        break;
      case RateKind::PressureLog: {
        auto entries = r.plog;
        std::stable_sort(entries.begin(), entries.end(), [](auto& a, auto& b) { return a.first < b.first; });
        for (auto& [P, rate] : entries) {
          if (!(P > 0)) return error(DiagCode::InvalidRateParameter, "PLOG pressure must be > 0");
          if (k.plog.empty() || std::fabs(P - k.plog.back().P_Pa) > 1e-9 * P) k.plog.push_back({P, {}});
          k.plog.back().rates.push_back(conv(rate, nf));
        }
        for (const auto& lvl : k.plog) {
          double sum = 0;
          for (const auto& a : lvl.rates) sum += a.k(1000.0);
          if (!(sum > 0)) return error(DiagCode::InvalidRateParameter, "PLOG rate at one pressure must sum to > 0");
        }
        break;
      }
    }
    if (r.troe) k.troe = r.troe;
    if (r.rev.set) k.reverse = conv(r.rev, nr + (k.kind == RateKind::ThirdBody ? 1 : 0));
    k.third_body = k.kind == RateKind::ThirdBody || k.kind == RateKind::FallOff || k.kind == RateKind::ChemicallyActivated;
    if (k.third_body) {
      if (r.collider >= 0) {
        k.efficiency.assign(m_.species_count(), 0.0);
        k.efficiency[static_cast<std::size_t>(r.collider)] = 1.0;
      } else {
        k.efficiency.assign(m_.species_count(), r.default_eff);
        for (auto& [s, v] : r.eff) k.efficiency[s] = v;
      }
    }

    // Element and charge balance from the species compositions.
    std::array<double, kElementSlots> delta{};
    double dq = 0;
    for (auto& [s, nu] : k.reactants) {
      for (std::size_t z = 1; z < kElementSlots; ++z) delta[z] -= nu * m_.state.species[s].composition[z];
      dq -= nu * m_.state.species[s].net_charge;
    }
    for (auto& [s, nu] : k.products) {
      for (std::size_t z = 1; z < kElementSlots; ++z) delta[z] += nu * m_.state.species[s].composition[z];
      dq += nu * m_.state.species[s].net_charge;
    }
    bool balanced = true;
    for (std::size_t z = 1; z < kElementSlots; ++z) {
      if (std::fabs(delta[z]) < 1e-9) continue;
      balanced = false;
      char buf[48];
      std::snprintf(buf, sizeof buf, "%+g", delta[z]);
      d_.error(DiagCode::ElementImbalance, at_line(r.line) + std::string(element_symbol(static_cast<std::uint8_t>(z))) +
                                               " is not conserved in " + r.equation + " (products - reactants = " + buf + ")")
          .at(r.pos)
          .with_element(static_cast<std::uint8_t>(z), delta[z]);
    }
    if (std::fabs(dq) > 1e-9) {
      balanced = false;
      d_.error(DiagCode::ChargeImbalance, at_line(r.line) + "charge is not conserved in " + r.equation).at(r.pos).with_element(0, dq);
    }
    if (balanced) m_.reactions.push_back(std::move(k));
  }

  // CHEMKIN semantics: identical reactions must both carry DUPLICATE.
  void check_duplicates() {
    const auto& R = m_.reactions;
    std::vector<char> partner(R.size(), 0);
    std::size_t reported = 0;
    for (std::size_t i = 0; i < R.size(); ++i) {
      for (std::size_t j = i + 1; j < R.size(); ++j) {
        const bool same = R[i].kind == R[j].kind && R[i].third_body == R[j].third_body &&
                          ((R[i].reactants == R[j].reactants && R[i].products == R[j].products) ||
                           (R[i].reversible && R[j].reversible && R[i].reactants == R[j].products &&
                            R[i].products == R[j].reactants));
        if (!same) continue;
        partner[i] = partner[j] = 1;
        if (!(R[i].duplicate && R[j].duplicate) && reported++ < 20)
          d_.warn(DiagCode::DuplicateReaction, "line " + std::to_string(R[j].line) + ": " + R[j].equation +
                                                   " repeats line " + std::to_string(R[i].line) + " without DUPLICATE on both");
      }
    }
    for (std::size_t i = 0; i < R.size(); ++i)
      if (R[i].duplicate && !partner[i] && reported++ < 20)
        d_.warn(DiagCode::DuplicateReaction, "line " + std::to_string(R[i].line) + ": DUPLICATE given but no matching reaction");
  }

 private:
  bool strip_falloff(std::string& side, std::string& who) {
    const std::size_t a = side.find("(+");
    if (a == std::string::npos) return false;
    const std::size_t b = side.find(')', a);
    if (b == std::string::npos) {
      error(DiagCode::MechanismSyntax, "unterminated '(+'");
      return false;
    }
    who = side.substr(a + 2, b - a - 2);
    side.erase(a, b - a + 1);
    if (side.find("(+") != std::string::npos) error(DiagCode::MechanismSyntax, "more than one '(+M)' on one side");
    return true;
  }

  void parse_side(const std::string& side, std::vector<std::pair<std::uint32_t, double>>& out, int& third) {
    if (side.empty()) {
      error(DiagCode::EmptyReactionSide, "reaction side is empty");
      return;
    }
    std::size_t s = 0;
    while (s <= side.size()) {
      std::size_t e = side.find('+', s);
      if (e == std::string::npos) e = side.size();
      const std::string term = side.substr(s, e - s);
      s = e + 1;
      if (term.empty()) {
        error(DiagCode::MechanismSyntax, "empty term (stray '+')");
      } else {
        std::size_t k = 0;
        while (k < term.size() && (std::isdigit(static_cast<unsigned char>(term[k])) || term[k] == '.')) ++k;
        double coef = 1.0;
        std::string name = term.substr(k);
        if (k > 0 && (!parse_number(term.substr(0, k), coef) || coef <= 0)) {
          error(DiagCode::InvalidCoefficient, "invalid coefficient in '" + term + "'");
          coef = 1.0;
        }
        if (k > 0 && m_.find(name) < 0 && m_.find(term) >= 0) {  // names that start with a digit, e.g. 1-C4H8
          name = term;
          coef = 1.0;
        }
        if (upper(name) == "M") {
          ++third;
        } else {
          const std::int32_t sp = m_.find(name);
          if (sp < 0) error(DiagCode::UnknownSpecies, "unknown species '" + name + "'");
          else out.push_back({static_cast<std::uint32_t>(sp), coef});
        }
      }
      if (e == side.size()) break;
    }
  }

  void keyword(const std::string& name, const std::string& vals, bool slash, std::size_t pos, std::size_t line,
               const Units& u) {
    const std::string K = upper(name);
    auto fail = [&](DiagCode c, const std::string& msg) {
      d_.error(c, at_line(line) + msg).at(pos);
      r_.ok = false;
    };
    std::vector<double> num;
    for (const auto& tk : split_ws(vals, 0)) {
      double v;
      if (!parse_number(tk.text, v)) return fail(DiagCode::InvalidRateParameter, "'" + tk.text + "' is not a number in " + name);
      num.push_back(v);
    }
    if (K == "DUPLICATE" || K == "DUP") { r_.duplicate = true; return; }
    static const char* unsupported[] = {"SRI", "FORD", "RORD", "LT", "RLT", "JAN", "FIT1", "HV", "TDEP", "EXCI",
                                        "MOME", "XSMI", "UNITS", "CHEB", "PCHEB", "TCHEB", "USRPROG"};
    for (const char* u2 : unsupported)
      if (K == u2) return fail(DiagCode::UnsupportedReactionType, name + " is not supported yet");
    auto rate3 = [&](RawRate& out) {
      if (num.size() != 3) return fail(DiagCode::InvalidRateParameter, name + " needs 3 values: A b Ea");
      out = {num[0], num[1], num[2] * u.ea, true};
    };
    if (K == "LOW") return rate3(r_.low);
    if (K == "HIGH") return rate3(r_.high);
    if (K == "REV") return rate3(r_.rev);
    if (K == "TROE") {
      if (num.size() != 3 && num.size() != 4) return fail(DiagCode::InvalidRateParameter, "TROE needs 3 or 4 values: a T3 T1 [T2]");
      TroeParameters t;
      t.a = num[0];
      t.T3 = num[1];
      t.T1 = num[2];
      t.has_T2 = num.size() == 4;
      if (t.has_T2) t.T2 = num[3];
      r_.troe = t;
      return;
    }
    if (K == "PLOG") {
      if (num.size() != 4) return fail(DiagCode::InvalidRateParameter, "PLOG needs 4 values: P A b Ea");
      r_.plog.push_back({num[0] * u.pressure, RawRate{num[1], num[2], num[3] * u.ea, true}});
      return;
    }
    // Otherwise: third-body efficiency NAME/value/.
    const std::int32_t sp = m_.find(name);
    if (sp < 0) return fail(DiagCode::UnknownSpecies, "unknown keyword or species '" + name + "'");
    if (!slash || num.size() != 1 || num[0] < 0)
      return fail(DiagCode::InvalidRateParameter, "efficiency for " + name + " must be one number >= 0 between slashes");
    r_.eff.push_back({static_cast<std::uint32_t>(sp), num[0]});
    r_.have_eff = true;
  }

  Mechanism& m_;
  Diagnostics& d_;
  Draft r_;
};

// Reaction line: EQUATION A b Ea  (the last three tokens must be numbers).
void reaction_line(ReactionBuilder& b, const std::vector<Tok>& toks, const Units& u, std::size_t line, Diagnostics& d) {
  std::size_t numeric_tail = 0;
  for (std::size_t i = toks.size(); i-- > 0;) {
    double v;
    if (!parse_number(toks[i].text, v)) break;
    ++numeric_tail;
  }
  if (toks.size() < 4 || numeric_tail < 3) {
    b.finish();
    d.error(DiagCode::MechanismSyntax, at_line(line) + "reaction lines are 'EQUATION A b Ea' (three rate parameters)").at(toks[0].pos);
    return;
  }
  double A, beta, Ea;
  parse_number(toks[toks.size() - 3].text, A);
  parse_number(toks[toks.size() - 2].text, beta);
  parse_number(toks[toks.size() - 1].text, Ea);
  std::string eq;
  for (std::size_t i = 0; i + 3 < toks.size(); ++i) eq += (i ? " " : "") + toks[i].text;
  b.begin(eq, RawRate{A, beta, Ea * u.ea, true}, u.conc, toks[0].pos, line);
}

// Thermo for every species: mechanism cards > built-in database (canonical SMILES) > Joback.
void assign_thermo(Mechanism& m, const std::map<std::string, Card>& cards, Diagnostics& d) {
  m.thermo.assign(m.species_count(), std::nullopt);
  std::vector<char> needs(m.species_count(), 0);
  for (const auto& r : m.reactions) {
    if (!r.reversible || r.reverse) continue;
    for (auto& [s, nu] : r.reactants) needs[s] = 1;
    for (auto& [s, nu] : r.products) needs[s] = 1;
  }
  for (std::uint32_t i = 0; i < m.species_count(); ++i) {
    const Species& sp = m.state.species[i];
    auto card = cards.find(upper(m.names[i]));
    if (card != cards.end()) {
      if (sp.atom_count > 0 && sp.composition != card->second.comp) {
        d.error(DiagCode::CompositionMismatch, m.names[i] + ": SMILES gives " + sp.formula + " but its thermo card (line " +
                                                   std::to_string(card->second.line) + ") gives " +
                                                   hill_formula(card->second.comp, 0))
            .on_species(i);
      }
      SpeciesThermo t;
      t.source = ThermoSource::Nasa7Database;
      t.nasa = card->second.nasa;
      t.label = "thermo card " + m.names[i];
      m.thermo[i] = t;
      continue;
    }
    if (sp.atom_count > 0) {
      Diagnostics local;
      m.thermo[i] = resolve_thermo(m.state, i, kTref, local);
      if (!needs[i]) continue;
      for (const auto& x : local.items()) {
        if (x.code == DiagCode::NoJobackGroup || x.code == DiagCode::TemperatureOutOfRange) continue;
        (x.severity == Severity::Error ? d.error(x.code, m.names[i] + ": " + x.message)
                                       : d.warn(x.code, m.names[i] + ": " + x.message))
            .on_species(i);
      }
    } else if (needs[i]) {
      d.error(DiagCode::NoThermoData, m.names[i] + " has no thermo data, needed for reverse rates").on_species(i);
    }
  }
}

// ---- Cantera YAML helpers ------------------------------------------------------

bool yaml_number(const yaml::Node* n, double& v, std::string* unit = nullptr) {
  if (!n || !n->is_scalar()) return false;
  auto toks = split_ws(n->scalar, 0);
  if (toks.empty() || !parse_number(toks[0].text, v)) return false;
  if (unit) {
    unit->clear();
    for (std::size_t i = 1; i < toks.size(); ++i) *unit += (i > 1 ? " " : "") + toks[i].text;
  } else if (toks.size() > 1) {
    return false;
  }
  return true;
}

std::optional<double> energy_unit(const std::string& u) {
  const std::string t = lower(u);
  if (t == "cal/mol") return 4.184;
  if (t == "kcal/mol") return 4184.0;
  if (t == "j/mol") return 1.0;
  if (t == "kj/mol") return 1000.0;
  if (t == "j/kmol") return 1e-3;
  if (t == "kj/kmol") return 1.0;
  if (t == "cal/kmol") return 4.184e-3;
  if (t == "kcal/kmol") return 4.184;
  if (t == "k") return kGasConstant;
  if (t == "ev") return 96485.33212;
  return std::nullopt;
}

std::optional<double> pressure_unit(const std::string& u) {
  const std::string t = lower(u);
  if (t == "pa") return 1.0;
  if (t == "kpa") return 1e3;
  if (t == "mpa") return 1e6;
  if (t == "bar") return 1e5;
  if (t == "atm") return kAtm;
  if (t == "torr") return 133.32236842;
  return std::nullopt;
}

}  // namespace

// ---- .ocm --------------------------------------------------------------------

std::optional<Mechanism> parse_mechanism(std::string_view text, Diagnostics& diag) {
  Mechanism m;
  const std::size_t errors_before = diag.error_count();
  enum class Section { None, Species, Reactions, Thermo } section = Section::None;
  Units units;
  ReactionBuilder builder(m, diag);
  std::vector<Line> thermo_lines;
  std::map<std::string, Card> cards;

  for (const Line& ln : split_lines(text)) {
    if (section == Section::Thermo) {
      const std::string t = upper(trim(ln.text));
      if (t == "END") { section = Section::None; continue; }
      thermo_lines.push_back(ln);
      continue;
    }
    const std::string body = strip_comment(ln.text, true);
    auto toks = split_ws(body, ln.pos);
    if (toks.empty()) continue;
    const std::string head = lower(toks[0].text);
    const std::string where = at_line(ln.no);
    if (head == "end") {
      if (section == Section::Reactions) builder.finish();
      section = Section::None;
      continue;
    }
    if (section == Section::None) {
      if (head == "species") {
        section = Section::Species;
        if (toks.size() > 1) diag.error(DiagCode::MechanismSyntax, where + "put species on their own lines after 'species'").at(toks[1].pos);
      } else if (head == "thermo") {
        section = Section::Thermo;
      } else if (head == "reactions") {
        section = Section::Reactions;
        for (std::size_t i = 1; i < toks.size(); ++i) {
          const UnitResult res = apply_unit(toks[i].text, units);
          if (res == UnitResult::Unknown) diag.error(DiagCode::MechanismSyntax, where + "unknown unit '" + toks[i].text + "'").at(toks[i].pos);
          if (res == UnitResult::Unsupported) diag.error(DiagCode::UnsupportedReactionType, where + "unit '" + toks[i].text + "' is not supported").at(toks[i].pos);
        }
      } else {
        diag.error(DiagCode::MechanismSyntax, where + "expected 'species', 'thermo' or 'reactions' section").at(toks[0].pos);
      }
      continue;
    }
    if (section == Section::Species) {
      if (toks.size() != 2) {
        diag.error(DiagCode::MechanismSyntax, where + "species lines are 'NAME SMILES'").at(toks[0].pos);
        continue;
      }
      if (m.find(toks[0].text) >= 0) {
        diag.error(DiagCode::DuplicateSpecies, where + "species '" + toks[0].text + "' is declared twice").at(toks[0].pos);
        continue;
      }
      add_smiles_species(m, toks[0].text, toks[1].text, toks[1].pos, where, diag);
      continue;
    }
    // reactions section
    if (body.find('=') != std::string::npos) reaction_line(builder, toks, units, ln.no, diag);
    else builder.aux(body, ln.pos, ln.no, units);
  }
  builder.finish();
  if (section != Section::None) diag.error(DiagCode::MechanismSyntax, "missing 'end' after the last section").at(text.size());
  if (m.names.empty()) diag.error(DiagCode::MechanismSyntax, "mechanism declares no species").at(0);
  parse_cards(thermo_lines, cards, diag);
  builder.check_duplicates();
  assign_thermo(m, cards, diag);
  if (diag.error_count() != errors_before) return std::nullopt;
  return m;
}

// ---- CHEMKIN-II --------------------------------------------------------------

std::optional<Mechanism> import_chemkin(std::string_view chem, std::string_view thermo, std::string_view smiles_map,
                                        Diagnostics& diag) {
  Mechanism m;
  const std::size_t errors_before = diag.error_count();
  enum class Section { None, Elements, Species, Thermo, Reactions } section = Section::None;
  std::vector<Tok> species_names;
  std::vector<Line> inline_thermo, reaction_lines;
  bool thermo_all = false;
  Units units;
  bool saw_elements = false;

  for (const Line& ln : split_lines(chem)) {
    if (section == Section::Thermo) {
      const auto toks = split_ws(strip_comment(ln.text, false), 0);
      if (!toks.empty() && upper(toks[0].text) == "END") { section = Section::None; continue; }
      inline_thermo.push_back(ln);
      continue;
    }
    const std::string body = strip_comment(ln.text, false);
    auto toks = split_ws(body, ln.pos);
    if (toks.empty()) continue;
    std::size_t i = 0;
    if (section == Section::None) {
      const std::string K = upper(toks[0].text);
      i = 1;
      if (K.rfind("ELEM", 0) == 0) { section = Section::Elements; saw_elements = true; }
      else if (K.rfind("SPEC", 0) == 0) section = Section::Species;
      else if (K.rfind("THER", 0) == 0) {
        section = Section::Thermo;
        thermo_all = toks.size() > 1 && upper(toks[1].text) == "ALL";
        continue;
      } else if (K.rfind("REAC", 0) == 0) {
        section = Section::Reactions;
        for (std::size_t k = 1; k < toks.size(); ++k) {
          const UnitResult res = apply_unit(toks[k].text, units);
          if (res == UnitResult::Unknown) diag.error(DiagCode::MechanismSyntax, at_line(ln.no) + "unknown unit '" + toks[k].text + "'").at(toks[k].pos);
          if (res == UnitResult::Unsupported) diag.error(DiagCode::UnsupportedReactionType, at_line(ln.no) + "unit '" + toks[k].text + "' is not supported").at(toks[k].pos);
        }
        continue;
      } else {
        diag.error(DiagCode::MechanismSyntax, at_line(ln.no) + "expected ELEMENTS, SPECIES, THERMO or REACTIONS").at(toks[0].pos);
        continue;
      }
    }
    if (section == Section::Reactions) {
      if (upper(toks[0].text) == "END") { section = Section::None; continue; }
      reaction_lines.push_back(ln);
      continue;
    }
    for (; i < toks.size(); ++i) {
      if (upper(toks[i].text) == "END") { section = Section::None; break; }
      if (section == Section::Elements) {
        std::string s = toks[i].text;
        const std::size_t slash = s.find('/');
        if (slash != std::string::npos) s = s.substr(0, slash);
        std::string norm = lower(s);
        if (!norm.empty()) norm[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(norm[0])));
        if (norm == "E") diag.error(DiagCode::UnsupportedReactionType, at_line(ln.no) + "electrons / ionic species are not supported").at(toks[i].pos);
        else if (!find_element_by_symbol(norm)) diag.error(DiagCode::UnknownElement, at_line(ln.no) + "unknown element '" + s + "'").at(toks[i].pos);
      } else if (section == Section::Species) {
        species_names.push_back(toks[i]);
      }
    }
  }
  if (!saw_elements) diag.warn(DiagCode::MechanismSyntax, "no ELEMENTS section; elements are taken from species compositions");
  if (section != Section::None) diag.error(DiagCode::MechanismSyntax, "missing END after the last section").at(chem.size());

  std::map<std::string, Card> cards;
  if (!thermo_all && !thermo.empty()) parse_cards(split_lines(thermo), cards, diag);
  {
    std::map<std::string, Card> inline_cards;
    parse_cards(inline_thermo, inline_cards, diag);
    for (auto& [k, v] : inline_cards) cards[k] = v;  // inline THERMO overrides the external file
  }
  const SmilesMap smap = parse_smiles_map(smiles_map, diag);

  for (const Tok& t : species_names) {
    if (m.find(t.text) >= 0) {
      diag.warn(DiagCode::DuplicateSpecies, "species '" + t.text + "' listed twice; second ignored").at(t.pos);
      continue;
    }
    const std::string key = upper(t.text);
    auto sm = smap.find(key);
    if (sm != smap.end()) {
      add_smiles_species(m, t.text, sm->second.first, sm->second.second, "SMILES map: ", diag);
    } else if (auto c = cards.find(key); c != cards.end()) {
      add_formula_species(m, t.text, c->second.comp);
    } else {
      diag.error(DiagCode::UnknownSpecies, "species '" + t.text + "' has neither a thermo card nor a SMILES entry; composition unknown").at(t.pos);
    }
  }

  ReactionBuilder builder(m, diag);
  for (const Line& ln : reaction_lines) {
    const std::string body = strip_comment(ln.text, false);
    auto toks = split_ws(body, ln.pos);
    if (toks.empty()) continue;
    if (body.find('=') != std::string::npos) reaction_line(builder, toks, units, ln.no, diag);
    else builder.aux(body, ln.pos, ln.no, units);
  }
  builder.finish();
  builder.check_duplicates();
  assign_thermo(m, cards, diag);
  if (m.names.empty()) diag.error(DiagCode::MechanismSyntax, "no species").at(0);
  if (diag.error_count() != errors_before) return std::nullopt;
  return m;
}

// ---- Cantera YAML ------------------------------------------------------------

std::optional<Mechanism> import_cantera_yaml(std::string_view text, std::string_view smiles_map, Diagnostics& diag) {
  const std::size_t errors_before = diag.error_count();
  yaml::Node root;
  yaml::ParseError perr;
  if (!yaml::parse(text, root, perr)) {
    diag.error(DiagCode::MechanismSyntax, "YAML line " + std::to_string(perr.line) + ": " + perr.message);
    return std::nullopt;
  }
  auto yerr = [&](DiagCode c, std::size_t line, const std::string& msg) { diag.error(c, "YAML line " + std::to_string(line) + ": " + msg); };
  if (!root.is_map()) {
    yerr(DiagCode::MechanismSyntax, 1, "top level must be a mapping");
    return std::nullopt;
  }

  // Units (Cantera defaults: m, kmol, J/kmol, Pa).
  double length = 1.0, quantity = 1000.0, ea_unit = 1e-3, p_unit = 1.0;
  if (const yaml::Node* u = root.get("units"); u && u->is_map()) {
    for (const auto& [k, v] : u->map) {
      const std::string val = lower(v.scalar);
      if (k == "length") {
        if (val == "cm") length = 0.01;
        else if (val == "m") length = 1.0;
        else if (val == "mm") length = 1e-3;
        else yerr(DiagCode::MechanismSyntax, v.line, "unsupported length unit '" + v.scalar + "'");
      } else if (k == "quantity") {
        if (val == "mol") quantity = 1.0;
        else if (val == "kmol") quantity = 1000.0;
        else yerr(DiagCode::UnsupportedReactionType, v.line, "unsupported quantity unit '" + v.scalar + "'");
      } else if (k == "activation-energy") {
        auto f = energy_unit(v.scalar);
        if (f) ea_unit = *f;
        else yerr(DiagCode::MechanismSyntax, v.line, "unsupported activation-energy unit '" + v.scalar + "'");
      } else if (k == "pressure") {
        auto f = pressure_unit(v.scalar);
        if (f) p_unit = *f;
        else yerr(DiagCode::MechanismSyntax, v.line, "unsupported pressure unit '" + v.scalar + "'");
      }
    }
  }
  const double conc = length * length * length / quantity;

  // Species restricted to the first phase's list, if one is given.
  std::vector<std::string> phase_species;
  if (const yaml::Node* ph = root.get("phases"); ph && ph->is_seq() && !ph->seq.empty())
    if (const yaml::Node* sl = ph->seq[0].get("species"); sl && sl->is_seq())
      for (const auto& s : sl->seq)
        if (s.is_scalar()) phase_species.push_back(s.scalar);

  Mechanism m;
  std::map<std::string, Card> cards;
  const SmilesMap smap = parse_smiles_map(smiles_map, diag);
  const yaml::Node* species = root.get("species");
  if (!species || !species->is_seq()) {
    yerr(DiagCode::MechanismSyntax, root.line, "missing 'species' list");
    return std::nullopt;
  }
  for (const yaml::Node& s : species->seq) {
    const yaml::Node* name = s.get("name");
    if (!name || !name->is_scalar()) { yerr(DiagCode::MechanismSyntax, s.line, "species without a name"); continue; }
    if (!phase_species.empty() && std::find(phase_species.begin(), phase_species.end(), name->scalar) == phase_species.end()) continue;
    Card card;
    card.name = name->scalar;
    card.line = s.line;
    bool ok = true;
    if (const yaml::Node* comp = s.get("composition"); comp && comp->is_map()) {
      for (const auto& [sym, cnt] : comp->map) {
        double n = 0;
        const ElementInfo* e = find_element_by_symbol(sym);
        if (sym == "E") { yerr(DiagCode::UnsupportedReactionType, s.line, card.name + ": electrons are not supported"); ok = false; }
        else if (!e) { yerr(DiagCode::UnknownElement, s.line, "unknown element '" + sym + "' in " + card.name); ok = false; }
        else if (!yaml_number(&cnt, n) || std::fabs(n - std::round(n)) > 1e-9) { yerr(DiagCode::MechanismSyntax, s.line, card.name + ": bad count for " + sym); ok = false; }
        else card.comp[e->z] += static_cast<std::int32_t>(std::lround(n));
      }
    } else {
      yerr(DiagCode::MechanismSyntax, s.line, card.name + ": missing composition");
      ok = false;
    }
    const yaml::Node* th = s.get("thermo");
    bool have_thermo = false;
    if (th && th->is_map()) {
      const yaml::Node* model = th->get("model");
      const yaml::Node* ranges = th->get("temperature-ranges");
      const yaml::Node* data = th->get("data");
      if (!model || model->scalar != "NASA7") {
        yerr(DiagCode::UnsupportedReactionType, th->line, card.name + ": only NASA7 thermo is supported");
        ok = false;
      } else if (!ranges || !ranges->is_seq() || !data || !data->is_seq() || (ranges->seq.size() != data->seq.size() + 1) ||
                 data->seq.empty() || data->seq.size() > 2) {
        yerr(DiagCode::MechanismSyntax, th->line, card.name + ": NASA7 needs 1-2 data arrays and matching temperature-ranges");
        ok = false;
      } else {
        std::vector<double> T;
        for (const auto& x : ranges->seq) { double v = 0; if (!yaml_number(&x, v)) ok = false; T.push_back(v); }
        std::array<std::array<double, 7>, 2> a{};
        for (std::size_t r = 0; r < data->seq.size(); ++r) {
          if (!data->seq[r].is_seq() || data->seq[r].seq.size() != 7) { ok = false; break; }
          for (std::size_t c = 0; c < 7; ++c) if (!yaml_number(&data->seq[r].seq[c], a[r][c])) ok = false;
        }
        if (!ok) {
          yerr(DiagCode::MechanismSyntax, th->line, card.name + ": malformed NASA7 numbers");
        } else {
          card.nasa.t_low = T.front();
          card.nasa.t_high = T.back();
          card.nasa.t_mid = data->seq.size() == 2 ? T[1] : T.back();
          card.nasa.low = a[0];
          card.nasa.high = data->seq.size() == 2 ? a[1] : a[0];
          have_thermo = true;
        }
      }
    }
    if (!ok) continue;
    if (m.find(card.name) >= 0) { yerr(DiagCode::DuplicateSpecies, s.line, "species '" + card.name + "' defined twice"); continue; }
    auto sm = smap.find(upper(card.name));
    if (sm != smap.end()) add_smiles_species(m, card.name, sm->second.first, sm->second.second, "SMILES map: ", diag);
    else add_formula_species(m, card.name, card.comp);
    if (have_thermo || sm != smap.end()) cards.emplace(upper(card.name), card);
    if (sm != smap.end() && !have_thermo) cards.erase(upper(card.name));
    if (sm != smap.end() && m.find(card.name) >= 0) {
      const Species& last = m.state.species[static_cast<std::size_t>(m.find(card.name))];
      if (last.composition != card.comp)
        yerr(DiagCode::CompositionMismatch, s.line, card.name + ": SMILES gives " + last.formula + " but composition is " + hill_formula(card.comp, 0));
    }
  }

  ReactionBuilder b(m, diag);
  const yaml::Node* reactions = root.get("reactions");
  if (reactions && reactions->is_seq()) {
    for (const yaml::Node& rx : reactions->seq) {
      const yaml::Node* eq = rx.get("equation");
      if (!eq || !eq->is_scalar()) { yerr(DiagCode::MechanismSyntax, rx.line, "reaction without 'equation'"); continue; }
      std::string type = "elementary";
      if (const yaml::Node* t = rx.get("type")) type = t->scalar;
      auto rate = [&](const char* key, RawRate& out) -> bool {
        const yaml::Node* n = rx.get(key);
        if (!n) return false;
        if (!n->is_map()) { yerr(DiagCode::MechanismSyntax, n->line, std::string(key) + " must be a mapping {A, b, Ea}"); return false; }
        double A = 0, beta = 0, Ea = 0;
        std::string unit;
        if (!yaml_number(n->get("A"), A) || !yaml_number(n->get("b"), beta) || !yaml_number(n->get("Ea"), Ea, &unit)) {
          yerr(DiagCode::InvalidRateParameter, n->line, std::string(key) + ": A and b must be plain numbers, Ea a number with optional unit");
          return false;
        }
        double f = ea_unit;
        if (!unit.empty()) {
          auto u2 = energy_unit(unit);
          if (!u2) { yerr(DiagCode::MechanismSyntax, n->line, "unsupported energy unit '" + unit + "'"); return false; }
          f = *u2;
        }
        out = {A, beta, Ea * f, true};
        return true;
      };
      RawRate main;
      if (type == "elementary" || type == "three-body") rate("rate-constant", main);
      else if (type == "falloff") rate("high-P-rate-constant", main);
      else if (type == "chemically-activated") rate("low-P-rate-constant", main);
      else if (type != "pressure-dependent-Arrhenius") {
        yerr(DiagCode::UnsupportedReactionType, rx.line, "reaction type '" + type + "' is not supported");
        continue;
      }
      b.begin(eq->scalar, main, conc, Diagnostic::npos, rx.line);
      Draft& d = b.draft();
      if (type == "falloff") { RawRate lo; if (rate("low-P-rate-constant", lo)) d.low = lo; else b.error(DiagCode::MechanismSyntax, "falloff needs low-P-rate-constant"); }
      if (type == "chemically-activated") { RawRate hi; if (rate("high-P-rate-constant", hi)) d.high = hi; else b.error(DiagCode::MechanismSyntax, "chemically-activated needs high-P-rate-constant"); }
      if (const yaml::Node* troe = rx.get("Troe"); troe && troe->is_map()) {
        TroeParameters t;
        if (!yaml_number(troe->get("A"), t.a) || !yaml_number(troe->get("T3"), t.T3) || !yaml_number(troe->get("T1"), t.T1))
          b.error(DiagCode::InvalidRateParameter, "Troe needs A, T3, T1");
        if (troe->get("T2")) { t.has_T2 = yaml_number(troe->get("T2"), t.T2); }
        d.troe = t;
      }
      if (rx.get("SRI")) b.error(DiagCode::UnsupportedReactionType, "SRI fall-off is not supported yet");
      if (rx.get("orders")) b.error(DiagCode::UnsupportedReactionType, "non-elementary reaction orders are not supported yet");
      if (type == "pressure-dependent-Arrhenius") {
        const yaml::Node* rcs = rx.get("rate-constants");
        if (!rcs || !rcs->is_seq() || rcs->seq.empty()) b.error(DiagCode::MechanismSyntax, "pressure-dependent-Arrhenius needs rate-constants");
        else
          for (const auto& rc : rcs->seq) {
            double P = 0, A = 0, beta = 0, Ea = 0;
            std::string pu, eu;
            if (!yaml_number(rc.get("P"), P, &pu) || !yaml_number(rc.get("A"), A) || !yaml_number(rc.get("b"), beta) ||
                !yaml_number(rc.get("Ea"), Ea, &eu)) {
              b.error(DiagCode::InvalidRateParameter, "rate-constants entries need P, A, b, Ea");
              break;
            }
            auto pf = pu.empty() ? std::optional<double>(p_unit) : pressure_unit(pu);
            auto ef = eu.empty() ? std::optional<double>(ea_unit) : energy_unit(eu);
            if (!pf || !ef) { b.error(DiagCode::MechanismSyntax, "unsupported unit in rate-constants"); break; }
            d.plog.push_back({P * *pf, RawRate{A, beta, Ea * *ef, true}});
          }
      }
      if (const yaml::Node* de = rx.get("default-efficiency")) {
        double v = 1;
        if (yaml_number(de, v)) { d.default_eff = v; d.have_eff = true; }
      }
      if (const yaml::Node* eff = rx.get("efficiencies"); eff && eff->is_map()) {
        for (const auto& [sp, v] : eff->map) {
          double val = 0;
          const std::int32_t idx = m.find(sp);
          if (idx < 0) continue;  // Cantera ignores efficiencies of species not in the phase
          if (!yaml_number(&v, val) || val < 0) { b.error(DiagCode::InvalidRateParameter, "bad efficiency for " + sp); continue; }
          d.eff.push_back({static_cast<std::uint32_t>(idx), val});
          d.have_eff = true;
        }
      }
      if (const yaml::Node* dup = rx.get("duplicate")) d.duplicate = lower(dup->scalar) == "true";
      if (const yaml::Node* neg = rx.get("negative-A")) d.allow_negative_A = lower(neg->scalar) == "true";
      b.finish();
    }
  }
  b.check_duplicates();
  assign_thermo(m, cards, diag);
  if (m.names.empty()) diag.error(DiagCode::MechanismSyntax, "no species");
  if (diag.error_count() != errors_before) return std::nullopt;
  return m;
}

// ---- dispatch ----------------------------------------------------------------

MechanismFormat detect_mechanism_format(std::string_view text) noexcept {
  for (const Line& ln : split_lines(text)) {
    const std::string& s = ln.text;
    if (s.empty() || s[0] == ' ' || s[0] == '\t' || s[0] == '#' || s[0] == '!') continue;
    if (s.rfind("species:", 0) == 0 || s.rfind("phases:", 0) == 0 || s.rfind("reactions:", 0) == 0 ||
        s.rfind("units:", 0) == 0 || s.rfind("description:", 0) == 0 || s.rfind("generator:", 0) == 0)
      return MechanismFormat::CanteraYaml;
    const auto toks = split_ws(strip_comment(s, false), 0);
    if (toks.empty()) continue;
    const std::string K = upper(toks[0].text);
    if (K.rfind("ELEM", 0) == 0) return MechanismFormat::Chemkin;
    if (K == "SPECIES" || K == "SPEC") return toks.size() > 1 ? MechanismFormat::Chemkin : MechanismFormat::OmniChem;
  }
  return MechanismFormat::OmniChem;
}

std::optional<Mechanism> parse_any_mechanism(std::string_view text, std::string_view thermo, std::string_view smiles_map,
                                             MechanismFormat format, Diagnostics& diag) {
  if (format == MechanismFormat::Auto) format = detect_mechanism_format(text);
  switch (format) {
    case MechanismFormat::Chemkin: return import_chemkin(text, thermo, smiles_map, diag);
    case MechanismFormat::CanteraYaml: return import_cantera_yaml(text, smiles_map, diag);
    default:
      if (!thermo.empty() || !smiles_map.empty())
        diag.warn(DiagCode::InvalidInput, ".ocm mechanisms carry their own SMILES and thermo; extra thermo / SMILES inputs ignored");
      return parse_mechanism(text, diag);
  }
}

namespace {
bool read_file(const std::string& path, std::string& out, Diagnostics& diag) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    diag.error(DiagCode::InvalidInput, "cannot open '" + path + "'");
    return false;
  }
  std::stringstream ss;
  ss << in.rdbuf();
  out = ss.str();
  return true;
}
}  // namespace

std::optional<Mechanism> load_mechanism(const std::string& path, Diagnostics& diag, const std::string& thermo_path,
                                        const std::string& smiles_path) {
  std::string text, thermo, smiles;
  if (!read_file(path, text, diag)) return std::nullopt;
  if (!thermo_path.empty() && !read_file(thermo_path, thermo, diag)) return std::nullopt;
  if (!smiles_path.empty() && !read_file(smiles_path, smiles, diag)) return std::nullopt;
  return parse_any_mechanism(text, thermo, smiles, MechanismFormat::Auto, diag);
}

}  // namespace omnichem
