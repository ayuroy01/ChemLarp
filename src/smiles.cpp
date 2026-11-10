#include "omnichem/smiles.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "omnichem/canonical.hpp"
#include "omnichem/elements.hpp"
#include "omnichem/formula.hpp"
#include "omnichem/rings.hpp"

namespace omnichem {
namespace {

constexpr std::uint8_t kOrderImplicit = 0;
constexpr std::uint8_t kOrderAromatic = 5;
constexpr std::size_t kKekuleStepLimit = 200000;
constexpr std::int64_t kRefPending = -1;  // ring-closure slot not yet closed
constexpr std::int64_t kRefH = -2;        // implicit H of a bracket atom
constexpr std::int64_t kRefLP = -3;       // lone pair of a 3-coordinate stereocentre

struct PAtom {
  std::uint8_t z = 0;  // 0 = placeholder for an atom that failed to parse
  std::int8_t charge = 0;
  int hcount = 0;  // explicit for bracket atoms, computed for organic-subset atoms
  bool bracket = false;
  bool aromatic = false;
  std::uint16_t isotope = 0;
  Chirality chiral = Chirality::None;
  std::uint16_t map = 0;
  std::size_t pos = 0;
  bool has_prev = false;
  std::vector<std::int64_t> refs;  // neighbours in order of appearance (OpenSMILES chirality order)
};

struct PBond {
  std::uint32_t a = 0, b = 0;
  std::uint8_t order = kOrderImplicit;
  std::uint8_t flags = 0;
  std::size_t pos = 0;
};

struct Pending {
  bool present = false;
  std::uint8_t order = kOrderImplicit;
  std::uint8_t flags = 0;
  std::size_t pos = 0;
};

struct RingOpen {
  bool open = false;
  std::uint32_t atom = 0;
  Pending bond;
  std::size_t pos = 0;
  std::size_t slot = 0;  // index in atoms[atom].refs reserved for the partner
};

struct LocalTetra {
  std::uint32_t atom;
  std::array<std::int64_t, 4> refs;
  Chirality parity;
};

struct LocalDoubleBond {
  std::uint32_t bond, ref_a, ref_b;
  bool cis;
};

bool is_digit(char c) { return c >= '0' && c <= '9'; }
bool is_lower(char c) { return c >= 'a' && c <= 'z'; }
bool is_upper(char c) { return c >= 'A' && c <= 'Z'; }

std::string quoted(std::string_view s) { return "'" + std::string(s) + "'"; }

std::string describe_atom(const PAtom& a) {
  std::string s(element_symbol(a.z));
  if (a.charge != 0) s += (a.charge > 0 ? "+" : "") + std::to_string(a.charge);
  return s;
}

class Parser {
 public:
  Parser(std::string_view s, Diagnostics& d, std::size_t base) : s_(s), d_(d), base_(base) {}

  std::vector<PAtom> atoms;
  std::vector<PBond> bonds;

  std::vector<LocalTetra> tetra;
  std::vector<LocalDoubleBond> dbstereo;

  bool parse();     // syntax; returns true when no syntax error was found
  bool finalize();  // aromaticity, implicit H, kekulization, valence, stereo

 private:
  Diagnostic& error(DiagCode code, std::size_t pos, std::string msg) {
    ++errors_;
    return d_.error(code, std::move(msg)).at(base_ + pos);
  }

  std::size_t parse_bracket(std::size_t i);
  std::size_t parse_organic(std::size_t i);
  void ring_closure(int r, std::size_t at);
  void attach(std::uint32_t idx);
  void derive_stereo();
  int ring_size_through(std::uint32_t bond_index) const;
  bool has_bond(std::uint32_t a, std::uint32_t b) const {
    for (const auto& bd : bonds)
      if ((bd.a == a && bd.b == b) || (bd.a == b && bd.b == a)) return true;
    return false;
  }

  std::string_view s_;
  Diagnostics& d_;
  std::size_t base_;
  std::size_t errors_ = 0;

  std::int64_t prev_ = -1;
  Pending pending_;
  std::array<RingOpen, 100> rings_{};
};

void Parser::attach(std::uint32_t idx) {
  if (prev_ >= 0) {
    PBond b;
    b.a = static_cast<std::uint32_t>(prev_);
    b.b = idx;
    b.order = pending_.present ? pending_.order : kOrderImplicit;
    b.flags = pending_.flags;
    b.pos = pending_.present ? pending_.pos : atoms[idx].pos;
    bonds.push_back(b);
    auto& r = atoms[idx].refs;
    r.insert(r.begin(), prev_);  // the preceding atom comes first, before any [C@H] hydrogen
    atoms[idx].has_prev = true;
    atoms[static_cast<std::size_t>(prev_)].refs.push_back(idx);
  }
  pending_ = {};
  prev_ = idx;
}

void Parser::ring_closure(int r, std::size_t at) {
  if (prev_ < 0) {
    error(DiagCode::UnexpectedCharacter, at, "ring-closure digit must follow an atom");
    pending_ = {};
    return;
  }
  RingOpen& ro = rings_[static_cast<std::size_t>(r)];
  if (!ro.open) {
    ro.open = true;
    ro.atom = static_cast<std::uint32_t>(prev_);
    ro.bond = pending_;
    ro.pos = at;
    ro.slot = atoms[static_cast<std::size_t>(prev_)].refs.size();
    atoms[static_cast<std::size_t>(prev_)].refs.push_back(kRefPending);
    pending_ = {};
    return;
  }
  const std::uint32_t a = ro.atom;
  const std::uint32_t b = static_cast<std::uint32_t>(prev_);
  Pending chosen = pending_.present ? pending_ : ro.bond;
  // A '/' or '\\' written at the closing digit reads from the closing atom to
  // the opening one; directions are stored in order of appearance, so flip it.
  if (pending_.present && (chosen.flags & (kBondUp | kBondDown)))
    chosen.flags = static_cast<std::uint8_t>(chosen.flags ^ (kBondUp | kBondDown));
  if (ro.bond.present && pending_.present && ro.bond.order != pending_.order) {
    error(DiagCode::RingBondConflict, at,
          "ring bond " + std::to_string(r) + " has conflicting bond symbols at its two ends");
  }
  if (a == b) {
    error(DiagCode::RingSelfBond, at, "ring bond " + std::to_string(r) + " closes onto the same atom");
  } else if (has_bond(a, b)) {
    error(DiagCode::DuplicateBond, at,
          "ring bond " + std::to_string(r) + " duplicates an existing bond between the same two atoms");
  } else {
    PBond bd;
    bd.a = a;
    bd.b = b;
    bd.order = chosen.present ? chosen.order : kOrderImplicit;
    bd.flags = static_cast<std::uint8_t>(chosen.flags | kBondRingClosure);
    bd.pos = at;
    bonds.push_back(bd);
    atoms[a].refs[ro.slot] = b;
    atoms[b].refs.push_back(a);
  }
  ro = {};
  pending_ = {};
}

std::size_t Parser::parse_organic(std::size_t i) {
  const char c = s_[i];
  const char nx = i + 1 < s_.size() ? s_[i + 1] : '\0';
  PAtom a;
  a.pos = i;
  std::string_view sym;
  std::size_t len = 1;
  switch (c) {
    case 'B': if (nx == 'r') { sym = "Br"; len = 2; } else { sym = "B"; } break;
    case 'C': if (nx == 'l') { sym = "Cl"; len = 2; } else { sym = "C"; } break;
    case 'N': sym = "N"; break;
    case 'O': sym = "O"; break;
    case 'P': sym = "P"; break;
    case 'S': sym = "S"; break;
    case 'F': sym = "F"; break;
    case 'I': sym = "I"; break;
    case 'b': sym = "B"; a.aromatic = true; break;
    case 'c': sym = "C"; a.aromatic = true; break;
    case 'n': sym = "N"; a.aromatic = true; break;
    case 'o': sym = "O"; a.aromatic = true; break;
    case 'p': sym = "P"; a.aromatic = true; break;
    case 's': sym = "S"; a.aromatic = true; break;
    default: {
      if (c == 'H') {
        error(DiagCode::UnknownElement, i, "hydrogen outside brackets; write it as [H] or as an H count, e.g. [CH4]");
      } else if (is_lower(c) && i > 0 && is_upper(s_[i - 1]) &&
                 find_element_by_symbol(std::string{s_[i - 1], c})) {
        // e.g. "NaCl": 'N' was read as nitrogen, then 'a' is left over.
        const std::string tok{s_[i - 1], c};
        error(DiagCode::UnknownElement, i - 1,
              "element " + quoted(tok) + " is outside the organic subset; write it as [" + tok + "]");
      } else {
        if (is_upper(c) && is_lower(nx)) len = 2;
        std::string tok(s_.substr(i, len));
        if (len == 2 && !find_element_by_symbol(tok)) { len = 1; tok = std::string(1, c); }
        if (find_element_by_symbol(tok)) {
          error(DiagCode::UnknownElement, i,
                "element " + quoted(tok) + " is outside the organic subset; write it as [" + tok + "]");
        } else {
          error(DiagCode::UnknownElement, i, "unknown element symbol " + quoted(tok));
        }
      }
      atoms.push_back(a);  // placeholder keeps the graph shape for further syntax checking
      return i + len;
    }
  }
  a.z = find_element_by_symbol(sym)->z;
  atoms.push_back(a);
  return i + len;
}

std::size_t Parser::parse_bracket(std::size_t i) {
  const std::size_t close = s_.find(']', i);
  if (close == std::string_view::npos) {
    error(DiagCode::UnclosedBracketAtom, i, "'[' is never closed by ']'");
    return s_.size();
  }
  const std::string_view body = s_.substr(i + 1, close - i - 1);
  const std::size_t off = i + 1;  // position of body[0] in s_
  PAtom a;
  a.bracket = true;
  a.pos = i;
  std::size_t k = 0;

  // isotope
  if (k < body.size() && is_digit(body[k])) {
    unsigned iso = 0;
    std::size_t start = k;
    while (k < body.size() && is_digit(body[k])) iso = iso * 10 + static_cast<unsigned>(body[k++] - '0');
    if (k - start > 3) error(DiagCode::InvalidBracketAtom, off + start, "isotope number has more than 3 digits");
    a.isotope = static_cast<std::uint16_t>(std::min(iso, 999u));
  }

  // element symbol
  if (k >= body.size()) {
    error(DiagCode::InvalidBracketAtom, i, "bracket atom has no element symbol");
    atoms.push_back(a);
    return close + 1;
  }
  {
    const char c = body[k];
    const char nx = k + 1 < body.size() ? body[k + 1] : '\0';
    if (c == '*') {
      error(DiagCode::UnsupportedFeature, off + k, "wildcard atom '*' is not supported");
      ++k;
    } else if (is_upper(c)) {
      std::string two{c, nx};
      if (is_lower(nx) && find_element_by_symbol(two)) {
        a.z = find_element_by_symbol(two)->z;
        k += 2;
      } else if (const ElementInfo* e = find_element_by_symbol(std::string(1, c))) {
        a.z = e->z;
        k += 1;
      } else {
        error(DiagCode::UnknownElement, off + k,
              "unknown element symbol " + quoted(is_lower(nx) ? two : std::string(1, c)));
        k += is_lower(nx) ? 2 : 1;
      }
    } else if (is_lower(c)) {
      std::string two{c, nx};
      std::string upper2 = two;
      upper2[0] = static_cast<char>(upper2[0] - 'a' + 'A');
      if ((two == "se" || two == "as" || two == "te") && find_element_by_symbol(upper2)) {
        a.z = find_element_by_symbol(upper2)->z;
        a.aromatic = true;
        k += 2;
      } else if (c == 'b' || c == 'c' || c == 'n' || c == 'o' || c == 'p' || c == 's') {
        a.z = find_element_by_symbol(std::string(1, static_cast<char>(c - 'a' + 'A')))->z;
        a.aromatic = true;
        k += 1;
      } else {
        error(DiagCode::InvalidAromaticAtom, off + k, quoted(std::string(1, c)) + " is not a valid aromatic element");
        k += 1;
      }
    } else {
      error(DiagCode::InvalidBracketAtom, off + k,
            "expected element symbol in bracket atom, found " + quoted(std::string(1, c)));
      atoms.push_back(a);
      return close + 1;
    }
  }

  // chirality
  if (k < body.size() && body[k] == '@') {
    ++k;
    a.chiral = Chirality::CounterClockwise;
    if (k < body.size() && body[k] == '@') { ++k; a.chiral = Chirality::Clockwise; }
    if (k + 1 < body.size() && is_upper(body[k]) && is_upper(body[k + 1]) && body[k] != 'H') {
      error(DiagCode::UnsupportedFeature, off + k,
            "extended chirality classes (@TH, @AL, @SP, @TB, @OH) are not supported");
      k += 2;
      while (k < body.size() && is_digit(body[k])) ++k;
    }
  }

  // hydrogen count
  if (k < body.size() && body[k] == 'H') {
    ++k;
    int h = 1;
    if (k < body.size() && is_digit(body[k])) {
      h = 0;
      while (k < body.size() && is_digit(body[k])) h = h * 10 + (body[k++] - '0');
    }
    if (h > 99) error(DiagCode::InvalidBracketAtom, off + k, "hydrogen count is unreasonably large");
    a.hcount = std::min(h, 99);
    for (int k2 = 0; k2 < a.hcount; ++k2) a.refs.push_back(kRefH);
  }

  // charge
  if (k < body.size() && (body[k] == '+' || body[k] == '-')) {
    const char sign = body[k];
    const std::size_t start = k++;
    int mag = 1;
    if (k < body.size() && is_digit(body[k])) {
      mag = 0;
      while (k < body.size() && is_digit(body[k])) mag = mag * 10 + (body[k++] - '0');
    } else {
      while (k < body.size() && body[k] == sign) { ++mag; ++k; }
    }
    if (mag > 15) {
      error(DiagCode::InvalidBracketAtom, off + start, "formal charge magnitude exceeds 15");
      mag = 15;
    }
    a.charge = static_cast<std::int8_t>(sign == '+' ? mag : -mag);
  }

  // atom class
  if (k < body.size() && body[k] == ':') {
    ++k;
    unsigned m = 0;
    const std::size_t start = k;
    while (k < body.size() && is_digit(body[k])) m = m * 10 + static_cast<unsigned>(body[k++] - '0');
    if (k == start) error(DiagCode::InvalidBracketAtom, off + k, "atom class ':' must be followed by digits");
    a.map = static_cast<std::uint16_t>(std::min(m, 65535u));
  }

  if (k < body.size()) {
    error(DiagCode::InvalidBracketAtom, off + k,
          "unexpected " + quoted(std::string(1, body[k])) + " inside bracket atom");
  }
  if (a.aromatic && a.z != 0) {
    const ElementInfo* e = find_element_by_z(a.z);
    if (!e || !e->aromatic_allowed)
      error(DiagCode::InvalidAromaticAtom, i, std::string(element_symbol(a.z)) + " cannot be aromatic");
  }
  atoms.push_back(a);
  return close + 1;
}

bool Parser::parse() {
  const std::size_t n = s_.size();
  if (n == 0) {
    error(DiagCode::EmptyInput, 0, "empty SMILES string");
    return false;
  }
  struct Branch { std::int64_t atom; std::size_t pos; std::size_t atoms_at_open; };
  std::vector<Branch> branches;

  std::size_t i = 0;
  while (i < n) {
    const char c = s_[i];
    switch (c) {
      case '(':
        if (prev_ < 0) error(DiagCode::UnexpectedCharacter, i, "branch '(' must follow an atom");
        if (pending_.present) {
          error(DiagCode::BondWithoutAtom, pending_.pos, "bond symbol before '(' has no atom to bond to");
          pending_ = {};
        }
        branches.push_back({prev_, i, atoms.size()});
        ++i;
        continue;
      case ')':
        if (branches.empty()) {
          error(DiagCode::UnmatchedCloseParen, i, "')' has no matching '('");
          ++i;
          continue;
        }
        if (pending_.present) {
          error(DiagCode::BondWithoutAtom, pending_.pos, "bond symbol at end of branch has no atom to bond to");
          pending_ = {};
        }
        if (atoms.size() == branches.back().atoms_at_open)
          error(DiagCode::EmptyBranch, branches.back().pos, "empty branch '()'");
        prev_ = branches.back().atom;
        branches.pop_back();
        ++i;
        continue;
      case '-': case '=': case '#': case '$': case ':': case '/': case '\\': {
        if (pending_.present) {
          error(DiagCode::UnexpectedCharacter, i, "two consecutive bond symbols");
        } else if (prev_ < 0) {
          error(DiagCode::BondWithoutAtom, i, "bond symbol " + quoted(std::string(1, c)) + " has no preceding atom");
        } else {
          pending_.present = true;
          pending_.pos = i;
          pending_.flags = 0;
          switch (c) {
            case '-': pending_.order = 1; break;
            case '=': pending_.order = 2; break;
            case '#': pending_.order = 3; break;
            case '$': pending_.order = 4; break;
            case ':': pending_.order = kOrderAromatic; break;
            case '/': pending_.order = 1; pending_.flags = kBondUp; break;
            default: pending_.order = 1; pending_.flags = kBondDown; break;
          }
        }
        ++i;
        continue;
      }
      case '.':
        if (pending_.present) {
          error(DiagCode::BondWithoutAtom, pending_.pos, "bond symbol before '.' has no atom to bond to");
          pending_ = {};
        }
        if (prev_ < 0) error(DiagCode::UnexpectedCharacter, i, "'.' must follow an atom");
        prev_ = -1;
        ++i;
        continue;
      case '%': {
        if (i + 2 < n && is_digit(s_[i + 1]) && is_digit(s_[i + 2])) {
          ring_closure((s_[i + 1] - '0') * 10 + (s_[i + 2] - '0'), i);
          i += 3;
        } else {
          error(DiagCode::UnexpectedCharacter, i, "'%' must be followed by exactly two digits");
          ++i;
        }
        continue;
      }
      case '[': {
        const std::size_t before = atoms.size();
        const std::size_t next = parse_bracket(i);
        if (atoms.size() > before) attach(static_cast<std::uint32_t>(atoms.size() - 1));
        i = next;
        continue;
      }
      case '*': {
        error(DiagCode::UnsupportedFeature, i, "wildcard atom '*' is not supported");
        PAtom a;
        a.pos = i;
        atoms.push_back(a);
        attach(static_cast<std::uint32_t>(atoms.size() - 1));
        ++i;
        continue;
      }
      default:
        if (is_digit(c)) {
          ring_closure(c - '0', i);
          ++i;
        } else if (is_upper(c) || is_lower(c)) {
          i = parse_organic(i);
          attach(static_cast<std::uint32_t>(atoms.size() - 1));
        } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
          error(DiagCode::UnexpectedCharacter, i, "whitespace inside SMILES");
          ++i;
        } else {
          error(DiagCode::UnexpectedCharacter, i, "unexpected character " + quoted(std::string(1, c)));
          ++i;
        }
        continue;
    }
  }

  if (pending_.present) error(DiagCode::BondWithoutAtom, pending_.pos, "trailing bond symbol has no atom to bond to");
  for (const auto& br : branches) error(DiagCode::UnmatchedOpenParen, br.pos, "'(' is never closed");
  for (std::size_t r = 0; r < rings_.size(); ++r)
    if (rings_[r].open)
      error(DiagCode::UnclosedRing, rings_[r].pos, "ring bond " + std::to_string(r) + " is opened but never closed");
  if (atoms.empty() && errors_ == 0) error(DiagCode::EmptyInput, 0, "SMILES contains no atoms");
  return errors_ == 0;
}

bool Parser::finalize() {
  const std::size_t n = atoms.size();
  std::vector<int> sum(n, 0);         // bond-order sum with aromatic bonds counted as 1
  std::vector<int> arom_bonds(n, 0);
  const std::size_t errors_before = errors_;

  // 1. resolve implicit / aromatic bond orders
  for (auto& b : bonds) {
    const PAtom& A = atoms[b.a];
    const PAtom& B = atoms[b.b];
    if (b.order == kOrderImplicit) b.order = (A.aromatic && B.aromatic) ? kOrderAromatic : 1;
    if (b.order == kOrderAromatic) {
      if (!A.aromatic || !B.aromatic) {
        error(DiagCode::InvalidAromaticAtom, b.pos, "aromatic bond ':' requires both atoms to be aromatic");
      }
      b.flags |= kBondAromatic;
      ++arom_bonds[b.a];
      ++arom_bonds[b.b];
      sum[b.a] += 1;
      sum[b.b] += 1;
    } else {
      sum[b.a] += b.order;
      sum[b.b] += b.order;
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (atoms[i].aromatic && arom_bonds[i] == 0) {
      error(DiagCode::InvalidAromaticAtom, atoms[i].pos,
            "aromatic atom " + quoted(describe_atom(atoms[i])) + " is not part of any aromatic bond (is it in a ring?)")
          .on_atom(static_cast<std::int64_t>(i));
    }
  }
  if (errors_ != errors_before) return false;

  // 2. implicit hydrogens and which aromatic atoms still need a double bond
  std::vector<char> needs(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    PAtom& a = atoms[i];
    if (a.bracket) {
      if (a.aromatic) {
        const int v = lowest_default_valence(a.z, a.charge);
        if (v > 0 && v - sum[i] - a.hcount >= 1) needs[i] = 1;
      }
      continue;
    }
    if (a.aromatic) {
      const int v = lowest_default_valence(a.z, 0);
      const int free = v - sum[i];
      if (free >= 1) {
        needs[i] = 1;
        a.hcount = free - 1;
      } else {
        a.hcount = 0;
      }
    } else {
      a.hcount = 0;
      for (std::uint8_t v : default_valences(a.z, 0)) {
        if (v >= sum[i]) { a.hcount = v - sum[i]; break; }
      }
    }
  }

  // 3. kekulize: perfect matching of `needs` atoms over aromatic bonds
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> nb(n);
  for (std::uint32_t k = 0; k < bonds.size(); ++k) {
    const auto& b = bonds[k];
    if ((b.flags & kBondAromatic) && needs[b.a] && needs[b.b]) {
      nb[b.a].push_back({b.b, k});
      nb[b.b].push_back({b.a, k});
    }
  }
  std::vector<std::int64_t> mate(n, -1);
  std::size_t steps = 0;
  bool limit_hit = false;
  std::function<bool()> solve = [&]() -> bool {
    if (++steps > kKekuleStepLimit) { limit_hit = true; return false; }
    std::int64_t best = -1;
    int best_opts = INT_MAX;
    for (std::size_t v = 0; v < n; ++v) {
      if (!needs[v] || mate[v] >= 0) continue;
      int opts = 0;
      for (auto [w, bk] : nb[v]) opts += (mate[w] < 0);
      if (opts < best_opts) { best = static_cast<std::int64_t>(v); best_opts = opts; if (opts == 0) break; }
    }
    if (best < 0) return true;
    if (best_opts == 0) return false;
    for (auto [w, bk] : nb[static_cast<std::size_t>(best)]) {
      if (mate[w] >= 0) continue;
      mate[static_cast<std::size_t>(best)] = bk;
      mate[w] = bk;
      if (solve()) return true;
      mate[static_cast<std::size_t>(best)] = -1;
      mate[w] = -1;
      if (limit_hit) return false;
    }
    return false;
  };
  const bool kekule_ok = solve();
  if (!kekule_ok) {
    std::size_t first = 0;
    std::size_t unmatched = 0;
    for (std::size_t v = 0; v < n; ++v)
      if (needs[v]) { if (unmatched == 0) first = v; ++unmatched; }
    if (limit_hit) {
      error(DiagCode::KekulizationLimit, atoms[first].pos, "aromatic system too large to kekulize within the search limit");
    } else {
      error(DiagCode::KekulizationFailed, atoms[first].pos,
            "cannot assign alternating double bonds to the aromatic system; "
            "pyrrole-type nitrogens must be written [nH]")
          .on_atom(static_cast<std::int64_t>(first));
    }
  }
  for (std::uint32_t k = 0; k < bonds.size(); ++k) {
    auto& b = bonds[k];
    if (b.flags & kBondAromatic) b.order = (kekule_ok && mate[b.a] == static_cast<std::int64_t>(k)) ? 2 : 1;
  }

  // 4. valence validation on the Kekulé structure (aromatic atoms are skipped if
  //    kekulization failed, so the remaining atoms still get checked and reported)
  std::vector<int> total(n, 0);
  for (const auto& b : bonds) { total[b.a] += b.order; total[b.b] += b.order; }
  for (std::size_t i = 0; i < n; ++i) {
    const PAtom& a = atoms[i];
    if (!kekule_ok && a.aromatic) continue;
    const int v = total[i] + a.hcount;
    const int mx = max_valence(a.z, a.charge);
    if (mx >= 0 && v > mx) {
      error(DiagCode::ValenceExceeded, a.pos,
            "atom " + quoted(describe_atom(a)) + " has valence " + std::to_string(v) +
                "; the maximum allowed is " + std::to_string(mx))
          .on_atom(static_cast<std::int64_t>(i));
    }
  }
  if (errors_ != errors_before) return false;
  derive_stereo();
  return errors_ == errors_before;
}

// Shortest ring containing the bond (0 if none): BFS between its ends without it.
int Parser::ring_size_through(std::uint32_t bond_index) const {
  const std::uint32_t src = bonds[bond_index].a, dst = bonds[bond_index].b;
  std::vector<int> dist(atoms.size(), -1);
  std::vector<std::uint32_t> queue{src};
  dist[src] = 0;
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::uint32_t x = queue[head];
    for (std::uint32_t k = 0; k < bonds.size(); ++k) {
      if (k == bond_index) continue;
      std::uint32_t y;
      if (bonds[k].a == x) y = bonds[k].b;
      else if (bonds[k].b == x) y = bonds[k].a;
      else continue;
      if (dist[y] >= 0) continue;
      dist[y] = dist[x] + 1;
      if (y == dst) return dist[y] + 1;
      queue.push_back(y);
    }
  }
  return 0;
}

// Tetrahedral centres from @/@@ and cis/trans from '/' '\'.
//
// Directional bonds: for a double-bond end e and a neighbour x joined by a
// directional single bond, define t(e,x) = s * (x written before e ? +1 : -1)
// with s = +1 for '/' and -1 for '\'. Two neighbours on the same end lie on
// opposite sides iff their t differ; neighbours a (of u) and b (of v) are cis
// iff t(u,a) == t(v,b). Atom index order equals order of appearance.
void Parser::derive_stereo() {
  for (std::uint32_t i = 0; i < atoms.size(); ++i) {
    PAtom& a = atoms[i];
    if (a.chiral == Chirality::None) continue;
    auto refs = a.refs;
    if (refs.size() == 3 && a.hcount == 0) refs.insert(refs.begin() + (a.has_prev ? 1 : 0), kRefLP);
    if (refs.size() != 4) {
      d_.warn(DiagCode::StereoIgnored, "chirality on " + describe_atom(a) + " with " + std::to_string(refs.size()) +
                                           " neighbours is not tetrahedral; ignored")
          .at(base_ + a.pos)
          .on_atom(i);
      a.chiral = Chirality::None;
      continue;
    }
    tetra.push_back({i, {refs[0], refs[1], refs[2], refs[3]}, a.chiral});
  }

  auto t_of = [&](std::uint32_t e, std::uint32_t k) {
    const PBond& b = bonds[k];
    const std::uint32_t x = b.a == e ? b.b : b.a;
    const int s = (b.flags & kBondUp) ? 1 : -1;
    return std::pair<std::uint32_t, int>{x, s * (x < e ? 1 : -1)};
  };
  for (std::uint32_t k = 0; k < bonds.size(); ++k) {
    const PBond& db = bonds[k];
    if (db.order != 2 || (db.flags & kBondAromatic)) continue;
    std::array<std::vector<std::pair<std::uint32_t, int>>, 2> marks;
    const std::array<std::uint32_t, 2> ends{db.a, db.b};
    for (int side = 0; side < 2; ++side) {
      for (std::uint32_t j = 0; j < bonds.size(); ++j) {
        if (j == k || bonds[j].order != 1 || !(bonds[j].flags & (kBondUp | kBondDown))) continue;
        if (bonds[j].a != ends[side] && bonds[j].b != ends[side]) continue;
        marks[side].push_back(t_of(ends[side], j));
      }
    }
    if (marks[0].empty() || marks[1].empty()) continue;
    bool conflict = false;
    for (int side = 0; side < 2; ++side) {
      if (marks[side].size() == 2 && marks[side][0].second == marks[side][1].second) {
        error(DiagCode::ConflictingBondStereo, atoms[ends[side]].pos,
              "conflicting '/' '\\' marks: both neighbours of this double-bond atom are placed on the same side");
        conflict = true;
      }
    }
    if (conflict) continue;
    const int ring = ring_size_through(k);
    if (ring > 0 && ring < 8) {
      d_.warn(DiagCode::StereoIgnored, "cis/trans marks on a double bond in a " + std::to_string(ring) +
                                           "-membered ring are implied by the ring; ignored")
          .at(base_ + db.pos);
      continue;
    }
    dbstereo.push_back({k, marks[0][0].first, marks[1][0].first, marks[0][0].second == marks[1][0].second});
  }
}

// Tarjan bridge-finding: a bond is in a ring iff it is not a bridge.
void perceive_rings(ChemicalSystemState& st, const Species& s) {
  const std::uint32_t base = s.atom_begin, n = s.atom_count;
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> adj(n);
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    adj[st.bonds.begin[b] - base].push_back({st.bonds.end[b] - base, b});
    adj[st.bonds.end[b] - base].push_back({st.bonds.begin[b] - base, b});
  }
  std::vector<int> disc(n, -1), low(n, 0);
  std::vector<char> bridge(s.bond_count, 0);
  struct Frame { std::uint32_t v; std::uint32_t parent_bond; std::size_t it; };
  std::vector<Frame> stack;
  int timer = 0;
  for (std::uint32_t root = 0; root < n; ++root) {
    if (disc[root] >= 0) continue;
    disc[root] = low[root] = timer++;
    stack.push_back({root, UINT32_MAX, 0});
    while (!stack.empty()) {
      Frame& f = stack.back();
      if (f.it < adj[f.v].size()) {
        auto [w, b] = adj[f.v][f.it++];
        if (b == f.parent_bond) continue;
        if (disc[w] < 0) {
          disc[w] = low[w] = timer++;
          stack.push_back({w, b, 0});
        } else {
          low[f.v] = std::min(low[f.v], disc[w]);
        }
      } else {
        const Frame done = f;
        stack.pop_back();
        if (!stack.empty()) {
          Frame& p = stack.back();
          low[p.v] = std::min(low[p.v], low[done.v]);
          if (low[done.v] > disc[p.v]) bridge[done.parent_bond - s.bond_begin] = 1;
        }
      }
    }
  }
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    if (bridge[b - s.bond_begin]) continue;
    st.bonds.flags[b] |= kBondInRing;
    st.atoms.flags[st.bonds.begin[b]] |= kAtomInRing;
    st.atoms.flags[st.bonds.end[b]] |= kAtomInRing;
  }
}

}  // namespace

std::int32_t parse_smiles(ChemicalSystemState& state, std::string_view smiles, Diagnostics& diag,
                          std::size_t base_offset) {
  Parser p(smiles, diag, base_offset);
  if (!p.parse()) return -1;
  if (!p.finalize()) return -1;

  const std::uint32_t atom_base = static_cast<std::uint32_t>(state.atoms.size());
  const std::uint32_t bond_base = static_cast<std::uint32_t>(state.bonds.size());
  const std::uint32_t sp = static_cast<std::uint32_t>(state.species.size());

  state.atoms.reserve(atom_base + p.atoms.size());
  for (const PAtom& a : p.atoms) {
    Atom out;
    out.element = a.z;
    out.formal_charge = a.charge;
    out.implicit_h = static_cast<std::uint8_t>(a.hcount);
    out.flags = static_cast<std::uint8_t>((a.aromatic ? kAtomAromatic : 0) | (a.bracket ? kAtomBracket : 0));
    out.isotope = a.isotope;
    out.chirality = a.chiral;
    out.atom_map = a.map;
    out.species = sp;
    state.atoms.push(out);
  }
  state.bonds.reserve(bond_base + p.bonds.size());
  for (const PBond& b : p.bonds) {
    state.bonds.push(Bond{atom_base + b.a, atom_base + b.b, b.order, b.flags, sp});
  }
  auto global_ref = [&](std::int64_t r) -> std::uint32_t {
    if (r == kRefH) return kImplicitHRef;
    if (r == kRefLP) return kLonePairRef;
    return atom_base + static_cast<std::uint32_t>(r);
  };
  for (const LocalTetra& t : p.tetra) {
    TetrahedralStereo g;
    g.atom = atom_base + t.atom;
    for (int k = 0; k < 4; ++k) g.refs[static_cast<std::size_t>(k)] = global_ref(t.refs[static_cast<std::size_t>(k)]);
    g.parity = t.parity;
    state.tetrahedral.push_back(g);
  }
  for (const LocalDoubleBond& d : p.dbstereo)
    state.double_bonds.push_back({bond_base + d.bond, atom_base + d.ref_a, atom_base + d.ref_b, d.cis});

  Species s;
  s.smiles = std::string(smiles);
  s.atom_begin = atom_base;
  s.atom_count = static_cast<std::uint32_t>(p.atoms.size());
  s.bond_begin = bond_base;
  s.bond_count = static_cast<std::uint32_t>(p.bonds.size());
  for (const PAtom& a : p.atoms) s.net_charge += a.charge;
  perceive_rings(state, s);
  state.append_species(std::move(s));

  Species& added = state.species[sp];
  added.composition = compute_composition(state, sp);
  added.formula = hill_formula(added.composition, added.net_charge);
  added.molar_mass = molar_mass(added.composition);
  state.rebuild_adjacency();
  perceive_aromaticity(state, sp);
  added.canonical_smiles = canonical_smiles(state, sp, diag);
  return static_cast<std::int32_t>(sp);
}

}  // namespace omnichem
