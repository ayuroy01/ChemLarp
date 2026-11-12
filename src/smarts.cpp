#include "omnichem/smarts.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "omnichem/elements.hpp"
#include "omnichem/rings.hpp"

namespace omnichem {

// ---- compiled query -------------------------------------------------------------------

enum class Op : std::uint8_t { Prim, Not, AndHi, Or, AndLo };

enum class AP : std::uint8_t {  // atom primitives
  True, Element, Aromatic, Aliphatic, AtomicNum, Degree, TotalH, ImplicitH, RingCount, RingSize, Valence,
  Connectivity, RingConn, Charge, Isotope, Recursive,
};
enum class BP : std::uint8_t { Default, Single, Double, Triple, Aromatic, Any, Ring };

struct Node {
  Op op = Op::Prim;
  std::uint8_t prim = 0;  // AP or BP
  int value = 0;          // primitive argument; -1 = "any positive" where the primitive allows it
  int aromatic = -1;      // Element: -1 either, 0 aliphatic, 1 aromatic
  int a = -1, b = -1;     // children
};

struct QAtom {
  int expr = -1;
  bool has_prev = false;   // written after another atom (SMILES neighbour order)
  int chiral = 0;          // 0 none, 1 '@', 2 '@@'
  std::vector<int> refs;   // neighbours in SMILES order; -1 = the hydrogen written in the bracket
  std::size_t pos = 0;
};

struct QBond {
  int a = 0, b = 0, expr = -1;
};

}  // namespace omnichem

struct omnichem::SmartsQueryData {
  std::vector<Node> nodes;
  std::vector<QAtom> atoms;
  std::vector<QBond> bonds;
  std::vector<std::shared_ptr<const SmartsQueryData>> recursive;
};

namespace omnichem {

std::size_t SmartsQuery::atom_count() const noexcept { return data_ ? data_->atoms.size() : 0; }

namespace {

// ---- parser -----------------------------------------------------------------------------

class Parser {
 public:
  Parser(std::string_view s, Diagnostics& d, std::size_t base) : s_(s), d_(d), base_(base) {}

  bool run(SmartsQueryData& q);

 private:
  bool fail(std::size_t at, const std::string& msg, DiagCode code = DiagCode::SmartsSyntax) {
    d_.error(code, msg).at(base_ + at);
    ok_ = false;
    return false;
  }
  int node(Node n) {
    q_->nodes.push_back(n);
    return static_cast<int>(q_->nodes.size()) - 1;
  }
  int prim(AP p, int v = 0, int arom = -1) {
    Node n;
    n.prim = static_cast<std::uint8_t>(p);
    n.value = v;
    n.aromatic = arom;
    return node(n);
  }
  int bprim(BP p) {
    Node n;
    n.prim = static_cast<std::uint8_t>(p);
    return node(n);
  }
  int binary(Op op, int a, int b) {
    if (a < 0) return b;
    if (b < 0) return a;
    Node n;
    n.op = op;
    n.a = a;
    n.b = b;
    return node(n);
  }

  bool digit(std::size_t i) const { return i < s_.size() && s_[i] >= '0' && s_[i] <= '9'; }
  int number(std::size_t& i, int fallback) {
    if (!digit(i)) return fallback;
    int v = 0;
    while (digit(i)) v = v * 10 + (s_[i++] - '0');
    return v;
  }

  // atom expressions inside [...]
  int atom_lo(std::size_t& i, std::size_t end, QAtom& qa, bool& first);
  int atom_or(std::size_t& i, std::size_t end, QAtom& qa, bool& first);
  int atom_hi(std::size_t& i, std::size_t end, QAtom& qa, bool& first);
  int atom_not(std::size_t& i, std::size_t end, QAtom& qa, bool& first);
  int atom_prim(std::size_t& i, std::size_t end, QAtom& qa, bool& first);

  // bond expressions
  static bool bond_char(char c) {
    return c == '-' || c == '=' || c == '#' || c == ':' || c == '~' || c == '@' || c == '/' || c == '\\' || c == '!' ||
           c == '&' || c == ',' || c == ';';
  }
  int bond_lo(std::size_t& i);
  int bond_or(std::size_t& i);
  int bond_hi(std::size_t& i);
  int bond_not(std::size_t& i);

  int add_atom(int expr, std::size_t pos, int chiral, bool chiral_h);
  void connect(int bond_expr);
  void ring_closure(int label, std::size_t at);

  std::string_view s_;
  Diagnostics& d_;
  std::size_t base_;
  SmartsQueryData* q_ = nullptr;
  bool ok_ = true;
  int prev_ = -1;
  int pending_bond_ = -1;
  bool has_pending_ = false;
  struct Open { int atom; int bond; std::size_t slot; std::size_t pos; };
  std::map<int, Open> rings_;
  int bracket_h_ = -1;
};

int Parser::add_atom(int expr, std::size_t pos, int chiral, bool chiral_h) {
  QAtom a;
  a.expr = expr;
  a.chiral = chiral;
  a.pos = pos;
  const int idx = static_cast<int>(q_->atoms.size());
  a.has_prev = prev_ >= 0;
  if (prev_ >= 0) a.refs.push_back(prev_);
  if (chiral && chiral_h) a.refs.push_back(-1);
  q_->atoms.push_back(a);
  if (prev_ >= 0) {
    q_->atoms[static_cast<std::size_t>(prev_)].refs.push_back(idx);
    QBond b;
    b.a = prev_;
    b.b = idx;
    b.expr = has_pending_ ? pending_bond_ : bprim(BP::Default);
    q_->bonds.push_back(b);
  }
  has_pending_ = false;
  pending_bond_ = -1;
  prev_ = idx;
  return idx;
}

void Parser::ring_closure(int label, std::size_t at) {
  if (prev_ < 0) { fail(at, "ring-closure digit must follow an atom"); return; }
  auto it = rings_.find(label);
  if (it == rings_.end()) {
    auto& refs = q_->atoms[static_cast<std::size_t>(prev_)].refs;
    refs.push_back(-2);  // slot for the partner
    rings_[label] = {prev_, has_pending_ ? pending_bond_ : -1, refs.size() - 1, at};
  } else {
    const Open o = it->second;
    rings_.erase(it);
    if (o.atom == prev_) { fail(at, "ring closure bonds an atom to itself"); return; }
    q_->atoms[static_cast<std::size_t>(o.atom)].refs[o.slot] = prev_;
    q_->atoms[static_cast<std::size_t>(prev_)].refs.push_back(o.atom);
    QBond b;
    b.a = o.atom;
    b.b = prev_;
    b.expr = has_pending_ ? pending_bond_ : (o.bond >= 0 ? o.bond : bprim(BP::Default));
    q_->bonds.push_back(b);
  }
  has_pending_ = false;
  pending_bond_ = -1;
}

int Parser::bond_lo(std::size_t& i) {
  int e = bond_or(i);
  while (ok_ && i < s_.size() && s_[i] == ';') {
    ++i;
    e = binary(Op::AndLo, e, bond_or(i));
  }
  return e;
}
int Parser::bond_or(std::size_t& i) {
  int e = bond_hi(i);
  while (ok_ && i < s_.size() && s_[i] == ',') {
    ++i;
    e = binary(Op::Or, e, bond_hi(i));
  }
  return e;
}
int Parser::bond_hi(std::size_t& i) {
  int e = bond_not(i);
  while (ok_ && i < s_.size() && bond_char(s_[i]) && s_[i] != ',' && s_[i] != ';') {
    if (s_[i] == '&') ++i;
    e = binary(Op::AndHi, e, bond_not(i));
  }
  return e;
}
int Parser::bond_not(std::size_t& i) {
  if (i < s_.size() && s_[i] == '!') {
    ++i;
    Node n;
    n.op = Op::Not;
    n.a = bond_not(i);
    return node(n);
  }
  if (i >= s_.size()) { fail(i, "bond expression ends without a bond"); return -1; }
  const char c = s_[i++];
  switch (c) {
    case '-': return bprim(BP::Single);
    case '/': case '\\': return bprim(BP::Single);
    case '=': return bprim(BP::Double);
    case '#': return bprim(BP::Triple);
    case ':': return bprim(BP::Aromatic);
    case '~': return bprim(BP::Any);
    case '@': return bprim(BP::Ring);
    default: fail(i - 1, std::string("'") + c + "' is not a bond primitive"); return -1;
  }
}

int Parser::atom_lo(std::size_t& i, std::size_t end, QAtom& qa, bool& first) {
  int e = atom_or(i, end, qa, first);
  while (ok_ && i < end && s_[i] == ';') {
    ++i;
    e = binary(Op::AndLo, e, atom_or(i, end, qa, first));
  }
  return e;
}
int Parser::atom_or(std::size_t& i, std::size_t end, QAtom& qa, bool& first) {
  int e = atom_hi(i, end, qa, first);
  while (ok_ && i < end && s_[i] == ',') {
    ++i;
    e = binary(Op::Or, e, atom_hi(i, end, qa, first));
  }
  return e;
}
int Parser::atom_hi(std::size_t& i, std::size_t end, QAtom& qa, bool& first) {
  int e = atom_not(i, end, qa, first);
  while (ok_ && i < end && s_[i] != ',' && s_[i] != ';') {
    if (s_[i] == '&') ++i;
    e = binary(Op::AndHi, e, atom_not(i, end, qa, first));
  }
  return e;
}
int Parser::atom_not(std::size_t& i, std::size_t end, QAtom& qa, bool& first) {
  if (i < end && s_[i] == '!') {
    ++i;
    Node n;
    n.op = Op::Not;
    n.a = atom_not(i, end, qa, first);
    return node(n);
  }
  return atom_prim(i, end, qa, first);
}

int Parser::atom_prim(std::size_t& i, std::size_t end, QAtom& qa, bool& first) {
  if (i >= end) { fail(i, "empty atom expression"); return -1; }
  const bool was_first = first;
  first = false;
  const char c = s_[i];
  auto upper = [](char x) { return x >= 'A' && x <= 'Z'; };
  auto lower = [](char x) { return x >= 'a' && x <= 'z'; };
  if (digit(i)) {
    first = was_first;  // [2H]: an element may still follow the mass
    return prim(AP::Isotope, number(i, 0));
  }
  switch (c) {
    case '*': ++i; return prim(AP::True);
    case '#': {
      ++i;
      if (!digit(i)) { fail(i, "'#' must be followed by an atomic number"); return -1; }
      return prim(AP::AtomicNum, number(i, 0));
    }
    case '+': case '-': {
      const int sign = c == '+' ? 1 : -1;
      ++i;
      int n = 1;
      if (digit(i)) n = number(i, 1);
      else
        while (i < end && s_[i] == c) { ++n; ++i; }
      return prim(AP::Charge, sign * n);
    }
    case '@': {
      ++i;
      int ch = 1;
      if (i < end && s_[i] == '@') { ch = 2; ++i; }
      if (i < end && s_[i] == '?') { fail(i, "'@?' (unspecified-or) chirality is not supported", DiagCode::SmartsUnsupported); return -1; }
      if (i + 1 < end && upper(s_[i]) && upper(s_[i + 1]) && (s_[i] == 'T' || s_[i] == 'A' || s_[i] == 'S' || s_[i] == 'O')) {
        fail(i, "chirality classes (@TH, @AL, @SP, @TB, @OH) are not supported", DiagCode::SmartsUnsupported);
        return -1;
      }
      qa.chiral = ch;
      return prim(AP::True);
    }
    case ':': {
      ++i;
      number(i, 0);  // atom map: not a constraint
      return prim(AP::True);
    }
    case '$': {
      if (i + 1 >= end || s_[i + 1] != '(') { fail(i, "'$' must be followed by '(' (recursive SMARTS)"); return -1; }
      std::size_t j = i + 2;
      int depth = 1;
      while (j < end && depth > 0) {
        if (s_[j] == '(') ++depth;
        if (s_[j] == ')') --depth;
        if (depth > 0) ++j;
      }
      if (depth != 0) { fail(i, "recursive SMARTS '$(' is not closed"); return -1; }
      Parser p(s_.substr(i + 2, j - i - 2), d_, base_ + i + 2);
      auto data = std::make_shared<SmartsQueryData>();
      if (!p.run(*data)) { ok_ = false; return -1; }
      if (data->atoms.empty()) { fail(i, "empty recursive SMARTS"); return -1; }
      q_->recursive.push_back(std::move(data));
      const int idx = static_cast<int>(q_->recursive.size()) - 1;
      i = j + 1;
      return prim(AP::Recursive, idx);
    }
    default: break;
  }
  // Hydrogen as an element: first primitive, e.g. [H], [H+], [2H] (isotope already consumed).
  if (c == 'H' && was_first && !(i + 1 < end && lower(s_[i + 1]))) {
    ++i;
    return prim(AP::AtomicNum, 1);
  }
  if (upper(c)) {
    if (i + 1 < end && lower(s_[i + 1])) {
      const std::string two{c, s_[i + 1]};
      if (const ElementInfo* e = find_element_by_symbol(two)) {
        i += 2;
        return prim(AP::Element, e->z, 0);
      }
    }
    switch (c) {
      case 'A': ++i; return prim(AP::Aliphatic);
      case 'D': ++i; return prim(AP::Degree, number(i, 1));
      case 'H': ++i; bracket_h_ = number(i, 1); return prim(AP::TotalH, bracket_h_);
      case 'R': ++i; return prim(AP::RingCount, number(i, -1));
      case 'X': ++i; return prim(AP::Connectivity, number(i, 1));
      default: break;
    }
    if (const ElementInfo* e = find_element_by_symbol(std::string(1, c))) {
      ++i;
      return prim(AP::Element, e->z, 0);
    }
    fail(i, std::string("unknown atom primitive '") + c + "'");
    return -1;
  }
  if (lower(c)) {
    if (i + 1 < end && lower(s_[i + 1])) {
      const std::string two{c, s_[i + 1]};
      if (two == "se" || two == "as" || two == "te") {
        std::string sym = two;
        sym[0] = static_cast<char>(sym[0] - 'a' + 'A');
        i += 2;
        return prim(AP::Element, find_element_by_symbol(sym)->z, 1);
      }
    }
    switch (c) {
      case 'a': ++i; return prim(AP::Aromatic);
      case 'c': ++i; return prim(AP::Element, 6, 1);
      case 'n': ++i; return prim(AP::Element, 7, 1);
      case 'o': ++i; return prim(AP::Element, 8, 1);
      case 's': ++i; return prim(AP::Element, 16, 1);
      case 'p': ++i; return prim(AP::Element, 15, 1);
      case 'b': ++i; return prim(AP::Element, 5, 1);
      case 'h': ++i; return prim(AP::ImplicitH, number(i, -1));
      case 'r': ++i; return prim(AP::RingSize, number(i, -1));
      case 'v': ++i; return prim(AP::Valence, number(i, 1));
      case 'x': ++i; return prim(AP::RingConn, number(i, -1));
      default: break;
    }
  }
  fail(i, std::string("unknown atom primitive '") + c + "'");
  return -1;
}

}  // namespace


namespace {

bool Parser::run(SmartsQueryData& q) {
  q_ = &q;
  const std::size_t n = s_.size();
  if (n == 0) return fail(0, "empty SMARTS");
  std::vector<int> branch;
  std::size_t i = 0;
  auto upper_organic = [&](std::size_t k, int& z, std::size_t& len) {
    if (s_.substr(k, 2) == "Cl") { z = 17; len = 2; return true; }
    if (s_.substr(k, 2) == "Br") { z = 35; len = 2; return true; }
    len = 1;
    switch (s_[k]) {
      case 'B': z = 5; return true;
      case 'C': z = 6; return true;
      case 'N': z = 7; return true;
      case 'O': z = 8; return true;
      case 'P': z = 15; return true;
      case 'S': z = 16; return true;
      case 'F': z = 9; return true;
      case 'I': z = 53; return true;
      default: return false;
    }
  };
  while (i < n && ok_) {
    const char c = s_[i];
    if (c == '(') {
      if (prev_ < 0) return fail(i, "component grouping '(...)' is not supported", DiagCode::SmartsUnsupported);
      if (has_pending_) return fail(i, "bond before '(' has no atom");
      branch.push_back(prev_);
      ++i;
      continue;
    }
    if (c == ')') {
      if (branch.empty()) return fail(i, "')' has no matching '('");
      if (has_pending_) return fail(i, "bond at the end of a branch has no atom");
      prev_ = branch.back();
      branch.pop_back();
      ++i;
      continue;
    }
    if (c == '.') {
      if (has_pending_) return fail(i, "bond before '.' has no atom");
      prev_ = -1;
      ++i;
      continue;
    }
    if (c == '>') return fail(i, "reaction SMARTS is not supported", DiagCode::SmartsUnsupported);
    if (digit(i)) {
      ring_closure(s_[i] - '0', i);
      ++i;
      continue;
    }
    if (c == '%') {
      if (!digit(i + 1) || !digit(i + 2)) return fail(i, "'%' must be followed by two digits");
      ring_closure((s_[i + 1] - '0') * 10 + (s_[i + 2] - '0'), i);
      i += 3;
      continue;
    }
    if (bond_char(c)) {
      if (prev_ < 0) return fail(i, "bond has no atom before it");
      if (has_pending_) return fail(i, "two bonds in a row");
      pending_bond_ = bond_lo(i);
      has_pending_ = true;
      continue;
    }
    if (c == '[') {
      std::size_t end = i + 1;
      int depth = 0;
      while (end < n && (depth > 0 || s_[end] != ']')) {
        if (s_[end] == '[' || s_[end] == '(') ++depth;
        if (s_[end] == ']' || s_[end] == ')') --depth;
        ++end;
      }
      if (end >= n) return fail(i, "'[' is not closed");
      QAtom tmp;
      bool first = true;
      bracket_h_ = -1;
      std::size_t k = i + 1;
      const int e = atom_lo(k, end, tmp, first);
      if (!ok_) return false;
      if (k != end) return fail(k, "unexpected character in atom expression");
      add_atom(e, i, tmp.chiral, bracket_h_ == 1);
      i = end + 1;
      continue;
    }
    int z = 0;
    std::size_t len = 1;
    if (c == '*') add_atom(prim(AP::True), i, 0, false);
    else if (c == 'A') add_atom(prim(AP::Aliphatic), i, 0, false);
    else if (c == 'a') add_atom(prim(AP::Aromatic), i, 0, false);
    else if (upper_organic(i, z, len)) add_atom(prim(AP::Element, z, 0), i, 0, false);
    else if (c == 'c') add_atom(prim(AP::Element, 6, 1), i, 0, false);
    else if (c == 'n') add_atom(prim(AP::Element, 7, 1), i, 0, false);
    else if (c == 'o') add_atom(prim(AP::Element, 8, 1), i, 0, false);
    else if (c == 's') add_atom(prim(AP::Element, 16, 1), i, 0, false);
    else if (c == 'p') add_atom(prim(AP::Element, 15, 1), i, 0, false);
    else if (c == 'b') add_atom(prim(AP::Element, 5, 1), i, 0, false);
    else return fail(i, std::string("unexpected character '") + c + "' (atoms outside the organic subset need brackets)");
    i += len;
  }
  if (!ok_) return false;
  if (!branch.empty()) return fail(n, "'(' is not closed");
  if (has_pending_) return fail(n, "SMARTS ends with a bond");
  if (!rings_.empty()) return fail(rings_.begin()->second.pos, "ring closure " + std::to_string(rings_.begin()->first) + " is never closed");
  return true;
}

// ---- target ------------------------------------------------------------------------------

struct Target {
  const ChemicalSystemState& st;
  std::uint32_t sp = 0, base = 0, n = 0;
  std::vector<char> hidden;  // plain H folded into its neighbour
  std::vector<int> z, arom, charge, iso, tot_h, imp_h, deg, conn, val, ring_count, min_ring, ring_conn;
  std::vector<std::vector<std::pair<int, std::uint32_t>>> adj;  // local neighbour, global bond

  Target(const ChemicalSystemState& s, std::uint32_t species) : st(s), sp(species) {
    const Species& S = st.species[sp];
    base = S.atom_begin;
    n = S.atom_count;
    hidden.assign(n, 0);
    for (std::uint32_t i = 0; i < n; ++i) {
      const std::uint32_t g = base + i;
      if (st.atoms.element[g] != 1 || st.atoms.formal_charge[g] || st.atoms.isotope[g] || st.atoms.atom_map[g] ||
          st.atoms.implicit_h[g])
        continue;
      const auto nb = st.adjacency.neighbors(g);
      if (nb.size() == 1 && st.atoms.element[nb[0]] != 1 && st.bonds.order[st.adjacency.bonds_of(g)[0]] == 1) hidden[i] = 1;
    }
    for (auto* v : {&z, &arom, &charge, &iso, &tot_h, &imp_h, &deg, &conn, &val, &ring_count, &min_ring, &ring_conn})
      v->assign(n, 0);
    adj.assign(n, {});
    for (std::uint32_t i = 0; i < n; ++i) {
      const std::uint32_t g = base + i;
      z[i] = st.atoms.element[g];
      arom[i] = (st.atoms.flags[g] & kAtomAromatic) ? 1 : 0;
      charge[i] = st.atoms.formal_charge[g];
      iso[i] = st.atoms.isotope[g];
      const int implicit = st.atoms.implicit_h[g];
      imp_h[i] = implicit;  // Daylight: hydrogens not written as atoms (bracket H included); folded [H] added below
      int hid = 0, explicit_h = 0, valence = implicit;
      const auto nb = st.adjacency.neighbors(g);
      const auto bs = st.adjacency.bonds_of(g);
      for (std::size_t k = 0; k < nb.size(); ++k) {
        valence += st.bonds.order[bs[k]];
        const std::uint32_t l = nb[k] - base;
        if (st.bonds.flags[bs[k]] & kBondInRing) ++ring_conn[i];
        if (hidden[l]) { ++hid; continue; }
        if (st.atoms.element[nb[k]] == 1) ++explicit_h;
        adj[i].push_back({static_cast<int>(l), bs[k]});
      }
      deg[i] = static_cast<int>(adj[i].size());
      tot_h[i] = implicit + hid + explicit_h;
      imp_h[i] += hid;
      conn[i] = deg[i] + implicit + hid;
      val[i] = valence;
    }
    for (const Ring& r : find_relevant_rings(st, sp)) {
      for (std::uint32_t g : r.atoms) {
        const std::uint32_t l = g - base;
        ++ring_count[l];
        const int size = static_cast<int>(r.atoms.size());
        if (min_ring[l] == 0 || size < min_ring[l]) min_ring[l] = size;
      }
    }
  }
};

// ---- matcher ------------------------------------------------------------------------------

class Matcher {
 public:
  Matcher(const Target& t, const SmartsQueryData& q, bool chirality)
      : t_(t), q_(q), chirality_(chirality), cache_(q.atoms.size() * t.n, 0) {
    plan();
  }

  // Calls `found` for each complete mapping (local target atoms); stops when it returns false.
  template <class F>
  void search(F&& found, int fixed_first = -1) {
    map_.assign(q_.atoms.size(), -1);
    used_.assign(t_.n, 0);
    stop_ = false;
    fixed_ = fixed_first;
    extend(0, found);
  }

 private:
  struct Step {
    int atom;
    int parent;      // earlier query atom bonded to it, or -1 (new fragment)
    int parent_bond; // query bond to the parent
    std::vector<std::pair<int, int>> closures;  // (earlier query atom, query bond)
  };

  void plan() {
    const std::size_t qn = q_.atoms.size();
    std::vector<std::vector<std::pair<int, int>>> qadj(qn);
    for (std::size_t b = 0; b < q_.bonds.size(); ++b) {
      qadj[static_cast<std::size_t>(q_.bonds[b].a)].push_back({q_.bonds[b].b, static_cast<int>(b)});
      qadj[static_cast<std::size_t>(q_.bonds[b].b)].push_back({q_.bonds[b].a, static_cast<int>(b)});
    }
    std::vector<int> pos(qn, -1);
    for (std::size_t root = 0; root < qn; ++root) {
      if (pos[root] >= 0) continue;
      std::vector<std::pair<int, int>> stack{{static_cast<int>(root), -1}};
      while (!stack.empty()) {
        auto [a, via] = stack.back();
        stack.pop_back();
        if (pos[static_cast<std::size_t>(a)] >= 0) continue;
        Step s;
        s.atom = a;
        s.parent = via < 0 ? -1 : (q_.bonds[static_cast<std::size_t>(via)].a == a ? q_.bonds[static_cast<std::size_t>(via)].b
                                                                                   : q_.bonds[static_cast<std::size_t>(via)].a);
        s.parent_bond = via;
        for (auto [nb, b] : qadj[static_cast<std::size_t>(a)])
          if (pos[static_cast<std::size_t>(nb)] >= 0 && b != via) s.closures.push_back({nb, b});
        pos[static_cast<std::size_t>(a)] = static_cast<int>(steps_.size());
        steps_.push_back(std::move(s));
        auto nbs = qadj[static_cast<std::size_t>(a)];
        for (auto it = nbs.rbegin(); it != nbs.rend(); ++it)
          if (pos[static_cast<std::size_t>(it->first)] < 0) stack.push_back({it->first, it->second});
      }
    }
  }

  bool eval_bond(int node, std::uint32_t tb) const {
    const Node& nd = q_.nodes[static_cast<std::size_t>(node)];
    switch (nd.op) {
      case Op::Not: return !eval_bond(nd.a, tb);
      case Op::AndHi: case Op::AndLo: return eval_bond(nd.a, tb) && eval_bond(nd.b, tb);
      case Op::Or: return eval_bond(nd.a, tb) || eval_bond(nd.b, tb);
      case Op::Prim: break;
    }
    const bool ar = (t_.st.bonds.flags[tb] & kBondAromatic) != 0;
    const int o = t_.st.bonds.order[tb];
    switch (static_cast<BP>(nd.prim)) {
      case BP::Default: return ar || o == 1;
      case BP::Single: return !ar && o == 1;
      case BP::Double: return !ar && o == 2;
      case BP::Triple: return o == 3;
      case BP::Aromatic: return ar;
      case BP::Any: return true;
      case BP::Ring: return (t_.st.bonds.flags[tb] & kBondInRing) != 0;
    }
    return false;
  }

  bool eval_atom(int node, int ta) {
    const Node& nd = q_.nodes[static_cast<std::size_t>(node)];
    switch (nd.op) {
      case Op::Not: return !eval_atom(nd.a, ta);
      case Op::AndHi: case Op::AndLo: return eval_atom(nd.a, ta) && eval_atom(nd.b, ta);
      case Op::Or: return eval_atom(nd.a, ta) || eval_atom(nd.b, ta);
      case Op::Prim: break;
    }
    const auto i = static_cast<std::size_t>(ta);
    const int v = nd.value;
    switch (static_cast<AP>(nd.prim)) {
      case AP::True: return true;
      case AP::Element: return t_.z[i] == v && (nd.aromatic < 0 || t_.arom[i] == nd.aromatic);
      case AP::Aromatic: return t_.arom[i] == 1;
      case AP::Aliphatic: return t_.arom[i] == 0;
      case AP::AtomicNum: return t_.z[i] == v;
      case AP::Degree: return t_.deg[i] == v;
      case AP::TotalH: return t_.tot_h[i] == v;
      case AP::ImplicitH: return v < 0 ? t_.imp_h[i] > 0 : t_.imp_h[i] == v;
      case AP::RingCount: return v < 0 ? t_.ring_count[i] > 0 : t_.ring_count[i] == v;
      case AP::RingSize: return v < 0 ? t_.min_ring[i] > 0 : t_.min_ring[i] == v;
      case AP::Valence: return t_.val[i] == v;
      case AP::Connectivity: return t_.conn[i] == v;
      case AP::RingConn: return v < 0 ? t_.ring_conn[i] > 0 : t_.ring_conn[i] == v;
      case AP::Charge: return t_.charge[i] == v;
      case AP::Isotope: return t_.iso[i] == v;
      case AP::Recursive: return recursive(v, ta);
    }
    return false;
  }

  bool recursive(int idx, int ta) {
    if (rec_.size() < q_.recursive.size()) {
      rec_.resize(q_.recursive.size());
      rec_cache_.resize(q_.recursive.size());
    }
    auto& cache = rec_cache_[static_cast<std::size_t>(idx)];
    if (cache.empty()) cache.assign(t_.n, 0);
    char& c = cache[static_cast<std::size_t>(ta)];
    if (c) return c == 1;
    auto& m = rec_[static_cast<std::size_t>(idx)];
    if (!m) m = std::make_unique<Matcher>(t_, *q_.recursive[static_cast<std::size_t>(idx)], false);
    bool any = false;
    m->search([&](const std::vector<int>&) { any = true; return false; }, ta);
    c = any ? 1 : 2;
    return any;
  }

  bool atom_ok(int qa, int ta) {
    if (t_.hidden[static_cast<std::size_t>(ta)]) return false;
    char& c = cache_[static_cast<std::size_t>(qa) * t_.n + static_cast<std::size_t>(ta)];
    if (!c) c = eval_atom(q_.atoms[static_cast<std::size_t>(qa)].expr, ta) ? 1 : 2;
    return c == 1;
  }

  std::int64_t target_bond(int a, int b) const {
    for (auto [nb, bond] : t_.adj[static_cast<std::size_t>(a)])
      if (nb == b) return bond;
    return -1;
  }

  // Parity check for '@'/'@@' query atoms against the target's stored stereo.
  bool chirality_ok() const {
    constexpr long long kH = -1, kLP = -3;
    for (std::size_t qa = 0; qa < q_.atoms.size(); ++qa) {
      const QAtom& a = q_.atoms[qa];
      if (!a.chiral) continue;
      const std::uint32_t g = t_.base + static_cast<std::uint32_t>(map_[qa]);
      const TetrahedralStereo* ts = nullptr;
      for (const auto& t : t_.st.tetrahedral)
        if (t.atom == g) { ts = &t; break; }
      if (!ts || ts->parity == Chirality::None) return false;
      std::vector<long long> T;
      for (std::uint32_t r : ts->refs) {
        if (r == kImplicitHRef) T.push_back(kH);
        else if (r == kLonePairRef) T.push_back(kLP);
        else if (t_.hidden[r - t_.base]) T.push_back(kH);
        else T.push_back(static_cast<long long>(r));
      }
      std::vector<long long> Q;
      for (int r : a.refs) Q.push_back(r < 0 ? kH : static_cast<long long>(t_.base + static_cast<std::uint32_t>(map_[static_cast<std::size_t>(r)])));
      if (Q.size() == 3) {  // the neighbour the query leaves out takes the implicit-H position
        long long missing = 0;
        int count = 0;
        for (long long x : T)
          if (std::find(Q.begin(), Q.end(), x) == Q.end()) { missing = x; ++count; }
        if (count != 1) return false;
        Q.insert(Q.begin() + (a.has_prev ? 1 : 0), missing);
      }
      if (Q.size() != 4) return false;
      std::vector<int> perm;
      for (long long x : Q) {
        const auto it = std::find(T.begin(), T.end(), x);
        if (it == T.end()) return false;
        perm.push_back(static_cast<int>(it - T.begin()));
      }
      int inversions = 0;
      for (int x = 0; x < 4; ++x)
        for (int y = x + 1; y < 4; ++y) inversions += perm[static_cast<std::size_t>(x)] > perm[static_cast<std::size_t>(y)];
      const bool target_ccw = ts->parity == Chirality::CounterClockwise;
      const bool ccw_in_query_order = (inversions % 2 == 0) ? target_ccw : !target_ccw;
      if (ccw_in_query_order != (a.chiral == 1)) return false;
    }
    return true;
  }

  template <class F>
  void extend(std::size_t k, F& found) {
    if (stop_) return;
    if (k == steps_.size()) {
      if (chirality_ && !chirality_ok()) return;
      if (!found(map_)) stop_ = true;
      return;
    }
    const Step& s = steps_[k];
    auto try_atom = [&](int ta, std::int64_t via_bond) {
      if (used_[static_cast<std::size_t>(ta)] || !atom_ok(s.atom, ta)) return;
      if (via_bond >= 0 && !eval_bond(q_.bonds[static_cast<std::size_t>(s.parent_bond)].expr, static_cast<std::uint32_t>(via_bond)))
        return;
      for (auto [qb_atom, qb] : s.closures) {
        const std::int64_t tb = target_bond(ta, map_[static_cast<std::size_t>(qb_atom)]);
        if (tb < 0 || !eval_bond(q_.bonds[static_cast<std::size_t>(qb)].expr, static_cast<std::uint32_t>(tb))) return;
      }
      map_[static_cast<std::size_t>(s.atom)] = ta;
      used_[static_cast<std::size_t>(ta)] = 1;
      extend(k + 1, found);
      used_[static_cast<std::size_t>(ta)] = 0;
      map_[static_cast<std::size_t>(s.atom)] = -1;
    };
    if (k == 0 && fixed_ >= 0) {
      try_atom(fixed_, -1);
    } else if (s.parent >= 0) {
      for (auto [nb, bond] : t_.adj[static_cast<std::size_t>(map_[static_cast<std::size_t>(s.parent)])]) {
        try_atom(nb, bond);
        if (stop_) return;
      }
    } else {
      for (std::uint32_t ta = 0; ta < t_.n && !stop_; ++ta) try_atom(static_cast<int>(ta), -1);
    }
  }

  const Target& t_;
  const SmartsQueryData& q_;
  bool chirality_;
  std::vector<char> cache_;
  std::vector<Step> steps_;
  std::vector<int> map_;
  std::vector<char> used_;
  bool stop_ = false;
  int fixed_ = -1;
  std::vector<std::unique_ptr<Matcher>> rec_;
  std::vector<std::vector<char>> rec_cache_;
};

}  // namespace

SmartsQuery parse_smarts(std::string_view smarts, Diagnostics& diag) {
  auto data = std::make_shared<SmartsQueryData>();
  Parser p(smarts, diag, 0);
  SmartsQuery q;
  if (p.run(*data)) q.data_ = std::move(data);
  return q;
}

std::vector<std::vector<std::uint32_t>> match_smarts(const ChemicalSystemState& state, std::uint32_t species,
                                                     const SmartsQuery& query, const SmartsMatchOptions& options) {
  std::vector<std::vector<std::uint32_t>> out;
  if (!query.valid() || query.atom_count() == 0) return out;
  const Target t(state, species);
  Matcher m(t, query.data(), options.use_chirality);
  std::set<std::vector<int>> seen;
  m.search([&](const std::vector<int>& map) {
    if (options.unique) {
      std::vector<int> key = map;
      std::sort(key.begin(), key.end());
      if (!seen.insert(key).second) return true;
    }
    std::vector<std::uint32_t> g;
    g.reserve(map.size());
    for (int a : map) g.push_back(t.base + static_cast<std::uint32_t>(a));
    out.push_back(std::move(g));
    return options.max_matches == 0 || out.size() < options.max_matches;
  });
  return out;
}

bool has_match(const ChemicalSystemState& state, std::uint32_t species, const SmartsQuery& query, bool use_chirality) {
  SmartsMatchOptions o;
  o.max_matches = 1;
  o.use_chirality = use_chirality;
  o.unique = false;
  return !match_smarts(state, species, query, o).empty();
}

}  // namespace omnichem
