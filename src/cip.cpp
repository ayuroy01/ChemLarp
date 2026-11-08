#include "omnichem/cip.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>

#include "omnichem/elements.hpp"

namespace omnichem {
namespace {

constexpr std::size_t kMaxNodes = 200000;
constexpr int kMaxDepth = 64;

struct Node {
  std::int64_t atom;  // global atom id; -1 implicit H; -2 lone pair / phantom
  std::uint8_t z;
  std::uint32_t mass;  // milli-daltons, for Rule 2
  bool duplicate;
  int parent;  // node index, -1 for root
};

class Digraph {
 public:
  Digraph(const ChemicalSystemState& st, std::uint32_t root, std::uint32_t exclude) : st_(st), exclude_(exclude) {
    nodes_.push_back(make_real(root, -1));
  }

  int branch(std::uint32_t ref) {
    Node n;
    if (ref == kImplicitHRef) n = {-1, 1, mass_of(1, 0), false, 0};
    else if (ref == kLonePairRef) n = {-2, 0, 0, false, 0};
    else n = make_real(ref, 0);
    return add(n);
  }

  // >0 if branch a outranks b
  int compare(int a, int b) {
    for (int pass = 1; pass <= 2; ++pass) {
      const int c = compare_branches(a, b, pass);
      if (c != 0) return c;
    }
    return 0;
  }

  bool overflow() const { return overflow_; }

 private:
  Node make_real(std::uint32_t atom, int parent) const {
    return {static_cast<std::int64_t>(atom), st_.atoms.element[atom], mass_of(st_.atoms.element[atom], st_.atoms.isotope[atom]),
            false, parent};
  }
  static std::uint32_t mass_of(std::uint8_t z, std::uint16_t isotope) {
    if (isotope) return static_cast<std::uint32_t>(isotope) * 1000u;
    const ElementInfo* e = find_element_by_z(z);
    return e ? static_cast<std::uint32_t>(std::lround(e->mass * 1000.0)) : 0u;
  }
  int add(Node n) {
    if (nodes_.size() >= kMaxNodes) overflow_ = true;
    nodes_.push_back(n);
    kids_.emplace_back();
    expanded_.push_back(0);
    return static_cast<int>(nodes_.size() - 1);
  }
  std::tuple<int, std::uint32_t> key(int n, int pass) const {
    return {nodes_[static_cast<std::size_t>(n)].z, pass >= 2 ? nodes_[static_cast<std::size_t>(n)].mass : 0u};
  }
  bool is_ancestor(int n, std::int64_t atom) const {
    for (int p = n; p >= 0; p = nodes_[static_cast<std::size_t>(p)].parent)
      if (nodes_[static_cast<std::size_t>(p)].atom == atom) return true;
    return false;
  }

  const std::vector<int>& children(int n) {
    // nodes_ / kids_ may reallocate while expanding: index, never hold references across add().
    if (static_cast<std::size_t>(n) >= kids_.size()) kids_.resize(nodes_.size());
    if (expanded_.size() < nodes_.size()) expanded_.resize(nodes_.size(), 0);
    if (expanded_[static_cast<std::size_t>(n)] || overflow_) return kids_[static_cast<std::size_t>(n)];
    expanded_[static_cast<std::size_t>(n)] = 1;
    const Node self = nodes_[static_cast<std::size_t>(n)];
    if (self.duplicate || self.atom < 0) return kids_[static_cast<std::size_t>(n)];
    const std::uint32_t a = static_cast<std::uint32_t>(self.atom);
    const std::int64_t parent_atom = self.parent >= 0 ? nodes_[static_cast<std::size_t>(self.parent)].atom : -9;
    std::vector<int> out;
    auto nb = st_.adjacency.neighbors(a);
    auto bd = st_.adjacency.bonds_of(a);
    for (std::size_t k = 0; k < nb.size(); ++k) {
      const std::uint32_t w = nb[k];
      const int order = st_.bonds.order[bd[k]];
      if (self.parent < 0 && w == exclude_) continue;  // root of an E/Z end: skip the double-bond partner
      Node dup = make_real(w, n);
      dup.duplicate = true;
      if (static_cast<std::int64_t>(w) == parent_atom) {
        for (int d = 1; d < order; ++d) out.push_back(add(dup));
        continue;
      }
      if (is_ancestor(n, w)) out.push_back(add(dup));  // ring closure
      else out.push_back(add(make_real(w, n)));
      for (int d = 1; d < order; ++d) out.push_back(add(dup));
    }
    for (int h = 0; h < st_.atoms.implicit_h[a]; ++h) out.push_back(add({-1, 1, mass_of(1, 0), false, n}));
    kids_[static_cast<std::size_t>(n)] = std::move(out);
    return kids_[static_cast<std::size_t>(n)];
  }

  std::vector<int> sorted_children(int n, int pass) {
    auto memo = sorted_.find({n, pass});
    if (memo != sorted_.end()) return memo->second;
    std::vector<int> c = children(n);
    std::stable_sort(c.begin(), c.end(), [&](int x, int y) {
      if (key(x, pass) != key(y, pass)) return key(x, pass) > key(y, pass);
      return compare_branches(x, y, pass) > 0;
    });
    sorted_[{n, pass}] = c;
    return c;
  }

  // Breadth-first, sphere by sphere: the sets of each sphere are ordered by the
  // ranking of their parents, atoms within a set by their own (recursive) rank.
  int compare_branches(int a, int b, int pass) {
    if (a == b) return 0;
    auto memo = cmp_.find({a, b, pass});
    if (memo != cmp_.end()) return memo->second;
    int result = 0;
    if (key(a, pass) != key(b, pass)) {
      result = key(a, pass) > key(b, pass) ? 1 : -1;
    } else {
      std::vector<std::vector<int>> la{{a}}, lb{{b}};
      for (int depth = 0; depth < kMaxDepth && result == 0 && !overflow_; ++depth) {
        std::vector<std::vector<int>> na, nbv;
        for (const auto& set : la)
          for (int x : set) na.push_back(sorted_children(x, pass));
        for (const auto& set : lb)
          for (int x : set) nbv.push_back(sorted_children(x, pass));
        bool any = false;
        const std::size_t m = std::max(na.size(), nbv.size());
        static const std::vector<int> empty;
        for (std::size_t i = 0; i < m && result == 0; ++i) {
          const auto& sa = i < na.size() ? na[i] : empty;
          const auto& sb = i < nbv.size() ? nbv[i] : empty;
          any = any || !sa.empty() || !sb.empty();
          const std::size_t len = std::max(sa.size(), sb.size());
          for (std::size_t j = 0; j < len; ++j) {
            const auto ka = j < sa.size() ? key(sa[j], pass) : std::tuple<int, std::uint32_t>{0, 0};
            const auto kb = j < sb.size() ? key(sb[j], pass) : std::tuple<int, std::uint32_t>{0, 0};
            if (ka != kb) { result = ka > kb ? 1 : -1; break; }
          }
        }
        if (!any) break;
        la = std::move(na);
        lb = std::move(nbv);
      }
    }
    cmp_[{a, b, pass}] = result;
    cmp_[{b, a, pass}] = -result;
    return result;
  }

  const ChemicalSystemState& st_;
  std::uint32_t exclude_;
  std::vector<Node> nodes_;
  std::vector<std::vector<int>> kids_{1};
  std::vector<char> expanded_{0};
  std::map<std::pair<int, int>, std::vector<int>> sorted_;
  std::map<std::tuple<int, int, int>, int> cmp_;
  bool overflow_ = false;
};

int permutation_parity(const std::array<std::size_t, 4>& p) {
  int inv = 0;
  for (std::size_t i = 0; i < 4; ++i)
    for (std::size_t j = i + 1; j < 4; ++j) inv += p[i] > p[j];
  return inv & 1;
}

// Ranks substituent refs (highest first). Returns false on a tie.
bool rank_refs(const ChemicalSystemState& st, std::uint32_t root, std::uint32_t exclude,
               const std::vector<std::uint32_t>& refs, std::vector<std::size_t>& order, bool& overflow) {
  Digraph g(st, root, exclude);
  std::vector<int> nodes;
  for (std::uint32_t r : refs) nodes.push_back(g.branch(r));
  order.resize(refs.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  bool tie = false;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
    return g.compare(nodes[x], nodes[y]) > 0;
  });
  for (std::size_t i = 0; i + 1 < order.size(); ++i)
    if (g.compare(nodes[order[i]], nodes[order[i + 1]]) == 0) tie = true;
  overflow = g.overflow();
  return !tie && !overflow;
}

}  // namespace

const char* to_string(CipLabel l) noexcept {
  switch (l) {
    case CipLabel::R: return "R";
    case CipLabel::S: return "S";
    case CipLabel::E: return "E";
    case CipLabel::Z: return "Z";
    default: return "";
  }
}

int cip_compare(const ChemicalSystemState& st, std::uint32_t root, std::uint32_t a, std::uint32_t b,
                std::uint32_t exclude) {
  Digraph g(st, root, exclude);
  const int na = g.branch(a), nb = g.branch(b);
  return g.compare(na, nb);
}

CipLabel cip_label(const ChemicalSystemState& st, const TetrahedralStereo& c, Diagnostics* diag) {
  if (c.parity == Chirality::None) return CipLabel::None;
  std::vector<std::uint32_t> refs(c.refs.begin(), c.refs.end());
  std::vector<std::size_t> order;
  bool overflow = false;
  if (!rank_refs(st, c.atom, UINT32_MAX, refs, order, overflow)) {
    if (diag)
      diag->warn(DiagCode::CipUndetermined,
                 overflow ? "CIP digraph too large; descriptor not assigned"
                          : "two substituents tie under CIP rules 1-2 (centre may not be stereogenic, or needs rules 3-5)")
          .on_atom(c.atom);
    return CipLabel::None;
  }
  // Reorder refs as [lowest, highest, 2nd, 3rd]. '@' in that order means
  // highest->2nd->3rd is counter-clockwise viewed from the lowest, i.e.
  // clockwise with the lowest pointing away: R.
  const std::array<std::size_t, 4> target{order[3], order[0], order[1], order[2]};
  const bool flip = permutation_parity(target) != 0;
  const bool ccw = (c.parity == Chirality::CounterClockwise) != flip;
  return ccw ? CipLabel::R : CipLabel::S;
}

CipLabel cip_label(const ChemicalSystemState& st, const DoubleBondStereo& d, Diagnostics* diag) {
  const std::uint32_t u = st.bonds.begin[d.bond], v = st.bonds.end[d.bond];
  auto highest = [&](std::uint32_t end, std::uint32_t other, std::uint32_t& best) -> bool {
    std::vector<std::uint32_t> refs;
    for (std::uint32_t w : st.adjacency.neighbors(end))
      if (w != other) refs.push_back(w);
    for (int h = 0; h < st.atoms.implicit_h[end]; ++h) refs.push_back(kImplicitHRef);
    if (refs.empty()) return false;
    if (refs.size() == 1) { best = refs[0]; return true; }
    std::vector<std::size_t> order;
    bool overflow = false;
    if (!rank_refs(st, end, other, refs, order, overflow)) return false;
    best = refs[order[0]];
    return true;
  };
  std::uint32_t hu = 0, hv = 0;
  if (!highest(u, v, hu) || !highest(v, u, hv)) {
    if (diag)
      diag->warn(DiagCode::CipUndetermined, "double bond has identical substituents on one end; not stereogenic")
          .on_atom(u);
    return CipLabel::None;
  }
  const bool cis = d.cis != (hu != d.ref_begin) != (hv != d.ref_end);
  return cis ? CipLabel::Z : CipLabel::E;
}

std::vector<CipAssignment> cip_centres(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics* diag) {
  std::vector<CipAssignment> out;
  for (const auto& t : st.tetrahedral)
    if (st.atoms.species[t.atom] == sp) out.push_back({t.atom, cip_label(st, t, diag)});
  return out;
}

std::vector<CipAssignment> cip_double_bonds(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics* diag) {
  std::vector<CipAssignment> out;
  for (const auto& d : st.double_bonds)
    if (st.bonds.species[d.bond] == sp) out.push_back({d.bond, cip_label(st, d, diag)});
  return out;
}

}  // namespace omnichem
