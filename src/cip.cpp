#include "chemlarp/cip.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <tuple>
#include <utility>

#include "chemlarp/elements.hpp"

namespace chemlarp {
namespace {

constexpr std::size_t kMaxNodes = 200000;
constexpr int kMaxSpheres = 1000;

// Sequence rules in the order they are applied.
enum Rule : int { k1a = 0, k1b, k2, k3, k4a, k4b, k4c, k5 };

struct Node {
  std::int64_t atom;   // global atom id; -1 implicit H; -2 lone pair / phantom
  std::uint32_t z6;    // 6 x atomic number (duplicates in mancude rings carry a fractional value)
  std::uint32_t mass;  // milli-daltons, for Rule 2
  std::uint8_t dup;    // 0 real, 1 multiple-bond duplicate, 2 ring-closure duplicate
  int parent;          // node index, -1 for the root
  int depth;           // sphere in the digraph as built
  int target;          // ring-closure duplicate: the node it duplicates
  CipLabel aux;        // auxiliary descriptor, set by compute_aux()
};

// A branch: the subtree at `to` seen from its neighbour `from`. The digraph is
// a tree, so it can be read from any node without rebuilding it.
struct Edge {
  int from;
  int to;
};

class Digraph {
 public:
  Digraph(const ChemicalSystemState& st, std::uint32_t root) : st_(st) {
    add(real(root, -1, 0));
    tetra_of_.assign(st.atoms.element.size(), -1);
    for (std::size_t i = 0; i < st.tetrahedral.size(); ++i) tetra_of_[st.tetrahedral[i].atom] = static_cast<int>(i);
    for (std::size_t i = 0; i < st.double_bonds.size(); ++i) db_of_[st.double_bonds[i].bond] = static_cast<int>(i);
    // A tetrahedral root can only hold a multiple bond through an expanded
    // octet (>S=O, >P=O). Those bonds are ranked as single: no duplicates.
    root_is_centre_ = tetra_of_[root] >= 0;
  }

  bool overflow() const { return overflow_; }

  // Descriptor of the tetrahedral centre at node n using rules up to `level`.
  CipLabel label_tetra(int n, int level) {
    const int idx = tetra_of_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(n)].atom)];
    if (idx < 0) return CipLabel::None;
    const TetrahedralStereo& t = st_.tetrahedral[static_cast<std::size_t>(idx)];
    if (t.parity == Chirality::None) return CipLabel::None;
    set_root(n);
    std::array<int, 4> nd{};
    for (std::size_t i = 0; i < 4; ++i) nd[i] = find_ligand(n, t.refs[i]);
    auto rel = [&](std::size_t x, std::size_t y, int lv) {
      if (nd[x] < 0 || nd[y] < 0) return (nd[x] >= 0) - (nd[y] >= 0);  // a lone pair ranks lowest
      return cmp({n, nd[x]}, {n, nd[y]}, lv);
    };
    std::array<std::size_t, 4> order{0, 1, 2, 3};
    std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return rel(x, y, level) > 0; });
    bool pseudo = false;
    if (level >= k5 && rel(order[0], order[1], level) == 0 && rel(order[1], order[2], level) > 0 &&
        rel(order[2], order[3], level) == 0) {
      // Two equivalent pairs, a = b > a' = b' (a spiro atom joining two alike
      // rings): take a, b, then the ring partner of a, then that of b. The
      // result does not depend on which of the equals is called a.
      if (!same_ring(n, nd[order[0]], nd[order[2]])) std::swap(order[2], order[3]);
      if (!same_ring(n, nd[order[0]], nd[order[2]]) || !same_ring(n, nd[order[1]], nd[order[3]]) ||
          same_ring(n, nd[order[0]], nd[order[1]]))
        return CipLabel::None;
      return descriptor(t, order, false);
    }
    for (std::size_t i = 0; i + 1 < 4; ++i) {
      if (rel(order[i], order[i + 1], level) == 0) return CipLabel::None;
      if (level >= k5 && rel(order[i], order[i + 1], k4c) == 0) pseudo = true;  // only Rule 5 separates them
    }
    if (overflow_) return CipLabel::None;
    return descriptor(t, order, pseudo);
  }

  // Descriptor of the double bond between node p and its child c.
  CipLabel label_db(int p, int c, int level) {
    const std::int64_t bond = bond_between(nodes_[static_cast<std::size_t>(p)].atom, nodes_[static_cast<std::size_t>(c)].atom);
    if (bond < 0) return CipLabel::None;
    const auto it = db_of_.find(static_cast<std::uint32_t>(bond));
    if (it == db_of_.end()) return CipLabel::None;
    const DoubleBondStereo& d = st_.double_bonds[static_cast<std::size_t>(it->second)];
    std::uint32_t top_p = 0, top_c = 0;
    if (!top_ligand(p, c, level, top_p) || !top_ligand(c, p, level, top_c) || overflow_) return CipLabel::None;
    const bool p_is_begin = static_cast<std::int64_t>(st_.bonds.begin[d.bond]) == nodes_[static_cast<std::size_t>(p)].atom;
    const std::uint32_t hu = p_is_begin ? top_p : top_c, hv = p_is_begin ? top_c : top_p;
    const bool cis = d.cis != (hu != d.ref_begin) != (hv != d.ref_end);
    return cis ? CipLabel::Z : CipLabel::E;
  }

  // The child of the root that stands for neighbour atom / pseudo-ref `ref`.
  int root_ligand(std::uint32_t ref) { return find_ligand(0, ref); }

  int compare_root_ligands(int a, int b, int level) {
    set_root(0);
    return cmp({0, a}, {0, b}, level);
  }

  // Auxiliary descriptors for every stereogenic unit in the digraph, working
  // from the outermost sphere towards the root (P-92.6): each unit is ranked
  // from its own node, using the descriptors already found farther out.
  void compute_aux() {
    for (std::size_t i = 0; i < nodes_.size() && !overflow_; ++i) kids(static_cast<int>(i));
    if (overflow_) return;
    int deepest = 0;
    for (const Node& n : nodes_) deepest = std::max(deepest, n.depth);
    std::vector<std::vector<int>> by_depth(static_cast<std::size_t>(deepest) + 1);
    for (std::size_t i = 0; i < nodes_.size(); ++i)
      if (!nodes_[i].dup && nodes_[i].atom >= 0) by_depth[static_cast<std::size_t>(nodes_[i].depth)].push_back(static_cast<int>(i));
    for (int depth = deepest; depth >= 1; --depth) {
      std::vector<std::pair<int, CipLabel>> found;
      for (int n : by_depth[static_cast<std::size_t>(depth)]) {
        const Node self = nodes_[static_cast<std::size_t>(n)];
        if (tetra_of_[static_cast<std::size_t>(self.atom)] >= 0) {
          const CipLabel l = label_tetra(n, k5);
          if (l != CipLabel::None) found.push_back({n, l});
        }
        if (!nodes_[static_cast<std::size_t>(self.parent)].dup) {
          const CipLabel l = label_db(self.parent, n, k5);
          if (l != CipLabel::None) {
            found.push_back({n, l});
            found.push_back({self.parent, l});
          }
        }
      }
      for (const auto& [n, l] : found) nodes_[static_cast<std::size_t>(n)].aux = l;
      if (!found.empty()) any_aux_ = true;
    }
    set_root(0);
  }

 private:
  // Reorder refs as [lowest, highest, 2nd, 3rd]. '@' in that order means
  // highest->2nd->3rd is counter-clockwise viewed from the lowest, i.e.
  // clockwise with the lowest pointing away: R.
  static CipLabel descriptor(const TetrahedralStereo& t, const std::array<std::size_t, 4>& order, bool pseudo) {
    const std::array<std::size_t, 4> p{order[3], order[0], order[1], order[2]};
    int inv = 0;
    for (std::size_t i = 0; i < 4; ++i)
      for (std::size_t j = i + 1; j < 4; ++j) inv += p[i] > p[j];
    const bool r = (t.parity == Chirality::CounterClockwise) != ((inv & 1) != 0);
    if (pseudo) return r ? CipLabel::r : CipLabel::s;
    return r ? CipLabel::R : CipLabel::S;
  }
  // True if ligand nodes x and y of node n are joined by a path that avoids n's atom.
  bool same_ring(int n, int x, int y) const {
    if (x < 0 || y < 0) return false;
    const std::int64_t c = nodes_[static_cast<std::size_t>(n)].atom, ax = nodes_[static_cast<std::size_t>(x)].atom,
                       ay = nodes_[static_cast<std::size_t>(y)].atom;
    if (ax < 0 || ay < 0) return false;
    std::vector<std::uint32_t> stack{static_cast<std::uint32_t>(ax)}, seen{static_cast<std::uint32_t>(ax)};
    while (!stack.empty()) {
      const std::uint32_t a = stack.back();
      stack.pop_back();
      if (static_cast<std::int64_t>(a) == ay) return true;
      for (std::uint32_t w : st_.adjacency.neighbors(a)) {
        if (static_cast<std::int64_t>(w) == c || std::find(seen.begin(), seen.end(), w) != seen.end()) continue;
        seen.push_back(w);
        stack.push_back(w);
      }
    }
    return false;
  }
  Node real(std::uint32_t atom, int parent, int depth) const {
    const std::uint8_t z = st_.atoms.element[atom];
    return {static_cast<std::int64_t>(atom), 6u * z, mass_of(z, st_.atoms.isotope[atom]), 0, parent, depth, -1, CipLabel::None};
  }
  static std::uint32_t mass_of(std::uint8_t z, std::uint16_t isotope) {
    if (isotope) return static_cast<std::uint32_t>(isotope) * 1000u;
    const ElementInfo* e = find_element_by_z(z);
    return e ? static_cast<std::uint32_t>(std::lround(e->mass * 1000.0)) : 0u;
  }
  int add(const Node& n) {
    if (nodes_.size() >= kMaxNodes) overflow_ = true;
    nodes_.push_back(n);
    kids_.emplace_back();
    expanded_.push_back(0);
    return static_cast<int>(nodes_.size() - 1);
  }
  int ancestor(int n, std::int64_t atom) const {
    for (int p = n; p >= 0; p = nodes_[static_cast<std::size_t>(p)].parent)
      if (nodes_[static_cast<std::size_t>(p)].atom == atom) return p;
    return -1;
  }
  std::int64_t bond_between(std::int64_t a, std::int64_t b) const {
    if (a < 0 || b < 0) return -1;
    auto nb = st_.adjacency.neighbors(static_cast<std::uint32_t>(a));
    auto bd = st_.adjacency.bonds_of(static_cast<std::uint32_t>(a));
    for (std::size_t k = 0; k < nb.size(); ++k)
      if (nb[k] == static_cast<std::uint32_t>(b)) return bd[k];
    return -1;
  }

  // --- mancude rings: the duplicate of an aromatic atom carries the mean
  // atomic number of the atoms it is double-bonded to over all Kekulé forms.
  std::int64_t arom_mate(std::uint32_t x) const {
    auto nb = st_.adjacency.neighbors(x);
    auto bd = st_.adjacency.bonds_of(x);
    for (std::size_t k = 0; k < nb.size(); ++k)
      if (st_.bonds.order[bd[k]] == 2 && (st_.bonds.flags[bd[k]] & kBondAromatic)) return nb[k];
    return -1;
  }
  bool alternating(std::uint32_t from, std::uint32_t goal, std::vector<std::uint32_t>& seen) const {
    auto nb = st_.adjacency.neighbors(from);
    auto bd = st_.adjacency.bonds_of(from);
    for (std::size_t k = 0; k < nb.size(); ++k) {
      if (!(st_.bonds.flags[bd[k]] & kBondAromatic) || st_.bonds.order[bd[k]] != 1) continue;
      const std::uint32_t x = nb[k];
      if (x == goal) return true;
      if (std::find(seen.begin(), seen.end(), x) != seen.end()) continue;
      const std::int64_t y = arom_mate(x);
      if (y < 0 || std::find(seen.begin(), seen.end(), static_cast<std::uint32_t>(y)) != seen.end()) continue;
      seen.push_back(x);
      seen.push_back(static_cast<std::uint32_t>(y));
      if (alternating(static_cast<std::uint32_t>(y), goal, seen)) return true;
    }
    return false;
  }
  std::pair<std::uint32_t, std::uint32_t> mancude(std::uint32_t a) {  // {6 x Z, mass}
    const auto memo = mancude_.find(a);
    if (memo != mancude_.end()) return memo->second;
    const std::int64_t m0 = arom_mate(a);
    std::uint32_t sum = 0, mass = 0, count = 0;
    auto nb = st_.adjacency.neighbors(a);
    auto bd = st_.adjacency.bonds_of(a);
    for (std::size_t k = 0; k < nb.size() && m0 >= 0; ++k) {
      if (!(st_.bonds.flags[bd[k]] & kBondAromatic)) continue;
      const std::uint32_t w = nb[k];
      bool ok = static_cast<std::int64_t>(w) == m0;
      if (!ok) {
        const std::int64_t m = arom_mate(w);
        std::vector<std::uint32_t> seen{a, w, static_cast<std::uint32_t>(m0)};
        ok = m >= 0 && alternating(static_cast<std::uint32_t>(m0), static_cast<std::uint32_t>(m), seen);
      }
      if (ok) {
        sum += st_.atoms.element[w];
        mass += mass_of(st_.atoms.element[w], 0);
        ++count;
      }
    }
    const std::pair<std::uint32_t, std::uint32_t> out{count ? 6u * sum / count : 0u, count ? mass / count : 0u};
    mancude_[a] = out;
    return out;
  }

  const std::vector<int>& kids(int n) {
    const std::size_t un = static_cast<std::size_t>(n);
    if (expanded_[un] || overflow_) return kids_[un];
    expanded_[un] = 1;
    const Node self = nodes_[un];  // copy: nodes_ reallocates while expanding
    if (self.dup || self.atom < 0) return kids_[un];
    const std::uint32_t a = static_cast<std::uint32_t>(self.atom);
    const std::int64_t parent_atom = self.parent >= 0 ? nodes_[static_cast<std::size_t>(self.parent)].atom : -9;
    std::vector<int> out;
    auto nb = st_.adjacency.neighbors(a);
    auto bd = st_.adjacency.bonds_of(a);
    for (std::size_t k = 0; k < nb.size(); ++k) {
      const std::uint32_t w = nb[k];
      int order = st_.bonds.order[bd[k]];
      if (root_is_centre_ && (n == 0 || (self.parent == 0 && static_cast<std::int64_t>(w) == parent_atom))) order = 1;
      Node dup = real(w, n, self.depth + 1);
      dup.dup = 1;
      dup.mass = mass_of(st_.atoms.element[w], 0);  // duplicates carry no isotope label
      if (order > 1 && (st_.bonds.flags[bd[k]] & kBondAromatic)) std::tie(dup.z6, dup.mass) = mancude(a);
      if (static_cast<std::int64_t>(w) == parent_atom) {
        for (int d = 1; d < order; ++d) out.push_back(add(dup));
        continue;
      }
      Node sub = real(w, n, self.depth + 1);
      const int anc = ancestor(n, w);
      if (anc >= 0) {
        sub.dup = 2;
        sub.target = anc;
        sub.mass = dup.mass;
      }
      out.push_back(add(sub));
      for (int d = 1; d < order; ++d) out.push_back(add(dup));
    }
    for (int h = 0; h < st_.atoms.implicit_h[a]; ++h)
      out.push_back(add({-1, 6u, mass_of(1, 0), 0, n, self.depth + 1, -1, CipLabel::None}));
    kids_[un] = std::move(out);
    return kids_[un];
  }

  // Ligand node of n for a stereo ref (atom id, implicit H or lone pair); -1 if absent.
  int find_ligand(int n, std::uint32_t ref) {
    if (ref == kLonePairRef) return -1;
    const std::vector<int> ks = kids(n);
    for (int k : ks) {
      const Node& x = nodes_[static_cast<std::size_t>(k)];
      if (ref == kImplicitHRef ? x.atom == -1 : (x.dup != 1 && x.atom == static_cast<std::int64_t>(ref))) return k;
    }
    const int p = nodes_[static_cast<std::size_t>(n)].parent;
    if (p >= 0 && ref != kImplicitHRef && nodes_[static_cast<std::size_t>(p)].atom == static_cast<std::int64_t>(ref)) return p;
    return -1;
  }

  // Highest-ranked substituent of double-bond end `e` (partner node `q`), as a stereo ref.
  bool top_ligand(int e, int q, int level, std::uint32_t& ref) {
    set_root(e);
    const std::int64_t partner = nodes_[static_cast<std::size_t>(q)].atom;
    std::vector<int> subs;
    for (int x : nbrs({q, e})) {
      const Node& nx = nodes_[static_cast<std::size_t>(x)];
      if (nx.dup == 1 && nx.atom == partner) continue;
      subs.push_back(x);
    }
    if (subs.empty() || subs.size() > 2) return false;
    int best = subs[0];
    if (subs.size() == 2) {
      const int c = cmp({e, subs[0]}, {e, subs[1]}, level);
      if (c == 0) return false;
      if (c < 0) best = subs[1];
    }
    const std::int64_t atom = nodes_[static_cast<std::size_t>(best)].atom;
    ref = atom >= 0 ? static_cast<std::uint32_t>(atom) : kImplicitHRef;
    return true;
  }

  void set_root(int n) {
    root_ = n;
    sorted_.clear();
    bfs_.clear();
    pairs_.clear();
  }

  int tree_distance(int a, int b) const {
    int d = 0;
    while (a != b) {
      if (nodes_[static_cast<std::size_t>(a)].depth >= nodes_[static_cast<std::size_t>(b)].depth) a = nodes_[static_cast<std::size_t>(a)].parent;
      else b = nodes_[static_cast<std::size_t>(b)].parent;
      ++d;
    }
    return d;
  }

  std::vector<int> nbrs(Edge e) {
    std::vector<int> out;
    for (int x : kids(e.to))
      if (x != e.from) out.push_back(x);
    const int p = nodes_[static_cast<std::size_t>(e.to)].parent;
    if (p >= 0 && p != e.from) out.push_back(p);
    return out;
  }

  // One rule applied to one node (larger ranks higher). Rule 4b is not a
  // per-node property and is handled by like_unlike().
  std::int64_t key(int n, int rule) const {
    if (n < 0) return 0;
    const Node& x = nodes_[static_cast<std::size_t>(n)];
    switch (rule) {
      case k1a: return x.z6;
      case k1b:  // a ring-closure duplicate takes the root distance of the atom it duplicates
        return x.dup == 2 ? tree_distance(root_, n) - tree_distance(root_, x.target) : 0;
      case k2: return x.mass;
      case k3: return x.aux == CipLabel::Z ? 2 : x.aux == CipLabel::E ? 1 : 0;
      case k4a: return (x.aux == CipLabel::R || x.aux == CipLabel::S) ? 2 : (x.aux == CipLabel::r || x.aux == CipLabel::s) ? 1 : 0;
      case k4c: return x.aux == CipLabel::r ? 2 : x.aux == CipLabel::s ? 1 : 0;
      case k5: return x.aux == CipLabel::R ? 2 : x.aux == CipLabel::S ? 1 : 0;
      default: return 0;
    }
  }
  int node_cmp(int x, int y, int level) const {
    for (int r = k1a; r <= level; ++r) {
      if (r == k4b) continue;
      const std::int64_t kx = key(x, r), ky = key(y, r);
      if (kx != ky) return kx > ky ? 1 : -1;
    }
    return 0;
  }

  std::vector<int> sorted(Edge e, int level) {
    const std::array<int, 3> id{e.from, e.to, level};
    const auto memo = sorted_.find(id);
    if (memo != sorted_.end()) return memo->second;
    std::vector<int> c = nbrs(e);
    std::stable_sort(c.begin(), c.end(), [&](int x, int y) { return cmp({e.to, x}, {e.to, y}, level) > 0; });
    sorted_[id] = c;
    return c;
  }

  // Rules applied hierarchically: each rule is exhausted over the whole
  // branches before the next one is consulted.
  int cmp(Edge a, Edge b, int level) {
    if (!any_aux_) level = std::min<int>(level, k2);
    for (int r = k1a; r <= level; ++r) {
      const int c = r == k4b ? like_unlike(a, b) : bfs(a, b, r);
      if (c != 0) return c;
    }
    return 0;
  }

  // Breadth-first, sphere by sphere: the sets of each sphere are ordered by the
  // ranking of their parents, atoms within a set by their own (recursive) rank.
  int bfs(Edge a, Edge b, int level) {
    if (a.to == b.to && a.from == b.from) return 0;
    const std::array<int, 5> id{a.from, a.to, b.from, b.to, level};
    const auto memo = bfs_.find(id);
    if (memo != bfs_.end()) return memo->second;
    int result = node_cmp(a.to, b.to, level);
    if (result == 0) {
      std::vector<Edge> la{a}, lb{b};
      for (int sphere = 0; sphere < kMaxSpheres && result == 0 && !overflow_; ++sphere) {
        std::vector<Edge> na, nbv;
        bool any = false;
        const std::size_t m = std::max(la.size(), lb.size());
        for (std::size_t i = 0; i < m && result == 0; ++i) {
          const std::vector<int> sa = i < la.size() ? sorted(la[i], level) : std::vector<int>{};
          const std::vector<int> sb = i < lb.size() ? sorted(lb[i], level) : std::vector<int>{};
          any = any || !sa.empty() || !sb.empty();
          const std::size_t len = std::max(sa.size(), sb.size());
          for (std::size_t j = 0; j < len && result == 0; ++j)
            result = node_cmp(j < sa.size() ? sa[j] : -1, j < sb.size() ? sb[j] : -1, level);
          for (int x : sa) na.push_back({la[i].to, x});
          for (int x : sb) nbv.push_back({lb[i].to, x});
        }
        if (!any) break;
        la = std::move(na);
        lb = std::move(nbv);
      }
    }
    bfs_[id] = result;
    bfs_[{b.from, b.to, a.from, a.to, level}] = -result;
    return result;
  }

  // Rule 4b. The chiral descriptors of a branch, read in rank order, become a
  // string of like (1) / unlike (0) against the reference descriptor `ref`.
  // Substituents that tie under rules 1-4a are ordered like-first so the
  // string does not depend on how the molecule was written.
  std::vector<char> pairs(Edge e, CipLabel ref) {
    const std::array<int, 3> id{e.from, e.to, static_cast<int>(ref)};
    const auto memo = pairs_.find(id);
    if (memo != pairs_.end()) return memo->second;
    std::vector<char> out;
    std::vector<Edge> queue{e};
    for (std::size_t h = 0; h < queue.size(); ++h) {
      const Edge x = queue[h];
      const CipLabel l = nodes_[static_cast<std::size_t>(x.to)].aux;
      if (l == CipLabel::R || l == CipLabel::S) out.push_back(l == ref ? 1 : 0);
      std::vector<int> ks = sorted(x, k4a);
      for (std::size_t i = 0; i < ks.size();) {
        std::size_t j = i + 1;
        while (j < ks.size() && cmp({x.to, ks[i]}, {x.to, ks[j]}, k4a) == 0) ++j;
        if (j - i > 1)
          std::stable_sort(ks.begin() + static_cast<std::ptrdiff_t>(i), ks.begin() + static_cast<std::ptrdiff_t>(j),
                           [&](int p, int q) { return pairs({x.to, p}, ref) > pairs({x.to, q}, ref); });
        i = j;
      }
      for (int k : ks) queue.push_back({x.to, k});
    }
    pairs_[id] = out;
    return out;
  }
  // The reference is the highest-ranked chiral descriptor of the branch; when R
  // and S tie for that place both are tried and the like-richer string kept.
  std::vector<char> like_string(Edge e) {
    std::vector<char> lr = pairs(e, CipLabel::R);
    if (lr.empty()) return lr;
    std::vector<char> ls = pairs(e, CipLabel::S);
    if (lr[0] && ls[0]) return std::max(lr, ls);
    return lr[0] ? lr : ls;
  }
  int like_unlike(Edge a, Edge b) {
    const std::vector<char> la = like_string(a), lb = like_string(b);
    if (la == lb) return 0;
    return la > lb ? 1 : -1;
  }

  const ChemicalSystemState& st_;
  std::vector<Node> nodes_;
  std::vector<std::vector<int>> kids_;
  std::vector<char> expanded_;
  std::vector<int> tetra_of_;
  std::map<std::uint32_t, int> db_of_;
  std::map<std::uint32_t, std::pair<std::uint32_t, std::uint32_t>> mancude_;
  std::map<std::array<int, 3>, std::vector<int>> sorted_;
  std::map<std::array<int, 5>, int> bfs_;
  std::map<std::array<int, 3>, std::vector<char>> pairs_;
  int root_ = 0;
  bool root_is_centre_ = false;
  bool any_aux_ = false;
  bool overflow_ = false;
};

}  // namespace

const char* to_string(CipLabel l) noexcept {
  switch (l) {
    case CipLabel::R: return "R";
    case CipLabel::S: return "S";
    case CipLabel::E: return "E";
    case CipLabel::Z: return "Z";
    case CipLabel::r: return "r";
    case CipLabel::s: return "s";
    case CipLabel::Ra: return "Ra";
    case CipLabel::Sa: return "Sa";
    default: return "";
  }
}

int cip_compare(const ChemicalSystemState& st, std::uint32_t root, std::uint32_t a, std::uint32_t b) {
  Digraph g(st, root);
  const int na = g.root_ligand(a), nb = g.root_ligand(b);
  if (na < 0 || nb < 0) return (na >= 0) - (nb >= 0);
  int c = g.compare_root_ligands(na, nb, k2);
  if (c == 0 && !g.overflow()) {
    g.compute_aux();
    c = g.compare_root_ligands(na, nb, k5);
  }
  return c;
}

CipLabel cip_label(const ChemicalSystemState& st, const TetrahedralStereo& c, Diagnostics* diag) {
  if (c.parity == Chirality::None) return CipLabel::None;
  Digraph g(st, c.atom);
  CipLabel l = g.label_tetra(0, k2);
  if (l == CipLabel::None && !g.overflow()) {
    g.compute_aux();
    l = g.label_tetra(0, k5);
  }
  if (l == CipLabel::None && diag)
    diag->warn(DiagCode::CipUndetermined, g.overflow() ? "CIP digraph too large; descriptor not assigned"
                                                       : "two substituents tie under CIP rules 1-5; the centre is not stereogenic")
        .on_atom(c.atom);
  return l;
}

CipLabel cip_label(const ChemicalSystemState& st, const DoubleBondStereo& d, Diagnostics* diag) {
  const std::uint32_t u = st.bonds.begin[d.bond], v = st.bonds.end[d.bond];
  Digraph g(st, u);
  const int nv = g.root_ligand(v);
  CipLabel l = nv < 0 ? CipLabel::None : g.label_db(0, nv, k2);
  if (l == CipLabel::None && nv >= 0 && !g.overflow()) {
    g.compute_aux();
    l = g.label_db(0, nv, k5);
  }
  if (l == CipLabel::None && diag)
    diag->warn(DiagCode::CipUndetermined, g.overflow() ? "CIP digraph too large; descriptor not assigned"
                                                       : "double bond has equivalent substituents on one end; not stereogenic")
        .on_atom(u);
  return l;
}

CipLabel cip_label(const ChemicalSystemState& st, const AxialStereo& a, Diagnostics* diag) {
  if (a.parity == Chirality::None) return CipLabel::None;
  const int near = cip_compare(st, a.begin, a.refs[0], a.refs[1]);
  const int far = cip_compare(st, a.end, a.refs[2], a.refs[3]);
  if (near == 0 || far == 0) {
    if (diag)
      diag->warn(DiagCode::CipUndetermined, "two ligands on one end of the chirality axis tie under CIP rules 1-5; not stereogenic")
          .on_atom(a.centre != kNoAtom ? a.centre : a.begin);
    return CipLabel::None;
  }
  // Priorities 1..4: near pair, then far pair. Read like a tetrahedral centre: the refs reordered as
  // (4, 1, 2, 3); looking from 4, 1 -> 2 -> 3 counterclockwise is R.
  std::array<int, 4> prio = {near > 0 ? 1 : 2, near > 0 ? 2 : 1, far > 0 ? 3 : 4, far > 0 ? 4 : 3};
  std::array<int, 4> order{};  // order[i] = index of the ref with priority target[i]
  const std::array<int, 4> target = {4, 1, 2, 3};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      if (prio[static_cast<std::size_t>(j)] == target[static_cast<std::size_t>(i)]) order[static_cast<std::size_t>(i)] = j;
  int inversions = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j) inversions += order[static_cast<std::size_t>(i)] > order[static_cast<std::size_t>(j)];
  const bool ccw = (a.parity == Chirality::CounterClockwise) != ((inversions & 1) != 0);
  return ccw ? CipLabel::Ra : CipLabel::Sa;
}

std::vector<CipAssignment> cip_axes(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics* diag) {
  std::vector<CipAssignment> out;
  for (const auto& a : st.axial) {
    if (st.atoms.species[a.begin] != sp) continue;
    std::uint32_t index = a.centre;
    if (index == kNoAtom)
      for (std::uint32_t b : st.adjacency.bonds_of(a.begin))
        if (st.bonds.begin[b] == a.end || st.bonds.end[b] == a.end) index = b;
    out.push_back({index, cip_label(st, a, diag)});
  }
  return out;
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

}  // namespace chemlarp
