#include "chemlarp/inchi.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "chemlarp/elements.hpp"
#include "chemlarp/smiles.hpp"

namespace chemlarp {
namespace {

constexpr std::size_t kMaxLeaves = 40000;   // minimal numberings kept per component
constexpr std::size_t kMaxNodes = 400000;   // search-tree nodes per component

struct MBond {
  int a = 0, b = 0;
  int order = 1;
  bool arom = false, ring = false, taut = false;
  int ring_size = 0;       // smallest ring through the bond, 0 if acyclic
  std::uint32_t global = 0;
};

struct Centre {
  int atom = 0;
  std::array<int, 4> refs{};  // heavy-atom indices; -1 = hydrogen or lone pair
  int parity = 0;             // 1: '@' in refs order, 2: '@@', 0: undefined candidate
  bool active = true;
};

struct DBond {
  int a = 0, b = 0;
  int ref_a = -1, ref_b = -1;  // heavy-atom indices; -1 = hydrogen
  bool cis = false;
  bool active = true;
  bool undefined = false;  // geometry not given in the input
};

struct Mol {
  std::vector<int> el, h, charge, iso;
  std::vector<std::array<int, 2>> hiso;  // deuterium and tritium counts (included in h)
  std::vector<std::vector<std::pair<int, int>>> adj;  // (neighbour, bond index)
  std::vector<MBond> bonds;
  std::vector<int> tgroup;  // mobile-H group id per atom, -1 if none
  std::vector<int> comp;    // connected-component id
  std::vector<Centre> centres;
  std::vector<DBond> dbonds;
  std::size_t size() const { return el.size(); }
};

bool is_metal(int z) {
  if (z == 1 || z == 2) return false;
  if (z >= 5 && z <= 10) return false;
  if (z >= 14 && z <= 18) return false;
  if (z >= 32 && z <= 36) return false;
  if (z == 52 || z == 53 || z == 54 || z == 85 || z == 86) return false;
  return true;
}

bool is_halogen(int z) { return z == 9 || z == 17 || z == 35 || z == 53; }
bool is_endpoint_el(int z) { return z == 7 || z == 8 || z == 16 || z == 34 || z == 52; }

// Position in Hill order: carbon first, the rest alphabetical by symbol.
int hill_index(int z) {
  static const std::array<int, kElementSlots> table = [] {
    std::vector<int> zs;
    for (int i = 1; i < static_cast<int>(kElementSlots); ++i) zs.push_back(i);
    std::sort(zs.begin(), zs.end(), [](int x, int y) {
      if ((x == 6) != (y == 6)) return x == 6;
      return element_symbol(static_cast<std::uint8_t>(x)) < element_symbol(static_cast<std::uint8_t>(y));
    });
    std::array<int, kElementSlots> t{};
    for (std::size_t i = 0; i < zs.size(); ++i) t[static_cast<std::size_t>(zs[i])] = static_cast<int>(i);
    return t;
  }();
  return table[static_cast<std::size_t>(z)];
}

int rounded_mass(int z) {
  const ElementInfo* e = find_element_by_z(static_cast<std::uint8_t>(z));
  return e ? static_cast<int>(e->mass + 0.5) : 0;
}

// ---- 1. normalisation ----------------------------------------------------------

Mol build(const ChemicalSystemState& st, const Species& s, bool& metal_cut) {
  Mol m;
  const std::uint32_t base = s.atom_begin;
  std::vector<int> local(s.atom_count, -1), host(s.atom_count, -1);
  // which hydrogens fold into a neighbour's count
  for (std::uint32_t g = base; g < base + s.atom_count; ++g) {
    if (st.atoms.element[g] != 1 || st.atoms.formal_charge[g] != 0 || st.atoms.implicit_h[g] != 0) continue;
    const auto nb = st.adjacency.neighbors(g);
    if (nb.size() != 1) continue;
    const std::uint32_t o = nb[0];
    if (st.atoms.element[o] == 1 && (o > g || host[o - base] >= 0)) continue;  // H-H: keep the first as a vertex
    if (is_metal(st.atoms.element[o])) continue;  // metal hydrides are disconnected below
    host[g - base] = static_cast<int>(o - base);
  }
  for (std::uint32_t g = base; g < base + s.atom_count; ++g) {
    if (host[g - base] >= 0) continue;
    local[g - base] = static_cast<int>(m.el.size());
    m.el.push_back(st.atoms.element[g]);
    m.h.push_back(st.atoms.implicit_h[g]);
    m.charge.push_back(st.atoms.formal_charge[g]);
    m.iso.push_back(st.atoms.isotope[g]);
    m.hiso.push_back({0, 0});
  }
  for (std::uint32_t g = base; g < base + s.atom_count; ++g) {
    if (host[g - base] < 0) continue;
    const auto v = static_cast<std::size_t>(local[static_cast<std::size_t>(host[g - base])]);
    ++m.h[v];
    if (st.atoms.isotope[g] == 2) ++m.hiso[v][0];
    if (st.atoms.isotope[g] == 3) ++m.hiso[v][1];
  }
  metal_cut = false;
  for (std::size_t v = 0, n0 = m.size(); v < n0; ++v) {
    if (!is_metal(m.el[v])) continue;
    for (; m.h[v] > 0; --m.h[v]) {  // M-H: the hydrogens become components of their own
      metal_cut = true;
      m.el.push_back(1);
      m.h.push_back(0);
      m.charge.push_back(0);
      m.iso.push_back(0);
      m.hiso.push_back({0, 0});
    }
  }
  m.adj.resize(m.size());
  const std::vector<int> charge0(m.charge.begin(), m.charge.end());  // a ligand that is already an anion keeps its charge
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    const int x = local[st.bonds.begin[b] - base], y = local[st.bonds.end[b] - base];
    if (x < 0 || y < 0) continue;
    const auto ux = static_cast<std::size_t>(x), uy = static_cast<std::size_t>(y);
    const int order = st.bonds.order[b];
    if (is_metal(m.el[ux]) || is_metal(m.el[uy])) {  // standard InChI disconnects metals
      metal_cut = true;
      const bool mx = is_metal(m.el[ux]), my = is_metal(m.el[uy]);
      if (mx != my) {
        const std::size_t metal = mx ? ux : uy, other = mx ? uy : ux;
        if (order == 1 && charge0[other] == 0 && (is_endpoint_el(m.el[other]) || is_halogen(m.el[other]))) {  // ionic: M-X -> M+ X-
          ++m.charge[metal];
          --m.charge[other];
        }
      }
      continue;
    }
    MBond mb;
    mb.a = x;
    mb.b = y;
    mb.order = order;
    mb.arom = (st.bonds.flags[b] & kBondAromatic) != 0;
    mb.ring = (st.bonds.flags[b] & kBondInRing) != 0;
    mb.global = b;
    m.adj[ux].push_back({y, static_cast<int>(m.bonds.size())});
    m.adj[uy].push_back({x, static_cast<int>(m.bonds.size())});
    m.bonds.push_back(mb);
  }
  // adjacent +/- pairs become a neutral multiple bond (nitro, N-oxide, diazo, ylide)
  for (MBond& b : m.bonds) {
    auto& ca = m.charge[static_cast<std::size_t>(b.a)];
    auto& cb = m.charge[static_cast<std::size_t>(b.b)];
    while ((ca > 0 && cb < 0) || (ca < 0 && cb > 0)) {
      ca += ca > 0 ? -1 : 1;
      cb += cb > 0 ? -1 : 1;
      ++b.order;
    }
  }
  // stereo
  auto ref_of = [&](std::uint32_t r) -> int {
    if (r == kImplicitHRef || r == kLonePairRef) return -1;
    return local[r - base];  // -1 for a folded explicit hydrogen
  };
  for (const auto& t : st.tetrahedral) {
    if (st.atoms.species[t.atom] != static_cast<std::uint32_t>(&s - st.species.data())) continue;
    if (t.parity == Chirality::None || local[t.atom - base] < 0) continue;
    Centre c;
    c.atom = local[t.atom - base];
    int blanks = 0;
    for (std::size_t i = 0; i < 4; ++i) blanks += (c.refs[i] = ref_of(t.refs[i])) < 0;
    if (blanks > 1) continue;  // two hydrogens: not a stereocentre
    bool cut = false;  // a ligand lost to metal disconnection
    for (int r : c.refs) {
      if (r < 0) continue;
      bool found = false;
      for (auto [nb, bi] : m.adj[static_cast<std::size_t>(c.atom)]) found |= nb == r;
      cut |= !found;
    }
    if (cut) continue;
    c.parity = t.parity == Chirality::CounterClockwise ? 1 : 2;
    m.centres.push_back(c);
  }
  for (const auto& d : st.double_bonds) {
    if (d.bond < s.bond_begin || d.bond >= s.bond_begin + s.bond_count) continue;
    DBond e;
    e.a = local[st.bonds.begin[d.bond] - base];
    e.b = local[st.bonds.end[d.bond] - base];
    if (e.a < 0 || e.b < 0) continue;
    e.ref_a = ref_of(d.ref_begin);
    e.ref_b = ref_of(d.ref_end);
    e.cis = d.cis;
    m.dbonds.push_back(e);
  }
  return m;
}

void find_components(Mol& m) {
  m.comp.assign(m.size(), -1);
  int n = 0;
  for (std::size_t i = 0; i < m.size(); ++i) {
    if (m.comp[i] >= 0) continue;
    std::vector<std::size_t> stack{i};
    m.comp[i] = n;
    while (!stack.empty()) {
      const std::size_t v = stack.back();
      stack.pop_back();
      for (auto [nb, bi] : m.adj[v])
        if (m.comp[static_cast<std::size_t>(nb)] < 0) {
          m.comp[static_cast<std::size_t>(nb)] = n;
          stack.push_back(static_cast<std::size_t>(nb));
        }
    }
    ++n;
  }
}

void ring_sizes(Mol& m) {
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi) {
    MBond& b = m.bonds[bi];
    if (!b.ring) continue;
    std::vector<int> dist(m.size(), -1);
    std::queue<int> q;
    dist[static_cast<std::size_t>(b.a)] = 0;
    q.push(b.a);
    while (!q.empty() && dist[static_cast<std::size_t>(b.b)] < 0) {
      const int v = q.front();
      q.pop();
      for (auto [nb, e] : m.adj[static_cast<std::size_t>(v)]) {
        if (static_cast<std::size_t>(e) == bi || dist[static_cast<std::size_t>(nb)] >= 0) continue;
        dist[static_cast<std::size_t>(nb)] = dist[static_cast<std::size_t>(v)] + 1;
        q.push(nb);
      }
    }
    b.ring_size = dist[static_cast<std::size_t>(b.b)] + 1;
  }
}

// Move protons so every component is as neutral as it can be; returns protons added per component.
std::vector<int> neutralise(Mol& m, int ncomp) {
  std::vector<int> protons(static_cast<std::size_t>(ncomp), 0);
  auto acid = [&](std::size_t v) {  // cation that can lose H+
    const int z = m.el[v];
    return m.charge[v] > 0 && m.h[v] > 0 && (z == 7 || z == 15 || z == 8 || z == 16);
  };
  auto unsaturated = [&](std::size_t v) {
    for (auto [o, bi] : m.adj[v])
      if (m.bonds[static_cast<std::size_t>(bi)].order > 1 || m.bonds[static_cast<std::size_t>(bi)].arom) return true;
    return false;
  };
  auto oxo = [&](std::size_t v, std::size_t except) {  // v carries =O, =S or =N
    for (auto [o, bi] : m.adj[v]) {
      const auto uo = static_cast<std::size_t>(o);
      if (uo != except && m.bonds[static_cast<std::size_t>(bi)].order == 2 && is_endpoint_el(m.el[uo])) return true;
    }
    return false;
  };
  // Sites whose conjugate acid standard InChI writes: acidic O, any S/Se/Te, amide-type N, free halide.
  auto acidic_site = [&](std::size_t v) {
    const int z = m.el[v];
    if (is_halogen(z)) return m.adj[v].empty();
    if (z == 16 || z == 34 || z == 52) return !(m.adj[v].empty() && m.h[v] == 0 && m.charge[v] < -1);
    if (z == 8) {
      if (m.adj[v].empty()) return m.h[v] > 0 || m.charge[v] == 0;
      if (m.adj[v].size() != 1) return false;
      const auto o = static_cast<std::size_t>(m.adj[v][0].first);
      const int zo = m.el[o];
      if (zo == 6) return unsaturated(o);
      if (zo == 7 && oxo(o, v)) return false;  // nitrate, nitro: stays an anion
      return !(zo == 5 || zo == 14 || is_halogen(zo) || is_metal(zo));
    }
    if (z == 7) {
      for (auto [o, bi] : m.adj[v])
        if (m.el[static_cast<std::size_t>(o)] == 6 && oxo(static_cast<std::size_t>(o), v)) return true;
    }
    return false;
  };
  auto base = [&](std::size_t v) { return m.charge[v] < 0 && acidic_site(v); };  // anion that can take H+
  auto oxoacid = [&](std::size_t v) {  // neutral X(=O)-OH
    if (m.charge[v] != 0 || m.h[v] == 0 || m.el[v] != 8 || m.adj[v].size() != 1) return false;
    return oxo(static_cast<std::size_t>(m.adj[v][0].first), v);
  };
  for (int c = 0; c < ncomp; ++c) {
    std::vector<std::size_t> acids, bases;
    int net = 0;
    for (std::size_t v = 0; v < m.size(); ++v) {
      if (m.comp[v] != c) continue;
      net += m.charge[v];
      if (acid(v)) acids.push_back(v);
      if (base(v)) bases.push_back(v);
    }
    // zwitterions: proton transfer inside the component
    std::size_t ai = 0, bi = 0;
    while (ai < acids.size() && bi < bases.size()) {
      --m.h[acids[ai]];
      --m.charge[acids[ai]];
      ++m.h[bases[bi]];
      ++m.charge[bases[bi]];
      if (!acid(acids[ai])) ++ai;
      if (!base(bases[bi])) ++bi;
    }
    for (; net < 0 && bi < bases.size(); ++net) {
      ++m.h[bases[bi]];
      ++m.charge[bases[bi]];
      ++protons[static_cast<std::size_t>(c)];
      if (!base(bases[bi])) ++bi;
    }
    for (; net > 0 && ai < acids.size(); --net) {
      --m.h[acids[ai]];
      --m.charge[acids[ai]];
      --protons[static_cast<std::size_t>(c)];
      if (!acid(acids[ai])) ++ai;
    }
    // a cation that cannot lose a proton is balanced by ionising an oxoacid (betaine form)
    for (std::size_t v = 0; v < m.size() && net > 0; ++v) {
      if (m.comp[v] != c || m.charge[v] != 0 || m.h[v] == 0 || !(oxoacid(v) || acidic_site(v))) continue;
      if (m.el[v] == 7) {  // amide N-H ionises, amidine and guanidine N-H do not
        bool amide = false;
        for (auto [o, bo] : m.adj[v])
          for (auto [x, bx] : m.adj[static_cast<std::size_t>(o)])
            amide |= m.el[static_cast<std::size_t>(o)] == 6 && m.bonds[static_cast<std::size_t>(bx)].order == 2 && m.el[static_cast<std::size_t>(x)] != 7 &&
                     is_endpoint_el(m.el[static_cast<std::size_t>(x)]);
        if (!amide) continue;
      }
      --m.h[v];
      --m.charge[v];
      --protons[static_cast<std::size_t>(c)];
      --net;
    }
  }
  return protons;
}

// ---- 2. mobile hydrogen ------------------------------------------------------------

void find_mobile_h(Mol& m) {
  const std::size_t n = m.size();
  std::vector<int> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> find = [&](int x) { return parent[static_cast<std::size_t>(x)] == x ? x : parent[static_cast<std::size_t>(x)] = find(parent[static_cast<std::size_t>(x)]); };
  std::vector<char> in_group(n, 0), has_double(n, 0);
  for (const MBond& b : m.bonds)
    if (b.order == 2) has_double[static_cast<std::size_t>(b.a)] = has_double[static_cast<std::size_t>(b.b)] = 1;
  auto endpoint = [&](int v) {  // a negative charge moves with the hydrogens
    const int q = m.charge[static_cast<std::size_t>(v)];
    if (m.el[static_cast<std::size_t>(v)] != 7 && m.adj[static_cast<std::size_t>(v)].size() != 1) return false;
    return is_endpoint_el(m.el[static_cast<std::size_t>(v)]) && (q == 0 || q == -1);
  };
  auto donor = [&](int v) { return endpoint(v) && (m.h[static_cast<std::size_t>(v)] > 0 || m.charge[static_cast<std::size_t>(v)] == -1 || in_group[static_cast<std::size_t>(v)]); };
  auto single_ok = [&](const MBond& b) { return b.order == 1 || b.arom || b.taut; };

  // A bond drawn single can be double in another tautomer or Kekule structure
  // when the atoms left without a partner can be re-paired along an alternating
  // path of ring bonds; a mobile-H endpoint ends the path (it takes the hydrogen).
  std::vector<int> mate(n, -1);
  for (const MBond& b : m.bonds)
    if (b.order == 2) mate[static_cast<std::size_t>(b.a)] = b.b, mate[static_cast<std::size_t>(b.b)] = b.a;
  std::vector<signed char> memo(m.bonds.size(), -1);
  auto flexible = [&](int v) { return v >= 0 && in_group[static_cast<std::size_t>(v)] && endpoint(v); };
  auto feasible = [&](const MBond& b) {
    const std::size_t bi = static_cast<std::size_t>(&b - m.bonds.data());
    if (memo[bi] >= 0) return memo[bi] == 1;
    bool ok = false;
    for (int turn = 0; turn < 2 && !ok; ++turn) {
      const int qa = turn ? b.b : b.a, ya = turn ? b.a : b.b;  // free qa's partner; ya's partner is the goal
      const int s0 = mate[static_cast<std::size_t>(qa)], goal = mate[static_cast<std::size_t>(ya)];
      if (s0 < 0 || (goal < 0 && !flexible(ya))) continue;
      if (flexible(s0) && (goal < 0 || flexible(goal))) {
        ok = true;
        break;
      }
      std::vector<char> seen(n, 0);
      seen[static_cast<std::size_t>(b.a)] = seen[static_cast<std::size_t>(b.b)] = seen[static_cast<std::size_t>(s0)] = 1;
      int budget = 20000;
      std::function<bool(int, bool)> walk = [&](int u, bool need_goal) {
        if (--budget < 0) return true;  // give up: assume it can
        for (auto [v, e] : m.adj[static_cast<std::size_t>(u)]) {
          const MBond& eb = m.bonds[static_cast<std::size_t>(e)];
          if (eb.order == 2 || !eb.ring || seen[static_cast<std::size_t>(v)]) continue;
          if (v == goal) return true;
          const int w = mate[static_cast<std::size_t>(v)];
          if (w < 0 || seen[static_cast<std::size_t>(w)] || w == goal) continue;
          if (flexible(w) && !need_goal) return true;
          seen[static_cast<std::size_t>(v)] = seen[static_cast<std::size_t>(w)] = 1;
          if (walk(w, need_goal)) return true;
          seen[static_cast<std::size_t>(v)] = seen[static_cast<std::size_t>(w)] = 0;
        }
        return false;
      };
      // with a goal the freed partners must pair up with each other, unless the goal itself is an endpoint
      ok = walk(s0, goal >= 0 && !flexible(goal));
    }
    memo[bi] = ok;
    return ok;
  };
  auto double_ok = [&](const MBond& b) { return b.order == 2 || feasible(b); };
  auto acceptor = [&](int v, const MBond& b) { return endpoint(v) && (b.order == 2 || feasible(b)); };
  auto centre_ok = [&](int v) {
    const int z = m.el[static_cast<std::size_t>(v)];
    return z == 6 || z == 7 || z == 16 || z == 15 || z == 51 || z == 33 || z == 34 || z == 52 || z == 35 || z == 17 || z == 53;
  };
  bool changed = true;
  auto join = [&](int x, int y, std::initializer_list<int> path) {
    bool fresh = find(x) != find(y) || !in_group[static_cast<std::size_t>(x)] || !in_group[static_cast<std::size_t>(y)];
    for (int bi : path)
      if (!m.bonds[static_cast<std::size_t>(bi)].taut) m.bonds[static_cast<std::size_t>(bi)].taut = fresh = true;
    parent[static_cast<std::size_t>(find(x))] = find(y);
    in_group[static_cast<std::size_t>(x)] = in_group[static_cast<std::size_t>(y)] = 1;
    changed |= fresh;
  };
  while (changed) {
    changed = false;
    std::fill(memo.begin(), memo.end(), static_cast<signed char>(-1));
    for (int q = 0; q < static_cast<int>(n); ++q) {
      if (!centre_ok(q)) continue;
      const auto& nq = m.adj[static_cast<std::size_t>(q)];
      for (auto [x, bx] : nq) {
        if (!donor(x) || !single_ok(m.bonds[static_cast<std::size_t>(bx)])) continue;
        // 1,3 shift: X(H)-Q=Y
        const int zq0 = m.el[static_cast<std::size_t>(q)];
        const bool hetero_centre = zq0 == 16 || zq0 == 34 || zq0 == 52;  // S, Se, Te: only terminal atoms take part
        auto terminal = [&](int v) { return m.adj[static_cast<std::size_t>(v)].size() == 1; };
        for (auto [y, by] : nq)
          if (y != x && acceptor(y, m.bonds[static_cast<std::size_t>(by)]) && (!hetero_centre || (terminal(x) && terminal(y)))) join(x, y, {bx, by});
        // 1,5 shift X(H)-Q=B-C=Y, in the two ring patterns standard InChI knows:
        //   pyrazole:    X and Y are bonded nitrogens of one five-membered ring;
        //   4-pyridinol: a ring nitrogen and an exocyclic N/O/S para to it in a six-membered carbon ring.
        const int zq = m.el[static_cast<std::size_t>(q)];
        if (zq != 6 && zq != 7) continue;
        const bool x_ring = m.bonds[static_cast<std::size_t>(bx)].ring;
        for (auto [b2, bqb] : nq) {
          const MBond& e1 = m.bonds[static_cast<std::size_t>(bqb)];
          const int zb = m.el[static_cast<std::size_t>(b2)];
          if (b2 == x || !e1.ring || !double_ok(e1) || (zb != 6 && zb != 7)) continue;
          for (auto [c2, bbc] : m.adj[static_cast<std::size_t>(b2)]) {
            const MBond& e2 = m.bonds[static_cast<std::size_t>(bbc)];
            const int zc = m.el[static_cast<std::size_t>(c2)];
            if (c2 == q || c2 == x || !e2.ring || !single_ok(e2) || (zc != 6 && zc != 7)) continue;
            const bool carbons = zq == 6 && zb == 6 && zc == 6;
            for (auto [y, bcy] : m.adj[static_cast<std::size_t>(c2)]) {
              const MBond& e3 = m.bonds[static_cast<std::size_t>(bcy)];
              if (y == b2 || y == x || y == q || !acceptor(y, e3)) continue;
              bool ok = false;
              const bool both_n = m.el[static_cast<std::size_t>(x)] == 7 && m.el[static_cast<std::size_t>(y)] == 7;
              if (both_n && x_ring && e3.ring)  // pyrazole
                for (auto [o, bo] : m.adj[static_cast<std::size_t>(x)]) ok |= o == y && single_ok(m.bonds[static_cast<std::size_t>(bo)]);
              if (!ok && carbons && x_ring != e3.ring && m.el[static_cast<std::size_t>(x_ring ? x : y)] == 7) {  // 4-pyridinol
                // the ring (five or six atoms) must be unsaturated all the way round
                const int from = x_ring ? c2 : q, to = x_ring ? x : y;
                auto in_path = [&](int v) { return v == x || v == y || v == q || v == b2 || v == c2; };
                for (auto [r, br] : m.adj[static_cast<std::size_t>(from)]) {
                  if (in_path(r) || !m.bonds[static_cast<std::size_t>(br)].ring || !has_double[static_cast<std::size_t>(r)]) continue;
                  for (auto [t, bt] : m.adj[static_cast<std::size_t>(r)]) {
                    if (!m.bonds[static_cast<std::size_t>(bt)].ring) continue;
                    if (t == to) ok = true;
                    if (in_path(t) || t == from || !has_double[static_cast<std::size_t>(t)]) continue;
                    for (auto [u, bu] : m.adj[static_cast<std::size_t>(t)]) ok |= u == to && m.bonds[static_cast<std::size_t>(bu)].ring;
                  }
                }
              }
              if (ok) join(x, y, {bx, bqb, bbc, bcy});
            }
          }
        }
      }
    }
  }
  m.tgroup.assign(n, -1);
  for (std::size_t v = 0; v < n; ++v)
    if (in_group[v]) m.tgroup[v] = find(static_cast<int>(v));
}

// ---- 3. canonical numbering ----------------------------------------------------------
// Partition refinement plus individualisation; every discrete partition (leaf)
// is scored and the minimal ones are kept, so the stereo layers can later be
// minimised over the whole automorphism class.

struct Canon {
  const Mol& m;
  std::vector<int> verts;                 // molecule atom per local index
  std::vector<int> loc;                   // local index per molecule atom, -1 outside
  std::vector<std::vector<int>> nb;       // local adjacency
  std::vector<int> best;                  // minimal leaf key
  std::vector<std::vector<int>> leaves;   // 0-based canonical position per local index
  std::size_t nodes = 0;
  bool truncated = false;

  Canon(const Mol& mol, int comp) : m(mol), loc(mol.size(), -1) {
    for (std::size_t v = 0; v < m.size(); ++v)
      if (m.comp[v] == comp) {
        loc[v] = static_cast<int>(verts.size());
        verts.push_back(static_cast<int>(v));
      }
    nb.resize(verts.size());
    for (std::size_t i = 0; i < verts.size(); ++i)
      for (auto [o, bi] : m.adj[static_cast<std::size_t>(verts[i])]) nb[i].push_back(loc[static_cast<std::size_t>(o)]);
  }

  std::size_t n() const { return verts.size(); }
  std::size_t atom(std::size_t i) const { return static_cast<std::size_t>(verts[i]); }
  int fixed_h(std::size_t i) const { return m.tgroup[atom(i)] >= 0 ? 0 : m.h[atom(i)]; }

  void refine(std::vector<int>& r) const {
    const std::size_t cnt = n();
    std::vector<std::vector<int>> key(cnt);
    std::vector<std::size_t> order(cnt);
    std::size_t cells = 0;
    for (;;) {
      for (std::size_t v = 0; v < cnt; ++v) {
        key[v].assign(1, r[v]);
        for (int o : nb[v]) key[v].push_back(r[static_cast<std::size_t>(o)]);
        std::sort(key[v].begin() + 1, key[v].end());
      }
      std::iota(order.begin(), order.end(), std::size_t{0});
      std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return key[x] < key[y]; });
      std::size_t nc = 0, first = 0;
      for (std::size_t i = 0; i < cnt; ++i) {
        if (i == 0 || key[order[i]] != key[order[i - 1]]) {
          first = i;
          ++nc;
        }
        r[order[i]] = static_cast<int>(first);
      }
      if (nc == cells) break;
      cells = nc;
    }
  }

  std::vector<int> leaf_key(const std::vector<int>& r) const {
    const std::size_t cnt = n();
    std::vector<std::size_t> inv(cnt);
    for (std::size_t v = 0; v < cnt; ++v) inv[static_cast<std::size_t>(r[v])] = v;
    std::vector<int> key, tmp;
    for (std::size_t i = 0; i < cnt; ++i) {  // connection table
      key.push_back(static_cast<int>(i));
      tmp.clear();
      for (int o : nb[inv[i]])
        if (r[static_cast<std::size_t>(o)] < static_cast<int>(i)) tmp.push_back(r[static_cast<std::size_t>(o)]);
      std::sort(tmp.begin(), tmp.end());
      key.insert(key.end(), tmp.begin(), tmp.end());
    }
    std::map<int, int> gid;  // mobile-H groups, labelled by first appearance
    for (std::size_t i = 0; i < cnt; ++i) {
      const int g = m.tgroup[atom(inv[i])];
      key.push_back(g < 0 ? 1000 : gid.try_emplace(g, static_cast<int>(gid.size()) + 1).first->second);
    }
    for (std::size_t i = 0; i < cnt; ++i) key.push_back(fixed_h(inv[i]) ? fixed_h(inv[i]) : 1000);  // atoms with H first
    for (std::size_t i = 0; i < cnt; ++i) {
      const std::size_t a = atom(inv[i]);
      key.push_back(m.iso[a] ? m.iso[a] - rounded_mass(m.el[a]) + 1000 : 1 << 20);  // labelled atoms first
      key.push_back(-m.hiso[a][1]);
      key.push_back(-m.hiso[a][0]);
    }
    return key;
  }

  // Interchangeable terminal atoms on one neighbour: branching over them only permutes labels.
  bool twins(const std::vector<std::size_t>& cell) const {
    const std::size_t a0 = atom(cell[0]);
    if (nb[cell[0]].size() != 1) return false;
    for (std::size_t v : cell) {
      const std::size_t a = atom(v);
      if (nb[v].size() != 1 || nb[v][0] != nb[cell[0]][0]) return false;
      if (m.h[a] != m.h[a0] || m.iso[a] != m.iso[a0] || m.hiso[a] != m.hiso[a0] || m.tgroup[a] != m.tgroup[a0]) return false;
      if (m.adj[a][0].second != m.adj[a0][0].second &&
          m.bonds[static_cast<std::size_t>(m.adj[a][0].second)].order != m.bonds[static_cast<std::size_t>(m.adj[a0][0].second)].order &&
          m.tgroup[a] < 0)
        return false;
    }
    return true;
  }

  void search(std::vector<int> r) {
    if (++nodes > kMaxNodes) {
      truncated = true;
      if (!leaves.empty()) return;
    }
    refine(r);
    const std::size_t cnt = n();
    std::vector<int> size(cnt, 0);
    for (int x : r) ++size[static_cast<std::size_t>(x)];
    std::size_t target = cnt;
    for (std::size_t i = 0; i < cnt; ++i)
      if (size[i] > 1) {
        target = i;
        break;
      }
    if (target == cnt) {
      std::vector<int> key = leaf_key(r);
      if (leaves.empty() || key < best) {
        best = std::move(key);
        leaves.assign(1, r);
      } else if (key == best && leaves.size() < kMaxLeaves) {
        leaves.push_back(r);
      } else if (key == best) {
        truncated = true;
      }
      return;
    }
    std::vector<std::size_t> cell;
    for (std::size_t v = 0; v < cnt; ++v)
      if (r[v] == static_cast<int>(target)) cell.push_back(v);
    const bool one = twins(cell) || (truncated && !leaves.empty());
    for (std::size_t v : cell) {
      std::vector<int> child = r;
      for (std::size_t o : cell)
        if (o != v) ++child[o];
      search(std::move(child));
      if (one) break;
    }
  }

  void run() {
    const std::size_t cnt = n();
    std::vector<std::size_t> order(cnt);
    std::iota(order.begin(), order.end(), std::size_t{0});
    auto colour = [&](std::size_t v) { return std::pair<int, std::size_t>(hill_index(m.el[atom(v)]), nb[v].size()); };
    std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return colour(x) < colour(y); });
    std::vector<int> r(cnt, 0);
    std::size_t first = 0;
    for (std::size_t i = 0; i < cnt; ++i) {
      if (i > 0 && colour(order[i]) != colour(order[i - 1])) first = i;
      r[order[i]] = static_cast<int>(first);
    }
    search(std::move(r));
  }
};

// ---- 4. layer writers ------------------------------------------------------------------

std::string hill_formula(const Canon& c) {
  std::map<std::string, int> count;
  int hydrogens = 0, carbons = 0;
  for (std::size_t i = 0; i < c.n(); ++i) {
    const std::size_t a = c.atom(i);
    hydrogens += c.m.h[a];
    if (c.m.el[a] == 6) ++carbons;
    else if (c.m.el[a] == 1) ++hydrogens;
    else ++count[std::string(element_symbol(static_cast<std::uint8_t>(c.m.el[a])))];
  }
  auto item = [](const std::string& sym, int k) { return k == 0 ? std::string() : k == 1 ? sym : sym + std::to_string(k); };
  std::string out;
  if (carbons > 0) {
    out = item("C", carbons) + item("H", hydrogens);
  } else if (hydrogens > 0) {
    count["H"] = hydrogens;
  }
  for (const auto& [sym, k] : count) out += item(sym, k);
  return out;
}

// Connection layer for the numbering `pos` (0-based canonical position per local atom).
std::string write_connections(const Canon& c, const std::vector<int>& pos) {
  const std::size_t n = c.n();
  if (n < 2) return {};
  std::vector<std::vector<int>> adj(n);
  for (std::size_t v = 0; v < n; ++v)
    for (int o : c.nb[v]) adj[static_cast<std::size_t>(pos[v])].push_back(pos[static_cast<std::size_t>(o)]);
  for (auto& a : adj) std::sort(a.begin(), a.end());
  std::size_t start = 0;
  for (std::size_t i = 1; i < n; ++i)
    if (adj[i].size() < adj[start].size()) start = i;
  std::vector<int> parent(n, -2), size(n, 1), order(n, -1);
  std::vector<std::vector<int>> kids(n), closures(n);
  int clock = 0;
  std::function<void(int, int)> dfs = [&](int v, int from) {
    const auto uv = static_cast<std::size_t>(v);
    parent[uv] = from;
    order[uv] = clock++;
    for (int o : adj[uv]) {
      const auto uo = static_cast<std::size_t>(o);
      if (o == from) continue;
      if (parent[uo] == -2) {
        dfs(o, v);
        kids[uv].push_back(o);
        size[uv] += size[uo];
      } else if (order[uo] < order[uv]) {
        closures[uv].push_back(o);  // back edge, written at the later atom
        ++size[uv];                 // branches are ordered by how many numbers they print
      }
    }
  };
  dfs(static_cast<int>(start), -1);
  std::function<std::string(int)> write = [&](int v) {
    const auto uv = static_cast<std::size_t>(v);
    std::string out = std::to_string(v + 1);
    std::vector<std::string> items;
    for (int o : closures[uv]) items.push_back(std::to_string(o + 1));
    std::stable_sort(kids[uv].begin(), kids[uv].end(),
                     [&](int x, int y) { return size[static_cast<std::size_t>(x)] < size[static_cast<std::size_t>(y)]; });
    for (int k : kids[uv]) items.push_back(write(k));
    if (items.empty()) return out;
    if (items.size() > 1) {
      out += '(';
      for (std::size_t i = 0; i + 1 < items.size(); ++i) out += (i ? "," : "") + items[i];
      out += ')';
      out += items.back();
    } else {
      out += '-' + items.back();
    }
    return out;
  };
  return write(static_cast<int>(start));
}

// "1-3,5" style list of 1-based numbers
std::string number_list(const std::vector<int>& xs) {
  std::string out;
  for (std::size_t i = 0; i < xs.size();) {
    std::size_t j = i;
    while (j + 1 < xs.size() && xs[j + 1] == xs[j] + 1) ++j;
    if (!out.empty()) out += ',';
    out += std::to_string(xs[i]);
    if (j > i) out += '-' + std::to_string(xs[j]);
    i = j + 1;
  }
  return out;
}

std::string write_hydrogens(const Canon& c, const std::vector<int>& pos) {
  std::map<int, std::vector<int>> by_count;
  std::map<int, std::pair<int, std::vector<int>>> groups;  // group id -> (H, members)
  std::map<int, int> minus;                                // group id -> negative charges
  for (std::size_t v = 0; v < c.n(); ++v) {
    const std::size_t a = c.atom(v);
    if (c.m.tgroup[a] >= 0) {
      auto& g = groups[c.m.tgroup[a]];
      g.first += c.m.h[a];
      if (c.m.charge[a] < 0) minus[c.m.tgroup[a]] -= c.m.charge[a];
      g.second.push_back(pos[v] + 1);
    } else if (c.m.h[a] > 0) {
      by_count[c.m.h[a]].push_back(pos[v] + 1);
    }
  }
  std::string out;
  for (auto& [k, xs] : by_count) {
    std::sort(xs.begin(), xs.end());
    if (!out.empty()) out += ',';
    out += number_list(xs) + 'H' + (k > 1 ? std::to_string(k) : "");
  }
  std::vector<std::vector<int>> mobile;
  for (auto& [id, g] : groups) {
    if (g.first == 0) continue;
    std::sort(g.second.begin(), g.second.end());
    std::vector<int> row{g.first * 16 + std::min(minus[id], 15)};
    row.insert(row.end(), g.second.begin(), g.second.end());
    mobile.push_back(row);
  }
  std::sort(mobile.begin(), mobile.end(), [](const auto& x, const auto& y) {
    if (x.size() != y.size()) return x.size() < y.size();
    return std::vector<int>(x.begin() + 1, x.end()) < std::vector<int>(y.begin() + 1, y.end());
  });
  if (!mobile.empty() && !out.empty()) out += ',';
  for (const auto& row : mobile) {
    const int hs = row[0] / 16, neg = row[0] % 16;
    out += "(H" + (hs > 1 ? std::to_string(hs) : std::string());
    if (neg > 0) out += (neg > 1 ? std::to_string(neg) : std::string()) + "-";
    for (std::size_t i = 1; i < row.size(); ++i) out += ',' + std::to_string(row[i]);
    out += ')';
  }
  return out;
}

// ---- 5. stereo -----------------------------------------------------------------------------

struct Stereo {
  const Canon& c;
  std::vector<Centre> centres;  // local indices
  std::vector<DBond> dbonds;    // local indices
  std::vector<int> t_entries, b_entries;  // flattened (number, parity) / (high, low, parity)
  int mirror = -1;                        // /m value, -1 when absent
  std::size_t leaf = 0;                   // numbering that realises the printed layers

  explicit Stereo(const Canon& canon) : c(canon) {}

  bool twin_pair(int x, int y) const {
    if (x < 0 || y < 0) return false;
    const auto ux = static_cast<std::size_t>(x), uy = static_cast<std::size_t>(y);
    const std::size_t a = c.atom(ux), b = c.atom(uy);
    return c.nb[ux].size() == 1 && c.nb[uy].size() == 1 && c.m.el[a] == c.m.el[b] && c.m.h[a] == c.m.h[b] &&
           c.m.iso[a] == c.m.iso[b] && c.m.hiso[a] == c.m.hiso[b] &&
           (c.m.tgroup[a] >= 0 ? c.m.tgroup[a] == c.m.tgroup[b]
                               : c.m.tgroup[b] < 0 && c.m.bonds[static_cast<std::size_t>(c.m.adj[a][0].second)].order ==
                                                          c.m.bonds[static_cast<std::size_t>(c.m.adj[b][0].second)].order);
  }

  void collect() {
    const Mol& m = c.m;
    std::vector<char> has_centre(c.n(), 0);
    auto centre_element = [&](std::size_t v) {
      const int z = m.el[c.atom(v)];
      const std::size_t links = c.nb[v].size() + static_cast<std::size_t>(m.h[c.atom(v)]);
      if (z == 6 || z == 14 || z == 32 || z == 50 || z == 5) return links == 4;
      if (z == 7 || z == 15 || z == 33) return links == 4 || (z != 7 && links == 3);
      if (z == 16 || z == 34) return links == 3 || links == 4;
      return false;
    };
    auto usable = [&](const Centre& x) {
      const auto v = static_cast<std::size_t>(x.atom);
      if (!centre_element(v) || m.tgroup[c.atom(v)] >= 0) return false;
      for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = i + 1; j < 4; ++j)
          if (twin_pair(x.refs[i], x.refs[j])) return false;
      return true;
    };
    for (Centre x : m.centres) {
      if (c.loc[static_cast<std::size_t>(x.atom)] < 0) continue;
      x.atom = c.loc[static_cast<std::size_t>(x.atom)];
      for (int& r : x.refs)
        if (r >= 0) r = c.loc[static_cast<std::size_t>(r)];
      if (!usable(x)) continue;
      has_centre[static_cast<std::size_t>(x.atom)] = 1;
      centres.push_back(x);
    }
    undefined_begin_ = centres.size();
    // undefined candidates: saturated carbon-group atoms with at most one hydrogen
    for (std::size_t v = 0; v < c.n(); ++v) {
      const int z = m.el[c.atom(v)];
      if (has_centre[v]) continue;
      const std::size_t deg = c.nb[v].size();
      const int hv = m.h[c.atom(v)];
      bool saturated = true, aromatic = false, tautomeric = false;
      for (auto [o, bi] : m.adj[c.atom(v)]) {
        const MBond& eb = m.bonds[static_cast<std::size_t>(bi)];
        saturated &= eb.order == 1;
        aromatic |= eb.arom;
        tautomeric |= eb.taut;
      }
      if (aromatic) continue;
      const bool carbon_like = (z == 6 || z == 14) && saturated && deg >= 3 && deg + static_cast<std::size_t>(hv) == 4;
      const bool onium = (z == 7 || z == 15) && m.charge[c.atom(v)] == 1 && saturated && deg == 4;
      const bool sulfur = (z == 16 || z == 34) && hv == 0 && !tautomeric && (deg == 3 || deg == 4) && m.charge[c.atom(v)] <= 1;
      const bool phosphorus = z == 15 && hv == 0 && !tautomeric && deg == 4 && m.charge[c.atom(v)] == 0;
      if (!carbon_like && !onium && !sulfur && !phosphorus) continue;
      Centre x;
      x.atom = static_cast<int>(v);
      x.refs = {-1, -1, -1, -1};
      for (std::size_t i = 0; i < c.nb[v].size(); ++i) x.refs[i + (c.nb[v].size() == 3 ? 1 : 0)] = c.nb[v][i];
      x.parity = 1;
      x.active = true;
      if (!usable(x)) continue;
      x.parity = -1;  // marks "undefined"; given a trial value during the symmetry test
      centres.push_back(x);
    }
    for (DBond d : m.dbonds) {
      if (c.loc[static_cast<std::size_t>(d.a)] < 0) continue;
      const MBond* bond = nullptr;
      for (auto [o, bi] : m.adj[static_cast<std::size_t>(d.a)])
        if (o == d.b) bond = &m.bonds[static_cast<std::size_t>(bi)];
      if (!bond || bond->order != 2 || bond->arom || bond->taut || (bond->ring && bond->ring_size < 8)) continue;
      d.a = c.loc[static_cast<std::size_t>(d.a)];
      d.b = c.loc[static_cast<std::size_t>(d.b)];
      if (d.ref_a >= 0) d.ref_a = c.loc[static_cast<std::size_t>(d.ref_a)];
      if (d.ref_b >= 0) d.ref_b = c.loc[static_cast<std::size_t>(d.ref_b)];
      bool ok = true;
      for (int end : {d.a, d.b}) {
        const int other = end == d.a ? d.b : d.a;
        std::vector<int> subs;
        for (int o : c.nb[static_cast<std::size_t>(end)])
          if (o != other) subs.push_back(o);
        if (subs.size() == 2 && twin_pair(subs[0], subs[1])) ok = false;
        if (subs.empty() && m.h[c.atom(static_cast<std::size_t>(end))] != 1) ok = false;
        if (subs.size() == 1 && m.h[c.atom(static_cast<std::size_t>(end))] > 1) ok = false;
      }
      if (ok) dbonds.push_back(d);
    }
    // double bonds whose geometry was not given
    for (const MBond& mb : m.bonds) {
      const int la = c.loc[static_cast<std::size_t>(mb.a)], lb = c.loc[static_cast<std::size_t>(mb.b)];
      if (la < 0 || lb < 0 || mb.order != 2 || mb.arom || mb.taut || (mb.ring && mb.ring_size < 8)) continue;
      bool ok = true;
      for (const DBond& d : dbonds) ok &= !((d.a == la && d.b == lb) || (d.a == lb && d.b == la));
      DBond d;
      d.a = la, d.b = lb, d.undefined = true;
      for (int end : {la, lb}) {
        const int other = end == la ? lb : la;
        const std::size_t ga = c.atom(static_cast<std::size_t>(end));
        const int z = m.el[ga];
        if (z != 6 && z != 7 && z != 14) ok = false;
        std::vector<int> subs;
        for (int o : c.nb[static_cast<std::size_t>(end)])
          if (o != other) subs.push_back(o);
        for (auto [o, bi] : m.adj[ga])  // cumulenes are not handled here
          if (c.loc[static_cast<std::size_t>(o)] != other && m.bonds[static_cast<std::size_t>(bi)].order != 1) ok = false;
        if (subs.size() == 2 && twin_pair(subs[0], subs[1])) ok = false;
        if (subs.empty() && (m.h[ga] != 1 || z == 7)) ok = false;
        if (subs.size() == 1 && m.h[ga] > 1) ok = false;
        if (subs.size() > 2) ok = false;
        (end == la ? d.ref_a : d.ref_b) = subs.empty() ? -1 : subs[0];
      }
      if (ok) dbonds.push_back(d);
    }
  }

  // parity value of one centre in a numbering: 0 '-', 1 '+', 2 '?'
  int centre_value(const Centre& x, const std::vector<int>& pos, bool invert) const {
    if (x.parity == 0) return 2;
    std::array<int, 4> k{};
    for (std::size_t i = 0; i < 4; ++i) k[i] = x.refs[i] < 0 ? -1 : pos[static_cast<std::size_t>(x.refs[i])];
    int inversions = 0;
    for (std::size_t i = 0; i < 4; ++i)
      for (std::size_t j = i + 1; j < 4; ++j) inversions += k[i] > k[j];
    const bool ccw = (x.parity == 1) != (inversions % 2 == 1);
    return (ccw ? 0 : 1) ^ (invert ? 1 : 0);
  }

  int bond_value(const DBond& d, const std::vector<int>& pos) const {
    bool cis = d.cis;
    auto top = [&](int end, int other) {
      int best = -1;
      for (int o : c.nb[static_cast<std::size_t>(end)])
        if (o != other && (best < 0 || pos[static_cast<std::size_t>(o)] > pos[static_cast<std::size_t>(best)])) best = o;
      return best;
    };
    if (top(d.a, d.b) != d.ref_a) cis = !cis;
    if (top(d.b, d.a) != d.ref_b) cis = !cis;
    return cis ? 0 : 1;
  }

  std::vector<int> vec(const std::vector<int>& pos, bool invert, std::vector<int>* b_out = nullptr, std::vector<int>* t_out = nullptr) const {
    std::vector<std::array<int, 3>> bs;
    std::vector<std::array<int, 2>> ts;
    for (const DBond& d : dbonds) {
      if (!d.active) continue;
      const int pa = pos[static_cast<std::size_t>(d.a)], pb = pos[static_cast<std::size_t>(d.b)];
      bs.push_back({std::max(pa, pb), std::min(pa, pb), unknown_bonds_ && d.undefined ? 2 : bond_value(d, pos)});
    }
    for (const Centre& x : centres)
      if (x.active) ts.push_back({pos[static_cast<std::size_t>(x.atom)], centre_value(x, pos, invert)});
    std::sort(bs.begin(), bs.end());
    std::sort(ts.begin(), ts.end());
    std::vector<int> out, b, t;
    for (const auto& e : bs) b.insert(b.end(), e.begin(), e.end());
    for (const auto& e : ts) t.insert(t.end(), e.begin(), e.end());
    out = b;
    out.push_back(-1);
    out.insert(out.end(), t.begin(), t.end());
    if (b_out) *b_out = b;
    if (t_out) *t_out = t;
    return out;
  }

  std::vector<int> best(bool invert, std::size_t* arg = nullptr) const {
    std::vector<int> out;
    for (std::size_t i = 0; i < c.leaves.size(); ++i) {
      std::vector<int> v = vec(c.leaves[i], invert);
      if (i == 0 || v < out) {
        out = std::move(v);
        if (arg) *arg = i;
      }
    }
    return out;
  }

  void solve() {
    collect();
    for (Centre& x : centres)
      if (x.parity < 0) x.parity = 1;  // trial value; restored to "undefined" below
    std::vector<char> undefined(centres.size(), 0);
    for (std::size_t i = undefined_begin_; i < centres.size(); ++i) undefined[i] = 1;
    // an element is stereogenic only if reversing it gives a different structure
    for (bool changed = true; changed;) {
      changed = false;
      const std::vector<int> base = best(false);
      for (Centre& x : centres) {
        if (!x.active) continue;
        x.parity = 3 - x.parity;
        const bool same = best(false) == base;
        x.parity = 3 - x.parity;
        if (same) {
          x.active = false;
          changed = true;
        }
      }
      for (DBond& d : dbonds) {
        if (!d.active) continue;
        d.cis = !d.cis;
        const bool same = best(false) == base;
        d.cis = !d.cis;
        if (same) {
          d.active = false;
          changed = true;
        }
      }
    }
    bool any_defined = false;
    for (std::size_t i = 0; i < centres.size(); ++i) {
      if (undefined[i]) centres[i].parity = 0;
      any_defined |= centres[i].active && !undefined[i];
    }
    // undefined elements are listed only next to a defined one of the same kind
    bool bond_defined = false;
    for (const DBond& d : dbonds) bond_defined |= d.active && !d.undefined;
    if (!any_defined)
      for (Centre& x : centres) x.active = false;
    if (!bond_defined)
      for (DBond& d : dbonds) d.active = false;
    unknown_bonds_ = true;
    if (!any_defined && !bond_defined) return;
    std::size_t la = 0, li = 0;
    const std::vector<int> actual = best(false, &la), inverted = best(true, &li);
    bool has_centres = false;
    for (const Centre& x : centres) has_centres |= x.active && x.parity != 0;
    if (!has_centres || actual == inverted) {
      leaf = la;
      vec(c.leaves[la], false, &b_entries, &t_entries);
    } else if (actual < inverted) {
      leaf = la;
      mirror = 0;
      vec(c.leaves[la], false, &b_entries, &t_entries);
    } else {
      leaf = li;
      mirror = 1;
      vec(c.leaves[li], true, &b_entries, &t_entries);
    }
  }

  bool unknown_bonds_ = false;  // false while testing which bonds are stereogenic
  std::size_t undefined_begin_ = 0;  // centres from here on have no parity in the input
};

// ---- 6. assembly ---------------------------------------------------------------------------

struct Layers {
  std::string formula, c, h, q, b, t, m, i;
  int heavy = 0, carbons = 0, hydrogens = 0;
  std::vector<std::pair<std::string, int>> hill;  // without hydrogen; carbon first
  std::vector<int> hvec;                           // fixed hydrogens in canonical order
  std::vector<int> ct;                             // connection table used to order components
  bool has_t = false;
};

std::string signed_number(int x) { return (x > 0 ? "+" : "") + std::to_string(x); }

std::string write_isotopes(const Canon& c, const std::vector<int>& pos) {
  std::map<int, std::string> rows;
  for (std::size_t v = 0; v < c.n(); ++v) {
    const std::size_t a = c.atom(v);
    std::string s;
    if (c.m.iso[a] && c.m.el[a] != 1) s += signed_number(c.m.iso[a] - rounded_mass(c.m.el[a]));
    if (c.m.el[a] == 1 && c.m.iso[a] > 1) s += signed_number(c.m.iso[a] - 1);
    auto hs = [&](int k, char sym) {
      if (k > 0) s += sym + (k > 1 ? std::to_string(k) : std::string());
    };
    hs(c.m.hiso[a][1], 'T');
    hs(c.m.hiso[a][0], 'D');
    if (!s.empty()) rows[pos[v] + 1] = s;
  }
  std::string out;
  for (const auto& [num, s] : rows) out += (out.empty() ? "" : ",") + std::to_string(num) + s;
  return out;
}

Layers describe(const Mol& m, int comp, Diagnostics* diag) {
  Canon c(m, comp);
  c.run();
  if (c.truncated && diag)
    diag->warn(DiagCode::InchiSymmetryLimit, "symmetry search truncated; the atom numbering may not be canonical");
  Stereo st(c);
  st.solve();
  const std::vector<int>& pos = c.leaves[st.leaf];
  Layers out;
  out.formula = hill_formula(c);
  out.c = write_connections(c, pos);
  out.h = write_hydrogens(c, pos);
  int charge = 0;
  for (std::size_t v = 0; v < c.n(); ++v) {
    charge += m.charge[c.atom(v)];
    out.heavy += m.el[c.atom(v)] != 1;
    out.carbons += m.el[c.atom(v)] == 6;
    out.hydrogens += m.h[c.atom(v)] + (m.el[c.atom(v)] == 1);
  }
  {
    std::map<std::string, int> counts;
    for (std::size_t v = 0; v < c.n(); ++v)
      if (const int z = m.el[c.atom(v)]; z != 1 && z != 6) ++counts[std::string(element_symbol(static_cast<std::uint8_t>(z)))];
    out.hvec.assign(c.n(), 0);
    for (std::size_t v = 0; v < c.n(); ++v) out.hvec[static_cast<std::size_t>(pos[v])] = m.h[c.atom(v)];
    if (out.carbons) out.hill.emplace_back("", out.carbons);  // "" sorts before every symbol
    for (const auto& kv : counts) out.hill.push_back(kv);
    // connection table: for each atom in canonical order, the atom and its lower-numbered neighbours
    std::vector<std::vector<int>> lower(c.n());
    for (std::size_t v = 0; v < c.n(); ++v)
      for (int o : c.nb[v])
        if (pos[static_cast<std::size_t>(o)] < pos[v]) lower[static_cast<std::size_t>(pos[v])].push_back(pos[static_cast<std::size_t>(o)]);
    for (std::size_t r = 1; r < c.n(); ++r) {
      std::sort(lower[r].begin(), lower[r].end());
      out.ct.push_back(static_cast<int>(r));
      out.ct.insert(out.ct.end(), lower[r].begin(), lower[r].end());
    }
  }
  if (charge != 0) out.q = signed_number(charge);
  for (std::size_t k = 0; k + 2 < st.b_entries.size(); k += 3)
    out.b += (out.b.empty() ? "" : ",") + std::to_string(st.b_entries[k] + 1) + '-' + std::to_string(st.b_entries[k + 1] + 1) +
             "-+?"[st.b_entries[k + 2]];
  for (std::size_t k = 0; k + 1 < st.t_entries.size(); k += 2)
    out.t += (out.t.empty() ? "" : ",") + std::to_string(st.t_entries[k] + 1) + "-+?"[st.t_entries[k + 1]];
  out.has_t = !out.t.empty();
  if (st.mirror >= 0) out.m = std::to_string(st.mirror);
  out.i = write_isotopes(c, pos);
  return out;
}

// "a;b;b" -> "a;2*b"; empty when no component contributes
std::string join_layer(const std::vector<Layers>& parts, std::string Layers::*field, char sep = ';', bool multiply = true) {
  bool any = false;
  for (const Layers& p : parts) any |= !(p.*field).empty();
  if (!any) return {};
  std::string out;
  for (std::size_t i = 0; i < parts.size();) {
    std::size_t j = i;
    const std::string& s = parts[i].*field;
    while (multiply && !s.empty() && j + 1 < parts.size() && parts[j + 1].*field == s) ++j;
    if (i > 0) out += sep;
    if (j > i) out += std::to_string(j - i + 1) + '*';
    out += s;
    i = j + 1;
  }
  return out;
}

}  // namespace

InchiResult inchi(const ChemicalSystemState& state, std::uint32_t species, Diagnostics* diag) {
  InchiResult res;
  if (species >= state.species.size()) return res;
  bool metal_cut = false;
  Mol m = build(state, state.species[species], metal_cut);
  if (metal_cut && diag) diag->warn(DiagCode::InchiMetalDisconnected, "bonds to metal atoms were disconnected (standard InChI)");
  if (diag)
    for (const auto& ax : state.axial)
      if (state.atoms.species[ax.begin] == species) {
        diag->warn(DiagCode::StereoIgnored, "axial chirality (allene or atropisomer) is not encoded in this InChI").on_atom(ax.begin);
        break;
      }
  if (diag)
    for (const auto& u : state.planar)
      if (!u.atoms.empty() && state.atoms.species[u.atoms[0]] == species) {
        diag->warn(DiagCode::StereoIgnored, "planar or helical chirality is not encoded in InChI").on_atom(u.atoms[0]);
        break;
      }
  find_components(m);
  int ncomp = 0;
  for (int x : m.comp) ncomp = std::max(ncomp, x + 1);
  ring_sizes(m);
  const std::vector<int> added = neutralise(m, ncomp);
  int protons = 0;
  for (int k : added) protons -= k;
  find_mobile_h(m);
  std::vector<Layers> parts;
  for (int k = 0; k < ncomp; ++k) {
    std::size_t first = 0;
    int count = 0;
    for (std::size_t v = 0; v < m.size(); ++v)
      if (m.comp[v] == k && count++ == 0) first = v;
    if (count == 1 && m.el[first] == 1 && m.h[first] == 0 && m.charge[first] == 1) {  // bare proton
      ++protons;
      continue;
    }
    parts.push_back(describe(m, k, diag));
  }
  // Component order: more carbons; then the Hill formula without hydrogen, element by
  // element; then more hydrogens; then the longer and the larger connection table.
  std::stable_sort(parts.begin(), parts.end(), [](const Layers& x, const Layers& y) {
    if (x.carbons != y.carbons) return x.carbons > y.carbons;
    for (std::size_t i = 0; i < x.hill.size() && i < y.hill.size(); ++i) {
      if (x.hill[i].first != y.hill[i].first) return x.hill[i].first < y.hill[i].first;
      if (x.hill[i].second != y.hill[i].second) return x.hill[i].second > y.hill[i].second;
    }
    if (x.hill.size() != y.hill.size()) return x.hill.size() > y.hill.size();
    if (x.hydrogens != y.hydrogens) return x.hydrogens > y.hydrogens;
    if (x.ct.size() != y.ct.size()) return x.ct.size() > y.ct.size();
    if (x.ct != y.ct) return x.ct > y.ct;
    if (x.hvec != y.hvec) return x.hvec > y.hvec;
    if (x.h != y.h) return x.h > y.h;
    if (x.q != y.q) return x.q < y.q;
    if (x.b.empty() != y.b.empty()) return !x.b.empty();
    if (x.t.empty() != y.t.empty()) return !x.t.empty();
    if (x.i.empty() != y.i.empty()) return !x.i.empty();
    return std::tie(x.b, x.t, y.m, x.i) < std::tie(y.b, y.t, x.m, y.i);
  });
  std::string out = "InChI=1S/";
  for (std::size_t i = 0; i < parts.size();) {
    std::size_t j = i;
    while (j + 1 < parts.size() && parts[j + 1].formula == parts[i].formula) ++j;
    if (i > 0) out += '.';
    if (j > i) out += std::to_string(j - i + 1);
    out += parts[i].formula;
    i = j + 1;
  }
  auto add = [&](const char* tag, const std::string& body) {
    if (!body.empty()) out += std::string("/") + tag + body;
  };
  add("c", join_layer(parts, &Layers::c));
  add("h", join_layer(parts, &Layers::h));
  add("q", join_layer(parts, &Layers::q));
  if (protons != 0) {
    if (out.back() == '/') out.pop_back();
    out += "/p" + signed_number(protons);
  }
  add("b", join_layer(parts, &Layers::b));
  const std::string t = join_layer(parts, &Layers::t);
  add("t", t);
  std::string mm;  // one character per component, '.' where a component has no /m
  bool any_m = false;
  for (const Layers& p : parts) {
    any_m |= !p.m.empty();
    mm += p.m.empty() ? std::string(".") : p.m;
  }
  if (!any_m) mm.clear();
  if (!mm.empty()) {
    add("m", mm);
    out += "/s1";
  }
  add("i", join_layer(parts, &Layers::i));
  if (out == "InChI=1S/") out.pop_back(), out += "/";
  res.ok = true;
  res.inchi = out;
  res.key = inchi_key(out);
  return res;
}

InchiResult inchi_from_smiles(std::string_view smiles, Diagnostics& diag) {
  ChemicalSystemState st;
  const auto sp = parse_smiles(st, smiles, diag);
  if (sp < 0) return {};
  return inchi(st, static_cast<std::uint32_t>(sp), &diag);
}

}  // namespace chemlarp
