#include "chemlarp/depict.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <deque>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "chemlarp/cip.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/rings.hpp"
#include "chemlarp/smiles.hpp"

namespace chemlarp {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

struct V2 {
  double x = 0.0, y = 0.0;
};
inline V2 operator+(V2 a, V2 b) { return {a.x + b.x, a.y + b.y}; }
inline V2 operator-(V2 a, V2 b) { return {a.x - b.x, a.y - b.y}; }
inline V2 operator*(V2 a, double k) { return {a.x * k, a.y * k}; }
inline double dot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
inline double cross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
inline double norm(V2 a) { return std::sqrt(dot(a, a)); }
inline V2 unit(V2 a) {
  const double n = norm(a);
  return n > 1e-12 ? a * (1.0 / n) : V2{1.0, 0.0};
}
inline V2 rot(V2 v, double a) {
  const double c = std::cos(a), s = std::sin(a);
  return {v.x * c - v.y * s, v.x * s + v.y * c};
}
inline V2 polar(double a) { return {std::cos(a), std::sin(a)}; }
inline double angle_of(V2 v) { return std::atan2(v.y, v.x); }
inline V2 perp(V2 v) { return {-v.y, v.x}; }
double wrap(double a) {
  while (a < 0.0) a += 2.0 * kPi;
  while (a >= 2.0 * kPi) a -= 2.0 * kPi;
  return a;
}

// ---- molecular graph of one species (local atom ids) ---------------------------------------

struct Graph {
  std::uint32_t begin = 0, n = 0;
  std::vector<std::vector<std::uint32_t>> nb;    // neighbours
  std::vector<std::vector<std::uint32_t>> bond;  // global bond ids, parallel to nb
  std::vector<std::vector<std::uint32_t>> rings; // atoms in ring order
  std::vector<std::vector<int>> atom_rings;
  std::vector<int> system;  // ring system of each atom, -1 when acyclic
  std::vector<std::vector<int>> system_rings;
  std::vector<std::vector<std::uint32_t>> system_atoms;
  std::vector<int> component;
  int components = 0;
  std::vector<char> linear;  // sp centre: drawn straight

  bool bonded(std::uint32_t u, std::uint32_t v) const {
    return std::find(nb[u].begin(), nb[u].end(), v) != nb[u].end();
  }
  // n1-a-n2 is an angle inside one ring
  bool ring_interior(std::uint32_t a, std::uint32_t n1, std::uint32_t n2) const {
    for (int r : atom_rings[a]) {
      const auto& R = rings[static_cast<std::size_t>(r)];
      const std::size_t m = R.size();
      for (std::size_t i = 0; i < m; ++i) {
        if (R[i] != a) continue;
        const std::uint32_t p = R[(i + m - 1) % m], q = R[(i + 1) % m];
        if ((p == n1 && q == n2) || (p == n2 && q == n1)) return true;
      }
    }
    return false;
  }
};

Graph build_graph(const ChemicalSystemState& st, std::uint32_t sp) {
  const Species& s = st.species[sp];
  Graph g;
  g.begin = s.atom_begin;
  g.n = s.atom_count;
  g.nb.resize(g.n);
  g.bond.resize(g.n);
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    const std::uint32_t u = st.bonds.begin[b] - g.begin, v = st.bonds.end[b] - g.begin;
    g.nb[u].push_back(v);
    g.bond[u].push_back(b);
    g.nb[v].push_back(u);
    g.bond[v].push_back(b);
  }
  g.atom_rings.resize(g.n);
  for (const Ring& r : find_sssr(st, sp)) {
    std::vector<std::uint32_t> loc;
    for (auto a : r.atoms) loc.push_back(a - g.begin);
    for (auto a : loc) g.atom_rings[a].push_back(static_cast<int>(g.rings.size()));
    g.rings.push_back(std::move(loc));
  }
  // ring systems: rings sharing at least one atom
  std::vector<int> parent(g.rings.size());
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](int x) {
    while (parent[static_cast<std::size_t>(x)] != x) x = parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
    return x;
  };
  for (std::uint32_t a = 0; a < g.n; ++a)
    for (std::size_t k = 1; k < g.atom_rings[a].size(); ++k) parent[static_cast<std::size_t>(find(g.atom_rings[a][k]))] = find(g.atom_rings[a][0]);
  std::vector<int> sys_of_root(g.rings.size(), -1);
  g.system.assign(g.n, -1);
  for (std::size_t r = 0; r < g.rings.size(); ++r) {
    const int root = find(static_cast<int>(r));
    if (sys_of_root[static_cast<std::size_t>(root)] < 0) {
      sys_of_root[static_cast<std::size_t>(root)] = static_cast<int>(g.system_rings.size());
      g.system_rings.emplace_back();
      g.system_atoms.emplace_back();
    }
    const int sys = sys_of_root[static_cast<std::size_t>(root)];
    g.system_rings[static_cast<std::size_t>(sys)].push_back(static_cast<int>(r));
    for (auto a : g.rings[r])
      if (g.system[a] < 0) {
        g.system[a] = sys;
        g.system_atoms[static_cast<std::size_t>(sys)].push_back(a);
      }
  }
  g.component.assign(g.n, -1);
  for (std::uint32_t s0 = 0; s0 < g.n; ++s0) {
    if (g.component[s0] >= 0) continue;
    std::deque<std::uint32_t> q{s0};
    g.component[s0] = g.components;
    while (!q.empty()) {
      const auto a = q.front();
      q.pop_front();
      for (auto w : g.nb[a])
        if (g.component[w] < 0) {
          g.component[w] = g.components;
          q.push_back(w);
        }
    }
    ++g.components;
  }
  g.linear.assign(g.n, 0);
  for (std::uint32_t a = 0; a < g.n; ++a) {
    if (g.nb[a].size() != 2) continue;
    int doubles = 0, triples = 0;
    for (auto b : g.bond[a]) {
      if (st.bonds.order[b] == 2) ++doubles;
      if (st.bonds.order[b] == 3) ++triples;
    }
    g.linear[a] = (triples > 0 || doubles == 2) ? 1 : 0;
  }
  return g;
}

// Direction into the largest angular gap around `a` among the neighbours for which
// `use` is true. Gaps inside a ring are avoided when another gap exists.
struct Gap {
  double start = 0.0, size = 2.0 * kPi;
  bool any = false;
};
Gap largest_gap(const Graph& g, std::uint32_t a, const std::vector<V2>& p, const std::vector<char>& use,
                bool avoid_ring_interior = true, const std::vector<V2>* extra = nullptr) {
  struct Dir {
    double ang;
    std::int64_t atom;
  };
  std::vector<Dir> d;
  for (auto w : g.nb[a])
    if (use[w]) d.push_back({wrap(angle_of(p[w] - p[a])), w});
  if (extra)
    for (V2 e : *extra) d.push_back({wrap(angle_of(e)), -1});
  Gap best;
  if (d.empty()) return best;
  best.any = true;
  if (d.size() == 1) {
    best.start = d[0].ang;
    best.size = 2.0 * kPi;
    return best;
  }
  std::sort(d.begin(), d.end(), [](const Dir& x, const Dir& y) { return x.ang < y.ang; });
  double best_open = -1.0, best_any = -1.0;
  Gap open, anyg;
  for (std::size_t i = 0; i < d.size(); ++i) {
    const Dir& x = d[i];
    const Dir& y = d[(i + 1) % d.size()];
    const double size = i + 1 < d.size() ? y.ang - x.ang : y.ang + 2.0 * kPi - x.ang;
    const bool interior = avoid_ring_interior && x.atom >= 0 && y.atom >= 0 &&
                          g.ring_interior(a, static_cast<std::uint32_t>(x.atom), static_cast<std::uint32_t>(y.atom));
    if (!interior && size > best_open + 1e-9) {
      best_open = size;
      open = {x.ang, size, true};
    }
    if (size > best_any + 1e-9) {
      best_any = size;
      anyg = {x.ang, size, true};
    }
  }
  return best_open > 0.0 ? open : anyg;
}

double crowding(const std::vector<V2>& pts, const std::vector<V2>& pos, const std::vector<char>& placed,
                std::uint32_t skip = UINT32_MAX) {
  double e = 0.0;
  for (V2 p : pts)
    for (std::size_t q = 0; q < pos.size(); ++q) {
      if (!placed[q] || q == skip) continue;
      const V2 d = p - pos[q];
      const double r2 = dot(d, d);
      if (r2 < 3.24) e += 1.0 / (r2 + 0.02);
    }
  return e;
}

// ---- ring systems in their own frame -----------------------------------------------------------

void close_arc(std::uint32_t A, std::uint32_t B, const std::vector<std::uint32_t>& run, std::vector<V2>& pos,
               std::vector<char>& placed) {
  const V2 pa = pos[A], pb = pos[B];
  const double d = norm(pb - pa);
  const double k = static_cast<double>(run.size());
  double theta = 0.0;
  if (d < k + 1.0 - 1e-9) {
    double lo = 1e-9, hi = 2.0 * kPi / (k + 1.0);
    for (int it = 0; it < 100; ++it) {
      const double mid = 0.5 * (lo + hi);
      const double r = std::sin((k + 1.0) * mid / 2.0) / std::sin(mid / 2.0);
      (r > d ? lo : hi) = mid;
    }
    theta = 0.5 * (lo + hi);
  }
  const V2 c = d > 1e-9 ? unit(pb - pa) : V2{1.0, 0.0};
  std::vector<V2> best;
  double best_e = std::numeric_limits<double>::infinity();
  for (int s : {1, -1}) {
    std::vector<V2> pts;
    V2 p = pa;
    for (std::size_t i = 0; i < run.size(); ++i) {
      p = p + rot(c, s * (k * theta / 2.0 - static_cast<double>(i) * theta));
      pts.push_back(p);
    }
    const double e = crowding(pts, pos, placed);
    if (e < best_e - 1e-9) {
      best_e = e;
      best = pts;
    }
  }
  for (std::size_t i = 0; i < run.size(); ++i) {
    pos[run[i]] = best[i];
    placed[run[i]] = 1;
  }
}

void close_ring(const Graph& g, const std::vector<std::uint32_t>& R, std::vector<V2>& pos, std::vector<char>& placed) {
  const std::size_t n = R.size();
  std::size_t m = 0;
  for (auto a : R) m += placed[a] ? 1 : 0;
  if (m == n) return;
  if (m == 1) {  // spiro: a whole polygon hanging off one atom
    std::size_t ia = 0;
    while (!placed[R[ia]]) ++ia;
    const std::uint32_t a = R[ia];
    const Gap gp = largest_gap(g, a, pos, placed);
    const V2 dir = gp.any ? polar(gp.start + gp.size / 2.0) : V2{1.0, 0.0};
    const double rad = 1.0 / (2.0 * std::sin(kPi / static_cast<double>(n)));
    const V2 c = pos[a] + dir * rad;
    const double a0 = angle_of(pos[a] - c);
    for (std::size_t k = 1; k < n; ++k) {
      pos[R[(ia + k) % n]] = c + polar(a0 + 2.0 * kPi * static_cast<double>(k) / static_cast<double>(n)) * rad;
      placed[R[(ia + k) % n]] = 1;
    }
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!placed[R[i]] || placed[R[(i + 1) % n]]) continue;
    std::vector<std::uint32_t> run;
    std::size_t j = (i + 1) % n;
    while (!placed[R[j]]) {
      run.push_back(R[j]);
      j = (j + 1) % n;
    }
    close_arc(R[i], R[j], run, pos, placed);
  }
}

void layout_system(const Graph& g, int sys, std::vector<V2>& pos, std::vector<char>& placed) {
  const auto& rs = g.system_rings[static_cast<std::size_t>(sys)];
  auto shared = [&](int r1, int r2) {
    int c = 0;
    for (auto a : g.rings[static_cast<std::size_t>(r1)])
      for (auto b : g.rings[static_cast<std::size_t>(r2)]) c += a == b;
    return c;
  };
  int first = rs[0], first_score = -1;
  for (int r : rs) {
    int fused = 0;
    for (int q : rs)
      if (q != r && shared(r, q) >= 2) ++fused;
    const int score = fused * 100 + static_cast<int>(g.rings[static_cast<std::size_t>(r)].size() == 6 ? 10 : 0);
    if (score > first_score) {
      first_score = score;
      first = r;
    }
  }
  const auto& R0 = g.rings[static_cast<std::size_t>(first)];
  pos[R0[0]] = {0.0, 0.0};
  pos[R0[1]] = {0.0, 1.0};
  placed[R0[0]] = placed[R0[1]] = 1;
  close_ring(g, R0, pos, placed);
  std::vector<char> done(g.rings.size(), 0);
  done[static_cast<std::size_t>(first)] = 1;
  for (;;) {
    int next = -1, best = 0;
    for (int r : rs) {
      if (done[static_cast<std::size_t>(r)]) continue;
      int m = 0;
      for (auto a : g.rings[static_cast<std::size_t>(r)]) m += placed[a] ? 1 : 0;
      if (m > best) {
        best = m;
        next = r;
      }
    }
    if (next < 0) break;
    done[static_cast<std::size_t>(next)] = 1;
    close_ring(g, g.rings[static_cast<std::size_t>(next)], pos, placed);
  }
}

// ---- whole-species layout -----------------------------------------------------------------------

struct DoubleStereoLocal {
  std::uint32_t a, b, ra, rb;
  bool cis;
};

bool same_side(V2 a, V2 b, V2 ra, V2 rb) {
  const V2 ax = b - a;
  return cross(ax, ra - a) * cross(ax, rb - a) > 0.0;
}

class Layout {
 public:
  Layout(const ChemicalSystemState& st, std::uint32_t sp, const Graph& g) : st_(st), g_(g) {
    pos.assign(g.n, {});
    placed_.assign(g.n, 0);
    turn_.assign(g.n, 0);
    for (const auto& d : st.double_bonds) {
      if (st.bonds.species[d.bond] != sp) continue;
      DoubleStereoLocal l{st.bonds.begin[d.bond] - g.begin, st.bonds.end[d.bond] - g.begin, d.ref_begin - g.begin,
                          d.ref_end - g.begin, d.cis};
      ez_.push_back(l);
    }
    // ring systems in their own frames
    local_.assign(g.n, {});
    std::vector<char> lp(g.n, 0);
    for (std::size_t s = 0; s < g.system_rings.size(); ++s) layout_system(g, static_cast<int>(s), local_, lp);
  }

  void run() {
    for (int c = 0; c < g_.components; ++c) grow_component(c);
    relieve();
    orient_and_arrange();
  }

  std::vector<V2> pos;

 private:
  const ChemicalSystemState& st_;
  const Graph& g_;
  std::vector<char> placed_;
  std::vector<int> turn_;
  std::vector<V2> local_;
  std::vector<DoubleStereoLocal> ez_;
  std::deque<std::uint32_t> queue_;
  std::vector<std::uint32_t> fresh_;

  void put(std::uint32_t a, V2 p) {
    pos[a] = p;
    placed_[a] = 1;
    fresh_.push_back(a);
  }

  // Places ring system `sys` so that atom v sits at p with its exocyclic bond toward `toward`.
  void place_system(int sys, std::uint32_t v, V2 p, const V2* toward) {
    const auto& atoms = g_.system_atoms[static_cast<std::size_t>(sys)];
    std::vector<char> in(g_.n, 0);
    for (auto a : atoms) in[a] = 1;
    if (!toward) {
      V2 c{};
      for (auto a : atoms) c = c + local_[a];
      c = c * (1.0 / static_cast<double>(atoms.size()));
      for (auto a : atoms) put(a, local_[a] - c + p);
      return;
    }
    const Gap gp = largest_gap(g_, v, local_, in);
    const V2 e = polar(gp.start + gp.size / 2.0);
    const V2 t = unit(*toward - p);
    std::vector<V2> best;
    double best_e = std::numeric_limits<double>::infinity();
    for (int mirror = 0; mirror < 2; ++mirror) {
      auto mir = [&](V2 x) { return mirror ? V2{x.x, -x.y} : x; };
      const double phi = angle_of(t) - angle_of(mir(e));
      std::vector<V2> pts;
      for (auto a : atoms) pts.push_back(rot(mir(local_[a] - local_[v]), phi) + p);
      const double en = crowding(pts, pos, placed_);
      if (en < best_e - 1e-9) {
        best_e = en;
        best = pts;
      }
    }
    for (std::size_t i = 0; i < atoms.size(); ++i) put(atoms[i], best[i]);
  }

  std::size_t subtree(std::uint32_t u, std::uint32_t from) const {
    std::vector<char> seen(g_.n, 0);
    seen[from] = 1;
    std::vector<std::uint32_t> stack{u};
    seen[u] = 1;
    std::size_t c = 0;
    while (!stack.empty()) {
      const auto a = stack.back();
      stack.pop_back();
      ++c;
      for (auto w : g_.nb[a])
        if (!seen[w] && !placed_[w]) {
          seen[w] = 1;
          stack.push_back(w);
        }
    }
    return c;
  }

  void grow_component(int comp) {
    std::vector<std::uint32_t> atoms;
    for (std::uint32_t a = 0; a < g_.n; ++a)
      if (g_.component[a] == comp) atoms.push_back(a);
    int big = -1;
    std::size_t big_n = 0;
    for (auto a : atoms)
      if (g_.system[a] >= 0 && g_.system_atoms[static_cast<std::size_t>(g_.system[a])].size() > big_n) {
        big = g_.system[a];
        big_n = g_.system_atoms[static_cast<std::size_t>(big)].size();
      }
    fresh_.clear();
    if (big >= 0) {
      place_system(big, g_.system_atoms[static_cast<std::size_t>(big)][0], {0.0, 0.0}, nullptr);
    } else {
      // start from one end of the longest chain
      auto farthest = [&](std::uint32_t s) {
        std::vector<int> dist(g_.n, -1);
        std::deque<std::uint32_t> q{s};
        dist[s] = 0;
        std::uint32_t last = s;
        while (!q.empty()) {
          const auto a = q.front();
          q.pop_front();
          last = a;
          for (auto w : g_.nb[a])
            if (dist[w] < 0) {
              dist[w] = dist[a] + 1;
              q.push_back(w);
            }
        }
        return last;
      };
      put(farthest(atoms[0]), {0.0, 0.0});
    }
    for (auto a : fresh_) queue_.push_back(a);
    while (!queue_.empty()) {
      const auto a = queue_.front();
      queue_.pop_front();
      expand(a);
    }
  }

  void expand(std::uint32_t a) {
    std::vector<std::uint32_t> P, U;
    for (auto w : g_.nb[a]) (placed_[w] ? P : U).push_back(w);
    if (U.empty()) return;
    std::stable_sort(U.begin(), U.end(), [&](std::uint32_t x, std::uint32_t y) { return subtree(x, a) > subtree(y, a); });
    std::vector<double> dirs;
    const V2 pa = pos[a];
    const std::size_t k = U.size();
    if (P.empty()) {
      if (k == 1) dirs = {-30 * kDeg};
      else if (k == 2) dirs = {-30 * kDeg, 210 * kDeg};
      else
        for (std::size_t i = 0; i < k; ++i) dirs.push_back(90 * kDeg + 2.0 * kPi * static_cast<double>(i) / static_cast<double>(k));
      turn_[a] = 1;
    } else if (P.size() == 1) {
      const std::uint32_t p = P[0];
      const double back = angle_of(pos[p] - pa), h = back + kPi;
      if (k == 1 && g_.linear[a]) {
        dirs = {h};
        turn_[a] = turn_[p];
      } else if (k <= 2) {
        int t = turn_[p] != 0 ? -turn_[p] : 0;
        if (t == 0) {
          const double ep = crowding({pa + polar(h + 60 * kDeg)}, pos, placed_, a);
          const double em = crowding({pa + polar(h - 60 * kDeg)}, pos, placed_, a);
          t = em < ep - 1e-9 ? -1 : 1;
        }
        dirs = {h + t * 60 * kDeg};
        if (k == 2) dirs.push_back(h - t * 60 * kDeg);
        turn_[a] = t;
      } else if (k == 3) {
        dirs = {back + kPi, back + 90 * kDeg, back + 270 * kDeg};
      } else {
        for (std::size_t i = 0; i < k; ++i) dirs.push_back(back + 2.0 * kPi * static_cast<double>(i + 1) / static_cast<double>(k + 1));
      }
    } else {
      const Gap gp = largest_gap(g_, a, pos, placed_);
      for (std::size_t i = 0; i < k; ++i)
        dirs.push_back(gp.start + gp.size * static_cast<double>(i + 1) / static_cast<double>(k + 1));
      // put the biggest branch in the middle of the gap
      if (k == 3) std::swap(dirs[0], dirs[1]);
    }
    fresh_.clear();
    for (std::size_t i = 0; i < k; ++i) {
      const std::uint32_t u = U[i];
      if (placed_[u]) continue;  // reached through a ring system placed a moment ago
      const V2 target = pa + polar(dirs[i]);
      if (g_.system[u] >= 0) place_system(g_.system[u], u, target, &pos[a]);
      else put(u, target);
    }
    enforce_ez(a);
    for (auto w : fresh_) queue_.push_back(w);
  }

  void enforce_ez(std::uint32_t a) {
    for (const auto& d : ez_) {
      if (d.a != a && d.b != a) continue;
      const std::uint32_t o = d.a == a ? d.b : d.a;
      if (!placed_[o] || !placed_[d.ra] || !placed_[d.rb] || !placed_[d.a] || !placed_[d.b]) continue;
      if (same_side(pos[d.a], pos[d.b], pos[d.ra], pos[d.rb]) == d.cis) continue;
      if (fresh_.empty()) continue;
      const V2 p0 = pos[a], dir = unit(pos[o] - pos[a]);
      for (auto w : fresh_) {
        const V2 r = pos[w] - p0;
        pos[w] = p0 + dir * (2.0 * dot(r, dir)) - r;
      }
      turn_[a] = -turn_[a];
    }
  }

  // ---- clean-up: reflect subtrees across acyclic single bonds ----
  static bool segments_cross(V2 p1, V2 p2, V2 q1, V2 q2) {
    const double d1 = cross(p2 - p1, q1 - p1), d2 = cross(p2 - p1, q2 - p1);
    const double d3 = cross(q2 - q1, p1 - q1), d4 = cross(q2 - q1, p2 - q1);
    return ((d1 > 1e-9 && d2 < -1e-9) || (d1 < -1e-9 && d2 > 1e-9)) && ((d3 > 1e-9 && d4 < -1e-9) || (d3 < -1e-9 && d4 > 1e-9));
  }

  double clash(const std::vector<char>& side, int comp) const {
    double e = 0.0;
    for (std::uint32_t i = 0; i < g_.n; ++i) {
      if (!side[i]) continue;
      for (std::uint32_t j = 0; j < g_.n; ++j) {
        if (side[j] || g_.component[j] != comp || g_.bonded(i, j)) continue;
        const V2 d = pos[i] - pos[j];
        const double r2 = dot(d, d);
        if (r2 < 1.96) e += 1.0 / (r2 + 0.01) - 1.0 / 1.97;
      }
    }
    for (std::uint32_t i = 0; i < g_.n; ++i) {
      if (!side[i]) continue;
      for (auto w : g_.nb[i]) {
        if (!side[w] || w < i) continue;
        for (std::uint32_t x = 0; x < g_.n; ++x) {
          if (side[x] || g_.component[x] != comp) continue;
          for (auto y : g_.nb[x]) {
            if (side[y] || y < x) continue;
            if (segments_cross(pos[i], pos[w], pos[x], pos[y])) e += 3.0;
          }
        }
      }
    }
    return e;
  }

  void relieve() {
    const Species& s = st_.species[st_.atoms.species[g_.begin]];
    for (int pass = 0; pass < 6; ++pass) {
      bool improved = false;
      for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
        if (st_.bonds.order[b] != 1 || (st_.bonds.flags[b] & kBondInRing)) continue;
        const std::uint32_t u = st_.bonds.begin[b] - g_.begin, v = st_.bonds.end[b] - g_.begin;
        if (g_.nb[u].size() < 2 || g_.nb[v].size() < 2) continue;
        std::vector<char> side(g_.n, 0);
        std::vector<std::uint32_t> stack{v};
        side[v] = 1;
        bool ring = false;
        std::size_t count = 0;
        while (!stack.empty()) {
          const auto a = stack.back();
          stack.pop_back();
          ++count;
          for (auto w : g_.nb[a]) {
            if (a == v && w == u) continue;
            if (w == u) ring = true;
            if (!side[w]) {
              side[w] = 1;
              stack.push_back(w);
            }
          }
        }
        if (ring) continue;
        std::size_t comp_n = 0;
        for (std::uint32_t i = 0; i < g_.n; ++i) comp_n += g_.component[i] == g_.component[u];
        if (2 * count > comp_n)
          for (std::uint32_t i = 0; i < g_.n; ++i)
            if (g_.component[i] == g_.component[u]) side[i] = !side[i];
        const double before = clash(side, g_.component[u]);
        if (before <= 1e-9) continue;
        const V2 p0 = pos[u], dir = unit(pos[v] - pos[u]);
        const std::vector<V2> saved = pos;
        for (std::uint32_t i = 0; i < g_.n; ++i)
          if (side[i]) {
            const V2 r = pos[i] - p0;
            pos[i] = p0 + dir * (2.0 * dot(r, dir)) - r;
          }
        const double after = clash(side, g_.component[u]);
        if (after < before - 1e-6) improved = true;
        else pos = saved;
      }
      if (!improved) break;
    }
  }

  void orient_and_arrange() {
    double cursor = 0.0;
    for (int c = 0; c < g_.components; ++c) {
      std::vector<std::uint32_t> atoms;
      for (std::uint32_t a = 0; a < g_.n; ++a)
        if (g_.component[a] == c) atoms.push_back(a);
      V2 m{};
      for (auto a : atoms) m = m + pos[a];
      m = m * (1.0 / static_cast<double>(atoms.size()));
      double sxx = 0, syy = 0, sxy = 0;
      for (auto a : atoms) {
        const V2 d = pos[a] - m;
        sxx += d.x * d.x;
        syy += d.y * d.y;
        sxy += d.x * d.y;
      }
      double turn = 0.0;
      bool mirror = false;
      if (atoms.size() >= 3 && (std::fabs(sxx - syy) > 1e-6 || std::fabs(sxy) > 1e-6)) {
        const double axis = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
        turn = -std::round(axis / (60 * kDeg)) * 60 * kDeg;  // multiples of 60° keep the hexagonal lattice
        mirror = axis + turn < -5 * kDeg;  // draw tilted skeletons rising to the right, as steroids are drawn
      }
      double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300;
      for (auto a : atoms) {
        pos[a] = rot(pos[a] - m, turn);
        if (mirror) pos[a].x = -pos[a].x;
        minx = std::min(minx, pos[a].x);
        maxx = std::max(maxx, pos[a].x);
        miny = std::min(miny, pos[a].y);
        maxy = std::max(maxy, pos[a].y);
      }
      const V2 shift{cursor - minx, -(miny + maxy) / 2.0};
      for (auto a : atoms) pos[a] = pos[a] + shift;
      cursor += (maxx - minx) + 1.6;
    }
  }
};

// ---- labels and drawing ------------------------------------------------------------------------

double glyph_width(char c) {
  static const double upper[26] = {.667, .667, .722, .722, .667, .611, .778, .722, .278, .5,   .667, .556, .833,
                                   .722, .778, .667, .778, .722, .667, .611, .722, .667, .944, .667, .667, .611};
  static const double lower[26] = {.556, .556, .5,   .556, .556, .278, .556, .556, .222, .222, .5, .222, .833,
                                   .556, .556, .556, .556, .333, .5,   .278, .556, .5,   .722, .5, .5, .5};
  if (c >= 'A' && c <= 'Z') return upper[c - 'A'];
  if (c >= 'a' && c <= 'z') return lower[c - 'a'];
  if (c >= '0' && c <= '9') return 0.556;
  if (c == '+' || c == '-') return 0.584;
  if (c == '(' || c == ')') return 0.333;
  return 0.6;
}

struct Run {
  std::string text;
  int kind = 0;  // 0 normal, 1 subscript, 2 superscript
};
double run_width(const Run& r, double fs) {
  double w = 0.0;
  for (char c : r.text) w += glyph_width(c);
  return w * fs * (r.kind ? 0.7 : 1.0);
}

struct Box {
  double x0, y0, x1, y1;
  bool inside(V2 p) const { return p.x > x0 && p.x < x1 && p.y > y0 && p.y < y1; }
};

struct TextLine {
  std::vector<Run> runs;
  double x = 0.0, y = 0.0;  // start of text, baseline
};

struct Label {
  bool present = false;
  std::string color = "currentColor";
  std::vector<TextLine> lines;
  std::vector<Box> boxes;
  std::string plain;
};

const char* element_color(std::uint8_t z) {
  switch (z) {
    case 7: return "#2f5fd0";
    case 8: return "#e0332b";
    case 9: case 17: return "#2a9d46";
    case 15: return "#e07b00";
    case 16: return "#c39a00";
    case 35: return "#a0302a";
    case 53: return "#8e2ba0";
    case 5: return "#d9822b";
    default: return "currentColor";
  }
}

std::string charge_text(int q) {
  if (q == 0) return "";
  std::string s = std::abs(q) > 1 ? std::to_string(std::abs(q)) : "";
  return s + (q > 0 ? "+" : "-");
}

std::string num(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%.2f", v);
  std::string s = b;
  if (s == "-0.00") s = "0.00";
  return s;
}

std::string escape_minus(const std::string& t) {
  std::string o;
  for (char c : t) o += c == '-' ? std::string("&#8722;") : std::string(1, c);
  return o;
}

std::string emit_line(const TextLine& L, double fs, const std::string& fill, const char* extra = "") {
  std::string o = "<text x=\"" + num(L.x) + "\" y=\"" + num(L.y) + "\" font-size=\"" + num(fs) + "\" fill=\"" + fill + "\"" + extra + ">";
  double shift = 0.0;
  for (const Run& r : L.runs) {
    const double want = r.kind == 1 ? 0.28 * fs : r.kind == 2 ? -0.42 * fs : 0.0;
    o += "<tspan";
    if (std::fabs(want - shift) > 1e-9) o += " dy=\"" + num(want - shift) + "\"";
    if (r.kind) o += " font-size=\"" + num(0.7 * fs) + "\"";
    o += ">" + escape_minus(r.text) + "</tspan>";
    shift = want;
  }
  return o + "</text>";
}

V2 clip_to_boxes(V2 p, V2 q, const std::vector<Box>& boxes) {
  for (int pass = 0; pass < 2; ++pass)
    for (const Box& b : boxes) {
      if (!b.inside(p)) continue;
      const V2 d = q - p;
      double t = 1e300;
      if (d.x > 1e-12) t = std::min(t, (b.x1 - p.x) / d.x);
      if (d.x < -1e-12) t = std::min(t, (b.x0 - p.x) / d.x);
      if (d.y > 1e-12) t = std::min(t, (b.y1 - p.y) / d.y);
      if (d.y < -1e-12) t = std::min(t, (b.y0 - p.y) / d.y);
      if (t > 1.0) t = 1.0;
      p = p + d * t;
    }
  return p;
}

}  // namespace

std::vector<std::array<double, 2>> layout_2d(const ChemicalSystemState& st, std::uint32_t sp) {
  const Graph g = build_graph(st, sp);
  Layout lay(st, sp, g);
  lay.run();
  std::vector<std::array<double, 2>> out;
  for (const V2& p : lay.pos) out.push_back({p.x, p.y});
  return out;
}

Depiction depict(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag, const DepictOptions& opt) {
  Depiction D;
  const Species& s = st.species[sp];
  const Graph g = build_graph(st, sp);
  Layout lay(st, sp, g);
  lay.run();
  std::vector<V2> P = lay.pos;  // math coordinates, unit bonds, y up

  // drawn graph: real atoms, then hydrogens added to carry wedges
  struct DBond {
    std::uint32_t a, b;
    std::uint8_t order;
    BondStereo stereo = BondStereo::None;
    std::int64_t global = -1;
  };
  std::vector<DBond> bonds;
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b)
    bonds.push_back({st.bonds.begin[b] - g.begin, st.bonds.end[b] - g.begin, st.bonds.order[b], BondStereo::None, b});
  std::vector<std::int32_t> atom_id;
  for (std::uint32_t i = 0; i < g.n; ++i) atom_id.push_back(static_cast<std::int32_t>(g.begin + i));
  std::vector<std::vector<std::uint32_t>> adj(g.n);
  auto rebuild_adj = [&] {
    adj.assign(atom_id.size(), {});
    for (std::uint32_t k = 0; k < bonds.size(); ++k) {
      adj[bonds[k].a].push_back(k);
      adj[bonds[k].b].push_back(k);
    }
  };
  rebuild_adj();
  auto other = [&](std::uint32_t k, std::uint32_t a) { return bonds[k].a == a ? bonds[k].b : bonds[k].a; };

  std::vector<const TetrahedralStereo*> centres;
  for (const auto& t : st.tetrahedral)
    if (st.atoms.species[t.atom] == sp) centres.push_back(&t);
  std::vector<char> is_centre(g.n, 0);
  for (auto* t : centres) is_centre[t->atom - g.begin] = 1;

  // ring-fusion stereocentres get an explicit hydrogen to carry the wedge
  std::vector<std::int64_t> drawn_h(g.n, -1);
  std::vector<char> all(g.n, 1);
  for (auto* t : centres) {
    const std::uint32_t c = t->atom - g.begin;
    const bool has_h = std::find(t->refs.begin(), t->refs.end(), kImplicitHRef) != t->refs.end();
    if (!has_h || g.nb[c].size() != 3) continue;
    bool all_ring = true;
    for (auto b : g.bond[c]) all_ring = all_ring && (st.bonds.flags[b] & kBondInRing);
    if (!all_ring) continue;
    const Gap gp = largest_gap(g, c, P, all);
    // the middle of the gap unless that crowds another label or wedged hydrogen
    V2 hp{};
    double best = std::numeric_limits<double>::infinity();
    for (double f : {0.5, 0.4, 0.6, 0.3, 0.7}) {
      const V2 cand = P[c] + polar(gp.start + gp.size * f) * 0.75;
      double e = std::fabs(f - 0.5) * 0.5;
      for (std::size_t q = 0; q < P.size(); ++q) {
        if (q == c) continue;
        const double r2 = dot(cand - P[q], cand - P[q]);
        const bool labelled = q >= g.n || st.atoms.element[g.begin + q] != 6;
        if (r2 < 1.44) e += (labelled ? 3.0 : 1.0) / (r2 + 0.05);
      }
      if (e < best) {
        best = e;
        hp = cand;
      }
    }
    drawn_h[c] = static_cast<std::int64_t>(atom_id.size());
    atom_id.push_back(-1);
    P.push_back(hp);
    bonds.push_back({c, static_cast<std::uint32_t>(atom_id.size() - 1), 1, BondStereo::None, -1});
  }
  rebuild_adj();
  const std::size_t N = atom_id.size();

  // ---- wedges ----
  auto local_of = [&](std::uint32_t c, std::uint32_t ref) -> std::int64_t {
    if (ref == kImplicitHRef) return drawn_h[c];
    if (ref == kLonePairRef) return -1;
    return static_cast<std::int64_t>(ref - g.begin);
  };
  auto volume = [&](std::uint32_t c, const std::array<std::int64_t, 4>& refs, std::uint32_t wedge_to, double z) {
    std::array<std::array<double, 3>, 4> v{};
    int missing = -1;
    for (int i = 0; i < 4; ++i) {
      if (refs[static_cast<std::size_t>(i)] < 0) {
        if (missing >= 0) return 0.0;
        missing = i;
        continue;
      }
      const auto r = static_cast<std::uint32_t>(refs[static_cast<std::size_t>(i)]);
      const V2 d = unit(P[r] - P[c]);
      v[static_cast<std::size_t>(i)] = {d.x, d.y, r == wedge_to ? z : 0.0};
    }
    if (missing >= 0) {
      std::array<double, 3> sum{0, 0, 0};
      for (int i = 0; i < 4; ++i)
        if (i != missing)
          for (int k = 0; k < 3; ++k) sum[static_cast<std::size_t>(k)] += v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
      const double n = std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
      if (n < 1e-9) return 0.0;
      v[static_cast<std::size_t>(missing)] = {-sum[0] / n, -sum[1] / n, -sum[2] / n};
    }
    auto sub = [](const std::array<double, 3>& x, const std::array<double, 3>& y) {
      return std::array<double, 3>{x[0] - y[0], x[1] - y[1], x[2] - y[2]};
    };
    const auto b = sub(v[1], v[0]), cc = sub(v[2], v[0]), d = sub(v[3], v[0]);
    const std::array<double, 3> x{cc[1] * d[2] - cc[2] * d[1], cc[2] * d[0] - cc[0] * d[2], cc[0] * d[1] - cc[1] * d[0]};
    return b[0] * x[0] + b[1] * x[1] + b[2] * x[2];
  };
  for (auto* t : centres) {
    ++D.stereo_centres;
    const std::uint32_t c = t->atom - g.begin;
    std::array<std::int64_t, 4> refs{};
    for (int i = 0; i < 4; ++i) refs[static_cast<std::size_t>(i)] = local_of(c, t->refs[static_cast<std::size_t>(i)]);
    std::vector<std::pair<int, std::uint32_t>> cand;
    for (auto k : adj[c]) {
      if (bonds[k].order != 1 || bonds[k].stereo != BondStereo::None) continue;
      const std::uint32_t w = other(k, c);
      int score = 0;
      if (atom_id[w] < 0) score += 16;
      if (bonds[k].global < 0 || !(st.bonds.flags[static_cast<std::size_t>(bonds[k].global)] & kBondInRing)) score += 8;
      if (w < g.n && !is_centre[w]) score += 4;
      if (adj[w].size() == 1) score += 2;
      if (atom_id[w] >= 0 && st.atoms.element[static_cast<std::size_t>(atom_id[w])] == 1) score += 1;
      cand.push_back({score, k});
    }
    std::stable_sort(cand.begin(), cand.end(), [](auto& x, auto& y) { return x.first > y.first; });
    bool ok = false;
    for (auto [score, k] : cand) {
      (void)score;
      const std::uint32_t w = other(k, c);
      const double up = volume(c, refs, w, 0.8), down = volume(c, refs, w, -0.8);
      if (std::fabs(up) < 1e-6 || std::fabs(down) < 1e-6 || (up > 0) == (down > 0)) continue;
      const bool want_negative = t->parity == Chirality::CounterClockwise;
      bonds[k].stereo = ((up < 0) == want_negative) ? BondStereo::Wedge : BondStereo::Hash;
      if (bonds[k].a != c) std::swap(bonds[k].a, bonds[k].b);
      ok = true;
      break;
    }
    if (ok) ++D.stereo_centres_ok;
    else {
      auto& d = diag.warn(DiagCode::StereoLost, "2D drawing cannot show the configuration of this stereocentre");
      d.atom = static_cast<std::int32_t>(t->atom);
      d.species = static_cast<std::int32_t>(sp);
    }
  }
  for (const auto& d : st.double_bonds) {
    if (st.bonds.species[d.bond] != sp) continue;
    ++D.stereo_bonds;
    const V2 a = P[st.bonds.begin[d.bond] - g.begin], b = P[st.bonds.end[d.bond] - g.begin];
    if (same_side(a, b, P[d.ref_begin - g.begin], P[d.ref_end - g.begin]) == d.cis) ++D.stereo_bonds_ok;
    else {
      auto& w = diag.warn(DiagCode::StereoLost, "2D drawing cannot show the configuration of this double bond");
      w.atom = static_cast<std::int32_t>(st.bonds.begin[d.bond]);
      w.species = static_cast<std::int32_t>(sp);
    }
  }

  // ---- px space ----
  const double L = opt.bond_length > 1.0 ? opt.bond_length : 30.0;
  const double fs = 0.62 * L, lw = 0.042 * L, sp2 = 0.18 * L;
  std::vector<V2> X(N);
  for (std::size_t i = 0; i < N; ++i) X[i] = {P[i].x * L, -P[i].y * L};

  std::vector<Label> labels(N);
  for (std::size_t i = 0; i < N; ++i) {
    Label& lab = labels[i];
    std::uint8_t z = 1;
    int q = 0, h = 0, iso = 0;
    if (atom_id[i] >= 0) {
      const auto ga = static_cast<std::size_t>(atom_id[i]);
      z = st.atoms.element[ga];
      q = st.atoms.formal_charge[ga];
      h = st.atoms.implicit_h[ga] - (i < g.n && drawn_h[i] >= 0 ? 1 : 0);
      iso = st.atoms.isotope[ga];
    }
    const std::size_t degree = adj[i].size();
    if (z == 6 && q == 0 && iso == 0 && degree > 0) continue;
    lab.present = true;
    if (opt.color) lab.color = element_color(z);
    const std::string sym(element_symbol(z));
    // which side may hold the hydrogens
    double dx_sum = 0.0;
    bool left_blocked = false, right_blocked = false, up_blocked = false, down_blocked = false;
    for (auto k : adj[i]) {
      const V2 d = unit(X[other(k, static_cast<std::uint32_t>(i))] - X[i]);
      dx_sum += d.x;
      if (d.x > 0.5) right_blocked = true;
      if (d.x < -0.5) left_blocked = true;
      if (d.y > 0.5) down_blocked = true;
      if (d.y < -0.5) up_blocked = true;
    }
    enum { Right, Left, Below, Above } hside = Right;
    if (degree == 0) {
      hside = (z == 8 || z == 16 || z == 34 || z == 9 || z == 17 || z == 35 || z == 53) ? Left : Right;
    } else if (h > 0) {
      if (dx_sum > 0.05) hside = Left;
      if (hside == Right && right_blocked) hside = left_blocked ? (down_blocked ? (up_blocked ? Right : Above) : Below) : Left;
      else if (hside == Left && left_blocked) hside = right_blocked ? (down_blocked ? (up_blocked ? Left : Above) : Below) : Right;
    }
    std::vector<Run> pre, post, hruns;
    if (h > 0) {
      hruns.push_back({"H", 0});
      if (h > 1) hruns.push_back({std::to_string(h), 1});
    }
    if (iso) pre.push_back({std::to_string(iso), 2});
    const std::string ch = charge_text(q);
    if (hside == Left) pre.insert(pre.begin(), hruns.begin(), hruns.end());
    if (hside == Right) post.insert(post.end(), hruns.begin(), hruns.end());
    if (!ch.empty()) post.push_back({ch, 2});
    TextLine line;
    double wpre = 0.0;
    for (auto& r : pre) wpre += run_width(r, fs);
    const Run main{sym, 0};
    const double wmain = run_width(main, fs);
    double total = wpre + wmain;
    for (auto& r : post) total += run_width(r, fs);
    line.runs = pre;
    line.runs.push_back(main);
    line.runs.insert(line.runs.end(), post.begin(), post.end());
    line.x = X[i].x - wpre - wmain / 2.0;
    line.y = X[i].y + 0.36 * fs;
    const double m = 0.12 * fs;
    lab.boxes.push_back({line.x - m, X[i].y - 0.5 * fs - m, line.x + total + m, X[i].y + 0.5 * fs + m});
    lab.lines.push_back(line);
    if ((hside == Below || hside == Above) && !hruns.empty()) {
      TextLine hl;
      hl.runs = hruns;
      double wh = 0.0;
      for (auto& r : hruns) wh += run_width(r, fs);
      const double cy = X[i].y + (hside == Below ? 1.0 : -1.0) * fs;
      hl.x = X[i].x - glyph_width('H') * fs / 2.0;
      hl.y = cy + 0.36 * fs;
      lab.boxes.push_back({hl.x - m, cy - 0.5 * fs - m, hl.x + wh + m, cy + 0.5 * fs + m});
      lab.lines.push_back(hl);
    }
    for (auto& l : lab.lines)
      for (auto& r : l.runs) lab.plain += r.text;
  }

  // ring centres for placing inner double-bond lines
  auto bond_between = [&](std::uint32_t u, std::uint32_t v) -> std::int64_t {
    for (std::size_t k = 0; k < g.nb[u].size(); ++k)
      if (g.nb[u][k] == v) return g.bond[u][k];
    return -1;
  };
  struct RingInfo {
    V2 centre;
    int doubles;
    std::size_t size;
  };
  std::vector<RingInfo> rinfo;
  for (const auto& R : g.rings) {
    V2 c{};
    int dbl = 0;
    for (std::size_t i = 0; i < R.size(); ++i) {
      c = c + X[R[i]];
      const auto b = bond_between(R[i], R[(i + 1) % R.size()]);
      if (b >= 0 && st.bonds.order[static_cast<std::size_t>(b)] == 2) ++dbl;
    }
    rinfo.push_back({c * (1.0 / static_cast<double>(R.size())), dbl, R.size()});
  }
  auto ring_for_bond = [&](std::uint32_t u, std::uint32_t v) -> int {
    int best = -1;
    for (int r : g.atom_rings[u]) {
      const auto& R = g.rings[static_cast<std::size_t>(r)];
      bool adjacent = false;
      for (std::size_t i = 0; i < R.size(); ++i) {
        const auto a = R[i], b = R[(i + 1) % R.size()];
        if ((a == u && b == v) || (a == v && b == u)) adjacent = true;
      }
      if (!adjacent) continue;
      if (best < 0) {
        best = r;
        continue;
      }
      const auto& x = rinfo[static_cast<std::size_t>(r)];
      const auto& y = rinfo[static_cast<std::size_t>(best)];
      if (x.doubles > y.doubles || (x.doubles == y.doubles && x.size < y.size)) best = r;
    }
    return best;
  };

  std::string lines_svg, fills_svg, text_svg;
  auto line_el = [&](V2 a, V2 b) {
    if (norm(b - a) < 0.05 * L) return;
    lines_svg += "<line x1=\"" + num(a.x) + "\" y1=\"" + num(a.y) + "\" x2=\"" + num(b.x) + "\" y2=\"" + num(b.y) + "\"/>";
  };
  auto has_label = [&](std::uint32_t i) { return labels[i].present; };
  auto clip_pair = [&](std::uint32_t i, std::uint32_t j, V2 a, V2 b) {
    const V2 ca = has_label(i) ? clip_to_boxes(a, b, labels[i].boxes) : a;
    const V2 cb = has_label(j) ? clip_to_boxes(b, a, labels[j].boxes) : b;
    return std::make_pair(ca, cb);
  };
  // side (+1 left of a->b, -1 right) with more substituents
  auto busy_side = [&](std::uint32_t i, std::uint32_t j) {
    double sum = 0.0;
    const V2 d = X[j] - X[i];
    for (std::uint32_t e : {i, j})
      for (auto k : adj[e]) {
        const std::uint32_t w = other(k, e);
        if (w == i || w == j) continue;
        const double c = cross(d, X[w] - X[i]);
        sum += c > 1e-9 ? 1.0 : c < -1e-9 ? -1.0 : 0.0;
      }
    return sum;
  };

  for (const DBond& bd : bonds) {
    const std::uint32_t i = bd.a, j = bd.b;
    const V2 a = X[i], b = X[j];
    const V2 dir = unit(b - a), n = perp(dir);
    if (bd.stereo == BondStereo::Wedge || bd.stereo == BondStereo::Hash) {
      auto [ca, cb] = clip_pair(i, j, a, b);
      const double frac = norm(cb - a) / std::max(1e-9, norm(b - a));
      const double half = 0.085 * L * frac;
      if (bd.stereo == BondStereo::Wedge) {
        const V2 p1 = cb + n * half, p2 = cb - n * half;
        fills_svg += "<polygon points=\"" + num(ca.x) + "," + num(ca.y) + " " + num(p1.x) + "," + num(p1.y) + " " + num(p2.x) + "," +
                     num(p2.y) + "\"/>";
      } else {
        const double len = norm(cb - ca);
        const int count = std::max(4, static_cast<int>(std::lround(len / (0.13 * L))));
        for (int k = 1; k <= count; ++k) {
          const double t = static_cast<double>(k) / static_cast<double>(count);
          const V2 c = ca + (cb - ca) * t;
          const double w = std::max(lw, half * (norm(c - a) / std::max(1e-9, norm(cb - a))));
          line_el(c + n * w, c - n * w);
        }
      }
      continue;
    }
    if (bd.order == 1 || bd.order > 3) {
      auto [ca, cb] = clip_pair(i, j, a, b);
      line_el(ca, cb);
      continue;
    }
    if (bd.order == 3) {
      for (double o : {-1.0, 0.0, 1.0}) {
        auto [ca, cb] = clip_pair(i, j, a + n * (o * sp2 * 0.85), b + n * (o * sp2 * 0.85));
        line_el(ca, cb);
      }
      continue;
    }
    // double bond
    const bool in_ring = bd.global >= 0 && (st.bonds.flags[static_cast<std::size_t>(bd.global)] & kBondInRing);
    const std::size_t di = adj[i].size(), dj = adj[j].size();
    const bool centred = !in_ring && ((di == 1 || dj == 1) ? (has_label(i) || has_label(j) || std::max(di, dj) >= 3)
                                                            : (di >= 3 && dj >= 3) || g.linear[i] || g.linear[j]);
    if (centred) {
      for (double o : {-0.5, 0.5}) {
        auto [ca, cb] = clip_pair(i, j, a + n * (o * sp2), b + n * (o * sp2));
        line_el(ca, cb);
      }
      continue;
    }
    double side = 1.0;
    const int r = in_ring ? ring_for_bond(i, j) : -1;
    if (r >= 0) side = dot(n, rinfo[static_cast<std::size_t>(r)].centre - a) >= 0.0 ? 1.0 : -1.0;
    else side = busy_side(i, j) >= 0.0 ? 1.0 : -1.0;
    {
      auto [ca, cb] = clip_pair(i, j, a, b);
      line_el(ca, cb);
    }
    const double trim = 0.8 * sp2;
    V2 ia = a + n * (side * sp2) + dir * (has_label(i) ? 0.0 : trim);
    V2 ib = b + n * (side * sp2) - dir * (has_label(j) ? 0.0 : trim);
    auto [ca, cb] = clip_pair(i, j, ia, ib);
    line_el(ca, cb);
  }

  // atom labels
  for (std::size_t i = 0; i < N; ++i)
    for (const auto& l : labels[i].lines) text_svg += emit_line(l, fs, labels[i].color);

  // annotations: atom numbers and CIP descriptors
  std::vector<Box> extra_boxes;
  std::vector<std::string> note(g.n);
  if (opt.atom_numbers)
    for (std::uint32_t i = 0; i < g.n; ++i) note[i] = std::to_string(i + 1);
  std::vector<std::pair<std::uint32_t, CipLabel>> ez_notes;
  if (opt.cip_labels) {
    for (const auto& c : cip_centres(st, sp))
      if (c.label != CipLabel::None) {
        auto& t = note[c.index - g.begin];
        t += std::string(t.empty() ? "" : " ") + "(" + to_string(c.label) + ")";
      }
    for (const auto& c : cip_double_bonds(st, sp))
      if (c.label != CipLabel::None) ez_notes.push_back({c.index, c.label});
    for (const auto& a : st.axial) {  // Ra/Sa on the allene centre, or on the first atom of a biaryl axis
      if (st.atoms.species[a.begin] != sp) continue;
      const CipLabel l = cip_label(st, a);
      if (l == CipLabel::None) continue;
      auto& t = note[(a.centre != kNoAtom ? a.centre : a.begin) - g.begin];
      t += std::string(t.empty() ? "" : " ") + "(" + to_string(l) + ")";
    }
  }
  const double nfs = 0.42 * L;
  auto note_box = [&](V2 c, const std::string& t) {
    double w = 0.0;
    for (char ch : t) w += glyph_width(ch) * nfs;
    extra_boxes.push_back({c.x - w / 2.0, c.y - 0.5 * nfs, c.x + w / 2.0, c.y + 0.5 * nfs});
    text_svg += "<text x=\"" + num(c.x) + "\" y=\"" + num(c.y + 0.36 * nfs) + "\" font-size=\"" + num(nfs) +
                "\" text-anchor=\"middle\" fill=\"currentColor\" fill-opacity=\"0.6\" font-style=\"italic\">" + t + "</text>";
  };
  std::vector<char> every(N, 1);
  for (std::uint32_t i = 0; i < g.n; ++i) {
    if (note[i].empty()) continue;
    std::vector<V2> dirs;
    for (auto k : adj[i]) dirs.push_back(X[other(k, i)] - X[i]);
    double ang = -kPi / 2.0;  // above (px y down)
    if (!dirs.empty()) {
      std::vector<double> a;
      for (V2 d : dirs) a.push_back(wrap(angle_of(d)));
      std::sort(a.begin(), a.end());
      double best = -1.0;
      for (std::size_t k = 0; k < a.size(); ++k) {
        const double gap = k + 1 < a.size() ? a[k + 1] - a[k] : a[0] + 2.0 * kPi - a[k];
        if (a.size() == 1 || gap > best) {
          best = gap;
          ang = a[k] + gap / 2.0;
        }
      }
    }
    const double dist = (labels[i].present ? 0.62 * fs : 0.0) + 0.5 * L;
    note_box(X[i] + polar(ang) * dist, note[i]);
  }
  for (auto [b, lab] : ez_notes) {
    const std::uint32_t i = st.bonds.begin[b] - g.begin, j = st.bonds.end[b] - g.begin;
    const V2 mid = (X[i] + X[j]) * 0.5, n = perp(unit(X[j] - X[i]));
    const double side = busy_side(i, j) >= 0.0 ? -1.0 : 1.0;
    note_box(mid + n * (side * 0.55 * L), std::string("(") + to_string(lab) + ")");
  }

  // bounds
  double minx = 1e300, miny = 1e300, maxx = -1e300, maxy = -1e300;
  auto grow = [&](double x0, double y0, double x1, double y1) {
    minx = std::min(minx, x0);
    miny = std::min(miny, y0);
    maxx = std::max(maxx, x1);
    maxy = std::max(maxy, y1);
  };
  for (std::size_t i = 0; i < N; ++i) {
    grow(X[i].x, X[i].y, X[i].x, X[i].y);
    for (const Box& b : labels[i].boxes) grow(b.x0, b.y0, b.x1, b.y1);
  }
  for (const Box& b : extra_boxes) grow(b.x0, b.y0, b.x1, b.y1);
  const double pad = 0.45 * L;
  const double ox = pad - minx, oy = pad - miny;
  D.width = std::ceil(maxx - minx + 2.0 * pad);
  D.height = std::ceil(maxy - miny + 2.0 * pad);

  D.svg = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " + num(D.width) + " " + num(D.height) + "\" width=\"" +
          num(D.width) + "\" height=\"" + num(D.height) +
          "\" font-family=\"Arial, Helvetica, sans-serif\" role=\"img\"><title>" + s.canonical_smiles + "</title>"
          "<g transform=\"translate(" + num(ox) + " " + num(oy) + ")\">"
          "<g stroke=\"currentColor\" stroke-width=\"" + num(lw) + "\" stroke-linecap=\"round\" fill=\"none\">" + lines_svg +
          "</g><g fill=\"currentColor\" stroke=\"currentColor\" stroke-width=\"" + num(lw * 0.6) + "\" stroke-linejoin=\"round\">" +
          fills_svg + "</g>" + text_svg + "</g></svg>";
  // the title must not break the markup
  {
    std::string safe;
    for (char c : s.canonical_smiles) {
      if (c == '&') safe += "&amp;";
      else if (c == '<') safe += "&lt;";
      else if (c == '>') safe += "&gt;";
      else safe += c;
    }
    const std::string raw = "<title>" + s.canonical_smiles + "</title>";
    const auto at = D.svg.find(raw);
    if (at != std::string::npos) D.svg.replace(at, raw.size(), "<title>" + safe + "</title>");
  }

  for (std::size_t i = 0; i < N; ++i)
    D.atoms.push_back({atom_id[i], X[i].x + ox, X[i].y + oy, labels[i].plain});
  for (const DBond& b : bonds) D.bonds.push_back({b.a, b.b, b.order, b.stereo});

  // ---- 2D molfile ----
  {
    std::string out = s.canonical_smiles + "\n  ChemLarp          2D\n\n";
    char buf[160];
    std::snprintf(buf, sizeof buf, "%3zu%3zu  0  0%3d  0  0  0  0  0999 V2000\n", N, bonds.size(), centres.empty() ? 0 : 1);
    out += buf;
    std::vector<std::pair<std::size_t, int>> charges, isotopes;
    std::vector<int> order_sum(N, 0), drawn_hydrogens(N, 0);
    for (const DBond& b : bonds) {
      order_sum[b.a] += b.order;
      order_sum[b.b] += b.order;
      if (atom_id[b.a] < 0) ++drawn_hydrogens[b.b];
      if (atom_id[b.b] < 0) ++drawn_hydrogens[b.a];
    }
    for (std::size_t i = 0; i < N; ++i) {
      std::uint8_t z = 1;
      int vvv = 0;
      if (atom_id[i] >= 0) {
        const auto ga = static_cast<std::size_t>(atom_id[i]);
        z = st.atoms.element[ga];
        if (st.atoms.formal_charge[ga]) charges.push_back({i + 1, st.atoms.formal_charge[ga]});
        if (st.atoms.isotope[ga]) isotopes.push_back({i + 1, st.atoms.isotope[ga]});
        vvv = mdl_valence_field(z, st.atoms.formal_charge[ga], order_sum[i],
                                std::max(0, static_cast<int>(st.atoms.implicit_h[ga]) - drawn_hydrogens[i]));
      }
      std::snprintf(buf, sizeof buf, "%10.4f%10.4f%10.4f %-3.3s 0  0  0  0  0%3d  0  0  0  0  0  0\n", P[i].x * 1.5, P[i].y * 1.5, 0.0,
                    std::string(element_symbol(z)).c_str(), vvv);
      out += buf;
    }
    // A drawn double bond fixes E/Z for any reader, so a stereogenic double bond
    // whose configuration the input leaves open is written as "either" (3).
    auto open_double_bond = [&](const DBond& b) {
      if (b.order != 2 || b.global < 0) return false;
      const auto gb = static_cast<std::uint32_t>(b.global);
      if (st.bonds.flags[gb] & kBondAromatic) return false;
      for (const auto& d : st.double_bonds)
        if (d.bond == gb) return false;
      const std::uint32_t u = st.bonds.begin[gb], v = st.bonds.end[gb];
      if (st.bonds.flags[gb] & kBondInRing) {  // E/Z only in rings of 8 or more
        std::vector<std::pair<std::uint32_t, int>> queue{{u, 0}};
        std::vector<std::uint32_t> seen{u};
        for (std::size_t h = 0; h < queue.size(); ++h) {
          const auto [x, d] = queue[h];
          if (d >= 6) continue;
          const auto nb = st.adjacency.neighbors(x);
          const auto bs = st.adjacency.bonds_of(x);
          for (std::size_t k = 0; k < nb.size(); ++k) {
            if (bs[k] == gb || std::find(seen.begin(), seen.end(), nb[k]) != seen.end()) continue;
            if (nb[k] == v) return false;
            seen.push_back(nb[k]);
            queue.push_back({nb[k], d + 1});
          }
        }
      }
      auto other = [&](std::uint32_t x, std::uint32_t not_this) -> std::int64_t {
        for (std::uint32_t n : st.adjacency.neighbors(x))
          if (n != not_this) return n;
        return -1;
      };
      return other(u, v) >= 0 && other(v, u) >= 0;
    };
    // The CIP labeller only ranks double bonds listed in the state, so test the
    // candidates on a copy that lists them.
    std::vector<char> either(bonds.size(), 0);
    {
      std::vector<std::size_t> cand;
      for (std::size_t k = 0; k < bonds.size(); ++k)
        if (open_double_bond(bonds[k])) cand.push_back(k);
      if (!cand.empty()) {
        ChemicalSystemState probe = st;
        for (std::size_t k : cand) {
          const auto gb = static_cast<std::uint32_t>(bonds[k].global);
          const std::uint32_t u = st.bonds.begin[gb], v = st.bonds.end[gb];
          DoubleBondStereo t;
          t.bond = gb;
          for (std::uint32_t n : st.adjacency.neighbors(u)) if (n != v) { t.ref_begin = n; break; }
          for (std::uint32_t n : st.adjacency.neighbors(v)) if (n != u) { t.ref_end = n; break; }
          probe.double_bonds.push_back(t);
        }
        for (std::size_t k : cand) {
          const auto gb = static_cast<std::uint32_t>(bonds[k].global);
          for (const auto& t : probe.double_bonds)
            if (t.bond == gb && cip_label(probe, t) != CipLabel::None) either[k] = 1;
        }
      }
    }
    for (std::size_t k = 0; k < bonds.size(); ++k) {
      const DBond& b = bonds[k];
      const unsigned stereo = either[k] ? 3u : static_cast<unsigned>(b.stereo);
      std::snprintf(buf, sizeof buf, "%3u%3u%3u%3u\n", b.a + 1, b.b + 1, static_cast<unsigned>(b.order), stereo);
      out += buf;
    }
    auto emit = [&](const char* tag, const std::vector<std::pair<std::size_t, int>>& v) {
      for (std::size_t k = 0; k < v.size(); k += 8) {
        const std::size_t m = std::min<std::size_t>(8, v.size() - k);
        std::snprintf(buf, sizeof buf, "M  %s%3zu", tag, m);
        out += buf;
        for (std::size_t q = 0; q < m; ++q) {
          std::snprintf(buf, sizeof buf, " %3zu %3d", v[k + q].first, v[k + q].second);
          out += buf;
        }
        out += "\n";
      }
    };
    emit("CHG", charges);
    emit("ISO", isotopes);
    out += "M  END\n";
    D.molfile = out;
  }
  D.ok = true;
  return D;
}

Depiction depict_smiles(std::string_view smiles, Diagnostics& diag, const DepictOptions& options) {
  ChemicalSystemState st;
  const std::int32_t sp = parse_smiles(st, smiles, diag);
  if (sp < 0 || diag.has_errors()) return {};
  return depict(st, static_cast<std::uint32_t>(sp), diag, options);
}

}  // namespace chemlarp
