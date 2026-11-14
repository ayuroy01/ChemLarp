#include "chemlarp/planar.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "chemlarp/cip.hpp"
#include "chemlarp/rings.hpp"

namespace chemlarp {

namespace {

struct P3 {
  double x = 0, y = 0, z = 0;
};
P3 operator-(P3 a, P3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
P3 operator+(P3 a, P3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
P3 cross(P3 a, P3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
double dot(P3 a, P3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double norm(P3 a) { return std::sqrt(dot(a, a)); }

bool contains(const std::vector<std::uint32_t>& v, std::uint32_t q) { return std::find(v.begin(), v.end(), q) != v.end(); }

bool bonded(const ChemicalSystemState& st, std::uint32_t a, std::uint32_t b) {
  for (std::uint32_t x : st.adjacency.neighbors(a))
    if (x == b) return true;
  return false;
}

// IUPAC torsion angle a-b-c-d in degrees (positive: clockwise rotation of the front bond, seen along b->c).
double torsion(P3 a, P3 b, P3 c, P3 d) {
  const P3 b1 = b - a, b2 = c - b, b3 = d - c;
  const P3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  return std::atan2(norm(b2) * dot(b1, n2), dot(n1, n2)) * 180.0 / 3.14159265358979323846;
}

}  // namespace

bool bridges_ring(const ChemicalSystemState& st, const std::vector<std::uint32_t>& ring, std::size_t ir, std::uint32_t x) {
  const std::size_t n = ring.size();
  const std::uint32_t r = ring[ir], prev = ring[(ir + n - 1) % n], next = ring[(ir + 1) % n];
  std::vector<std::uint32_t> stack = {x}, seen = {r, x};
  while (!stack.empty()) {
    const std::uint32_t u = stack.back();
    stack.pop_back();
    for (std::uint32_t w : st.adjacency.neighbors(u)) {
      if (contains(seen, w)) continue;
      if (contains(ring, w)) {
        if (w != prev && w != next) return true;
        continue;  // do not walk through the ring itself
      }
      seen.push_back(w);
      stack.push_back(w);
    }
  }
  return false;
}

void perceive_planar_stereo(ChemicalSystemState& st, std::uint32_t sp) {
  auto pos = [&](std::uint32_t i) { return P3{st.atoms.x[i], st.atoms.y[i], st.atoms.z[i]}; };
  auto aromatic = [&](std::uint32_t i) { return (st.atoms.flags[i] & kAtomAromatic) != 0; };
  const std::vector<Ring> rings = find_sssr(st, sp);

  // --- rings with a metal or a bridge ---------------------------------------------------------------------
  for (const Ring& rg : rings) {
    const auto& ring = rg.atoms;
    const std::size_t n = ring.size();
    if (n < 5 || n > 6) continue;
    P3 cen;
    for (std::uint32_t a : ring) cen = cen + pos(a);
    cen = {cen.x / n, cen.y / n, cen.z / n};
    P3 nrm;  // atoms[] run counterclockwise seen from the side nrm points to
    for (std::size_t i = 0; i < n; ++i) nrm = nrm + cross(pos(ring[i]) - cen, pos(ring[(i + 1) % n]) - cen);
    if (norm(nrm) < 1e-6) continue;
    nrm = {nrm.x / norm(nrm), nrm.y / norm(nrm), nrm.z / norm(nrm)};

    // Metal over the ring: bonded to every ring atom, or (hapto bonds are cut on reading) a metal atom of the
    // species within 2.6 Å of every ring atom.
    std::uint32_t metal = kNoAtom;
    for (std::uint32_t x : st.adjacency.neighbors(ring[0])) {
      if (contains(ring, x)) continue;
      bool all = true;
      for (std::uint32_t a : ring) all = all && bonded(st, a, x);
      if (all) metal = x;
    }
    if (metal == kNoAtom) {
      const Species& s = st.species[sp];
      for (std::uint32_t x = s.atom_begin; x < s.atom_begin + s.atom_count && metal == kNoAtom; ++x) {
        const int z = st.atoms.element[x];
        if (z < 21 || contains(ring, x)) continue;  // d- and f-block metals and heavier
        bool all = true;
        for (std::uint32_t a : ring) all = all && norm(pos(a) - pos(x)) < 2.6;
        if (all) metal = x;
      }
    }
    bool substituted = false;
    for (std::size_t i = 0; i < n; ++i) {
      if (st.atoms.element[ring[i]] != 6) substituted = true;
      for (std::uint32_t x : st.adjacency.neighbors(ring[i]))
        if (!contains(ring, x) && x != metal && st.atoms.element[x] > 1 && !bridges_ring(st, ring, i, x)) substituted = true;
    }
    if (!substituted) continue;
    PlanarStereo u;
    u.atoms = ring;
    if (metal != kNoAtom) {
      u.kind = PlanarKind::Metallocene;
      u.pilots = {metal};
      u.pilot_ccw = {static_cast<char>(dot(nrm, cen - pos(metal)) > 0.0)};
      st.planar.push_back(u);
      continue;
    }
    bool all_aromatic = true;
    for (std::uint32_t a : ring) all_aromatic = all_aromatic && aromatic(a);
    if (!all_aromatic) continue;
    std::vector<std::uint32_t> cand;
    for (std::size_t i = 0; i < n; ++i)
      for (std::uint32_t x : st.adjacency.neighbors(ring[i])) {
        if (contains(ring, x) || aromatic(x) || st.atoms.element[x] <= 1) continue;
        if (std::fabs(dot(pos(x) - cen, nrm)) < 0.25) continue;
        if (bridges_ring(st, ring, i, x) && !contains(cand, x)) cand.push_back(x);
      }
    if (cand.empty()) continue;
    // The highest-ranked candidates (all of them when tied).
    std::vector<std::uint32_t> top = {cand[0]};
    for (std::size_t i = 1; i < cand.size(); ++i) {
      const int c = cip_compare_atoms(st, cand[i], top[0]);
      if (c > 0) top = {cand[i]};
      else if (c == 0) top.push_back(cand[i]);
    }
    u.kind = PlanarKind::Cyclophane;
    u.pilots = top;
    for (std::uint32_t p : top) u.pilot_ccw.push_back(static_cast<char>(dot(nrm, pos(p) - cen) > 0.0));
    st.planar.push_back(u);
  }

  // --- helicenes ---------------------------------------------------------------------------------------
  std::vector<const Ring*> ar;
  for (const Ring& rg : rings) {
    if (rg.atoms.size() < 5 || rg.atoms.size() > 6) continue;
    bool ok = true;
    for (std::uint32_t a : rg.atoms) ok = ok && aromatic(a);
    if (ok) ar.push_back(&rg);
  }
  const std::size_t nr = ar.size();
  if (nr < 5) return;
  std::vector<int> ring_count(st.atoms.size(), 0);
  for (const Ring& rg : rings)
    for (std::uint32_t a : rg.atoms) ++ring_count[a];
  // fusion[i][j] = the shared bond (two atoms) of fused rings i and j.
  std::vector<std::vector<std::array<std::uint32_t, 2>>> fusion(nr, std::vector<std::array<std::uint32_t, 2>>(nr, {kNoAtom, kNoAtom}));
  for (std::size_t i = 0; i < nr; ++i)
    for (std::size_t j = i + 1; j < nr; ++j) {
      std::vector<std::uint32_t> shared;
      for (std::uint32_t a : ar[i]->atoms)
        if (contains(ar[j]->atoms, a)) shared.push_back(a);
      if (shared.size() == 2 && bonded(st, shared[0], shared[1])) fusion[i][j] = fusion[j][i] = {shared[0], shared[1]};
    }
  auto fused = [&](std::size_t i, std::size_t j) { return fusion[i][j][0] != kNoAtom; };
  // In ring m fused to rings l and r: the atoms of the two fusion bonds that are bonded to each other
  // (ortho-angular fusion), as {inner atom on the l side, inner atom on the r side}; none if the fusion
  // bonds share an atom or are not linked by one bond.
  auto inner = [&](std::size_t l, std::size_t m, std::size_t r, std::array<std::uint32_t, 2>& out) {
    const auto fl = fusion[l][m], fr = fusion[m][r];
    for (std::uint32_t a : fl)
      if (a == fr[0] || a == fr[1]) return false;
    int found = 0;
    for (std::uint32_t a : fl)
      for (std::uint32_t b : fr)
        if (bonded(st, a, b)) out = {a, b}, ++found;
    return found == 1;
  };
  struct Chain {
    std::vector<std::size_t> rings;
    std::vector<std::uint32_t> rim;  // u_1 .. u_{k-1}
  };
  std::vector<Chain> found;
  // Depth-first growth from every ordered pair of fused rings.
  std::vector<Chain> stack;
  for (std::size_t i = 0; i < nr; ++i)
    for (std::size_t j = 0; j < nr; ++j)
      if (fused(i, j)) stack.push_back({{i, j}, {}});
  while (!stack.empty()) {
    Chain c = stack.back();
    stack.pop_back();
    bool grown = false;
    const std::size_t l = c.rings[c.rings.size() - 2], m = c.rings.back();
    for (std::size_t r = 0; r < nr; ++r) {
      if (!fused(m, r) || std::find(c.rings.begin(), c.rings.end(), r) != c.rings.end()) continue;
      std::array<std::uint32_t, 2> in{};
      if (!inner(l, m, r, in)) continue;
      if (!c.rim.empty() && c.rim.back() != in[0]) continue;  // turns the other way (zigzag, as in chrysene)
      if (ring_count[in[0]] != 2 || ring_count[in[1]] != 2) continue;
      Chain d = c;
      d.rings.push_back(r);
      if (d.rim.empty()) d.rim.push_back(in[0]);
      d.rim.push_back(in[1]);
      stack.push_back(d);
      grown = true;
    }
    if (!grown && c.rings.size() >= 5) found.push_back(c);
  }
  // Keep maximal chains once (each is found from both ends and from inner starting pairs).
  std::vector<std::vector<std::size_t>> seen_sets;
  std::sort(found.begin(), found.end(), [](const Chain& a, const Chain& b) { return a.rings.size() > b.rings.size(); });
  for (const Chain& c : found) {
    std::vector<std::size_t> key = c.rings;
    std::sort(key.begin(), key.end());
    bool sub = false;
    for (const auto& s : seen_sets) sub = sub || std::includes(s.begin(), s.end(), key.begin(), key.end());
    if (sub) continue;
    seen_sets.push_back(key);
    // Inner rim with the terminal ring atoms: t0, u_1 .. u_{k-1}, t_end.
    std::vector<std::uint32_t> rim;
    const auto& first = ar[c.rings.front()]->atoms;
    const auto& last = ar[c.rings.back()]->atoms;
    const auto f0 = fusion[c.rings[0]][c.rings[1]];
    const auto fl = fusion[c.rings[c.rings.size() - 2]][c.rings.back()];
    std::uint32_t t0 = kNoAtom, t1 = kNoAtom;
    for (std::uint32_t a : first)
      if (a != f0[0] && a != f0[1] && bonded(st, a, c.rim.front())) t0 = a;
    for (std::uint32_t a : last)
      if (a != fl[0] && a != fl[1] && bonded(st, a, c.rim.back())) t1 = a;
    if (t0 == kNoAtom || t1 == kNoAtom) continue;
    rim.push_back(t0);
    rim.insert(rim.end(), c.rim.begin(), c.rim.end());
    rim.push_back(t1);
    double sum = 0.0;
    for (std::size_t k = 0; k + 3 < rim.size(); ++k) sum += torsion(pos(rim[k]), pos(rim[k + 1]), pos(rim[k + 2]), pos(rim[k + 3]));
    if (std::fabs(sum) < 15.0) continue;
    PlanarStereo u;
    u.kind = PlanarKind::Helicene;
    u.atoms = rim;
    u.plus = sum > 0.0;
    u.rings = static_cast<std::uint32_t>(c.rings.size());
    st.planar.push_back(u);
  }
}

}  // namespace chemlarp
