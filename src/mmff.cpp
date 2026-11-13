#include "chemlarp/mmff.hpp"

#include "chemlarp/conformer.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/json.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <numbers>
#include <set>
#include <tuple>

namespace chemlarp {
namespace {

struct MmffProp { int atno, crd, val, pilp, mltb, arom, lin, sbmb; };
struct MmffDef { const char* symbol; int level[4]; };
struct MmffPbci { double pbci, fcadj; bool present; };
struct MmffChg { int bt, i, j; double bci; };
struct MmffBond { int bt, i, j; double kb, r0; };
struct MmffBndk { int zi, zj; double r0, kb; };
struct MmffHL { int ri, rj; double a, d, dp; };
struct MmffCovRad { int z; double r, chi; };
struct MmffAngle { int at, i, j, k; double ka, theta0; };
struct MmffStbn { int sbt, i, j, k; double kijk, kkji; };
struct MmffDfsb { int ri, rj, rk; double f1, f2; };
struct MmffOop { int i, j, k, l; double koop; };
struct MmffTor { int tt, i, j, k, l; double V1, V2, V3; };
struct MmffVdw { int type; double alpha, N, A, G; int da; };

#include "mmff_data.inc"

constexpr double kDeg = 180.0 / std::numbers::pi;
constexpr double kAngleK = 143.9325 / (kDeg * kDeg);  // 0.043844 kcal/(mol deg^2) per md/Å

// ---------------------------------------------------------------------------------------------------------
// Lookup tables built once.

struct Tables {
  std::map<std::tuple<int, int, int>, double> chg;
  std::map<std::tuple<int, int, int>, std::pair<double, double>> bond;
  std::map<std::pair<int, int>, std::pair<double, double>> bndk;
  std::map<std::tuple<int, int, int, int>, std::pair<double, double>> angle;
  std::map<std::tuple<int, int, int, int>, std::pair<double, double>> stbn;
  std::map<std::tuple<int, int, int, int>, double> oop94, oop94s;
  std::map<std::tuple<int, int, int, int, int>, std::array<double, 3>> tor94, tor94s;
  std::map<int, MmffVdw> vdw;
  Tables() {
    for (const auto& c : kMmffChg) chg[{c.bt, c.i, c.j}] = c.bci;
    for (const auto& b : kMmffBond) bond[{b.bt, b.i, b.j}] = {b.kb, b.r0};
    for (const auto& b : kMmffBndk) bndk[{b.zi, b.zj}] = {b.r0, b.kb};
    for (const auto& a : kMmffAngle) angle[{a.at, a.i, a.j, a.k}] = {a.ka, a.theta0};
    for (const auto& s : kMmffStbn) stbn[{s.sbt, s.i, s.j, s.k}] = {s.kijk, s.kkji};
    for (const auto& o : kMmffOop) oop94[{o.i, o.j, o.k, o.l}] = o.koop;
    for (const auto& o : kMmffsOop) oop94s[{o.i, o.j, o.k, o.l}] = o.koop;
    for (const auto& t : kMmffTor) tor94[{t.tt, t.i, t.j, t.k, t.l}] = {t.V1, t.V2, t.V3};
    for (const auto& t : kMmffsTor) tor94s[{t.tt, t.i, t.j, t.k, t.l}] = {t.V1, t.V2, t.V3};
    for (const auto& v : kMmffVdw) vdw[v.type] = v;
  }
};
const Tables& tables() {
  static const Tables t;
  return t;
}

int level(int type, int lv) { return lv == 0 ? type : kMmffDef[type].level[lv]; }  // lv 0..3 (3: wildcard 0)

// ---------------------------------------------------------------------------------------------------------
// Molecule view: local atoms with neighbours and bond orders, rings up to six members.

struct Mol {
  std::size_t n = 0;
  std::vector<int> z, fc;
  std::vector<std::vector<int>> nb;      // neighbours
  std::vector<std::vector<int>> order;   // bond order to each neighbour
  std::vector<std::vector<int>> rings;   // simple rings of 3-6 atoms
  std::vector<bool> arom_ring;           // per ring, MMFF aromatic (5 or 6 members)
  int bo(int a, int b) const {
    for (std::size_t k = 0; k < nb[a].size(); ++k)
      if (nb[a][k] == b) return order[a][k];
    return 0;
  }
  bool bonded(int a, int b) const { return bo(a, b) > 0; }
  int nH(int a) const {
    int h = 0;
    for (int b : nb[a]) h += z[b] == 1;
    return h;
  }
  int deg(int a) const { return int(nb[a].size()); }
  int count_double(int a) const {
    int d = 0;
    for (int o : order[a]) d += o == 2;
    return d;
  }
  bool has_triple(int a) const {
    for (int o : order[a])
      if (o == 3) return true;
    return false;
  }
  int double_partner(int a) const {
    for (std::size_t k = 0; k < nb[a].size(); ++k)
      if (order[a][k] == 2) return nb[a][k];
    return -1;
  }
  bool in_ring_size(int a, std::size_t size) const {
    for (const auto& r : rings)
      if (r.size() == size && std::find(r.begin(), r.end(), a) != r.end()) return true;
    return false;
  }
  bool same_ring(std::size_t size, std::initializer_list<int> atoms) const {
    for (const auto& r : rings) {
      if (r.size() != size) continue;
      bool all = true;
      for (int a : atoms) all = all && std::find(r.begin(), r.end(), a) != r.end();
      if (all) return true;
    }
    return false;
  }
  bool in_arom_ring(int a, std::size_t size = 0) const {
    for (std::size_t r = 0; r < rings.size(); ++r)
      if (arom_ring[r] && (size == 0 || rings[r].size() == size) && std::find(rings[r].begin(), rings[r].end(), a) != rings[r].end())
        return true;
    return false;
  }
  // Number of terminal (one-connected) oxygen or sulfur neighbours.
  int terminal(int a, int el) const {
    int c = 0;
    for (int b : nb[a]) c += z[b] == el && deg(b) == 1;
    return c;
  }
};

void find_rings(Mol& m) {
  std::set<std::vector<int>> seen;
  std::vector<int> path;
  std::vector<bool> on(m.n, false);
  // Rings through their lowest-numbered atom, up to six members.
  std::function<void(int, int)> dfs = [&](int start, int a) {
    for (int b : m.nb[a]) {
      if (b == start && path.size() >= 3) {
        std::vector<int> key = path;
        std::sort(key.begin(), key.end());
        if (seen.insert(key).second) m.rings.push_back(path);
        continue;
      }
      if (on[b] || b < start || path.size() >= 6) continue;
      on[b] = true;
      path.push_back(b);
      dfs(start, b);
      path.pop_back();
      on[b] = false;
    }
  };
  for (std::size_t s = 0; s < m.n; ++s) {
    if (m.z[s] == 1) continue;
    path = {int(s)};
    on[s] = true;
    dfs(int(s), int(s));
    on[s] = false;
  }
}

// MMFF aromaticity on the Kekulé structure: a 6-ring is aromatic when each of its atoms has a double bond
// either inside the ring or to an atom of an already-aromatic ring sharing a bond with it; a 5-ring when
// four atoms qualify that way and the fifth is a π lone-pair donor (N with three neighbours, O, S, or an
// anionic N or C). Fused systems are resolved by iterating to a fixed point.
void find_aromatic(Mol& m) {
  m.arom_ring.assign(m.rings.size(), false);
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::size_t r = 0; r < m.rings.size(); ++r) {
      const auto& ring = m.rings[r];
      if (m.arom_ring[r] || (ring.size() != 5 && ring.size() != 6)) continue;
      const std::set<int> inring(ring.begin(), ring.end());
      int pi = 0, donors = 0;
      bool bad = false;
      for (int a : ring) {
        int partner = -1;
        int doubles = 0;
        for (std::size_t k = 0; k < m.nb[a].size(); ++k)
          if (m.order[a][k] == 2) partner = m.nb[a][k], ++doubles;
        if (m.has_triple(a) || doubles > 1) { bad = true; break; }
        bool ok = false;
        if (partner >= 0) {
          if (inring.count(partner)) ok = true;
          else
            for (std::size_t q = 0; q < m.rings.size() && !ok; ++q)
              if (q != r && m.arom_ring[q]) {
                const auto& o = m.rings[q];
                ok = std::find(o.begin(), o.end(), a) != o.end() && std::find(o.begin(), o.end(), partner) != o.end();
              }
        }
        if (ok) { ++pi; continue; }
        if (partner >= 0) { bad = true; break; }  // exocyclic double bond (C=O): not aromatic
        const int z = m.z[a];
        const bool donor = (z == 7 && m.deg(a) == 3 && (m.fc[a] == 0 || m.fc[a] == 1)) || (z == 8 && m.deg(a) == 2) ||
                           (z == 16 && m.deg(a) == 2) || (z == 7 && m.deg(a) == 2 && m.fc[a] == -1) ||
                           (z == 6 && m.deg(a) == 3 && m.fc[a] == -1);
        if (!donor) { bad = true; break; }
        ++donors;
      }
      if (bad) continue;
      if ((ring.size() == 6 && pi == 6) || (ring.size() == 5 && pi == 4 && donors == 1)) {
        m.arom_ring[r] = true;
        changed = true;
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------------------
// Atom typing.

struct Typer {
  const Mol& m;
  std::vector<int> t;
  explicit Typer(const Mol& mol) : m(mol), t(mol.n, 0) {}

  bool carbonyl_like(int c) const {  // carbon with a double bond to O or S (amide, thioamide, carbamate ...)
    if (m.z[c] != 6) return false;
    const int d = m.double_partner(c);
    return d >= 0 && (m.z[d] == 8 || m.z[d] == 16);
  }
  bool sulfonyl(int s) const {  // S with two or more terminal O (or =N) neighbours
    if (m.z[s] != 16) return false;
    int o = 0;
    for (int b : m.nb[s]) o += (m.z[b] == 8 && m.deg(b) == 1);
    return o >= 2 || (o >= 1 && m.deg(s) == 4);
  }
  // Central carbon of an amidinium or guanidinium ion: C=N+ with further N neighbours that are not in a
  // double bond; returns the number of such nitrogens (including the N+), 0 if not.
  int cationic_n_centre(int c) const {
    if (m.z[c] != 6 || m.deg(c) != 3 || m.in_arom_ring(c, 6)) return 0;
    const int d = m.double_partner(c);
    if (d < 0 || m.z[d] != 7 || m.deg(d) != 3) return 0;
    for (int b : m.nb[d])
      if (m.z[b] == 8 && m.deg(b) == 1) return 0;  // N-oxide, not an iminium
    // The =N and the single-bonded three-connected nitrogens share one positive charge, wherever the
    // structure places it.
    int n = 1, charge = m.fc[d];
    for (int b : m.nb[c])
      if (b != d && m.z[b] == 7 && m.deg(b) == 3 && m.count_double(b) == 0) ++n, charge += m.fc[b];
    return n >= 2 && charge == 1 ? n : 0;
  }

  void aromatic() {
    // Six-membered aromatic rings.
    for (std::size_t r = 0; r < m.rings.size(); ++r) {
      if (!m.arom_ring[r] || m.rings[r].size() != 6) continue;
      for (int a : m.rings[r]) {
        if (m.z[a] == 6) t[a] = 37;
        else if (m.z[a] == 7) {
          bool oxide = false;
          for (int b : m.nb[a]) oxide = oxide || (m.z[b] == 8 && m.deg(b) == 1);
          t[a] = m.deg(a) == 2 ? 38 : oxide ? 69 : 58;
        }
      }
    }
    // Five-membered rings (take precedence over the six-ring types at fusion atoms).
    for (std::size_t r = 0; r < m.rings.size(); ++r) {
      if (!m.arom_ring[r] || m.rings[r].size() != 5) continue;
      const auto& ring = m.rings[r];
      // The π donor: the ring atom without a double bond inside this or a fused aromatic ring.
      int donor = -1;
      for (int a : ring) {
        const int p = m.double_partner(a);
        if (p < 0) donor = a;
      }
      if (donor < 0) continue;
      const int zd = m.z[donor];
      // Cationic centre (imidazolium-like): a ring carbon double-bonded to a ring N+ (not an N-oxide) and
      // bonded to a further nitrogen with a lone pair, in the ring (the donor) or outside it.
      int cation_c = -1;
      for (int c : ring) {
        if (m.z[c] != 6) continue;
        const int p = m.double_partner(c);
        if (p < 0 || m.z[p] != 7 || m.deg(p) != 3) continue;
        bool oxide = false;
        for (int b : m.nb[p]) oxide = oxide || (m.z[b] == 8 && m.deg(b) == 1);
        if (oxide) continue;
        for (int b : m.nb[c])
          if (b != p && m.z[b] == 7 && m.deg(b) == 3 && m.count_double(b) == 0 && m.fc[p] + m.fc[b] == 1) cation_c = c;
      }
      const bool anion = zd == 7 && m.fc[donor] == -1;
      if (zd == 7 && m.fc[donor] == 1 && cation_c < 0) continue;  // a cationic donor needs a cationic centre
      const std::size_t di = std::find(ring.begin(), ring.end(), donor) - ring.begin();
      for (std::size_t k = 0; k < 5; ++k) {
        const int a = ring[k];
        const std::size_t dist = std::min((k + 5 - di) % 5, (di + 5 - k) % 5);
        const int z = m.z[a];
        if (cation_c >= 0) {
          if (a == cation_c) { t[a] = 80; continue; }
          if (z == 7 && m.bonded(a, cation_c)) { t[a] = 81; continue; }
          if (zd == 7 && m.bonded(donor, cation_c)) {  // imidazolium: the charge is shared by both N
            t[a] = z == 6 ? 78 : 79;
            continue;
          }
        }
        if (anion) {
          if (z == 7) t[a] = 76;
          else if (z == 6) t[a] = 78;
          continue;
        }
        if (a == donor) {
          t[a] = zd == 7 ? 39 : zd == 8 ? 59 : zd == 16 ? 44 : 78;
          continue;
        }
        bool oxide = false;
        for (int b : m.nb[a]) oxide = oxide || (m.z[b] == 8 && m.deg(b) == 1 && m.z[a] == 7);
        int ty = 0;
        if (z == 6) ty = dist == 1 ? 63 : 64;
        else if (z == 7) ty = oxide ? 82 : (m.fc[a] == 1 && m.deg(a) == 3) ? 81 : dist == 1 ? 65 : 66;
        // An atom alpha in one fused five-ring and beta in another takes the general type.
        if ((t[a] == 63 && ty == 64) || (t[a] == 64 && ty == 63)) ty = 78;
        if ((t[a] == 65 && ty == 66) || (t[a] == 66 && ty == 65)) ty = 79;
        if (t[a] == 78 || t[a] == 79) ty = t[a];
        t[a] = ty;
      }
    }
  }

  int carbon(int a) {
    const int d = m.deg(a);
    if (d == 4) return m.in_ring_size(a, 3) ? 22 : m.in_ring_size(a, 4) ? 20 : 1;
    if (d == 1 && m.fc[a] == -1 && m.has_triple(a)) return 60;  // isonitrile carbon
    if (m.has_triple(a) || m.count_double(a) == 2) return 4;
    if (d == 3) {
      const int p = m.double_partner(a);
      if (p < 0) return 2;
      const int zp = m.z[p];
      if (zp == 6) return m.in_ring_size(a, 4) ? 30 : 2;
      if (zp == 7) return cationic_n_centre(a) ? 57 : 3;
      if (zp == 15) return 3;
      if (zp == 8 || zp == 16) {
        // Carboxylate / thiocarboxylate: a second terminal O or S carrying a negative charge.
        for (int b : m.nb[a])
          if (b != p && (m.z[b] == 8 || m.z[b] == 16) && m.deg(b) == 1 && m.fc[b] == -1) return 41;
        return 3;
      }
      return 2;
    }
    if (d == 2) return 4;
    return 1;
  }

  int nitrogen(int a) {
    const int d = m.deg(a), q = m.fc[a];
    int oxides = 0, o_double = 0;
    for (std::size_t k = 0; k < m.nb[a].size(); ++k) {
      const int b = m.nb[a][k];
      if (m.z[b] == 8 && m.deg(b) == 1) {
        if (m.order[a][k] == 2) ++o_double;
        else if (m.fc[b] == -1) ++oxides;
      }
    }
    if (d == 4) return oxides ? 68 : 34;
    if (d == 3) {
      if (o_double + oxides >= 2) return 45;  // nitro, nitrate
      const int p = m.double_partner(a);
      if (q == 1 || p >= 0) {
        if (oxides || o_double) return 67;  // N-oxide of an imine (nitrone)
        if (p >= 0 && m.z[p] == 6) {
          const int n = cationic_n_centre(p);
          if (n == 3) return 56;
          if (n == 2) return 55;
        }
        return 54;
      }
      for (int b : m.nb[a]) {
        if (t[b] == 80) {  // exocyclic N on a cationic ring centre: amidinium- or guanidinium-like
          int nn = 0;
          for (int c : m.nb[b]) nn += m.z[c] == 7;
          return nn >= 3 ? 56 : 55;
        }
        const int n = cationic_n_centre(b);
        if (n == 3) return 56;
        if (n == 2) return 55;
      }
      for (int b : m.nb[a]) {
        if (sulfonyl(b)) return 43;
        if (m.z[b] == 15 && m.terminal(b, 8) >= 2) return 43;  // N on a PO2 group (phosphonamidate)
        if (m.z[b] == 6 && m.deg(b) == 2 && m.has_triple(b)) {   // cyanamide
          for (std::size_t k = 0; k < m.nb[b].size(); ++k)
            if (m.order[b][k] == 3 && m.z[m.nb[b][k]] == 7) return 43;
        }
      }
      for (int b : m.nb[a])
        if (carbonyl_like(b)) return 10;
      for (int b : m.nb[a]) {
        if (m.z[b] == 6 && m.deg(b) == 3 && (m.count_double(b) == 1 || m.in_arom_ring(b))) {
          const int pp = m.double_partner(b);
          if (pp < 0 || m.z[pp] == 6 || m.z[pp] == 7 || m.z[pp] == 15) return 40;
        }
        if (m.z[b] == 6 && m.has_triple(b)) return 40;  // ynamine
      }
      for (int b : m.nb[a])
        if (m.z[b] == 7 && m.count_double(b) == 1) {  // N-N=N (triazene); N-N=C stays an amine
          const int pp = m.double_partner(b);
          if (pp >= 0 && m.z[pp] == 7) return 10;
        }
      return 8;
    }
    if (d == 2) {
      if (q == -1) {
        // Anionic N on S: a sulfonamide anion when the S carries two terminal O, otherwise the S=N of a
        // sulfilimine or sulfoximine drawn as S+–N−.
        for (int b : m.nb[a])
          if (m.z[b] == 16 && m.terminal(b, 8) < 2 && m.fc[b] > 0) return 48;
        return 62;
      }
      if (m.count_double(a) == 2) return 53;
      if (m.has_triple(a)) return 61;  // isonitrile N+
      {
        const int p0 = m.double_partner(a);
        if (p0 >= 0 && m.z[p0] == 16)
          for (int b : m.nb[a])
            if (sulfonyl(b) && m.bo(a, b) == 1) return 43;  // N(=S)-SO2: sulfonamide-like
      }
      const int p = m.double_partner(a);
      if (p >= 0) {
        if (m.z[p] == 8) return 46;
        if (m.z[p] == 16) return 48;
        return 9;
      }
      return 9;
    }
    if (d == 1) {
      if (m.has_triple(a)) return 42;
      const int p = m.nb[a][0];
      if (m.z[p] == 7 && m.count_double(p) == 2) return 47;  // azide terminal
      return 42;
    }
    return 8;
  }

  int oxygen(int a) {
    const int d = m.deg(a);
    if (d == 3) return 49;
    if (d == 2) {
      if (m.fc[a] == 1) return 51;
      if (m.nH(a) == 2) return 70;
      return 6;
    }
    if (d == 1) {
      const int x = m.nb[a][0];
      const int bo = m.order[a][0];
      const int zx = m.z[x];
      if (zx == 6) {
        if (t[x] == 41) return 32;
        if (bo == 2) return 7;
        return 35;
      }
      if (zx == 7) {
        int ox = 0;
        for (int b : m.nb[x]) ox += (m.z[b] == 8 && m.deg(b) == 1);
        if (bo == 2 && ox == 1 && m.fc[x] == 0) return 7;  // nitroso
        if (bo == 1 && m.fc[a] == -1 && m.fc[x] == 0 && ox == 1) return 35;  // oximate O-
        return 32;
      }
      if (zx == 16) {
        const int ox = m.terminal(x, 8) + m.terminal(x, 16);
        if (ox == 1 && bo == 2 && m.deg(x) == 3) return 7;  // sulfoxide
        if (ox == 1 && bo == 2 && m.deg(x) == 2) return 7;  // sulfine C=S=O
        return 32;
      }
      if (zx == 15 || zx == 17) return 32;
      return bo == 2 ? 7 : 35;
    }
    return 6;
  }

  int sulfur(int a) {
    const int d = m.deg(a);
    if (d == 1) {
      const int x = m.nb[a][0];
      if (m.fc[a] == -1 || m.z[x] == 15) return 72;
      if (m.z[x] == 6 && t[x] == 41) return 72;  // thiocarboxylate
      if (m.z[x] == 16 && m.deg(x) == 3) return 72;  // thiosulfinate terminal S
      return 16;
    }
    if (d == 2) {
      if (m.count_double(a) == 2) return 74;
      if (m.count_double(a) == 1) {
        const int p = m.double_partner(a);
        if (p >= 0 && m.z[p] == 8) return 74;
      }
      return 15;
    }
    if (d == 3) {
      const int ot = m.terminal(a, 8);
      if (ot >= 2) {  // sulfene >C=SO2 (either drawing), else sulfinate
        const int p = m.double_partner(a);
        return (m.count_double(a) == 3 || (p >= 0 && m.z[p] == 6)) ? 18 : 73;
      }
      if (ot == 1 && m.terminal(a, 16) == 1) return 73;       // thiosulfinate
      int neg = 0;
      for (int b : m.nb[a]) neg += (m.z[b] == 8 && m.fc[b] == -1 && m.deg(b) == 1);
      if (ot == 1 && neg == 1 && m.fc[a] == 1) return 17;  // dative sulfoxide
      return 17;
    }
    return 18;
  }

  int phosphorus(int a) {
    const int d = m.deg(a);
    if (d <= 3) {
      const int p = m.double_partner(a);
      if (p >= 0 && m.z[p] == 6) return 75;  // phosphaalkene
    }
    if (d >= 4) return 25;
    if (d == 3) {
      const int p = m.double_partner(a);
      if (p >= 0 && m.z[p] == 6) return 75;
      return 26;
    }
    return 26;
  }

  int hydrogen(int a) {
    if (m.deg(a) != 1) return 5;
    const int x = m.nb[a][0];
    const int tx = t[x];
    switch (m.z[x]) {
      case 6:
      case 14: return 5;
      case 15: return 71;
      case 16: return 71;
      case 7:
        switch (tx) {
          case 8: case 39: case 62: case 79: return 23;
          case 9: return 27;
          case 10: case 40: case 43: case 48: case 42: return 28;  // includes H on S=N-H (HNSO)
          case 67: case 68: return 23;
          case 34: case 54: case 55: case 56: case 58: case 81: case 69: return 36;
          default: return 28;
        }
      case 8:
        switch (tx) {
          case 70: return 31;
          case 49: return 50;
          case 51: return 52;
          default: {
            for (int b : m.nb[x]) {
              if (b == a) continue;
              if (m.z[b] == 16) return 33;
              if (m.z[b] == 15) return 24;
              if (m.z[b] == 6 && carbonyl_like(b)) {
                const int p = m.double_partner(b);
                if (p >= 0 && m.z[p] == 8) return 24;  // carboxylic acid
                return 29;
              }
              if (m.z[b] == 6 && (m.count_double(b) == 1 || m.in_arom_ring(b))) return 29;  // enol, phenol, HO-C=N
            }
            return 21;
          }
        }
      default: return 5;
    }
  }

  bool run(Diagnostics& diag) {
    aromatic();
    // Carbon first (carboxylate recognition is needed for the oxygens), then the other heavy atoms.
    for (int pass = 0; pass < 2; ++pass)
      for (std::size_t a = 0; a < m.n; ++a) {
        if (t[a] != 0) continue;
        const int z = m.z[a];
        if (pass == 0 && z != 6) continue;
        if (pass == 1 && (z == 6 || z == 1)) continue;
        int ty = 0;
        switch (z) {
          case 6: ty = carbon(int(a)); break;
          case 7: ty = nitrogen(int(a)); break;
          case 8: ty = oxygen(int(a)); break;
          case 9: ty = m.deg(a) == 0 ? 89 : 11; break;
          case 17: ty = m.deg(a) == 0 ? 90 : m.terminal(int(a), 8) == 4 ? 77 : 12; break;
          case 35: ty = m.deg(a) == 0 ? 91 : 13; break;
          case 53: ty = 14; break;
          case 14: ty = 19; break;
          case 15: ty = phosphorus(int(a)); break;
          case 16: ty = sulfur(int(a)); break;
          case 3: ty = 92; break;
          case 11: ty = 93; break;
          case 19: ty = 94; break;
          case 30: ty = 95; break;
          case 20: ty = 96; break;
          case 29: ty = m.fc[a] == 1 ? 97 : 98; break;
          case 12: ty = 99; break;
          case 26: ty = m.fc[a] == 2 ? 87 : 88; break;
          default: break;
        }
        if (ty == 0) {
          diag.error(DiagCode::NoForceFieldType, "MMFF94 has no atom type for atom " + std::to_string(a + 1) +
                                                     " (atomic number " + std::to_string(z) + ")");
          return false;
        }
        t[a] = ty;
      }
    for (std::size_t a = 0; a < m.n; ++a)
      if (m.z[a] == 1) t[a] = hydrogen(int(a));
    return true;
  }
};

Mol make_mol(const ChemicalSystemState& st, std::uint32_t sp) {
  Mol m;
  const auto& s = st.species[sp];
  m.n = s.atom_count;
  m.z.resize(m.n);
  m.fc.resize(m.n);
  m.nb.resize(m.n);
  m.order.resize(m.n);
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    m.z[i] = st.atoms.element[s.atom_begin + i];
    m.fc[i] = st.atoms.formal_charge[s.atom_begin + i];
  }
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    const int i = int(st.bonds.begin[b] - s.atom_begin), j = int(st.bonds.end[b] - s.atom_begin);
    const int o = st.bonds.order[b];
    m.nb[i].push_back(j);
    m.order[i].push_back(o);
    m.nb[j].push_back(i);
    m.order[j].push_back(o);
  }
  // A charge the structure leaves implicit: four-valent N and three-valent O are cations (as MMFF's
  // own setup infers them).
  for (std::size_t a = 0; a < m.n; ++a) {
    int v = 0;
    for (int o : m.order[a]) v += o;
    if (m.fc[a] == 0 && ((m.z[a] == 7 && v == 4) || (m.z[a] == 8 && v == 3))) m.fc[a] = 1;
    // Dative representations: P with four single bonds is P+, S with three or four single bonds S+ / S2+.
    if (m.fc[a] == 0 && m.z[a] == 15 && v == 4 && m.deg(a) == 4) m.fc[a] = 1;
    if (m.fc[a] == 0 && m.z[a] == 16 && v == m.deg(a) && (v == 3 || v == 4)) m.fc[a] = v - 2;
  }
  find_rings(m);
  find_aromatic(m);
  return m;
}

// ---------------------------------------------------------------------------------------------------------
// Charges.

// MMFF formal charges: the net charge of a group (central atom and its terminal O, S or N atoms) is
// shared by its terminal atoms, so charge-separated neutral groups (nitro, N-oxide, azide, isonitrile,
// phosphoryl) carry none and a sulfonate O carries −1/3; the positive charge of an amidinium,
// guanidinium or imidazolium ion is shared by its nitrogens, and the negative charge of an azole anion by
// its ring nitrogens.
std::vector<double> shared_formal_charges(const Mol& m, const std::vector<int>& t) {
  // A sulfilimine / sulfoximine drawn S+–N− is the neutral S=N.
  std::vector<int> fc = m.fc;
  for (std::size_t a = 0; a < m.n; ++a)
    if (t[a] == 48 && fc[a] == -1)
      for (int b : m.nb[a])
        if (m.z[b] == 16 && fc[b] > 0) {
          fc[a] = 0;
          --fc[b];
          break;
        }
  std::vector<double> q0(m.n);
  for (std::size_t a = 0; a < m.n; ++a) q0[a] = fc[a];
  for (std::size_t c = 0; c < m.n; ++c) {
    std::vector<int> term;
    for (int b : m.nb[c])
      if (m.deg(b) == 1 && (t[b] == 32 || t[b] == 72 || t[b] == 47 || t[b] == 60)) term.push_back(b);
    if (term.empty()) continue;
    double total = fc[c];
    for (int b : term) total += fc[b];
    q0[c] = 0.0;
    for (int b : term) q0[b] = total / double(term.size());
  }
  for (std::size_t c = 0; c < m.n; ++c) {
    if (t[c] != 57 && t[c] != 80) continue;
    std::vector<int> g;
    for (int b : m.nb[c])
      if (t[b] == 55 || t[b] == 56 || t[b] == 81) g.push_back(b);
    if (g.empty()) continue;
    double total = 0.0;
    for (int b : g) total += fc[b];
    for (int b : g) q0[b] = total / double(g.size());
  }
  for (std::size_t r = 0; r < m.rings.size(); ++r) {
    if (!m.arom_ring[r]) continue;
    std::vector<int> g;
    for (int a : m.rings[r])
      if (t[a] == 76) g.push_back(a);
    if (g.size() < 2) continue;
    double total = 0.0;
    for (int a : g) total += fc[a];
    for (int a : g) q0[a] = total / double(g.size());
  }
  return q0;
}

int bond_type(const Mol& m, const std::vector<int>& t, int i, int j) {
  // Single bond between two atoms that can take part in single-multiple conjugation (sbmb), not within one
  // aromatic ring; or between two aromatic atoms of different rings (biphenyl).
  if (m.bo(i, j) != 1) return 0;
  const auto& pi = kMmffProp[t[i]];
  const auto& pj = kMmffProp[t[j]];
  bool same_arom = false;
  for (std::size_t r = 0; r < m.rings.size(); ++r)
    if (m.arom_ring[r]) {
      const auto& ring = m.rings[r];
      if (std::find(ring.begin(), ring.end(), i) != ring.end() && std::find(ring.begin(), ring.end(), j) != ring.end())
        same_arom = true;
    }
  if (same_arom) return 0;
  if (pi.arom && pj.arom) return 1;
  if (pi.sbmb && pj.sbmb) return 1;
  return 0;
}

std::vector<double> partial_charges(const Mol& m, const std::vector<int>& t, const std::vector<double>& q0,
                                    std::size_t& empirical) {
  const auto& T = tables();
  std::vector<double> q(m.n, 0.0);
  for (std::size_t i = 0; i < m.n; ++i) {
    const double u = kMmffPbci[t[i]].fcadj;
    const double M = kMmffProp[t[i]].crd;
    q[i] += (1.0 - M * u) * q0[i];
    for (int k : m.nb[i]) {
      const double uk = kMmffPbci[t[k]].fcadj;
      if (q0[k] < 0.0) q[i] += uk * q0[k];
    }
  }
  for (std::size_t i = 0; i < m.n; ++i)
    for (int j : m.nb[i]) {
      if (int(i) > j) continue;
      const int bt = bond_type(m, t, int(i), j);
      int a = t[i], b = t[j], lo = int(i), hi = j;
      if (a > b) std::swap(a, b), std::swap(lo, hi);
      double bci;
      auto it = T.chg.find({bt, a, b});
      if (it != T.chg.end()) bci = it->second;
      else {
        bci = kMmffPbci[b].pbci - kMmffPbci[a].pbci;
        ++empirical;
      }
      q[lo] -= bci;
      q[hi] += bci;
    }
  return q;
}

// ---------------------------------------------------------------------------------------------------------
// Parameters.

int period_row(int z) { return z <= 2 ? 0 : z <= 10 ? 1 : z <= 18 ? 2 : z <= 36 ? 3 : 4; }

const MmffCovRad* covrad(int z) {
  for (const auto& c : kMmffCovRad)
    if (c.z == z) return &c;
  return nullptr;
}

// Empirical bond: Schomaker–Stevenson r0 with Blom–Haaland corrections and Badger's rule for kb
// (Halgren, J. Comput. Chem. 17 (1996) 553, eqs. 18-20).
std::pair<double, double> empirical_bond(const Mol& m, const std::vector<int>& t, int i, int j) {
  const int zi = m.z[i], zj = m.z[j];
  const auto* ci = covrad(zi);
  const auto* cj = covrad(zj);
  double ri = ci ? ci->r : 1.0, rj = cj ? cj->r : 1.0;
  const double chi = ci && cj ? std::fabs(ci->chi - cj->chi) : 0.0;
  const int bo = m.bo(i, j);
  const auto& pi = kMmffProp[t[i]];
  const auto& pj = kMmffProp[t[j]];
  // Multiple-bond radius corrections.
  if (bo == 2 || (bo == 1 && bond_type(m, t, i, j) == 0 && pi.arom && pj.arom)) ri -= 0.1, rj -= 0.1;
  if (bo == 3) ri -= 0.17, rj -= 0.17;
  if (bo == 1 && pi.mltb == 1 && pj.mltb == 1) { ri -= 0.04; rj -= 0.04; }
  const double c = (zi == 1 || zj == 1) ? 0.05 : 0.085;
  const double r0 = ri + rj - c * std::pow(chi, 1.4);
  // Badger's rule from the reference pair, else Herschbach–Laurie.
  int a = std::min(zi, zj), b = std::max(zi, zj);
  const auto& T = tables();
  auto it = T.bndk.find({a, b});
  double kb;
  if (it != T.bndk.end()) kb = it->second.second * std::pow(it->second.first / r0, 6);
  else {
    int ra = period_row(a), rb = period_row(b);
    if (ra > rb) std::swap(ra, rb);
    double A = 1.0, D = 1.0;
    for (const auto& h : kMmffHL)
      if (h.ri == ra && h.rj == rb) A = h.a, D = h.d;
    kb = std::pow(10.0, -(r0 - A) / D);
  }
  return {kb, r0};
}

}  // namespace

std::optional<std::vector<int>> mmff_atom_types(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag) {
  const Mol m = make_mol(st, sp);
  Typer ty(m);
  if (!ty.run(diag)) return std::nullopt;
  return ty.t;
}

namespace {

// Forward-mode dual numbers over the ≤12 coordinates of one bonded term.
struct Dual {
  double v = 0.0;
  std::array<double, 12> d{};
};
inline Dual operator+(Dual a, const Dual& b) { a.v += b.v; for (int k = 0; k < 12; ++k) a.d[k] += b.d[k]; return a; }
inline Dual operator-(Dual a, const Dual& b) { a.v -= b.v; for (int k = 0; k < 12; ++k) a.d[k] -= b.d[k]; return a; }
inline Dual operator*(const Dual& a, const Dual& b) {
  Dual r; r.v = a.v * b.v; for (int k = 0; k < 12; ++k) r.d[k] = a.d[k] * b.v + a.v * b.d[k]; return r;
}
inline Dual operator*(Dual a, double s) { a.v *= s; for (double& x : a.d) x *= s; return a; }
inline Dual operator*(double s, Dual a) { return a * s; }
[[maybe_unused]] inline Dual operator+(Dual a, double s) { a.v += s; return a; }
[[maybe_unused]] inline Dual operator+(double s, Dual a) { a.v += s; return a; }
inline Dual operator-(Dual a, double s) { a.v -= s; return a; }
inline Dual operator-(double s, Dual a) { a.v = s - a.v; for (double& x : a.d) x = -x; return a; }
inline Dual operator/(const Dual& a, const Dual& b) {
  Dual r; r.v = a.v / b.v; for (int k = 0; k < 12; ++k) r.d[k] = (a.d[k] * b.v - a.v * b.d[k]) / (b.v * b.v); return r;
}
inline Dual chain(const Dual& a, double f, double df) { Dual r; r.v = f; for (int k = 0; k < 12; ++k) r.d[k] = df * a.d[k]; return r; }
inline Dual dsqrt(const Dual& a) { const double s = std::sqrt(a.v); return chain(a, s, s > 0 ? 0.5 / s : 0.0); }
inline Dual dcos(const Dual& a) { return chain(a, std::cos(a.v), -std::sin(a.v)); }
inline Dual dacos(const Dual& a) {
  const double x = std::clamp(a.v, -1.0, 1.0);
  const double s = std::sqrt(std::max(1.0 - x * x, 1e-30));
  return chain(a, std::acos(x), -1.0 / s);
}
inline Dual dasin(const Dual& a) {
  const double x = std::clamp(a.v, -1.0, 1.0);
  const double s = std::sqrt(std::max(1.0 - x * x, 1e-30));
  return chain(a, std::asin(x), 1.0 / s);
}
inline Dual datan2(const Dual& y, const Dual& x) {
  Dual r; r.v = std::atan2(y.v, x.v);
  const double den = x.v * x.v + y.v * y.v;
  for (int k = 0; k < 12; ++k) r.d[k] = den > 0 ? (x.v * y.d[k] - y.v * x.d[k]) / den : 0.0;
  return r;
}
using V3 = std::array<Dual, 3>;
V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Dual dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V3 cross(const V3& a, const V3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Dual norm(const V3& a) { return dsqrt(dot(a, a)); }

struct Seed {
  std::array<V3, 4> p;
  std::array<std::uint32_t, 4> atom{};
  int count = 0;
  Seed(const std::vector<double>& x, std::initializer_list<std::uint32_t> atoms) {
    for (std::uint32_t a : atoms) {
      for (int c = 0; c < 3; ++c) {
        Dual& q = p[count][c];
        q.v = x[3 * a + c];
        q.d[3 * count + c] = 1.0;
      }
      atom[count++] = a;
    }
  }
  void add(const Dual& e, std::vector<double>* g) const {
    if (!g)
      return;
    for (int i = 0; i < count; ++i)
      for (int c = 0; c < 3; ++c) (*g)[3 * atom[i] + c] += e.d[3 * i + c];
  }
};

Dual angle_deg(const V3& a, const V3& j, const V3& b) {
  const V3 u = sub(a, j), v = sub(b, j);
  return dacos(dot(u, v) / (norm(u) * norm(v))) * kDeg;
}

}  // namespace

std::optional<MmffForceField> MmffForceField::build(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag,
                                                    MmffVariant variant, double dielectric) {
  const Mol m = make_mol(st, sp);
  for (std::size_t a = 0; a < m.n; ++a)
    if (m.z[a] != 1) {
      // MMFF needs every hydrogen explicit.
      const std::uint32_t g = st.species[sp].atom_begin + std::uint32_t(a);
      if (st.atoms.implicit_h[g] != 0) {
        diag.error(DiagCode::ImplicitHydrogens, "MMFF94 needs explicit hydrogens (add them first)");
        return std::nullopt;
      }
    }
  Typer ty(m);
  if (!ty.run(diag)) return std::nullopt;
  MmffForceField ff;
  const auto& t = ty.t;
  ff.types_ = t;
  for (int x : t) ff.symbols_.push_back(kMmffDef[x].symbol);
  ff.q0_ = shared_formal_charges(m, t);
  ff.charges_ = partial_charges(m, t, ff.q0_, ff.empirical_);
  const auto& T = tables();
  const bool s94 = variant == MmffVariant::Mmff94s;

  // Bonds.
  std::map<std::pair<int, int>, double> r0_of;
  for (std::size_t i = 0; i < m.n; ++i)
    for (int j : m.nb[i]) {
      if (int(i) > j) continue;
      const int bt = bond_type(m, t, int(i), j);
      const int a = std::min(t[i], t[j]), b = std::max(t[i], t[j]);
      auto it = T.bond.find({bt, a, b});
      double kb, r0;
      if (it != T.bond.end()) kb = it->second.first, r0 = it->second.second;
      else {
        std::tie(kb, r0) = empirical_bond(m, t, int(i), j);
        ++ff.empirical_;
      }
      ff.bonds_.push_back({std::uint32_t(i), std::uint32_t(j), kb, r0});
      r0_of[{int(i), j}] = r0_of[{j, int(i)}] = r0;
    }

  // Angles and stretch-bends.
  for (std::size_t j = 0; j < m.n; ++j) {
    const auto& nbj = m.nb[j];
    for (std::size_t x = 0; x < nbj.size(); ++x)
      for (std::size_t y = x + 1; y < nbj.size(); ++y) {
        int i = nbj[x], k = nbj[y];
        if (t[i] > t[k]) std::swap(i, k);
        const int bti = bond_type(m, t, i, int(j)), btk = bond_type(m, t, int(j), k);
        const int sum = bti + btk;
        int at = sum;
        if (m.bonded(i, k)) at = sum == 0 ? 3 : sum == 1 ? 5 : 6;
        else {
          bool four = false;
          for (int c : m.nb[i])
            if (c != int(j) && c != k && m.bonded(c, k)) four = true;  // i-j-k-c ring (fused rings too)
          if (four) at = sum == 0 ? 4 : sum == 1 ? 7 : 8;
        }
        double ka = 0.0, th0 = 0.0;
        bool found = false;
        for (int lv = 0; lv < 4 && !found; ++lv) {
          int a = level(t[i], lv), c = level(t[k], lv);
          if (a > c) std::swap(a, c);
          auto it = T.angle.find({at, a, int(t[j]), c});
          if (it != T.angle.end()) {
            ka = it->second.first, th0 = it->second.second;
            found = true;
          }
        }
        if (!found || ka == 0.0) {
          // Empirical rule (Halgren 1996 III): ka from the bond lengths and θ0.
          ++ff.empirical_;
          if (!found) {
            const int crd = kMmffProp[t[j]].crd;
            th0 = kMmffProp[t[j]].lin ? 180.0 : crd == 3 ? 120.0 : crd == 2 && m.z[j] == 8 ? 105.0 : 109.47;
            if (at == 3 || at == 5 || at == 6) th0 = 60.0;
            if (at == 4 || at == 7 || at == 8) th0 = 90.0;
          }
          static const std::map<int, std::pair<double, double>> ZC = {
              {1, {1.395, 0.0}}, {6, {2.494, 1.016}}, {7, {2.711, 1.113}}, {8, {3.045, 1.337}}, {9, {2.847, 0.0}},
              {14, {2.350, 0.811}}, {15, {2.350, 1.068}}, {16, {2.980, 1.249}}, {17, {2.909, 1.078}},
              {35, {3.017, 0.0}}, {53, {3.086, 0.0}}};
          auto zi = ZC.find(m.z[i]), zj = ZC.find(m.z[j]), zk = ZC.find(m.z[k]);
          const double rij = r0_of[{i, int(j)}], rkj = r0_of[{k, int(j)}];
          const double D = (rij - rkj) * (rij - rkj) / ((rij + rkj) * (rij + rkj));
          double beta = 1.75;
          if (at == 4 || at == 7 || at == 8) beta *= 0.85;
          if (at == 3 || at == 5 || at == 6) beta *= 0.05;
          const double th = th0 / kDeg;
          if (zi != ZC.end() && zj != ZC.end() && zk != ZC.end())
            ka = beta * zi->second.first * zj->second.second * zk->second.first / ((rij + rkj) * th * th) * std::exp(-2.0 * D);
        }
        const bool lin = kMmffProp[t[j]].lin != 0;
        ff.angles_.push_back({std::uint32_t(i), std::uint32_t(j), std::uint32_t(k), ka, th0, lin});
        if (lin) continue;
        // Stretch-bend type from the angle type and which bond is a single bond between sp2 atoms.
        int sbt = 0;
        bool swap = false;
        switch (at) {
          case 1: sbt = bti ? 1 : 2; break;
          case 2: sbt = 3; break;
          case 4: sbt = 4; break;
          case 3: sbt = 5; break;
          case 5: sbt = bti ? 6 : 7; break;
          case 6: sbt = 8; break;
          case 7: sbt = bti ? 9 : 10; break;
          case 8: sbt = 11; break;
          default: sbt = 0;
        }
        int a = t[i], c = t[k];
        if (a == c && (sbt == 2 || sbt == 7 || sbt == 10)) {  // equal outer types: put the sp2 bond first
          swap = true;
          sbt = sbt == 2 ? 1 : sbt == 7 ? 6 : 9;
        }
        double kijk = 0.0, kkji = 0.0;
        auto it = T.stbn.find({sbt, a, int(t[j]), c});
        if (it != T.stbn.end()) kijk = it->second.first, kkji = it->second.second;
        else {
          int ri = period_row(m.z[i]), rj = period_row(m.z[j]), rk = period_row(m.z[k]);
          bool rev = ri > rk;
          if (rev) std::swap(ri, rk);
          for (const auto& f : kMmffDfsb)
            if (f.ri == ri && f.rj == rj && f.rk == rk) {
              kijk = rev ? f.f2 : f.f1;
              kkji = rev ? f.f1 : f.f2;
            }
        }
        int ii = i, kk = k;
        if (swap) std::swap(ii, kk);
        ff.stbn_.push_back({std::uint32_t(ii), std::uint32_t(j), std::uint32_t(kk), kijk, kkji, r0_of[{ii, int(j)}],
                            r0_of[{kk, int(j)}], th0});
      }
  }

  // Out-of-plane bending at tricoordinate centres.
  const auto& oopt = s94 ? T.oop94s : T.oop94;
  for (std::size_t j = 0; j < m.n; ++j) {
    if (m.nb[j].size() != 3) continue;
    double koop = 0.0;
    bool found = false;
    for (int lv = 0; lv < 4 && !found; ++lv) {
      std::array<int, 3> o = {level(t[m.nb[j][0]], lv), level(t[m.nb[j][1]], lv), level(t[m.nb[j][2]], lv)};
      std::sort(o.begin(), o.end());
      auto it = oopt.find({o[0], int(t[j]), o[1], o[2]});
      if (it != oopt.end()) koop = it->second, found = true;
    }
    if (!found || koop == 0.0) continue;
    const auto& n = m.nb[j];
    ff.oop_.push_back({std::uint32_t(n[0]), std::uint32_t(j), std::uint32_t(n[1]), std::uint32_t(n[2]), koop});
    ff.oop_.push_back({std::uint32_t(n[0]), std::uint32_t(j), std::uint32_t(n[2]), std::uint32_t(n[1]), koop});
    ff.oop_.push_back({std::uint32_t(n[1]), std::uint32_t(j), std::uint32_t(n[2]), std::uint32_t(n[0]), koop});
  }

  // Torsions.
  const auto& tort = s94 ? T.tor94s : T.tor94;
  for (std::size_t j = 0; j < m.n; ++j)
    for (int k : m.nb[j]) {
      if (int(j) > k) continue;
      if (kMmffProp[t[j]].lin || kMmffProp[t[k]].lin) continue;
      for (int i : m.nb[j]) {
        if (i == k) continue;
        for (int l : m.nb[k]) {
          if (l == int(j) || l == i) continue;
          // Torsion type.
          // Torsion type: 1 for a conjugated single central bond, 2 for a single central bond next to one;
          // otherwise 4 in a four-ring, 5 in a five-ring with a saturated carbon. (Bond-type classes take
          // precedence over the five-ring class, as in the validation suite.)
          const int btjk = bond_type(m, t, int(j), k);
          int tt = btjk ? 1 : 0;
          if (!btjk && m.bo(int(j), k) == 1 && (bond_type(m, t, i, int(j)) || bond_type(m, t, k, l))) tt = 2;
          if (m.same_ring(4, {i, int(j), k, l})) tt = 4;
          else if (tt == 0 && m.same_ring(5, {i, int(j), k, l})) {
            bool sat = false;
            for (int a : {i, int(j), k, l}) sat = sat || t[a] == 1;
            if (sat) tt = 5;
          }
          // Canonical order: j < k, or j == k and i <= l.
          int I = t[i], J = t[j], K = t[k], L = t[l];
          if (J > K || (J == K && I > L)) std::swap(I, L), std::swap(J, K);
          std::array<double, 3> V{0, 0, 0};
          bool found = false;
          // Step-down (MMFFDEF levels; 3 is the wildcard): exact; one outer atom at level 1 and the other a
          // wildcard; both wildcards. (Reproduces 41,127 of the 41,133 torsions of the MMFF94 validation
          // suite; the other six need the empirical rule.)
          const int seq[][2] = {{0, 0}, {1, 3}, {3, 1}, {3, 3}};
          for (int ttry : {tt, tt == 0 ? -1 : 0}) {
            if (ttry < 0) break;
            for (const auto& s : seq) {
              int a = level(I, s[0]), d = level(L, s[1]);
              int jj = J, kk = K;
              if (jj == kk && a > d) std::swap(a, d);
              auto it = tort.find({ttry, a, jj, kk, d});
              if (it != tort.end()) {
                V = it->second;
                found = true;
                break;
              }
            }
            if (found) break;
          }
          if (!found) {
            // Empirical rule (Halgren 1996 V) for the central atoms j, k: double and aromatic bonds get a
            // twofold barrier from U, saturated (or pyramidal, non-conjugated) atoms a threefold one from V.
            ++ff.empirical_;
            auto U = [&](int a) { return period_row(m.z[a]) == 1 ? 2.0 : 1.25; };
            auto Vz = [&](int a) {
              switch (m.z[a]) {
                case 6: return 2.12; case 7: return 1.5; case 8: return 0.2;
                case 14: return 1.22; case 15: return 2.4; case 16: return 0.49; default: return 0.0;
              }
            };
            const auto& pj = kMmffProp[t[j]];
            const auto& pk = kMmffProp[t[k]];
            if (pj.arom && pk.arom && !btjk) {
              const double beta = ((pj.val == 3 && pk.val == 4) || (pj.val == 4 && pk.val == 3)) ? 3.0 : 6.0;
              const double pi = (pj.pilp == 0 && pk.pilp == 0) ? 0.5 : 0.3;
              V[1] = beta * pi * std::sqrt(U(int(j)) * U(k));
            } else if (m.bo(int(j), k) == 2) {
              const double pi = (pj.mltb == 2 && pk.mltb == 2) ? 1.0 : 0.4;
              V[1] = 6.0 * pi * std::sqrt(U(int(j)) * U(k));
            } else if (!pj.mltb && !pk.mltb && !pj.arom && !pk.arom && pj.crd > 1 && pk.crd > 1) {
              V[2] = std::sqrt(Vz(int(j)) * Vz(k)) / double((pj.crd - 1) * (pk.crd - 1));
            }
          }
          if (V[0] == 0 && V[1] == 0 && V[2] == 0) continue;
          ff.torsions_.push_back({std::uint32_t(i), std::uint32_t(j), std::uint32_t(k), std::uint32_t(l), V[0], V[1], V[2]});
        }
      }
    }

  // Non-bonded pairs separated by three or more bonds (1-4 electrostatics scaled by 0.75).
  std::vector<std::vector<int>> dist(m.n);
  for (std::size_t a = 0; a < m.n; ++a) {
    std::vector<int> d(m.n, 99);
    d[a] = 0;
    std::vector<int> frontier = {int(a)};
    for (int step = 1; step <= 3; ++step) {
      std::vector<int> next;
      for (int u : frontier)
        for (int v : m.nb[u])
          if (d[v] > step) d[v] = step, next.push_back(v);
      frontier = next;
    }
    dist[a] = d;
  }
  for (std::size_t i = 0; i < m.n; ++i)
    for (std::size_t j = i + 1; j < m.n; ++j) {
      if (dist[i][j] < 3) continue;
      const auto vi = T.vdw.find(t[i]), vj = T.vdw.find(t[j]);
      double R = 0.0, eps = 0.0;
      if (vi != T.vdw.end() && vj != T.vdw.end()) {
        const auto &a = vi->second, &b = vj->second;
        const double Ri = a.A * std::pow(a.alpha, kMmffVdwPower), Rj = b.A * std::pow(b.alpha, kMmffVdwPower);
        const double g = (Ri - Rj) / (Ri + Rj);
        const bool donor = a.da == 1 || b.da == 1;
        R = 0.5 * (Ri + Rj) * (donor ? 1.0 : 1.0 + kMmffVdwB * (1.0 - std::exp(-kMmffVdwBeta * g * g)));
        eps = 181.16 * a.G * b.G * a.alpha * b.alpha / (std::sqrt(a.alpha / a.N) + std::sqrt(b.alpha / b.N)) / std::pow(R, 6);
        if ((a.da == 1 && b.da == 2) || (a.da == 2 && b.da == 1)) R *= kMmffDARAD, eps *= kMmffDAEPS;
      }
      const double scale = dist[i][j] == 3 ? 0.75 : 1.0;
      const double qq = 332.0716 * ff.charges_[i] * ff.charges_[j] * scale / dielectric;
      ff.pairs_.push_back({std::uint32_t(i), std::uint32_t(j), R, eps, qq});
    }
  return ff;
}

double MmffForceField::eval(const std::vector<double>& x, std::vector<double>* g, MmffEnergyBreakdown* parts) const {
  if (g) g->assign(x.size(), 0.0);
  MmffEnergyBreakdown e;
  for (const auto& b : bonds_) {
    Seed s(x, {b.i, b.j});
    const Dual dr = norm(sub(s.p[0], s.p[1])) - b.r0;
    const double cs = -2.0;
    const Dual E = 0.5 * 143.9325 * b.kb * dr * dr * (1.0 + cs * dr + (7.0 / 12.0) * cs * cs * dr * dr);
    e.bond += E.v;
    s.add(E, g);
  }
  for (const auto& a : angles_) {
    Seed s(x, {a.i, a.j, a.k});
    Dual E;
    if (a.linear) {
      const V3 u = sub(s.p[0], s.p[1]), v = sub(s.p[2], s.p[1]);
      E = 143.9325 * a.ka * (1.0 + dot(u, v) / (norm(u) * norm(v)));
    } else {
      const Dual dth = angle_deg(s.p[0], s.p[1], s.p[2]) - a.theta0;
      E = 0.5 * kAngleK * a.ka * dth * dth * (1.0 - 0.006981317 * dth);
    }
    e.angle += E.v;
    s.add(E, g);
  }
  for (const auto& b : stbn_) {
    Seed s(x, {b.i, b.j, b.k});
    const Dual dri = norm(sub(s.p[0], s.p[1])) - b.r0ij, drk = norm(sub(s.p[2], s.p[1])) - b.r0kj;
    const Dual dth = angle_deg(s.p[0], s.p[1], s.p[2]) - b.theta0;
    const Dual E = 2.51210 * (b.kijk * dri + b.kkji * drk) * dth;
    e.stretch_bend += E.v;
    s.add(E, g);
  }
  for (const auto& o : oop_) {
    Seed s(x, {o.i, o.j, o.k, o.l});
    const V3 ji = sub(s.p[0], s.p[1]), jk = sub(s.p[2], s.p[1]), jl = sub(s.p[3], s.p[1]);
    const V3 n = cross(ji, jk);
    const Dual chi = dasin(dot(n, jl) / (norm(n) * norm(jl))) * kDeg;
    const Dual E = 0.5 * kAngleK * o.koop * chi * chi;
    e.out_of_plane += E.v;
    s.add(E, g);
  }
  for (const auto& tq : torsions_) {
    Seed s(x, {tq.i, tq.j, tq.k, tq.l});
    const V3 b1 = sub(s.p[1], s.p[0]), b2 = sub(s.p[2], s.p[1]), b3 = sub(s.p[3], s.p[2]);
    const V3 n1 = cross(b1, b2), n2 = cross(b2, b3);
    const Dual y = dot(cross(n1, n2), b2) / norm(b2), xx = dot(n1, n2);
    const Dual phi = datan2(y, xx);
    const Dual E = 0.5 * (tq.V1 * (1.0 + dcos(phi)) + tq.V2 * (1.0 - dcos(2.0 * phi)) + tq.V3 * (1.0 + dcos(3.0 * phi)));
    e.torsion += E.v;
    s.add(E, g);
  }
  for (const auto& p : pairs_) {
    double d[3], r2 = 0.0;
    for (int c = 0; c < 3; ++c) d[c] = x[3 * p.i + c] - x[3 * p.j + c], r2 += d[c] * d[c];
    const double r = std::sqrt(r2);
    double dEdr = 0.0;
    if (p.eps != 0.0) {
      const double R = p.R, R7 = std::pow(R, 7), r7 = std::pow(r, 7);
      const double A = 1.07 * R / (r + 0.07 * R), A7 = std::pow(A, 7);
      const double Bt = 1.12 * R7 / (r7 + 0.12 * R7);
      const double E = p.eps * A7 * (Bt - 2.0);
      e.vdw += E;
      const double dA7 = -7.0 * A7 / (r + 0.07 * R);
      const double dB = -1.12 * R7 * 7.0 * std::pow(r, 6) / ((r7 + 0.12 * R7) * (r7 + 0.12 * R7));
      dEdr += p.eps * (dA7 * (Bt - 2.0) + A7 * dB);
    }
    if (p.qq != 0.0) {
      const double E = p.qq / (r + 0.05);
      e.electrostatic += E;
      dEdr += -p.qq / ((r + 0.05) * (r + 0.05));
    }
    if (g && r > 0)
      for (int c = 0; c < 3; ++c) {
        (*g)[3 * p.i + c] += dEdr * d[c] / r;
        (*g)[3 * p.j + c] -= dEdr * d[c] / r;
      }
  }
  if (parts) *parts = e;
  return e.total();
}

double MmffForceField::energy(const std::vector<double>& xyz) const { return eval(xyz, nullptr, nullptr); }
double MmffForceField::energy_gradient(const std::vector<double>& xyz, std::vector<double>& grad) const {
  return eval(xyz, &grad, nullptr);
}
MmffEnergyBreakdown MmffForceField::breakdown(const std::vector<double>& xyz) const {
  MmffEnergyBreakdown e;
  eval(xyz, nullptr, &e);
  return e;
}

MmffMinimizeResult mmff_minimize(ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag, MmffVariant variant,
                                 const LbfgsOptions& options) {
  MmffMinimizeResult res;
  auto ff = MmffForceField::build(st, sp, diag, variant);
  if (!ff) return res;
  const Species& s = st.species[sp];
  std::vector<double> x(3 * s.atom_count);
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    x[3 * i] = st.atoms.x[s.atom_begin + i];
    x[3 * i + 1] = st.atoms.y[s.atom_begin + i];
    x[3 * i + 2] = st.atoms.z[s.atom_begin + i];
  }
  res.energy_initial = ff->energy(x);
  auto r = minimize_lbfgs([&](const std::vector<double>& xx, std::vector<double>& gg) { return ff->energy_gradient(xx, gg); },
                          x, options);
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    st.atoms.x[s.atom_begin + i] = x[3 * i];
    st.atoms.y[s.atom_begin + i] = x[3 * i + 1];
    st.atoms.z[s.atom_begin + i] = x[3 * i + 2];
  }
  ++st.revision;
  res.ok = true;
  res.converged = r.converged;
  res.iterations = r.iterations;
  res.energy_final = r.f;
  res.grad_max = r.grad_max;
  if (!r.converged)
    diag.warn(DiagCode::MinimizationNotConverged, "MMFF94 minimisation stopped after " + std::to_string(r.iterations) + " iterations");
  return res;
}

std::string mmff_json(std::string_view smiles, MmffVariant variant, bool minimize) {
  Diagnostics d;
  ChemicalSystemState st;
  const std::int32_t sp = embed_smiles(st, smiles, d);
  std::optional<MmffForceField> ff;
  MmffMinimizeResult mr;
  if (sp >= 0) {
    const auto s = std::uint32_t(sp);
    if (minimize) mr = mmff_minimize(st, s, d, variant, LbfgsOptions{5000, 10, 1e-4, 1e-15, 0.3});
    ff = MmffForceField::build(st, s, d, variant);
  }
  auto jn = [](double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.10g", v);
    return std::isfinite(v) ? std::string(b) : std::string("null");
  };
  std::string o = std::string("{\"ok\":") + (ff && !d.has_errors() ? "true" : "false");
  if (ff) {
    const auto s = std::uint32_t(sp);
    const auto& spc = st.species[s];
    std::vector<double> x;
    o += ",\"variant\":" + json_escape(variant == MmffVariant::Mmff94s ? "MMFF94s" : "MMFF94") + ",\"elements\":[";
    for (std::uint32_t i = 0; i < spc.atom_count; ++i) {
      const auto* e = find_element_by_z(st.atoms.element[spc.atom_begin + i]);
      o += (i ? "," : "") + json_escape(e ? std::string(e->symbol) : "X");
      x.push_back(st.atoms.x[spc.atom_begin + i]);
      x.push_back(st.atoms.y[spc.atom_begin + i]);
      x.push_back(st.atoms.z[spc.atom_begin + i]);
    }
    o += "],\"bonds\":[";
    for (std::uint32_t b = spc.bond_begin; b < spc.bond_begin + spc.bond_count; ++b)
      o += (b > spc.bond_begin ? ",[" : "[") + std::to_string(st.bonds.begin[b] - spc.atom_begin) + "," +
           std::to_string(st.bonds.end[b] - spc.atom_begin) + "," + std::to_string(int(st.bonds.order[b])) + "]";
    auto ints = [&](const std::vector<int>& v) {
      std::string r = "[";
      for (std::size_t i = 0; i < v.size(); ++i) r += (i ? "," : "") + std::to_string(v[i]);
      return r + "]";
    };
    auto dbls = [&](const std::vector<double>& v) {
      std::string r = "[";
      for (std::size_t i = 0; i < v.size(); ++i) r += (i ? "," : "") + jn(v[i]);
      return r + "]";
    };
    std::string sy = "[";
    for (std::size_t i = 0; i < ff->symbols().size(); ++i) sy += (i ? "," : "") + json_escape(ff->symbols()[i]);
    sy += "]";
    const auto e = ff->breakdown(x);
    o += "],\"types\":" + ints(ff->types()) + ",\"symbols\":" + sy + ",\"charges\":" + dbls(ff->charges()) +
         ",\"formal_charges\":" + dbls(ff->formal_charges()) + ",\"energy\":{\"total\":" + jn(e.total()) +
         ",\"bond\":" + jn(e.bond) + ",\"angle\":" + jn(e.angle) + ",\"stretch_bend\":" + jn(e.stretch_bend) +
         ",\"out_of_plane\":" + jn(e.out_of_plane) + ",\"torsion\":" + jn(e.torsion) + ",\"vdw\":" + jn(e.vdw) +
         ",\"electrostatic\":" + jn(e.electrostatic) + "},\"minimized\":" + (minimize ? "true" : "false") +
         ",\"converged\":" + (mr.converged ? "true" : "false") + ",\"iterations\":" + std::to_string(mr.iterations) +
         ",\"empirical_parameters\":" + std::to_string(ff->empirical_parameters()) +
         ",\"sdf\":" + json_escape(write_sdf(st, s, std::string(smiles))) + ",\"xyz\":" + json_escape(write_xyz(st, s, std::string(smiles)));
  }
  return o + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp

