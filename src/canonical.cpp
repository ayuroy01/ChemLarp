#include "omnichem/canonical.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <utility>

#include "omnichem/cip.hpp"
#include "omnichem/elements.hpp"

namespace omnichem {
namespace {

constexpr int kTokH = -2;
constexpr int kTokLP = -3;

struct Edge {
  int to;
  std::uint32_t bond;
};

struct View {
  std::uint32_t base = 0;
  std::vector<std::uint32_t> atoms;  // local -> global
  std::vector<int> local;            // (global - base) -> local, -1 if folded
  std::vector<int> h;                // total hydrogens (implicit + folded)
  std::vector<std::vector<Edge>> adj;
};

bool foldable_h(const ChemicalSystemState& st, std::uint32_t g) {
  if (st.atoms.element[g] != 1 || st.atoms.formal_charge[g] != 0 || st.atoms.isotope[g] != 0 ||
      st.atoms.atom_map[g] != 0 || st.atoms.implicit_h[g] != 0)
    return false;
  auto nb = st.adjacency.neighbors(g);
  return nb.size() == 1 && st.atoms.element[nb[0]] != 1 && st.bonds.order[st.adjacency.bonds_of(g)[0]] == 1;
}

View build_view(const ChemicalSystemState& st, const Species& s) {
  View v;
  v.base = s.atom_begin;
  v.local.assign(s.atom_count, -1);
  for (std::uint32_t g = s.atom_begin; g < s.atom_begin + s.atom_count; ++g) {
    if (foldable_h(st, g)) continue;
    v.local[g - v.base] = static_cast<int>(v.atoms.size());
    v.atoms.push_back(g);
    v.h.push_back(st.atoms.implicit_h[g]);
  }
  for (std::uint32_t g = s.atom_begin; g < s.atom_begin + s.atom_count; ++g)
    if (v.local[g - v.base] < 0) ++v.h[static_cast<std::size_t>(v.local[st.adjacency.neighbors(g)[0] - v.base])];
  v.adj.resize(v.atoms.size());
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    const int x = v.local[st.bonds.begin[b] - v.base], y = v.local[st.bonds.end[b] - v.base];
    if (x < 0 || y < 0) continue;
    v.adj[static_cast<std::size_t>(x)].push_back({y, b});
    v.adj[static_cast<std::size_t>(y)].push_back({x, b});
  }
  return v;
}

int bond_code(const ChemicalSystemState& st, std::uint32_t b) {
  if (st.bonds.flags[b] & kBondAromatic) return 4;
  const int o = st.bonds.order[b];
  return o == 4 ? 5 : o;
}

template <class Key>
std::vector<long long> dense_rank(const std::vector<Key>& keys) {
  std::vector<std::size_t> idx(keys.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
  std::vector<long long> r(keys.size());
  long long cur = 0;
  for (std::size_t i = 0; i < idx.size(); ++i) {
    if (i > 0 && keys[idx[i - 1]] < keys[idx[i]]) ++cur;
    r[idx[i]] = cur;
  }
  return r;
}

std::size_t distinct(const std::vector<long long>& r) {
  return r.empty() ? 0 : static_cast<std::size_t>(*std::max_element(r.begin(), r.end()) + 1);
}

void refine(const ChemicalSystemState& st, const View& v, std::vector<long long>& rank) {
  for (;;) {
    std::vector<std::vector<long long>> keys(rank.size());
    for (std::size_t i = 0; i < rank.size(); ++i) {
      keys[i].push_back(rank[i]);
      std::vector<long long> nb;
      for (const Edge& e : v.adj[i]) nb.push_back(rank[static_cast<std::size_t>(e.to)] * 8 + bond_code(st, e.bond));
      std::sort(nb.begin(), nb.end());
      keys[i].insert(keys[i].end(), nb.begin(), nb.end());
    }
    auto next = dense_rank(keys);
    const bool stable = distinct(next) == distinct(rank);
    rank = std::move(next);
    if (stable) return;
  }
}

// Atom invariants refined to a stable partition; ties are not broken.
std::vector<long long> initial_ranks(const ChemicalSystemState& st, const View& v, const std::vector<int>& stereo_inv) {
  const std::size_t n = v.atoms.size();
  std::vector<std::array<long long, 9>> init(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::uint32_t g = v.atoms[i];
    init[i] = {static_cast<long long>(v.adj[i].size()),
               st.atoms.element[g],
               st.atoms.isotope[g],
               st.atoms.formal_charge[g] + 128,
               v.h[i],
               (st.atoms.flags[g] & kAtomInRing) ? 1 : 0,
               (st.atoms.flags[g] & kAtomAromatic) ? 1 : 0,
               stereo_inv[i],
               st.atoms.atom_map[g]};
  }
  auto rank = dense_rank(init);
  refine(st, v, rank);
  return rank;
}

// The lowest rank shared by two or more atoms (-1 when all ranks differ).
long long lowest_tie(const std::vector<long long>& rank) {
  std::vector<std::size_t> count(rank.size(), 0);
  for (long long r : rank) ++count[static_cast<std::size_t>(r)];
  for (std::size_t r = 0; r < count.size(); ++r)
    if (count[r] > 1) return static_cast<long long>(r);
  return -1;
}

// Gives atom `chosen` a rank of its own just below the rest of its class.
std::vector<long long> split(const std::vector<long long>& rank, std::size_t chosen) {
  std::vector<long long> r = rank;
  for (auto& x : r) x *= 2;
  r[chosen] -= 1;
  return dense_rank(r);
}

std::vector<long long> compute_ranks(const ChemicalSystemState& st, const View& v, const std::vector<int>& stereo_inv) {
  auto rank = initial_ranks(st, v, stereo_inv);
  for (long long tied; (tied = lowest_tie(rank)) >= 0;) {
    std::size_t chosen = 0;
    while (rank[chosen] != tied) ++chosen;
    rank = split(rank, chosen);
    refine(st, v, rank);
  }
  return rank;
}

struct ActiveDb {
  std::uint32_t bond;
  int u, v;      // local ends
  int ra, rb;    // local refs (-1 if a folded H)
  bool cis;
};

class Writer {
 public:
  Writer(const ChemicalSystemState& st, const Species& s, const SmilesWriteOptions& opt, Diagnostics* diag)
      : st_(st), sp_(s), opt_(opt), diag_(diag), v_(build_view(st, s)) {}

  std::string run(std::uint32_t species_index) {
    const std::size_t n = v_.atoms.size();
    collect_stereo(species_index);
    if (opt_.canonical) {
      std::vector<int> inv(n, 0);
      for (const auto& [loc, label] : centre_label_) inv[static_cast<std::size_t>(loc)] = static_cast<int>(label);
      for (const ActiveDb& d : dbs_) {
        inv[static_cast<std::size_t>(d.u)] += 8 * static_cast<int>(db_label_[d.bond]);
        inv[static_cast<std::size_t>(d.v)] += 8 * static_cast<int>(db_label_[d.bond]);
      }
      rank_ = stereo_ranks(inv);
    } else {
      rank_.resize(n);
      std::iota(rank_.begin(), rank_.end(), 0);
      std::mt19937_64 rng(opt_.random_seed);
      std::shuffle(rank_.begin(), rank_.end(), rng);
    }

    pre_.assign(n, -1);
    parent_.assign(n, -1);
    parent_bond_.assign(n, 0);
    children_.assign(n, {});
    ring_at_.assign(n, {});
    is_ring_.assign(sp_.bond_count, 0);
    std::vector<int> roots(n);
    std::iota(roots.begin(), roots.end(), 0);
    std::sort(roots.begin(), roots.end(), [&](int a, int b) { return rank_[static_cast<std::size_t>(a)] < rank_[static_cast<std::size_t>(b)]; });
    std::vector<int> starts;
    for (int r : roots)
      if (pre_[static_cast<std::size_t>(r)] < 0) { starts.push_back(r); dfs(r); }

    plan_double_bonds();
    digit_.assign(sp_.bond_count, -1);
    std::string out;
    for (std::size_t i = 0; i < starts.size(); ++i) {
      if (i) out += '.';
      write(starts[i], out);
    }
    return out;
  }

 private:
  bool out_aromatic(int x) const {
    return !opt_.kekule && (st_.atoms.flags[v_.atoms[static_cast<std::size_t>(x)]] & kAtomAromatic);
  }
  bool out_aromatic_bond(std::uint32_t b) const { return !opt_.kekule && (st_.bonds.flags[b] & kBondAromatic); }
  int other(std::uint32_t b, int x) const {
    const int a = v_.local[st_.bonds.begin[b] - v_.base];
    return a == x ? v_.local[st_.bonds.end[b] - v_.base] : a;
  }
  int to_token(std::uint32_t ref) const {
    if (ref == kImplicitHRef) return kTokH;
    if (ref == kLonePairRef) return kTokLP;
    const int l = v_.local[ref - v_.base];
    return l < 0 ? kTokH : l;
  }
  std::int64_t bond_between(int x, int y) const {
    for (const Edge& e : v_.adj[static_cast<std::size_t>(x)])
      if (e.to == y) return e.bond;
    return -1;
  }

  // Parity of each stereo element relative to the current ranks: 1 or 2 once
  // the neighbours that decide it have distinct ranks, else 0. It describes the
  // molecule, not the input order, so it can refine ranks like any invariant.
  std::vector<int> parity_invariants(const std::vector<long long>& rank, bool* undecided = nullptr) const {
    std::vector<int> out(rank.size(), 0);
    if (undecided) *undecided = false;
    auto key = [&](int tok) -> long long { return tok == kTokH ? -2 : tok == kTokLP ? -1 : rank[static_cast<std::size_t>(tok)]; };
    for (const Tetra& t : tetra_) {
      std::array<long long, 4> k{};
      for (std::size_t i = 0; i < 4; ++i) k[i] = key(t.refs[i]);
      int inversions = 0;
      bool tie = false;
      for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = i + 1; j < 4; ++j) {
          tie = tie || k[i] == k[j];
          inversions += k[i] > k[j];
        }
      if (tie) {
        if (undecided) *undecided = true;
        continue;
      }
      const bool cw = t.parity == Chirality::Clockwise;
      out[static_cast<std::size_t>(t.loc)] = ((inversions & 1) != 0) != cw ? 2 : 1;
    }
    for (const ActiveDb& d : dbs_) {
      // Is the stored reference the lowest-ranked substituent on its end? (-1 if undecided)
      auto ref_is_lowest = [&](int end, int partner, int ref) -> int {
        std::vector<long long> keys;
        for (const Edge& e : v_.adj[static_cast<std::size_t>(end)])
          if (e.to != partner) keys.push_back(rank[static_cast<std::size_t>(e.to)]);
        if (v_.h[static_cast<std::size_t>(end)] > 0) keys.push_back(-2);
        std::sort(keys.begin(), keys.end());
        if (keys.size() == 2 && keys[0] == keys[1]) return -1;
        const long long rk = ref < 0 ? -2 : rank[static_cast<std::size_t>(ref)];
        return rk == keys.front() ? 1 : 0;
      };
      const int a = ref_is_lowest(d.u, d.v, d.ra), b = ref_is_lowest(d.v, d.u, d.rb);
      if (a < 0 || b < 0) {
        if (undecided) *undecided = true;
        continue;
      }
      const bool cis = d.cis != (a == 0) != (b == 0);
      out[static_cast<std::size_t>(d.u)] += cis ? 4 : 8;
      out[static_cast<std::size_t>(d.v)] += cis ? 4 : 8;
    }
    return out;
  }

  // Refines with the parity invariants until the partition is stable.
  void settle(std::vector<long long>& rank) const {
    for (;;) {
      refine(st_, v_, rank);
      const std::vector<int> p = parity_invariants(rank);
      std::vector<std::pair<long long, int>> keys(rank.size());
      for (std::size_t i = 0; i < rank.size(); ++i) keys[i] = {rank[i], p[i]};
      auto next = dense_rank(keys);
      if (distinct(next) == distinct(rank)) return;
      rank = std::move(next);
    }
  }

  // Canonical ranks. Without stereo, ties between symmetric atoms are broken by
  // index, which cannot change the string. With stereo it can (mirror-image
  // spellings of an achiral molecule pick opposite branches), so each candidate
  // of the tied class is tried and the one whose stereo parities, read in rank
  // order, are largest is kept. This depends only on the molecule.
  std::vector<long long> stereo_ranks(const std::vector<int>& cip_inv) const {
    if (tetra_.empty() && dbs_.empty()) return compute_ranks(st_, v_, cip_inv);
    auto rank = initial_ranks(st_, v_, cip_inv);
    settle(rank);
    std::size_t budget = 4096;  // candidate evaluations before falling back to index order
    for (long long tied; (tied = lowest_tie(rank)) >= 0;) {
      // Once every stereo element is decided by distinct ranks, later splits only
      // refine and cannot reorder any element's neighbours: break ties by index.
      bool undecided = false;
      parity_invariants(rank, &undecided);
      if (!undecided) {
        std::size_t chosen = 0;
        while (rank[chosen] != tied) ++chosen;
        rank = split(rank, chosen);
        refine(st_, v_, rank);
        continue;
      }
      std::vector<long long> best;
      std::vector<std::pair<long long, int>> best_sig;
      for (std::size_t c = 0; c < rank.size(); ++c) {
        if (rank[c] != tied) continue;
        auto r = split(rank, c);
        settle(r);
        if (budget == 0) { best = std::move(r); break; }
        --budget;
        const std::vector<int> p = parity_invariants(r);
        std::vector<std::pair<long long, int>> sig(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) sig[i] = {r[i], p[i]};
        std::sort(sig.begin(), sig.end());
        if (best.empty() || std::lexicographical_compare(best_sig.begin(), best_sig.end(), sig.begin(), sig.end(),
                                                         [](const auto& x, const auto& y) { return x.second < y.second; })) {
          best = std::move(r);
          best_sig = std::move(sig);
        }
      }
      rank = std::move(best);
    }
    return rank;
  }

  void collect_stereo(std::uint32_t species_index) {
    if (!opt_.stereo) return;
    for (const auto& t : st_.tetrahedral) {
      if (st_.atoms.species[t.atom] != species_index) continue;
      const CipLabel l = cip_label(st_, t, opt_.canonical ? diag_ : nullptr);
      if (l == CipLabel::None) {
        if (diag_ && opt_.canonical)
          diag_->warn(DiagCode::StereoIgnored, "chirality mark on a non-stereogenic centre dropped from canonical SMILES")
              .on_atom(t.atom);
        continue;
      }
      const int loc = v_.local[t.atom - v_.base];
      centre_label_.push_back({loc, l});
      std::array<int, 4> refs{};
      for (int k = 0; k < 4; ++k) refs[static_cast<std::size_t>(k)] = to_token(t.refs[static_cast<std::size_t>(k)]);
      tetra_.push_back({loc, refs, t.parity});
    }
    for (const auto& d : st_.double_bonds) {
      if (st_.bonds.species[d.bond] != species_index) continue;
      const CipLabel l = cip_label(st_, d, opt_.canonical ? diag_ : nullptr);
      if (l == CipLabel::None) continue;
      db_label_[d.bond] = l;
      dbs_.push_back({d.bond, v_.local[st_.bonds.begin[d.bond] - v_.base], v_.local[st_.bonds.end[d.bond] - v_.base],
                      v_.local[d.ref_begin - v_.base], v_.local[d.ref_end - v_.base], d.cis});
    }
  }

  void dfs(int x) {
    pre_[static_cast<std::size_t>(x)] = counter_++;
    std::vector<Edge> nb = v_.adj[static_cast<std::size_t>(x)];
    std::sort(nb.begin(), nb.end(), [&](const Edge& a, const Edge& b) {
      return rank_[static_cast<std::size_t>(a.to)] < rank_[static_cast<std::size_t>(b.to)];
    });
    for (const Edge& e : nb) {
      if (parent_[static_cast<std::size_t>(x)] >= 0 && e.bond == parent_bond_[static_cast<std::size_t>(x)]) continue;
      if (pre_[static_cast<std::size_t>(e.to)] < 0) {
        parent_[static_cast<std::size_t>(e.to)] = x;
        parent_bond_[static_cast<std::size_t>(e.to)] = e.bond;
        children_[static_cast<std::size_t>(x)].push_back(e.to);
        dfs(e.to);
      } else if (!is_ring_[e.bond - sp_.bond_begin]) {
        is_ring_[e.bond - sp_.bond_begin] = 1;
        ring_at_[static_cast<std::size_t>(x)].push_back(e.bond);
        ring_at_[static_cast<std::size_t>(e.to)].push_back(e.bond);
      }
    }
  }

  // ---- '/' '\' assignment -------------------------------------------------
  int t_side(int e, int x, std::uint32_t b) const {
    return sign_[b - sp_.bond_begin] * (pre_[static_cast<std::size_t>(x)] < pre_[static_cast<std::size_t>(e)] ? 1 : -1);
  }
  bool markable(std::uint32_t b) const { return st_.bonds.order[b] == 1 && !(st_.bonds.flags[b] & kBondAromatic); }

  bool consistent(const ActiveDb& d) const {
    std::array<std::vector<std::pair<int, int>>, 2> marks;
    const std::array<int, 2> ends{d.u, d.v};
    for (int side = 0; side < 2; ++side) {
      for (const Edge& e : v_.adj[static_cast<std::size_t>(ends[side])]) {
        if (e.to == ends[1 - side] || sign_[e.bond - sp_.bond_begin] == 0) continue;
        marks[side].push_back({e.to, t_side(ends[side], e.to, e.bond)});
      }
      if (marks[side].size() == 2 && marks[side][0].second == marks[side][1].second) return false;
    }
    if (marks[0].empty() || marks[1].empty()) return true;
    const auto [x, tx] = marks[0][0];
    const auto [y, ty] = marks[1][0];
    const bool expected = d.cis != (x != d.ra) != (y != d.rb);
    return (tx == ty) == expected;
  }

  void plan_double_bonds() {
    sign_.assign(sp_.bond_count, 0);
    // Make the choice depend only on output order: orient each double bond from
    // its earlier-written end and process bonds in output order.
    for (ActiveDb& d : dbs_) {
      if (pre_[static_cast<std::size_t>(d.u)] > pre_[static_cast<std::size_t>(d.v)]) {
        std::swap(d.u, d.v);
        std::swap(d.ra, d.rb);
      }
    }
    std::sort(dbs_.begin(), dbs_.end(), [&](const ActiveDb& a, const ActiveDb& b) {
      return pre_[static_cast<std::size_t>(a.u)] < pre_[static_cast<std::size_t>(b.u)];
    });
    std::vector<ActiveDb> done;
    for (const ActiveDb& d : dbs_) {
      auto candidates = [&](int e, int far) {
        std::vector<int> c;
        for (const Edge& ed : v_.adj[static_cast<std::size_t>(e)])
          if (ed.to != far && markable(ed.bond)) c.push_back(ed.to);
        std::stable_sort(c.begin(), c.end(), [&](int a, int b) {
          const bool sa = sign_[static_cast<std::uint32_t>(bond_between(e, a)) - sp_.bond_begin] != 0;
          const bool sb = sign_[static_cast<std::uint32_t>(bond_between(e, b)) - sp_.bond_begin] != 0;
          if (sa != sb) return sa;
          return pre_[static_cast<std::size_t>(a)] < pre_[static_cast<std::size_t>(b)];
        });
        return c;
      };
      const auto cu = candidates(d.u, d.v), cv = candidates(d.v, d.u);
      bool placed = false;
      for (int a : cu) {
        for (int b : cv) {
          for (int su : {1, -1}) {
            const std::uint32_t ba = static_cast<std::uint32_t>(bond_between(d.u, a));
            const std::uint32_t bb = static_cast<std::uint32_t>(bond_between(d.v, b));
            const bool cis_ab = d.cis != (a != d.ra) != (b != d.rb);
            const int tv = cis_ab ? su : -su;
            const int need_a = su * (pre_[static_cast<std::size_t>(a)] < pre_[static_cast<std::size_t>(d.u)] ? 1 : -1);
            const int need_b = tv * (pre_[static_cast<std::size_t>(b)] < pre_[static_cast<std::size_t>(d.v)] ? 1 : -1);
            int& sa = sign_[ba - sp_.bond_begin];
            int& sb = sign_[bb - sp_.bond_begin];
            if ((sa != 0 && sa != need_a) || (sb != 0 && sb != need_b)) continue;
            const int old_a = sa, old_b = sb;
            sa = need_a;
            sb = need_b;
            bool ok = consistent(d);
            for (const ActiveDb& p : done) ok = ok && consistent(p);
            if (ok) { placed = true; break; }
            sa = old_a;
            sb = old_b;
          }
          if (placed) break;
        }
        if (placed) break;
      }
      if (placed) {
        done.push_back(d);
      } else if (diag_ && opt_.canonical) {
        diag_->warn(DiagCode::StereoIgnored, "cannot express this double-bond configuration with '/' '\\' marks")
            .on_atom(v_.atoms[static_cast<std::size_t>(d.u)]);
      }
    }
  }

  // ---- text ------------------------------------------------------------------
  std::string bond_symbol(std::uint32_t b, int from, int to) const {
    const int s = sign_.empty() ? 0 : sign_[b - sp_.bond_begin];
    if (s != 0) return s > 0 ? "/" : "\\";
    if (out_aromatic_bond(b)) return "";
    switch (st_.bonds.order[b]) {
      case 1: return (out_aromatic(from) && out_aromatic(to)) ? "-" : "";
      case 2: return "=";
      case 3: return "#";
      default: return "$";
    }
  }

  static std::string digit_text(int d) { return d < 10 ? std::to_string(d) : "%" + std::to_string(d); }

  int predicted_h(int x) const {
    const std::uint32_t g = v_.atoms[static_cast<std::size_t>(x)];
    int sum = 0;
    for (const Edge& e : v_.adj[static_cast<std::size_t>(x)]) sum += out_aromatic_bond(e.bond) ? 1 : st_.bonds.order[e.bond];
    if (out_aromatic(x)) {
      const int free = lowest_default_valence(st_.atoms.element[g], 0) - sum;
      return free >= 1 ? free - 1 : 0;
    }
    for (std::uint8_t val : default_valences(st_.atoms.element[g], 0))
      if (val >= sum) return val - sum;
    return 0;
  }

  std::string atom_text(int x, const std::string& chir) const {
    const std::uint32_t g = v_.atoms[static_cast<std::size_t>(x)];
    const std::uint8_t z = st_.atoms.element[g];
    const int q = st_.atoms.formal_charge[g];
    const bool arom = out_aromatic(x);
    std::string sym(element_symbol(z));
    if (arom) sym[0] = static_cast<char>(sym[0] - 'A' + 'a');
    const bool organic = z == 5 || z == 6 || z == 7 || z == 8 || z == 15 || z == 16 || z == 9 || z == 17 || z == 35 || z == 53;
    const bool organic_arom = z == 5 || z == 6 || z == 7 || z == 8 || z == 15 || z == 16;
    if (organic && (!arom || organic_arom) && q == 0 && st_.atoms.isotope[g] == 0 && st_.atoms.atom_map[g] == 0 &&
        chir.empty() && predicted_h(x) == v_.h[static_cast<std::size_t>(x)])
      return sym;
    std::string s = "[";
    if (st_.atoms.isotope[g]) s += std::to_string(st_.atoms.isotope[g]);
    s += sym + chir;
    const int h = v_.h[static_cast<std::size_t>(x)];
    if (h > 0) s += "H" + (h > 1 ? std::to_string(h) : std::string());
    if (q != 0) s += (q > 0 ? "+" : "-") + (std::abs(q) > 1 ? std::to_string(std::abs(q)) : std::string());
    if (st_.atoms.atom_map[g]) s += ":" + std::to_string(st_.atoms.atom_map[g]);
    return s + "]";
  }

  std::string chirality_text(int x, const std::vector<int>& out_tokens) const {
    for (const auto& t : tetra_) {
      if (t.loc != x) continue;
      if (out_tokens.size() != 4) return "";
      std::array<std::size_t, 4> perm{};
      std::array<char, 4> used{};
      for (std::size_t i = 0; i < 4; ++i) {
        bool found = false;
        for (std::size_t j = 0; j < 4; ++j) {
          if (!used[j] && t.refs[j] == out_tokens[i]) { perm[i] = j; used[j] = 1; found = true; break; }
        }
        if (!found) return "";
      }
      int inv = 0;
      for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = i + 1; j < 4; ++j) inv += perm[i] > perm[j];
      const bool ccw = (t.parity == Chirality::CounterClockwise) != ((inv & 1) != 0);
      return ccw ? "@" : "@@";
    }
    return "";
  }

  bool has_lone_pair_ref(int x) const {
    for (const auto& t : tetra_)
      if (t.loc == x) return std::find(t.refs.begin(), t.refs.end(), kTokLP) != t.refs.end();
    return false;
  }

  void write(int x, std::string& out) {
    const std::size_t xi = static_cast<std::size_t>(x);
    if (parent_[xi] >= 0) out += bond_symbol(parent_bond_[xi], parent_[xi], x);
    std::vector<std::uint32_t> closes, opens;
    for (std::uint32_t b : ring_at_[xi]) {
      const int o = other(b, x);
      (pre_[static_cast<std::size_t>(o)] < pre_[xi] ? closes : opens).push_back(b);
    }
    std::sort(closes.begin(), closes.end(), [&](std::uint32_t a, std::uint32_t b) {
      return digit_[a - sp_.bond_begin] < digit_[b - sp_.bond_begin];
    });
    std::sort(opens.begin(), opens.end(), [&](std::uint32_t a, std::uint32_t b) {
      return rank_[static_cast<std::size_t>(other(a, x))] < rank_[static_cast<std::size_t>(other(b, x))];
    });

    std::vector<int> tokens;
    if (parent_[xi] >= 0) tokens.push_back(parent_[xi]);
    if (v_.h[xi] > 0) tokens.push_back(kTokH);
    else if (has_lone_pair_ref(x)) tokens.push_back(kTokLP);
    for (std::uint32_t b : closes) tokens.push_back(other(b, x));
    for (std::uint32_t b : opens) tokens.push_back(other(b, x));
    for (int c : children_[xi]) tokens.push_back(c);

    out += atom_text(x, chirality_text(x, tokens));
    for (std::uint32_t b : closes) out += digit_text(digit_[b - sp_.bond_begin]);
    for (std::uint32_t b : opens) {
      int d = 1;
      while (std::find(in_use_.begin(), in_use_.end(), d) != in_use_.end()) ++d;
      in_use_.push_back(d);
      digit_[b - sp_.bond_begin] = d;
      out += bond_symbol(b, x, other(b, x)) + digit_text(d);
    }
    for (std::uint32_t b : closes) in_use_.erase(std::find(in_use_.begin(), in_use_.end(), digit_[b - sp_.bond_begin]));

    const auto& kids = children_[xi];
    for (std::size_t i = 0; i < kids.size(); ++i) {
      if (i + 1 < kids.size()) {
        out += '(';
        write(kids[i], out);
        out += ')';
      } else {
        write(kids[i], out);
      }
    }
  }

  struct Tetra {
    int loc;
    std::array<int, 4> refs;
    Chirality parity;
  };

  const ChemicalSystemState& st_;
  const Species& sp_;
  SmilesWriteOptions opt_;
  Diagnostics* diag_;
  View v_;
  std::vector<long long> rank_;
  std::vector<int> pre_, parent_;
  std::vector<std::uint32_t> parent_bond_;
  std::vector<std::vector<int>> children_;
  std::vector<std::vector<std::uint32_t>> ring_at_;
  std::vector<char> is_ring_;
  std::vector<int> sign_, digit_, in_use_;
  int counter_ = 0;
  std::vector<std::pair<int, CipLabel>> centre_label_;
  std::vector<Tetra> tetra_;
  std::vector<ActiveDb> dbs_;
  std::map<std::uint32_t, CipLabel> db_label_;
};

}  // namespace

std::string write_smiles(const ChemicalSystemState& st, std::uint32_t sp, const SmilesWriteOptions& opt, Diagnostics* diag) {
  Writer w(st, st.species[sp], opt, diag);
  return w.run(sp);
}

std::string canonical_smiles(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag) {
  return write_smiles(st, sp, SmilesWriteOptions{}, &diag);
}

std::vector<std::uint32_t> canonical_ranks(const ChemicalSystemState& st, std::uint32_t sp) {
  const Species& s = st.species[sp];
  View v = build_view(st, s);
  std::vector<int> inv(v.atoms.size(), 0);
  for (const auto& t : st.tetrahedral) {
    if (st.atoms.species[t.atom] != sp) continue;
    const int loc = v.local[t.atom - v.base];
    if (loc >= 0) inv[static_cast<std::size_t>(loc)] = static_cast<int>(cip_label(st, t));
  }
  auto rank = compute_ranks(st, v, inv);
  std::vector<std::uint32_t> out(s.atom_count, UINT32_MAX);
  for (std::size_t i = 0; i < v.atoms.size(); ++i) out[v.atoms[i] - s.atom_begin] = static_cast<std::uint32_t>(rank[i]);
  return out;
}

}  // namespace omnichem
