#include "omnichem/fingerprint.hpp"

#include <algorithm>
#include <bit>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include "omnichem/cip.hpp"

namespace omnichem {
namespace {

std::uint32_t fmix32(std::uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6bu;
  h ^= h >> 13;
  h *= 0xc2b2ae35u;
  h ^= h >> 16;
  return h;
}

// FNV-1a over 32-bit words (each word fed as four bytes, little-endian), then fmix32.
std::uint32_t hash_words(const std::vector<std::uint32_t>& w) {
  std::uint32_t h = 2166136261u;
  for (std::uint32_t x : w)
    for (int k = 0; k < 4; ++k) {
      h ^= (x >> (8 * k)) & 0xffu;
      h *= 16777619u;
    }
  return fmix32(h);
}

bool plain_h(const ChemicalSystemState& st, std::uint32_t g) {
  if (st.atoms.element[g] != 1 || st.atoms.formal_charge[g] || st.atoms.isotope[g] || st.atoms.implicit_h[g]) return false;
  const auto nb = st.adjacency.neighbors(g);
  return nb.size() == 1 && st.atoms.element[nb[0]] != 1;
}

}  // namespace

std::vector<std::pair<std::uint32_t, std::uint32_t>> morgan_environments(const ChemicalSystemState& st, std::uint32_t sp,
                                                                         const MorganOptions& opt) {
  const Species& S = st.species[sp];
  const std::uint32_t base = S.atom_begin, n = S.atom_count;
  std::vector<char> heavy(n, 0);
  for (std::uint32_t i = 0; i < n; ++i) heavy[i] = !plain_h(st, base + i);

  std::vector<int> cip(n, 0);
  if (opt.use_chirality)
    for (const auto& c : cip_centres(st, sp)) cip[c.index - base] = static_cast<int>(c.label);

  // Local graph over heavy atoms; bond type 1, 2, 3, or 4 for aromatic.
  struct E { std::uint32_t to; std::uint32_t bond; std::uint32_t type; };
  std::vector<std::vector<E>> adj(n);
  std::vector<std::uint32_t> id(n, 0);
  for (std::uint32_t i = 0; i < n; ++i) {
    if (!heavy[i]) continue;
    const std::uint32_t g = base + i;
    int hyd = st.atoms.implicit_h[g], heavy_nb = 0, valence = st.atoms.implicit_h[g];
    const auto nb = st.adjacency.neighbors(g);
    const auto bs = st.adjacency.bonds_of(g);
    for (std::size_t k = 0; k < nb.size(); ++k) {
      valence += st.bonds.order[bs[k]];
      const std::uint32_t l = nb[k] - base;
      if (!heavy[l]) { ++hyd; continue; }
      ++heavy_nb;
      const std::uint32_t type = (st.bonds.flags[bs[k]] & kBondAromatic) ? 4u : st.bonds.order[bs[k]];
      adj[i].push_back({l, bs[k] - S.bond_begin, type});
    }
    std::vector<std::uint32_t> inv{static_cast<std::uint32_t>(heavy_nb),
                                   static_cast<std::uint32_t>(valence - hyd),
                                   st.atoms.element[g],
                                   st.atoms.isotope[g],
                                   static_cast<std::uint32_t>(st.atoms.formal_charge[g] + 128),
                                   static_cast<std::uint32_t>(hyd),
                                   (st.atoms.flags[g] & kAtomInRing) ? 1u : 0u};
    if (opt.use_chirality && cip[i]) inv.push_back(0x1000u + static_cast<std::uint32_t>(cip[i]));
    id[i] = hash_words(inv);
  }

  std::map<std::uint32_t, std::uint32_t> counts;
  for (std::uint32_t i = 0; i < n; ++i)
    if (heavy[i]) ++counts[id[i]];

  // Bond set covered by each atom's environment, for structural duplicate removal.
  const std::size_t words = (S.bond_count + 63) / 64 + 1;
  std::vector<std::vector<std::uint64_t>> cover(n, std::vector<std::uint64_t>(words, 0));
  std::set<std::vector<std::uint64_t>> seen;
  std::vector<char> alive(n, 1);  // atoms whose environment can still grow
  for (int iter = 1; iter <= opt.radius; ++iter) {
    std::vector<std::uint32_t> next = id;
    std::vector<std::vector<std::uint64_t>> next_cover = cover;
    struct Candidate { std::vector<std::uint64_t> bonds; std::uint32_t id; std::uint32_t atom; };
    std::vector<Candidate> found;
    for (std::uint32_t i = 0; i < n; ++i) {
      if (!heavy[i] || !alive[i]) continue;
      std::vector<std::pair<std::uint32_t, std::uint32_t>> nb;
      std::vector<std::uint64_t> cov = cover[i];
      for (const E& e : adj[i]) {
        nb.push_back({e.type, id[e.to]});
        cov[e.bond / 64] |= std::uint64_t{1} << (e.bond % 64);
        for (std::size_t w = 0; w < words; ++w) cov[w] |= cover[e.to][w];
      }
      std::sort(nb.begin(), nb.end());
      std::vector<std::uint32_t> words_in{static_cast<std::uint32_t>(iter), id[i]};
      for (auto [t, x] : nb) { words_in.push_back(t); words_in.push_back(x); }
      next[i] = hash_words(words_in);
      next_cover[i] = cov;
      found.push_back({cov, next[i], i});
    }
    // Same-iteration duplicates: keep the one with the smaller identifier.
    std::sort(found.begin(), found.end(), [](const Candidate& a, const Candidate& b) {
      return a.bonds != b.bonds ? a.bonds < b.bonds : a.id < b.id;
    });
    for (std::size_t k = 0; k < found.size(); ++k) {
      const bool dup_same_iter = k > 0 && found[k].bonds == found[k - 1].bonds;
      const bool empty = std::all_of(found[k].bonds.begin(), found[k].bonds.end(), [](std::uint64_t w) { return w == 0; });
      // An environment that did not grow this round cannot grow later either.
      if (empty || found[k].bonds == cover[found[k].atom]) alive[found[k].atom] = 0;
      if (dup_same_iter || empty || seen.count(found[k].bonds)) continue;
      seen.insert(found[k].bonds);
      ++counts[found[k].id];
    }
    id = std::move(next);
    cover = std::move(next_cover);
  }
  return {counts.begin(), counts.end()};
}

int Fingerprint::count() const {
  int c = 0;
  for (std::uint64_t w : words) c += std::popcount(w);
  return c;
}

Fingerprint morgan_fingerprint(const ChemicalSystemState& st, std::uint32_t sp, const MorganOptions& opt) {
  if (opt.bits <= 0) throw std::invalid_argument("fingerprint length must be positive");
  Fingerprint f;
  f.bits = opt.bits;
  f.words.assign(static_cast<std::size_t>((opt.bits + 63) / 64), 0);
  for (const auto& [ident, count] : morgan_environments(st, sp, opt)) {
    (void)count;
    const std::uint32_t b = ident % static_cast<std::uint32_t>(opt.bits);
    f.words[b / 64] |= std::uint64_t{1} << (b % 64);
  }
  return f;
}

namespace {
std::pair<int, int> overlap(const Fingerprint& a, const Fingerprint& b) {
  if (a.bits != b.bits) throw std::invalid_argument("fingerprints have different lengths");
  int both = 0, total = 0;
  for (std::size_t i = 0; i < a.words.size(); ++i) {
    both += std::popcount(a.words[i] & b.words[i]);
    total += std::popcount(a.words[i]) + std::popcount(b.words[i]);
  }
  return {both, total};
}
}  // namespace

double tanimoto(const Fingerprint& a, const Fingerprint& b) {
  const auto [both, total] = overlap(a, b);
  const int uni = total - both;
  return uni == 0 ? 0.0 : static_cast<double>(both) / uni;
}

double dice(const Fingerprint& a, const Fingerprint& b) {
  const auto [both, total] = overlap(a, b);
  return total == 0 ? 0.0 : 2.0 * both / total;
}

}  // namespace omnichem
