#include "omnichem/rings.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <iterator>
#include <set>
#include <utility>

namespace omnichem {
namespace {

using BitRow = std::vector<std::uint64_t>;

BitRow to_bits(const std::vector<std::uint32_t>& bonds, std::uint32_t bond_begin, std::size_t nbonds) {
  BitRow r((nbonds + 63) / 64, 0);
  for (std::uint32_t b : bonds) {
    const std::size_t i = b - bond_begin;
    r[i / 64] |= std::uint64_t{1} << (i % 64);
  }
  return r;
}

// Orders a cycle's atoms by walking its bonds; returns false if the bond set is not one simple cycle.
bool cycle_atoms(const ChemicalSystemState& st, const std::vector<std::uint32_t>& bonds, std::vector<std::uint32_t>& out) {
  out.clear();
  if (bonds.size() < 3) return false;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
  for (std::uint32_t b : bonds) edges.push_back({st.bonds.begin[b], st.bonds.end[b]});
  std::vector<std::uint32_t> atoms;
  for (auto [u, v] : edges) { atoms.push_back(u); atoms.push_back(v); }
  std::sort(atoms.begin(), atoms.end());
  for (std::size_t i = 0; i < atoms.size(); i += 2)
    if (i + 1 >= atoms.size() || atoms[i] != atoms[i + 1] || (i + 2 < atoms.size() && atoms[i + 2] == atoms[i]))
      return false;  // every atom must have degree exactly 2
  std::vector<char> used(edges.size(), 0);
  std::uint32_t cur = edges[0].first;
  const std::uint32_t start = cur;
  for (std::size_t step = 0; step < edges.size(); ++step) {
    out.push_back(cur);
    bool moved = false;
    for (std::size_t e = 0; e < edges.size(); ++e) {
      if (used[e]) continue;
      if (edges[e].first == cur) { cur = edges[e].second; used[e] = 1; moved = true; break; }
      if (edges[e].second == cur) { cur = edges[e].first; used[e] = 1; moved = true; break; }
    }
    if (!moved) return false;
  }
  return cur == start && out.size() == bonds.size();
}

// Shortest cycle through bond `bond` using only ring bonds of the species.
std::vector<std::uint32_t> shortest_cycle(const ChemicalSystemState& st, std::uint32_t bond) {
  const std::uint32_t src = st.bonds.begin[bond], dst = st.bonds.end[bond];
  std::vector<std::int64_t> via_bond(st.atoms.size(), -1);
  std::vector<char> seen(st.atoms.size(), 0);
  std::vector<std::uint32_t> queue{src};
  seen[src] = 1;
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::uint32_t x = queue[head];
    auto nb = st.adjacency.neighbors(x);
    auto bd = st.adjacency.bonds_of(x);
    for (std::size_t k = 0; k < nb.size(); ++k) {
      if (bd[k] == bond || !(st.bonds.flags[bd[k]] & kBondInRing) || seen[nb[k]]) continue;
      seen[nb[k]] = 1;
      via_bond[nb[k]] = bd[k];
      if (nb[k] == dst) {
        std::vector<std::uint32_t> cyc{bond};
        std::uint32_t y = dst;
        while (y != src) {
          const std::uint32_t b = static_cast<std::uint32_t>(via_bond[y]);
          cyc.push_back(b);
          y = st.bonds.begin[b] == y ? st.bonds.end[b] : st.bonds.begin[b];
        }
        std::sort(cyc.begin(), cyc.end());
        return cyc;
      }
      queue.push_back(nb[k]);
    }
  }
  return {};
}

std::vector<std::vector<std::uint32_t>> candidate_cycles(const ChemicalSystemState& st, const Species& s) {
  std::set<std::vector<std::uint32_t>> uniq;
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    if (!(st.bonds.flags[b] & kBondInRing)) continue;
    auto c = shortest_cycle(st, b);
    if (!c.empty()) uniq.insert(std::move(c));
  }
  std::vector<std::vector<std::uint32_t>> out(uniq.begin(), uniq.end());
  std::stable_sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); });
  return out;
}

}  // namespace

std::vector<Ring> find_sssr(const ChemicalSystemState& st, std::uint32_t sp) {
  const Species& s = st.species[sp];
  std::vector<Ring> rings;
  if (s.bond_count == 0) return rings;
  auto cands = candidate_cycles(st, s);
  std::vector<BitRow> basis;  // reduced rows with their pivot bit
  std::vector<std::size_t> pivots;
  for (auto& c : cands) {
    BitRow row = to_bits(c, s.bond_begin, s.bond_count);
    for (std::size_t k = 0; k < basis.size(); ++k) {
      const std::size_t p = pivots[k];
      if (row[p / 64] >> (p % 64) & 1)
        for (std::size_t w = 0; w < row.size(); ++w) row[w] ^= basis[k][w];
    }
    std::size_t pivot = SIZE_MAX;
    for (std::size_t w = 0; w < row.size() && pivot == SIZE_MAX; ++w)
      if (row[w]) pivot = w * 64 + static_cast<std::size_t>(std::countr_zero(row[w]));
    if (pivot == SIZE_MAX) continue;  // dependent on smaller rings
    basis.push_back(row);
    pivots.push_back(pivot);
    Ring r;
    r.bonds = c;
    cycle_atoms(st, c, r.atoms);
    rings.push_back(std::move(r));
  }
  return rings;
}

int aromatic_pi_electrons(const ChemicalSystemState& st, std::uint32_t a) {
  if (!(st.atoms.flags[a] & kAtomInRing)) return -1;
  const std::uint8_t z = st.atoms.element[a];
  const int q = st.atoms.formal_charge[a];
  switch (z) {
    case 5: case 6: case 7: case 8: case 15: case 16: case 33: case 34: case 52: break;
    default: return -1;
  }
  int doubles = 0, degree = st.atoms.implicit_h[a];
  std::int64_t partner_bond = -1;
  auto nb = st.adjacency.neighbors(a);
  auto bd = st.adjacency.bonds_of(a);
  for (std::size_t k = 0; k < nb.size(); ++k) {
    const std::uint8_t o = st.bonds.order[bd[k]];
    ++degree;
    if (o >= 3) return -1;
    if (o == 2) { ++doubles; partner_bond = bd[k]; }
  }
  if (doubles > 1) return -1;
  if (doubles == 1) {
    const std::uint32_t pb = static_cast<std::uint32_t>(partner_bond);
    if (st.bonds.flags[pb] & kBondInRing) return 1;
    const std::uint32_t w = st.bonds.begin[pb] == a ? st.bonds.end[pb] : st.bonds.begin[pb];
    const std::uint8_t zw = st.atoms.element[w];
    return (zw == 8 || zw == 7 || zw == 16) ? 0 : -1;  // exocyclic C=O etc. contributes no pi electrons
  }
  switch (z) {
    case 6:
      if (q == -1 && degree == 3) return 2;
      if (q == 1 && degree == 3) return 0;
      return -1;
    case 7: case 15: case 33:
      if (q == 0 && degree == 3) return 2;
      if (q == -1 && degree == 2) return 2;
      return -1;
    case 8: case 16: case 34: case 52:
      return (q == 0 && degree == 2) ? 2 : -1;
    case 5:
      return (q == 0 && degree == 3) ? 0 : -1;
    default:
      return -1;
  }
}

void perceive_aromaticity(ChemicalSystemState& st, std::uint32_t sp) {
  const Species& s = st.species[sp];
  for (std::uint32_t a = s.atom_begin; a < s.atom_begin + s.atom_count; ++a)
    st.atoms.flags[a] = static_cast<std::uint8_t>(st.atoms.flags[a] & ~kAtomAromatic);
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b)
    st.bonds.flags[b] = static_cast<std::uint8_t>(st.bonds.flags[b] & ~kBondAromatic);
  if (s.bond_count == 0) return;

  std::vector<int> pi(s.atom_count);
  for (std::uint32_t a = 0; a < s.atom_count; ++a) pi[a] = aromatic_pi_electrons(st, s.atom_begin + a);

  const std::vector<Ring> sssr = find_sssr(st, sp);
  std::vector<std::vector<std::uint32_t>> cycles;
  for (const auto& r : sssr) cycles.push_back(r.bonds);
  // Envelopes of fused pairs (shared bonds removed), e.g. the 10-cycle of azulene.
  for (std::size_t i = 0; i < sssr.size(); ++i) {
    for (std::size_t j = i + 1; j < sssr.size(); ++j) {
      std::vector<std::uint32_t> shared, sym;
      std::set_intersection(sssr[i].bonds.begin(), sssr[i].bonds.end(), sssr[j].bonds.begin(), sssr[j].bonds.end(),
                            std::back_inserter(shared));
      if (shared.empty()) continue;
      std::set_symmetric_difference(sssr[i].bonds.begin(), sssr[i].bonds.end(), sssr[j].bonds.begin(),
                                    sssr[j].bonds.end(), std::back_inserter(sym));
      if (sym.size() <= 24) cycles.push_back(std::move(sym));
    }
  }

  std::vector<char> arom_atom(s.atom_count, 0);
  std::vector<std::uint32_t> atoms;
  std::vector<char> cycle_aromatic(cycles.size(), 0);
  for (std::size_t c = 0; c < cycles.size(); ++c) {
    if (!cycle_atoms(st, cycles[c], atoms)) continue;
    int total = 0;
    bool ok = true;
    for (std::uint32_t a : atoms) {
      const int e = pi[a - s.atom_begin];
      if (e < 0) { ok = false; break; }
      total += e;
    }
    if (!ok || total < 2 || total % 4 != 2) continue;
    cycle_aromatic[c] = 1;
    for (std::uint32_t a : atoms) arom_atom[a - s.atom_begin] = 1;
    for (std::uint32_t b : cycles[c]) st.bonds.flags[b] |= kBondAromatic;
  }
  // Bonds inside an SSSR ring (size >= 5) whose atoms all became aromatic (azulene's fusion bond).
  for (const auto& r : sssr) {
    if (r.atoms.size() < 5) continue;
    bool all = !r.atoms.empty();
    for (std::uint32_t a : r.atoms) all = all && arom_atom[a - s.atom_begin];
    if (!all) continue;
    for (std::uint32_t b : r.bonds) st.bonds.flags[b] |= kBondAromatic;
  }
  for (std::uint32_t a = 0; a < s.atom_count; ++a)
    if (arom_atom[a]) st.atoms.flags[s.atom_begin + a] |= kAtomAromatic;
}

}  // namespace omnichem
