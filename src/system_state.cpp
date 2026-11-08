#include "omnichem/system_state.hpp"

#include <limits>

namespace omnichem {

void AtomTable::reserve(std::size_t n) {
  x.reserve(n); y.reserve(n); z.reserve(n);
  element.reserve(n); formal_charge.reserve(n); implicit_h.reserve(n); flags.reserve(n);
  isotope.reserve(n); chirality.reserve(n); atom_map.reserve(n); species.reserve(n);
}

void AtomTable::push(const Atom& a) {
  x.push_back(a.x); y.push_back(a.y); z.push_back(a.z);
  element.push_back(a.element);
  formal_charge.push_back(a.formal_charge);
  implicit_h.push_back(a.implicit_h);
  flags.push_back(a.flags);
  isotope.push_back(a.isotope);
  chirality.push_back(static_cast<std::int8_t>(a.chirality));
  atom_map.push_back(a.atom_map);
  species.push_back(a.species);
}

Atom AtomTable::get(std::size_t i) const {
  Atom a;
  a.x = x[i]; a.y = y[i]; a.z = z[i];
  a.element = element[i];
  a.formal_charge = formal_charge[i];
  a.implicit_h = implicit_h[i];
  a.flags = flags[i];
  a.isotope = isotope[i];
  a.chirality = static_cast<Chirality>(chirality[i]);
  a.atom_map = atom_map[i];
  a.species = species[i];
  return a;
}

void AtomTable::truncate(std::size_t n) {
  x.resize(n); y.resize(n); z.resize(n);
  element.resize(n); formal_charge.resize(n); implicit_h.resize(n); flags.resize(n);
  isotope.resize(n); chirality.resize(n); atom_map.resize(n); species.resize(n);
}

void BondTable::reserve(std::size_t n) {
  begin.reserve(n); end.reserve(n); order.reserve(n); flags.reserve(n); species.reserve(n);
}

void BondTable::push(const Bond& b) {
  begin.push_back(b.begin);
  end.push_back(b.end);
  order.push_back(b.order);
  flags.push_back(b.flags);
  species.push_back(b.species);
}

Bond BondTable::get(std::size_t i) const {
  return Bond{begin[i], end[i], order[i], flags[i], species[i]};
}

void BondTable::truncate(std::size_t n) {
  begin.resize(n); end.resize(n); order.resize(n); flags.resize(n); species.resize(n);
}

std::uint32_t ChemicalSystemState::append_species(Species s) {
  species.push_back(std::move(s));
  moles.push_back(0.0);
  g_RT.push_back(std::numeric_limits<double>::quiet_NaN());
  ++revision;
  return static_cast<std::uint32_t>(species.size() - 1);
}

void ChemicalSystemState::rebuild_adjacency() {
  const std::size_t n = atoms.size();
  const std::size_t m = bonds.size();
  adjacency.offsets.assign(n + 1, 0);
  for (std::size_t b = 0; b < m; ++b) {
    ++adjacency.offsets[bonds.begin[b] + 1];
    ++adjacency.offsets[bonds.end[b] + 1];
  }
  for (std::size_t i = 0; i < n; ++i) adjacency.offsets[i + 1] += adjacency.offsets[i];
  adjacency.neighbor.assign(2 * m, 0);
  adjacency.bond.assign(2 * m, 0);
  AlignedVector<std::uint32_t> cursor(adjacency.offsets.begin(), adjacency.offsets.end() - 1);
  for (std::size_t b = 0; b < m; ++b) {
    const std::uint32_t u = bonds.begin[b], v = bonds.end[b];
    adjacency.neighbor[cursor[u]] = v;
    adjacency.bond[cursor[u]++] = static_cast<std::uint32_t>(b);
    adjacency.neighbor[cursor[v]] = u;
    adjacency.bond[cursor[v]++] = static_cast<std::uint32_t>(b);
  }
  ++revision;
}

}  // namespace omnichem
