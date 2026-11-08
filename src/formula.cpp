#include "omnichem/formula.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>

namespace omnichem {

Composition compute_composition(const ChemicalSystemState& state, std::uint32_t species) {
  Composition c{};
  const Species& s = state.species[species];
  for (std::uint32_t i = s.atom_begin; i < s.atom_begin + s.atom_count; ++i) {
    const std::uint8_t z = state.atoms.element[i];
    if (z > 0) ++c[z];
    c[1] += state.atoms.implicit_h[i];
  }
  return c;
}

std::string hill_formula(const Composition& c, int net_charge) {
  std::string out;
  auto emit = [&](std::size_t z) {
    if (c[z] == 0) return;
    out += element_symbol(static_cast<std::uint8_t>(z));
    if (c[z] != 1) out += std::to_string(c[z]);
  };
  std::vector<std::size_t> rest;
  const bool has_carbon = c[6] > 0;
  for (std::size_t z = 1; z < kElementSlots; ++z) {
    if (c[z] == 0) continue;
    if (has_carbon && (z == 6 || z == 1)) continue;
    rest.push_back(z);
  }
  std::sort(rest.begin(), rest.end(), [](std::size_t a, std::size_t b) {
    return element_symbol(static_cast<std::uint8_t>(a)) < element_symbol(static_cast<std::uint8_t>(b));
  });
  if (has_carbon) { emit(6); emit(1); }
  for (std::size_t z : rest) emit(z);
  if (net_charge != 0) {
    if (std::abs(net_charge) > 1) out += std::to_string(std::abs(net_charge));
    out += net_charge > 0 ? '+' : '-';
  }
  return out;
}

double molar_mass(const Composition& c) {
  double m = 0.0;
  for (std::size_t z = 1; z < kElementSlots; ++z) {
    if (c[z] == 0) continue;
    const ElementInfo* e = find_element_by_z(static_cast<std::uint8_t>(z));
    if (e) m += c[z] * e->mass;
  }
  return m;
}

}  // namespace omnichem
