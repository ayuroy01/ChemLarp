#include "omnichem/elements.hpp"

#include <array>

namespace omnichem {
namespace {

// Masses: IUPAC 2021 conventional standard atomic weights.
// default_valences: OpenSMILES organic-subset normal valences (implicit H rule).
// max_valence: hard upper bound used by the valence validator (-1 = not checked).
constexpr ElementInfo kTable[] = {
    {1, "H", 1.008, {1, 0, 0}, 1, false},
    {2, "He", 4.002602, {0, 0, 0}, 0, false},
    {3, "Li", 6.94, {0, 0, 0}, -1, false},
    {4, "Be", 9.0121831, {0, 0, 0}, -1, false},
    {5, "B", 10.81, {3, 0, 0}, 3, true},
    {6, "C", 12.011, {4, 0, 0}, 4, true},
    {7, "N", 14.007, {3, 5, 0}, 5, true},
    {8, "O", 15.999, {2, 0, 0}, 2, true},
    {9, "F", 18.998403163, {1, 0, 0}, 1, false},
    {10, "Ne", 20.1797, {0, 0, 0}, 0, false},
    {11, "Na", 22.98976928, {0, 0, 0}, -1, false},
    {12, "Mg", 24.305, {0, 0, 0}, -1, false},
    {13, "Al", 26.9815385, {0, 0, 0}, -1, false},
    {14, "Si", 28.085, {4, 0, 0}, 4, false},
    {15, "P", 30.973761998, {3, 5, 0}, 5, true},
    {16, "S", 32.06, {2, 4, 6}, 6, true},
    {17, "Cl", 35.45, {1, 0, 0}, 7, false},
    {18, "Ar", 39.948, {0, 0, 0}, 0, false},
    {19, "K", 39.0983, {0, 0, 0}, -1, false},
    {20, "Ca", 40.078, {0, 0, 0}, -1, false},
    {21, "Sc", 44.955908, {0, 0, 0}, -1, false},
    {22, "Ti", 47.867, {0, 0, 0}, -1, false},
    {23, "V", 50.9415, {0, 0, 0}, -1, false},
    {24, "Cr", 51.9961, {0, 0, 0}, -1, false},
    {25, "Mn", 54.938044, {0, 0, 0}, -1, false},
    {26, "Fe", 55.845, {0, 0, 0}, -1, false},
    {27, "Co", 58.933194, {0, 0, 0}, -1, false},
    {28, "Ni", 58.6934, {0, 0, 0}, -1, false},
    {29, "Cu", 63.546, {0, 0, 0}, -1, false},
    {30, "Zn", 65.38, {0, 0, 0}, -1, false},
    {31, "Ga", 69.723, {0, 0, 0}, -1, false},
    {32, "Ge", 72.630, {4, 0, 0}, 4, false},
    {33, "As", 74.921595, {3, 5, 0}, 5, true},
    {34, "Se", 78.971, {2, 4, 6}, 6, true},
    {35, "Br", 79.904, {1, 0, 0}, 7, false},
    {36, "Kr", 83.798, {0, 0, 0}, -1, false},
    {37, "Rb", 85.4678, {0, 0, 0}, -1, false},
    {38, "Sr", 87.62, {0, 0, 0}, -1, false},
    {39, "Y", 88.90584, {0, 0, 0}, -1, false},
    {40, "Zr", 91.224, {0, 0, 0}, -1, false},
    {41, "Nb", 92.90637, {0, 0, 0}, -1, false},
    {42, "Mo", 95.95, {0, 0, 0}, -1, false},
    {43, "Tc", 98.0, {0, 0, 0}, -1, false},
    {44, "Ru", 101.07, {0, 0, 0}, -1, false},
    {45, "Rh", 102.90550, {0, 0, 0}, -1, false},
    {46, "Pd", 106.42, {0, 0, 0}, -1, false},
    {47, "Ag", 107.8682, {0, 0, 0}, -1, false},
    {48, "Cd", 112.414, {0, 0, 0}, -1, false},
    {49, "In", 114.818, {0, 0, 0}, -1, false},
    {50, "Sn", 118.710, {0, 0, 0}, -1, false},
    {51, "Sb", 121.760, {0, 0, 0}, -1, false},
    {52, "Te", 127.60, {2, 4, 6}, 6, true},
    {53, "I", 126.90447, {1, 0, 0}, 7, false},
    {54, "Xe", 131.293, {0, 0, 0}, -1, false},
    {55, "Cs", 132.90545196, {0, 0, 0}, -1, false},
    {56, "Ba", 137.327, {0, 0, 0}, -1, false},
    {78, "Pt", 195.084, {0, 0, 0}, -1, false},
    {79, "Au", 196.966569, {0, 0, 0}, -1, false},
    {80, "Hg", 200.592, {0, 0, 0}, -1, false},
    {82, "Pb", 207.2, {0, 0, 0}, -1, false},
    {83, "Bi", 208.98040, {0, 0, 0}, -1, false},
};

const std::array<const ElementInfo*, kElementSlots>& by_z_index() {
  static const auto index = [] {
    std::array<const ElementInfo*, kElementSlots> a{};
    for (const auto& e : kTable) a[e.z] = &e;
    return a;
  }();
  return index;
}

const ElementInfo* shifted(std::uint8_t z, int formal_charge) noexcept {
  const ElementInfo* e = find_element_by_z(z);
  if (!e || e->max_valence < 0) return nullptr;
  if (formal_charge == 0) return e;
  const int zs = static_cast<int>(z) - formal_charge;  // isoelectronic neighbour
  if (zs < 1 || zs >= static_cast<int>(kElementSlots)) return nullptr;
  const ElementInfo* e2 = find_element_by_z(static_cast<std::uint8_t>(zs));
  if (!e2 || e2->max_valence < 0) return nullptr;
  return e2;
}

}  // namespace

const ElementInfo* find_element_by_symbol(std::string_view symbol) noexcept {
  for (const auto& e : kTable)
    if (e.symbol == symbol) return &e;
  return nullptr;
}

const ElementInfo* find_element_by_z(std::uint8_t z) noexcept {
  if (z >= kElementSlots) return nullptr;
  return by_z_index()[z];
}

int max_valence(std::uint8_t z, int formal_charge) noexcept {
  const ElementInfo* e = shifted(z, formal_charge);
  return e ? e->max_valence : -1;
}

std::span<const std::uint8_t> default_valences(std::uint8_t z, int formal_charge) noexcept {
  const ElementInfo* e = shifted(z, formal_charge);
  return e ? e->valences() : std::span<const std::uint8_t>{};
}

int lowest_default_valence(std::uint8_t z, int formal_charge) noexcept {
  auto v = default_valences(z, formal_charge);
  return v.empty() ? -1 : v.front();
}

std::string_view element_symbol(std::uint8_t z) noexcept {
  const ElementInfo* e = find_element_by_z(z);
  return e ? e->symbol : std::string_view("?");
}

}  // namespace omnichem
