#pragma once
// Periodic-table data used by the parser, valence validator and formula engine.

#include <array>
#include <cstdint>
#include <span>
#include <string_view>

namespace omnichem {

inline constexpr std::size_t kElementSlots = 119;  // indexed by atomic number; slot 0 unused

struct ElementInfo {
  std::uint8_t z;
  std::string_view symbol;
  double mass;                                    // standard atomic weight, g/mol
  std::array<std::uint8_t, 3> default_valences;   // OpenSMILES implicit-H valences, 0-terminated
  std::int8_t max_valence;                        // hard limit for the neutral atom; -1 = unchecked (metals)
  bool aromatic_allowed;

  std::span<const std::uint8_t> valences() const noexcept {
    std::size_t n = 0;
    while (n < default_valences.size() && default_valences[n] != 0) ++n;
    return {default_valences.data(), n};
  }
};

const ElementInfo* find_element_by_symbol(std::string_view symbol) noexcept;
const ElementInfo* find_element_by_z(std::uint8_t z) noexcept;

// Valence rules for charged main-group atoms use the isoelectronic neighbour:
// N+ behaves like C, O- like F, C- like N, B- like C, and so on.
// Returns -1 when the element is not valence-checked.
int max_valence(std::uint8_t z, int formal_charge) noexcept;
std::span<const std::uint8_t> default_valences(std::uint8_t z, int formal_charge) noexcept;
int lowest_default_valence(std::uint8_t z, int formal_charge) noexcept;

std::string_view element_symbol(std::uint8_t z) noexcept;  // "?" for unknown

}  // namespace omnichem
