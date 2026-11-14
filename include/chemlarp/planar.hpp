#pragma once
// Perception of planar (cyclophane, metallocene) and helical (helicene) chirality from 3D coordinates.
// Descriptors (Rp/Sp, P/M) come from cip_label(state, PlanarStereo); see cip.hpp for the rules.
//
// What is found:
//  - cyclophanes: an aromatic 5- or 6-membered ring with a bridge that leaves the ring plane (an atom outside
//    any aromatic ring, at least 0.25 Å off the ring's mean plane) and comes back to a non-adjacent ring atom;
//  - metallocenes and arene complexes: a 5- or 6-membered ring with one atom bonded to every ring atom;
//  - helicenes: five or more ortho-fused aromatic 5/6-rings turning the same way, whose inner-rim atoms belong
//    to two rings only (so circulenes such as coronene are not taken) and whose inner-rim torsion angles add
//    up to at least 15°.
// A cyclophane or metallocene ring is kept only when it carries a substituent other than the bridge or
// metal (or a ring heteroatom) and CIP gives it a descriptor.

#include <cstdint>
#include <vector>

#include "chemlarp/system_state.hpp"

namespace chemlarp {

// Appends the planar and helical units of `species` to state.planar, from state.atoms.x/y/z.
void perceive_planar_stereo(ChemicalSystemState& state, std::uint32_t species);

// True when neighbour x of ring[ir] (x outside the ring) leads back, without passing through ring[ir], to a
// ring atom that is neither ring[ir] nor one of its ring neighbours: x belongs to a bridge across the ring.
bool bridges_ring(const ChemicalSystemState& state, const std::vector<std::uint32_t>& ring, std::size_t ir, std::uint32_t x);

}  // namespace chemlarp
