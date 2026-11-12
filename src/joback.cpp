#include "chemlarp/joback.hpp"

#include <array>
#include <string>

#include "chemlarp/elements.hpp"

namespace chemlarp {
namespace {

struct GroupData {
  const char* name;
  double dHf;  // kJ/mol
  double dGf;  // kJ/mol
};

constexpr double kHfBase = 68.29;
constexpr double kGfBase = 53.88;

constexpr std::array<GroupData, static_cast<std::size_t>(JobackGroup::Count)> kGroups = {{
    {"-CH3", -76.45, -43.96},
    {"-CH2-", -20.64, 8.42},
    {">CH-", 29.89, 58.36},
    {">C<", 82.23, 116.02},
    {"=CH2", -9.63, 3.77},
    {"=CH-", 37.97, 48.53},
    {"=C<", 83.99, 92.36},
    {"=C=", 142.14, 136.70},
    {"#CH", 79.30, 77.71},
    {"#C-", 115.51, 109.82},
    {"ring -CH2-", -26.80, -3.68},
    {"ring >CH-", 8.67, 40.99},
    {"ring >C<", 79.72, 87.88},
    {"ring =CH-", 2.09, 11.30},
    {"ring =C<", 46.43, 54.05},
    {"-F", -251.92, -247.19},
    {"-Cl", -71.55, -64.31},
    {"-Br", -29.48, -38.06},
    {"-I", 21.06, 5.74},
    {"-OH (alcohol)", -208.04, -189.20},
    {"-OH (phenol)", -221.65, -197.37},
    {"-O- (nonring)", -132.22, -105.00},
    {"-O- (ring)", -138.16, -98.22},
    {">C=O (nonring)", -133.22, -120.50},
    {">C=O (ring)", -164.50, -126.27},
    {"O=CH- (aldehyde)", -162.03, -143.48},
    {"-COOH (acid)", -426.72, -387.87},
    {"-COO- (ester)", -337.92, -301.95},
    {"=O (other)", -247.61, -250.83},
    {"-NH2", -22.02, 14.07},
    {">NH (nonring)", 53.47, 89.39},
    {">NH (ring)", 31.65, 75.61},
    {">N- (nonring)", 123.34, 163.16},
    {"-N= (ring)", 55.52, 79.93},
    {"-CN", 88.43, 89.22},
    {"-NO2", -66.57, -16.83},
}};

struct Nb {
  std::uint32_t atom;
  std::uint8_t order;
};

}  // namespace

const char* to_string(JobackGroup g) noexcept { return kGroups[static_cast<std::size_t>(g)].name; }

JobackResult joback_estimate(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag) {
  JobackResult res;
  const Species& s = st.species[sp];
  const auto& A = st.atoms;
  const std::uint32_t base = s.atom_begin;
  std::vector<char> used(s.atom_count, 0);
  std::array<int, static_cast<std::size_t>(JobackGroup::Count)> counts{};

  auto el = [&](std::uint32_t g) { return A.element[g]; };
  auto is_used = [&](std::uint32_t g) { return used[g - base] != 0; };
  auto use = [&](std::uint32_t g) { used[g - base] = 1; };
  auto neighbors = [&](std::uint32_t g) {
    std::vector<Nb> out;
    auto nbr = st.adjacency.neighbors(g);
    auto bnd = st.adjacency.bonds_of(g);
    for (std::size_t k = 0; k < nbr.size(); ++k) out.push_back({nbr[k], st.bonds.order[bnd[k]]});
    return out;
  };
  auto total_h = [&](std::uint32_t g) {
    int h = A.implicit_h[g];
    for (const Nb& n : neighbors(g)) h += (el(n.atom) == 1);
    return h;
  };
  auto heavy_degree = [&](std::uint32_t g) {
    int d = 0;
    for (const Nb& n : neighbors(g)) d += (el(n.atom) != 1);
    return d;
  };
  auto in_ring = [&](std::uint32_t g) { return (A.flags[g] & kAtomInRing) != 0; };
  auto aromatic = [&](std::uint32_t g) { return (A.flags[g] & kAtomAromatic) != 0; };
  auto terminal_O = [&](std::uint32_t o) { return el(o) == 8 && heavy_degree(o) == 1 && total_h(o) == 0; };
  auto add = [&](JobackGroup grp) { ++counts[static_cast<std::size_t>(grp)]; };

  const std::uint32_t end = s.atom_begin + s.atom_count;

  // Pass 1: multi-atom groups (nitrile, carbonyl family) centred on carbon.
  for (std::uint32_t g = s.atom_begin; g < end; ++g) {
    if (el(g) != 6 || is_used(g) || aromatic(g)) continue;
    const auto nbs = neighbors(g);
    bool done = false;
    for (const Nb& n : nbs) {
      if (el(n.atom) == 7 && n.order == 3 && heavy_degree(n.atom) == 1 && total_h(n.atom) == 0 && !is_used(n.atom)) {
        add(JobackGroup::Nitrile);
        use(g);
        use(n.atom);
        done = true;
        break;
      }
    }
    if (done) continue;
    std::int64_t dO = -1;
    for (const Nb& n : nbs)
      if (n.order == 2 && !is_used(n.atom) && terminal_O(n.atom)) { dO = n.atom; break; }
    if (dO < 0) continue;
    std::int64_t oh = -1, oe = -1;
    for (const Nb& n : nbs) {
      if (el(n.atom) != 8 || n.order != 1 || is_used(n.atom) || A.formal_charge[n.atom] != 0) continue;
      if (total_h(n.atom) == 1 && heavy_degree(n.atom) == 1) oh = n.atom;
      else if (total_h(n.atom) == 0 && heavy_degree(n.atom) == 2) oe = n.atom;
    }
    const int h = total_h(g);
    if (oh >= 0) {
      add(JobackGroup::Acid);
      use(static_cast<std::uint32_t>(oh));
    } else if (oe >= 0) {
      add(JobackGroup::Ester);
      use(static_cast<std::uint32_t>(oe));
    } else if (h == 1) {
      add(JobackGroup::Aldehyde);
    } else if (h == 0) {
      add(in_ring(g) ? JobackGroup::RingKetone : JobackGroup::Ketone);
    } else {
      continue;  // formaldehyde: no Joback group; reported in pass 2
    }
    use(g);
    use(static_cast<std::uint32_t>(dO));
  }

  // Pass 1b: nitro groups (either N(=O)=O or charge-separated [N+](=O)[O-]).
  for (std::uint32_t g = s.atom_begin; g < end; ++g) {
    if (el(g) != 7 || is_used(g) || heavy_degree(g) != 3) continue;
    std::vector<std::uint32_t> os;
    for (const Nb& n : neighbors(g))
      if (terminal_O(n.atom) && !is_used(n.atom)) os.push_back(n.atom);
    if (os.size() != 2) continue;
    add(JobackGroup::Nitro);
    use(g);
    for (auto o : os) use(o);
  }

  // Pass 2: single-atom groups.
  res.ok = true;
  for (std::uint32_t g = s.atom_begin; g < end; ++g) {
    if (is_used(g) || el(g) == 1) continue;
    const auto nbs = neighbors(g);
    const int h = total_h(g);
    const int heavy = heavy_degree(g);
    int dbl = 0, tri = 0;
    for (const Nb& n : nbs) {
      if (aromatic(g) && aromatic(n.atom)) continue;  // Kekulé orders inside rings are not group-relevant
      dbl += (n.order == 2);
      tri += (n.order == 3);
    }
    const bool ring = in_ring(g);
    bool assigned = true;
    if (A.formal_charge[g] != 0) {
      assigned = false;
    } else {
      switch (el(g)) {
        case 6:
          if (aromatic(g)) {
            if (h == 1) add(JobackGroup::RingEqCH);
            else if (h == 0) add(JobackGroup::RingEqC);
            else assigned = false;
          } else if (tri == 1) {
            if (h == 1) add(JobackGroup::TrCH);
            else if (h == 0) add(JobackGroup::TrC);
            else assigned = false;
          } else if (dbl == 2 && h == 0) {
            add(JobackGroup::EqCEq);
          } else if (dbl == 1) {
            if (h == 2) add(JobackGroup::EqCH2);
            else if (h == 1) add(ring ? JobackGroup::RingEqCH : JobackGroup::EqCH);
            else if (h == 0) add(ring ? JobackGroup::RingEqC : JobackGroup::EqC);
            else assigned = false;
          } else if (dbl == 0 && tri == 0) {
            if (h == 3) add(JobackGroup::CH3);
            else if (h == 2) add(ring ? JobackGroup::RingCH2 : JobackGroup::CH2);
            else if (h == 1) add(ring ? JobackGroup::RingCH : JobackGroup::CH);
            else if (h == 0) add(ring ? JobackGroup::RingC : JobackGroup::C);
            else assigned = false;
          } else {
            assigned = false;
          }
          break;
        case 8:
          if (h == 1 && heavy == 1) {
            const bool phenol = aromatic(nbs.front().atom) && el(nbs.front().atom) == 6;
            add(phenol ? JobackGroup::OHPhenol : JobackGroup::OHAlcohol);
          } else if (h == 0 && heavy == 2 && dbl == 0) {
            add(ring ? JobackGroup::ORingEther : JobackGroup::OEther);
          } else if (h == 0 && heavy == 1 && dbl == 1) {
            add(JobackGroup::OxoOther);
          } else {
            assigned = false;
          }
          break;
        case 7:
          if (aromatic(g) && h == 0 && heavy == 2) add(JobackGroup::RingNEq);
          else if (aromatic(g) && h == 1) add(JobackGroup::RingNH);
          else if (!aromatic(g) && h == 2 && heavy == 1 && dbl == 0) add(JobackGroup::NH2);
          else if (!aromatic(g) && h == 1 && heavy == 2 && dbl == 0) add(ring ? JobackGroup::RingNH : JobackGroup::NH);
          else if (!aromatic(g) && h == 0 && heavy == 3 && dbl == 0 && !ring) add(JobackGroup::N);
          else assigned = false;
          break;
        case 9: case 17: case 35: case 53:
          if (heavy == 1 && h == 0) {
            add(el(g) == 9 ? JobackGroup::F : el(g) == 17 ? JobackGroup::Cl : el(g) == 35 ? JobackGroup::Br : JobackGroup::I);
          } else {
            assigned = false;
          }
          break;
        default:
          assigned = false;
      }
    }
    if (!assigned) {
      res.ok = false;
      diag.error(DiagCode::NoJobackGroup,
                 "no Joback group for atom " + std::to_string(g - base) + " (" +
                     std::string(element_symbol(el(g))) + " with " + std::to_string(h) + " H, " +
                     std::to_string(heavy) + " heavy neighbours) in " + s.smiles)
          .on_atom(g)
          .on_species(sp);
    }
  }
  int ngroups = 0;
  res.dHf_kJ = kHfBase;
  res.dGf_kJ = kGfBase;
  for (std::size_t k = 0; k < counts.size(); ++k) {
    if (counts[k] == 0) continue;
    ngroups += counts[k];
    res.dHf_kJ += counts[k] * kGroups[k].dHf;
    res.dGf_kJ += counts[k] * kGroups[k].dGf;
    res.groups.push_back({static_cast<JobackGroup>(k), counts[k]});
  }
  if (ngroups == 0 && res.ok) {
    res.ok = false;
    diag.error(DiagCode::NoJobackGroup, "species " + s.smiles + " contains no Joback groups").on_species(sp);
  }
  return res;
}

}  // namespace chemlarp
