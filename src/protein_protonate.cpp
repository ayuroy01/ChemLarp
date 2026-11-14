// Protonation states and hydrogens, with the hydrogen-bond network optimised over the movable groups.
#include <cstring>
#include <functional>
#include <numeric>
#include <set>

#include "protein_internal.hpp"

namespace chemlarp {

using namespace detail;

namespace {

// One atom as the network scoring sees it.
struct SAtom {
  Vec x{};
  char el = 'C';
  bool acceptor = false;  // O, or a ring N without hydrogen
  bool polar_h = false;   // H on N, O or S
  bool water = false;
  Vec parent{};           // of a polar H
  std::uint32_t res = 0;
};

struct State {
  std::vector<SAtom> atoms;                       // the atoms that differ between states
  std::vector<std::pair<std::string, Vec>> h;     // hydrogens to add
  std::vector<std::pair<std::string, std::string>> swaps;  // heavy atoms whose coordinates are exchanged
  int cls = 0;         // states of one class differ only by a rotation (not by flip or tautomer)
  std::string state;   // residue state name for histidine
  double prior = 0.0;
};

struct Group {
  std::size_t res = 0;
  std::string what;  // for the report
  std::vector<State> states;
  int choice = 0;
  Vec anchor{};
};

double ramp(double d, double d1, double d2) { return d <= d1 ? 1.0 : d >= d2 ? 0.0 : (d2 - d) / (d2 - d1); }

// Score of two atoms: hydrogen bonds negative (−1 for an ideal one), overlaps positive.
double pair_score(const SAtom& a, const SAtom& b) {
  const Vec dv = sub(a.x, b.x);
  const double d2 = dot(dv, dv);
  if (d2 > 3.3 * 3.3) return 0.0;
  const double d = std::sqrt(d2);
  if (a.polar_h && b.polar_h) return d < 2.0 ? 3.0 * (2.0 - d) : 0.0;
  if (a.polar_h != b.polar_h) {
    const SAtom& h = a.polar_h ? a : b;
    const SAtom& o = a.polar_h ? b : a;
    if (o.acceptor) {
      if (d < 1.5) return 10.0 * (1.5 - d);
      if (d < 2.7) {
        const double c = dot(unit(sub(h.x, h.parent)), unit(sub(o.x, h.x)));
        if (c > 0.3) return -ramp(d, 1.9, 2.7) * c * (o.water ? 0.7 : 1.0);
      }
      return 0.0;
    }
    const double lim = o.el == 'H' ? 1.9 : o.el == 'N' ? 2.3 : 2.4;
    return d < lim ? 3.0 * (lim - d) : 0.0;
  }
  if (a.el == 'H' || b.el == 'H') {
    const double lim = a.el == 'H' && b.el == 'H' ? 1.8 : 2.3;
    return d < lim ? 3.0 * (lim - d) : 0.0;
  }
  if (a.acceptor && b.acceptor) {
    if (a.water || b.water) return d < 3.2 ? -0.3 : 0.0;  // the water can donate
    return d < 3.1 ? 3.1 - d : 0.0;                        // two lone pairs facing each other
  }
  return 0.0;
}

bool is_metal(const std::string& e) {
  return e == "Zn" || e == "Fe" || e == "Cu" || e == "Mg" || e == "Ca" || e == "Mn" || e == "Co" || e == "Ni" || e == "Na" || e == "K" ||
         e == "Cd" || e == "Hg";
}

bool is_water(const std::string& n) { return n == "HOH" || n == "WAT" || n == "DOD" || n == "H2O"; }

double bond_len(const std::string& parent_el, bool aromatic) {
  if (parent_el == "N") return 1.01;
  if (parent_el == "O") return 0.96;
  if (parent_el == "S") return 1.34;
  return aromatic ? 1.08 : 1.09;
}

Vec bisect(Vec c, Vec n1, Vec n2, double d) {
  const Vec s = add(unit(sub(n1, c)), unit(sub(n2, c)));
  return add(c, mul(unit(mul(s, -1.0)), d));
}

}  // namespace

ProtonationReport protonate(Protein& p, const ProtonationOptions& opt, Diagnostics& diag) {
  ProtonationReport rep;
  const Pka K;
  const double pH = opt.pH;
  // Remove hydrogens.
  for (auto& r : p.residues) {
    if (r.hetero) continue;
    const std::size_t before = r.atoms.size();
    r.atoms.erase(std::remove_if(r.atoms.begin(), r.atoms.end(), [](const ProteinAtom& a) { return a.element == "H" || a.element == "D"; }),
                  r.atoms.end());
    rep.hydrogens_removed += before - r.atoms.size();
  }
  std::set<std::size_t> in_ss;
  for (const auto& [a, b] : p.disulfides) in_ss.insert(a), in_ss.insert(b);
  // pKa of each group.
  std::map<std::pair<std::size_t, std::string>, double> pk;
  if (opt.structure_pka) {
    rep.pka = predict_pka(p);
    for (const auto& s : rep.pka) pk[{s.residue, s.group}] = s.pka;
  }
  auto pka_of = [&](std::size_t ri, const char* group, double fallback) {
    auto it = pk.find({ri, group});
    return it != pk.end() ? it->second : fallback;
  };
  // Metal ions.
  Grid metals(4.0);
  bool any_metal = false;
  for (const auto& r : p.residues)
    for (const auto& a : r.atoms)
      if (is_metal(a.element)) metals.add(pos(a), 0), any_metal = true;
  auto near_metal = [&](Vec x, double cut) {
    bool hit = false;
    if (any_metal) metals.query(x, cut, [&](std::uint32_t, Vec, double) { hit = true; });
    return hit;
  };
  const auto& defs = residue_defs();
  std::vector<Group> groups;
  // Polar hydrogens placed in the fixed pass, with their parents (for the network scoring).
  struct FixedH {
    std::size_t res;
    std::size_t atom;
    Vec parent;
  };
  std::vector<FixedH> fixed_polar;
  // Heavy atoms that belong to a movable group (their identity depends on the state).
  std::set<std::pair<std::size_t, std::string>> variable;

  for (std::size_t ri = 0; ri < p.residues.size(); ++ri) {
    auto& r = p.residues[ri];
    const std::string base = canonical_name(r.name);
    auto dit = defs.find(base);
    if (dit == defs.end() || r.hetero) {
      if (!r.hetero && std::find(rep.unknown_residues.begin(), rep.unknown_residues.end(), r.name) == rep.unknown_residues.end())
        rep.unknown_residues.push_back(r.name);
      continue;
    }
    const ResidueDef& def = dit->second;
    // Missing heavy atoms.
    std::vector<std::string> need = {"N", "CA", "C", "O"};
    for (const auto& [a, b] : def.bonds) need.push_back(a), need.push_back(b);
    std::sort(need.begin(), need.end());
    need.erase(std::unique(need.begin(), need.end()), need.end());
    for (const auto& n : need)
      if (!r.atom(n)) rep.missing_atoms.push_back(residue_label(r) + " " + n);
    // Terminal OXT.
    if (r.c_terminal && opt.add_oxt && !r.atom("OXT") && r.atom("C") && r.atom("O") && r.atom("CA")) {
      const Vec c = pos(*r.atom("C"));
      const Vec dir = unit(add(unit(sub(pos(*r.atom("CA")), c)), unit(sub(pos(*r.atom("O")), c))));
      ProteinAtom o;
      o.name = "OXT";
      o.element = "O";
      set_pos(o, sub(c, mul(dir, 1.25)));
      o.added = true;
      r.atoms.push_back(o);
    }
    // Protonation state.
    int charge = 0;
    std::string state = base;
    bool his_charged = false;
    if (base == "ASP") {
      const bool prot = pH < pka_of(ri, "ASP", K.asp);
      state = prot ? "ASH" : "ASP";
      charge = prot ? 0 : -1;
    } else if (base == "GLU") {
      const bool prot = pH < pka_of(ri, "GLU", K.glu);
      state = prot ? "GLH" : "GLU";
      charge = prot ? 0 : -1;
    } else if (base == "LYS") {
      const bool prot = pH < pka_of(ri, "LYS", K.lys);
      state = prot ? "LYS" : "LYN";
      charge = prot ? 1 : 0;
    } else if (base == "ARG") {
      const bool prot = pH < pka_of(ri, "ARG", K.arg);
      state = prot ? "ARG" : "ARN";
      charge = prot ? 1 : 0;
    } else if (base == "TYR") {
      const bool prot = pH < pka_of(ri, "TYR", K.tyr);
      state = prot ? "TYR" : "TYM";
      charge = prot ? 0 : -1;
    } else if (base == "CYS") {
      if (in_ss.count(ri)) state = "CYX";
      else {
        const bool metal = r.atom("SG") && near_metal(pos(*r.atom("SG")), 2.7);
        const bool prot = !metal && pH < pka_of(ri, "CYS", K.cys);
        state = prot ? "CYS" : "CYM";
        charge = prot ? 0 : -1;
      }
    } else if (base == "HIS") {
      bool metal = false;
      for (const char* n : {"ND1", "NE2", "CD2", "CE1"})
        if (const ProteinAtom* a = r.atom(n)) metal = metal || near_metal(pos(*a), 2.6);
      his_charged = !metal && pH < pka_of(ri, "HIS", K.his);
      state = his_charged ? "HIP" : "HIE";  // the tautomer is settled by the network
      charge = his_charged ? 1 : 0;
    }
    const bool nh3 = r.n_terminal && pH < pka_of(ri, "N+", K.nterm);
    const bool cooh = r.c_terminal && pH < pka_of(ri, "C-", K.cterm);
    if (r.n_terminal) charge += nh3 ? 1 : 0;
    if (r.c_terminal) charge += cooh ? 0 : -1;
    r.state = state;
    r.charge = charge;
    rep.net_charge += charge;
    // Neighbours of an atom (inside the residue, plus the peptide links).
    auto neighbours = [&](const std::string& name) {
      std::vector<Vec> out;
      std::vector<std::string> names;
      auto addn = [&](const ProteinAtom* a, const std::string& n) {
        if (a) out.push_back(pos(*a)), names.push_back(n);
      };
      std::vector<std::pair<std::string, std::string>> bonds = {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"C", "OXT"}};
      bonds.insert(bonds.end(), def.bonds.begin(), def.bonds.end());
      for (const auto& [a, b] : bonds) {
        if (name == a) addn(r.atom(b), b);
        if (name == b) addn(r.atom(a), a);
      }
      if (name == "N" && !r.n_terminal)
        for (std::size_t j = ri; j-- > 0;)
          if (p.residues[j].chain == r.chain && defs.count(canonical_name(p.residues[j].name))) {
            addn(p.residues[j].atom("C"), "C-");
            break;
          }
      return std::make_pair(out, names);
    };
    auto push = [&](const std::string& name, Vec x, const ProteinAtom* parent) {
      ProteinAtom h;
      h.name = name;
      h.element = "H";
      set_pos(h, x);
      h.added = true;
      const bool polar = parent && (parent->element == "N" || parent->element == "O" || parent->element == "S");
      const Vec px = parent ? pos(*parent) : x;
      r.atoms.push_back(h);  // may move `parent`
      if (polar) fixed_polar.push_back({ri, r.atoms.size() - 1, px});
      ++rep.hydrogens_added;
    };
    // Hydrogens fixed by the heavy atoms.
    std::vector<HDef> hs;
    if (base == "GLY") hs.push_back({"CA", {"HA2", "HA3"}, 'T', "N", false});
    else hs.push_back({"CA", {"HA"}, 'T', nullptr, false});
    if (r.n_terminal) {
      if (base == "PRO") hs.push_back({"N", nh3 ? std::vector<const char*>{"H2", "H3"} : std::vector<const char*>{"H2"}, 'T', nullptr, false});
    } else if (base != "PRO") {
      hs.push_back({"N", {"H"}, 'P', nullptr, false});
    }
    for (const auto& h : def.h) {
      const std::string par = h.parent;
      if (h.geom == 'R') continue;                                             // Ser/Thr hydroxyl: network
      if ((base == "ASN" && par == "ND2") || (base == "GLN" && par == "NE2")) continue;  // amide: network
      if (base == "HIS" && (par == "CD2" || par == "CE1")) continue;           // ring: network
      hs.push_back(h);
    }
    if (base == "ARG" && state == "ARN")
      for (auto& x : hs)
        if (std::string(x.parent) == "NH2") x.names = {"HH21"};
    for (const auto& h : hs) {
      const ProteinAtom* P = r.atom(h.parent);
      if (!P) continue;
      const Vec c = pos(*P);
      const std::string pel = P->element;
      auto [nb, nbn] = neighbours(h.parent);
      if (nb.empty()) continue;
      const bool aromatic = h.geom == 'P' && pel == "C";
      const double d = bond_len(pel, aromatic);
      const std::size_t n = h.names.size();
      auto P2 = [&] { return r.atom(h.parent); };
      if (h.geom == 'T') {
        if (nb.size() >= 3 && n == 1) {
          Vec s{0, 0, 0};
          for (const auto& v : nb) s = add(s, unit(sub(v, c)));
          push(h.names[0], add(c, mul(unit(mul(s, -1.0)), d)), P2());
        } else if (nb.size() == 2 && n == 2) {
          // Prochirality: H2 has a negative sign of ((prev − C) × (next − C)) · (H − C).
          Vec prev = nb[0], next = nb[1];
          if (h.prev && nbn[1] == h.prev) std::swap(prev, next);
          const Vec u1 = unit(sub(prev, c)), u2 = unit(sub(next, c));
          const Vec b = unit(mul(add(u1, u2), -1.0));
          const Vec nn = unit(cross(u1, u2));
          const double a = 54.735 * std::numbers::pi / 180.0;
          const Vec ha = add(c, mul(add(mul(b, std::cos(a)), mul(nn, std::sin(a))), d));
          const Vec hb = add(c, mul(add(mul(b, std::cos(a)), mul(nn, -std::sin(a))), d));
          const bool a_is_2 = dot(cross(sub(prev, c), sub(next, c)), sub(ha, c)) < 0;
          push(h.names[0], a_is_2 ? ha : hb, P2());
          push(h.names[1], a_is_2 ? hb : ha, P2());
        } else if (nb.size() == 2 && n == 1) {
          push(h.names[0], bisect(c, nb[0], nb[1], d), P2());
        } else if (nb.size() == 1) {
          // Methyl: staggered against a reference atom bonded to the neighbour.
          const std::string nbname = nbn[0];
          auto [nb2, nb2n] = neighbours(nbname == "C-" ? "C" : nbname);
          Vec ref = c;
          for (std::size_t k = 0; k < nb2.size(); ++k)
            if (nb2n[k] != h.parent) {
              ref = nb2[k];
              break;
            }
          const double dihs[3] = {180.0, 60.0, -60.0};
          for (std::size_t k = 0; k < n && k < 3; ++k) push(h.names[k], place(ref, nb[0], c, d, 109.47, dihs[k]), P2());
        }
      } else if (h.geom == 'P') {
        if (nb.size() >= 2 && n == 1) {
          push(h.names[0], bisect(c, nb[0], nb[1], d), P2());
        } else if (nb.size() == 1) {
          // sp2 NH2 (or NH): in the plane of the neighbour X and the reference atom R bonded to X.
          const Vec X = nb[0];
          auto [nb2, nb2n] = neighbours(nbn[0]);
          Vec R = c;
          for (std::size_t k = 0; k < nb2.size(); ++k)
            if ((h.prev && nb2n[k] == h.prev) || (!h.prev && nb2n[k] != h.parent)) {
              R = nb2[k];
              break;
            }
          const Vec hc = place(R, X, c, d, 120.0, 0.0), ht = place(R, X, c, d, 120.0, 180.0);
          if (n == 2) {
            push(h.names[0], h.first_cis ? hc : ht, P2());
            push(h.names[1], h.first_cis ? ht : hc, P2());
          } else {
            push(h.names[0], ht, P2());
          }
        }
      }
    }

    // ---- Movable groups of this residue ---------------------------------------------------------------
    auto at = [&](const char* n) -> std::optional<Vec> {
      const ProteinAtom* a = r.atom(n);
      return a ? std::optional<Vec>(pos(*a)) : std::nullopt;
    };
    auto hatom = [&](Vec x, Vec parent) {
      SAtom s;
      s.x = x, s.el = 'H', s.polar_h = true, s.parent = parent, s.res = std::uint32_t(ri);
      return s;
    };
    auto heavy = [&](Vec x, char el, bool acceptor) {
      SAtom s;
      s.x = x, s.el = el, s.acceptor = acceptor, s.res = std::uint32_t(ri);
      return s;
    };
    auto new_group = [&](const std::string& what, Vec anchor) -> Group& {
      Group g;
      g.res = ri;
      g.what = what;
      g.anchor = anchor;
      groups.push_back(std::move(g));
      return groups.back();
    };
    // One hydrogen on O or S rotating about X–P.
    auto rotor = [&](const char* parent, const char* x, const char* ref, const char* hname, bool phenol) {
      auto P = at(parent), X = at(x), R = at(ref);
      if (!P || !X || !R) return;
      const bool sulfur = r.atom(parent)->element == "S";
      const double d = sulfur ? 1.34 : 0.96, angle = sulfur ? 96.0 : 109.47;
      Group& g = new_group(std::string(parent), *P);
      for (int t = 0; t < 360; t += phenol ? 180 : 30) {
        State s;
        const Vec h = place(*R, *X, *P, d, angle, double(t));
        s.atoms.push_back(hatom(h, *P));
        s.h.push_back({hname, h});
        if (!phenol) s.prior = 0.1 * (1.0 + std::cos(3.0 * t * std::numbers::pi / 180.0));  // staggered preferred
        g.states.push_back(std::move(s));
      }
      if (!phenol) std::rotate(g.states.begin(), g.states.begin() + 6, g.states.end());  // anti first
    };
    // NH3+ (three H) or NH2 (two of the three staggered positions).
    auto amine = [&](const char* parent, const char* x, const char* ref, std::vector<const char*> names, bool charged) {
      auto P = at(parent), X = at(x), R = at(ref);
      if (!P || !X || !R) return;
      Group& g = new_group(std::string(parent), *P);
      if (charged) {
        const double offs[4] = {0.0, 30.0, 60.0, 90.0}, pri[4] = {0.0, 0.15, 0.3, 0.15};
        for (int k = 0; k < 4; ++k) {
          State s;
          const double base_d[3] = {180.0, 60.0, -60.0};
          for (int j = 0; j < 3; ++j) {
            const Vec h = place(*R, *X, *P, 1.01, 109.47, base_d[j] + offs[k]);
            s.atoms.push_back(hatom(h, *P));
            s.h.push_back({names[j], h});
          }
          s.prior = pri[k];
          g.states.push_back(std::move(s));
        }
      } else {
        const double base_d[3] = {180.0, 60.0, -60.0};
        for (int skip = 2; skip >= 0; --skip) {
          State s;
          int j = 0;
          for (int q = 0; q < 3; ++q) {
            if (q == skip) continue;
            const Vec h = place(*R, *X, *P, 1.01, 109.47, base_d[q]);
            s.atoms.push_back(hatom(h, *P));
            s.h.push_back({names[j++], h});
          }
          g.states.push_back(std::move(s));
        }
      }
    };
    if (base == "SER") rotor("OG", "CB", "CA", "HG", false);
    if (base == "THR") rotor("OG1", "CB", "CA", "HG1", false);
    if (base == "CYS" && state == "CYS") rotor("SG", "CB", "CA", "HG", false);
    if (base == "TYR" && state == "TYR") rotor("OH", "CZ", "CE1", "HH", true);
    if (base == "LYS") amine("NZ", "CE", "CD", state == "LYS" ? std::vector<const char*>{"HZ1", "HZ2", "HZ3"} : std::vector<const char*>{"HZ1", "HZ2"}, state == "LYS");
    if (r.n_terminal && base != "PRO")
      amine("N", "CA", "C", nh3 ? std::vector<const char*>{"H1", "H2", "H3"} : std::vector<const char*>{"H1", "H2"}, nh3);
    // Amide of Asn/Gln: as deposited or flipped.
    auto amide = [&](const char* o, const char* n, const char* c, const char* h1, const char* h2) {
      auto O = at(o), N = at(n), C = at(c);
      if (!O || !N || !C) return;
      Group& g = new_group(base, *C);
      for (int flip = 0; flip < (opt.flips ? 2 : 1); ++flip) {
        State s;
        const Vec xo = flip ? *N : *O, xn = flip ? *O : *N;
        s.atoms.push_back(heavy(xo, 'O', true));
        s.atoms.push_back(heavy(xn, 'N', false));
        // h1 is trans to O, h2 cis (PDB names HD21/HE21 and HD22/HE22).
        const Vec ht = place(xo, *C, xn, 1.01, 120.0, 180.0), hc = place(xo, *C, xn, 1.01, 120.0, 0.0);
        s.atoms.push_back(hatom(ht, xn));
        s.atoms.push_back(hatom(hc, xn));
        s.h = {{h1, ht}, {h2, hc}};
        if (flip) s.swaps = {{o, n}}, s.prior = opt.flip_penalty;
        s.cls = flip;
        g.states.push_back(std::move(s));
      }
      variable.insert({ri, o}), variable.insert({ri, n});
    };
    if (base == "ASN") amide("OD1", "ND2", "CG", "HD21", "HD22");
    if (base == "GLN") amide("OE1", "NE2", "CD", "HE21", "HE22");
    if (base == "HIS") {
      auto CG = at("CG"), ND1 = at("ND1"), CD2 = at("CD2"), CE1 = at("CE1"), NE2 = at("NE2");
      if (CG && ND1 && CD2 && CE1 && NE2) {
        Group& g = new_group("HIS", *CG);
        for (int flip = 0; flip < (opt.flips ? 2 : 1); ++flip) {
          // Ring positions: a (next to CG), b (next to a) on one side; c, d on the other.
          const Vec nd = flip ? *CD2 : *ND1, cd = flip ? *ND1 : *CD2, ce = flip ? *NE2 : *CE1, ne = flip ? *CE1 : *NE2;
          const bool nd_metal = near_metal(nd, 2.6), ne_metal = near_metal(ne, 2.6);
          for (int taut = 0; taut < 3; ++taut) {  // 0 HIE, 1 HID, 2 HIP
            if ((taut == 2) != his_charged) continue;
            const bool hd = taut != 0, he = taut != 1;
            if ((hd && nd_metal) || (he && ne_metal)) continue;
            State s;
            s.atoms.push_back(heavy(nd, 'N', !hd));
            s.atoms.push_back(heavy(ne, 'N', !he));
            s.atoms.push_back(heavy(cd, 'C', false));
            s.atoms.push_back(heavy(ce, 'C', false));
            const Vec hcd = bisect(cd, *CG, ne, 1.08), hce = bisect(ce, nd, ne, 1.08);
            SAtom a1 = hatom(hcd, cd), a2 = hatom(hce, ce);
            a1.polar_h = a2.polar_h = false;
            s.atoms.push_back(a1), s.atoms.push_back(a2);
            s.h = {{"HD2", hcd}, {"HE1", hce}};
            if (hd) {
              const Vec h = bisect(nd, *CG, ce, 1.01);
              s.atoms.push_back(hatom(h, nd));
              s.h.push_back({"HD1", h});
            }
            if (he) {
              const Vec h = bisect(ne, cd, ce, 1.01);
              s.atoms.push_back(hatom(h, ne));
              s.h.push_back({"HE2", h});
            }
            if (flip) s.swaps = {{"ND1", "CD2"}, {"CE1", "NE2"}};
            s.prior = (flip ? opt.flip_penalty : 0.0) + (taut == 1 ? 0.2 : 0.0);  // HIE is the more common tautomer
            s.cls = flip * 4 + taut;
            s.state = taut == 0 ? "HIE" : taut == 1 ? "HID" : "HIP";
            g.states.push_back(std::move(s));
          }
        }
        if (g.states.empty()) {
          groups.pop_back();  // both nitrogens on metal ions: no N–H, ring hydrogens only
          push("HD2", bisect(*CD2, *CG, *NE2, 1.08), nullptr);
          push("HE1", bisect(*CE1, *ND1, *NE2, 1.08), nullptr);
          r.state = "HIM";
        } else {
          for (const char* n : {"ND1", "CD2", "CE1", "NE2"}) variable.insert({ri, n});
        }
      }
    }
    // Protonated carboxyl groups: the hydrogen on either oxygen, syn.
    auto acid = [&](const char* o1, const char* o2, const char* c, const char* hname) {
      auto O1 = at(o1), O2 = at(o2), C = at(c);
      if (!O1 || !O2 || !C) return;
      Group& g = new_group(std::string(c) + "OOH", *C);
      for (int flip = 0; flip < 2; ++flip) {
        State s;
        const Vec oh = flip ? *O1 : *O2, oc = flip ? *O2 : *O1;
        const Vec h = place(oc, *C, oh, 0.96, 108.0, 0.0);
        s.atoms.push_back(hatom(h, oh));
        s.h = {{hname, h}};
        if (flip) s.swaps = {{o1, o2}};
        g.states.push_back(std::move(s));
      }
    };
    if (state == "ASH") acid("OD1", "OD2", "CG", "HD2");
    if (state == "GLH") acid("OE1", "OE2", "CD", "HE2");
    if (cooh && r.atom("OXT")) acid("O", "OXT", "C", "HXT");
  }

  // ---- Network ----------------------------------------------------------------------------------------
  if (!groups.empty()) {
    // Environment: every atom that does not change with the states.
    std::vector<SAtom> env;
    Grid egrid(4.0);
    {
      std::map<std::pair<std::size_t, std::size_t>, Vec> polar;
      for (const auto& f : fixed_polar) polar[{f.res, f.atom}] = f.parent;
      for (std::size_t i = 0; i < p.residues.size(); ++i) {
        const auto& r = p.residues[i];
        const bool water = is_water(r.name);
        for (std::size_t k = 0; k < r.atoms.size(); ++k) {
          const auto& a = r.atoms[k];
          if (variable.count({i, a.name})) continue;
          if (water && (a.element == "H" || a.element == "D")) continue;  // waters are treated as free to turn
          SAtom s;
          s.x = pos(a);
          s.res = std::uint32_t(i);
          s.el = a.element == "D" ? 'H' : a.element.size() == 1 ? a.element[0] : 'X';
          s.acceptor = a.element == "O";
          s.water = water;
          auto it = polar.find({i, k});
          if (it != polar.end()) s.polar_h = true, s.parent = it->second;
          egrid.add(s.x, std::uint32_t(env.size()));
          env.push_back(s);
        }
      }
    }
    const std::size_t ng = groups.size();
    std::vector<std::vector<double>> self(ng);
    for (std::size_t g = 0; g < ng; ++g) {
      Group& G = groups[g];
      self[g].resize(G.states.size());
      for (std::size_t s = 0; s < G.states.size(); ++s) {
        double e = G.states[s].prior;
        for (const auto& a : G.states[s].atoms)
          egrid.query(a.x, 3.3, [&](std::uint32_t id, Vec, double) {
            if (env[id].res == a.res) return;
            e += pair_score(a, env[id]);
          });
        self[g][s] = e;
      }
      G.choice = int(std::min_element(self[g].begin(), self[g].end()) - self[g].begin());
    }
    if (opt.hydrogen_bonds) {
      // Pair tables between groups in contact.
      struct Edge {
        std::size_t a, b;
        std::vector<double> e;  // [sa * nb + sb]
      };
      std::vector<Edge> edges;
      std::vector<std::vector<std::size_t>> adj(ng);
      Grid ggrid(8.0);
      for (std::size_t g = 0; g < ng; ++g) ggrid.add(groups[g].anchor, std::uint32_t(g));
      for (std::size_t g = 0; g < ng; ++g)
        ggrid.query(groups[g].anchor, 9.0, [&](std::uint32_t h, Vec, double) {
          if (h <= g) return;
          const Group& A = groups[g];
          const Group& B = groups[h];
          Edge ed{g, h, std::vector<double>(A.states.size() * B.states.size(), 0.0)};
          bool any = false;
          for (std::size_t sa = 0; sa < A.states.size(); ++sa)
            for (std::size_t sb = 0; sb < B.states.size(); ++sb) {
              double e = 0.0;
              for (const auto& x : A.states[sa].atoms)
                for (const auto& y : B.states[sb].atoms) e += pair_score(x, y);
              ed.e[sa * B.states.size() + sb] = e;
              any = any || e != 0.0;
            }
          if (!any) return;
          adj[g].push_back(edges.size());
          adj[h].push_back(edges.size());
          edges.push_back(std::move(ed));
        });
      auto edge_energy = [&](const Edge& ed, std::size_t g, int sg, int other) {
        const std::size_t nb = groups[ed.b].states.size();
        return ed.a == g ? ed.e[std::size_t(sg) * nb + std::size_t(other)] : ed.e[std::size_t(other) * nb + std::size_t(sg)];
      };
      auto local_energy = [&](std::size_t g, int s) {
        double e = self[g][std::size_t(s)];
        for (std::size_t k : adj[g]) {
          const Edge& ed = edges[k];
          const std::size_t o = ed.a == g ? ed.b : ed.a;
          e += edge_energy(ed, g, s, groups[o].choice);
        }
        return e;
      };
      // Connected clusters.
      std::vector<int> comp(ng, -1);
      std::vector<std::vector<std::size_t>> clusters;
      for (std::size_t g = 0; g < ng; ++g) {
        if (comp[g] >= 0 || adj[g].empty()) continue;
        std::vector<std::size_t> stack = {g}, members;
        comp[g] = int(clusters.size());
        while (!stack.empty()) {
          const std::size_t u = stack.back();
          stack.pop_back();
          members.push_back(u);
          for (std::size_t k : adj[u]) {
            const std::size_t v = edges[k].a == u ? edges[k].b : edges[k].a;
            if (comp[v] < 0) comp[v] = comp[u], stack.push_back(v);
          }
        }
        std::sort(members.begin(), members.end());
        clusters.push_back(std::move(members));
      }
      for (const auto& cl : clusters) {
        // Iterated improvement from the individually best states.
        for (int sweep = 0; sweep < 30; ++sweep) {
          bool changed = false;
          for (std::size_t g : cl) {
            int best = groups[g].choice;
            double be = local_energy(g, best);
            for (int s = 0; s < int(groups[g].states.size()); ++s) {
              const double e = local_energy(g, s);
              if (e < be - 1e-9) be = e, best = s;
            }
            if (best != groups[g].choice) groups[g].choice = best, changed = true;
          }
          if (!changed) break;
        }
        double combos = 1.0;
        for (std::size_t g : cl) combos *= double(groups[g].states.size());
        if (combos > 1e5 || cl.size() < 2) continue;
        // Exhaustive search with a bound (energies of earlier members are complete, later ones bounded
        // below by their best self energy plus best edges).
        std::map<std::size_t, std::size_t> order;
        for (std::size_t i = 0; i < cl.size(); ++i) order[cl[i]] = i;
        auto total = [&] {
          double e = 0.0;
          for (std::size_t g : cl) e += self[g][std::size_t(groups[g].choice)];
          for (std::size_t g : cl)
            for (std::size_t k : adj[g])
              if (edges[k].a == g) e += edge_energy(edges[k], g, groups[g].choice, groups[edges[k].b].choice);
          return e;
        };
        double best_e = total();
        std::vector<int> best_choice, cur(cl.size(), 0);
        for (std::size_t g : cl) best_choice.push_back(groups[g].choice);
        std::vector<double> lower(cl.size() + 1, 0.0);  // lower bound of members i.. (self and edges to later/earlier)
        for (std::size_t i = cl.size(); i-- > 0;) {
          const std::size_t g = cl[i];
          double m = 1e300;
          for (std::size_t s = 0; s < groups[g].states.size(); ++s) {
            double e = self[g][s];
            for (std::size_t k : adj[g]) {
              const Edge& ed = edges[k];
              const std::size_t o = ed.a == g ? ed.b : ed.a;
              if (order[o] > i) continue;  // counted when the later member is placed
              double mn = 1e300;
              for (int t = 0; t < int(groups[o].states.size()); ++t) mn = std::min(mn, edge_energy(ed, g, int(s), t));
              e += mn;
            }
            m = std::min(m, e);
          }
          lower[i] = lower[i + 1] + m;
        }
        std::function<void(std::size_t, double)> dfs = [&](std::size_t i, double e) {
          if (e + lower[i] >= best_e - 1e-9) return;
          if (i == cl.size()) {
            best_e = e;
            best_choice = cur;
            return;
          }
          const std::size_t g = cl[i];
          for (int s = 0; s < int(groups[g].states.size()); ++s) {
            double de = self[g][std::size_t(s)];
            for (std::size_t k : adj[g]) {
              const Edge& ed = edges[k];
              const std::size_t o = ed.a == g ? ed.b : ed.a;
              if (order[o] < i) de += edge_energy(ed, g, s, cur[order[o]]);
            }
            cur[i] = s;
            dfs(i + 1, e + de);
          }
        };
        dfs(0, 0.0);
        for (std::size_t i = 0; i < cl.size(); ++i) groups[cl[i]].choice = best_choice[i];
      }
      // Decisions that the structure does not settle.
      for (std::size_t g = 0; g < ng; ++g) {
        const Group& G = groups[g];
        const State& chosen = G.states[std::size_t(G.choice)];
        const double e0 = local_energy(g, G.choice);
        double alt = 1e300;
        for (int s = 0; s < int(G.states.size()); ++s)
          if (G.states[std::size_t(s)].cls != chosen.cls) alt = std::min(alt, local_energy(g, s) - G.states[std::size_t(s)].prior + chosen.prior);
        if (alt - e0 < 0.3) {
          const auto& r = p.residues[G.res];
          rep.undetermined.push_back(residue_label(r) + (G.what == "HIS" ? ": tautomer or ring orientation" : ": amide orientation"));
        }
      }
    }
    // Apply.
    for (const Group& G : groups) {
      auto& r = p.residues[G.res];
      const State& s = G.states[std::size_t(G.choice)];
      for (const auto& [a, b] : s.swaps) {
        ProteinAtom* x = r.atom(a);
        ProteinAtom* y = r.atom(b);
        if (!x || !y) continue;
        std::swap(x->x, y->x), std::swap(x->y, y->y), std::swap(x->z, y->z);
        std::swap(x->bfactor, y->bfactor), std::swap(x->occupancy, y->occupancy);
      }
      if (!s.swaps.empty() && G.what != "CGOOH" && G.what != "CDOOH" && G.what != "COOH") rep.flipped.push_back(residue_label(r));
      if (!s.state.empty()) r.state = s.state;
      for (const auto& [name, x] : s.h) {
        ProteinAtom h;
        h.name = name;
        h.element = "H";
        set_pos(h, x);
        h.added = true;
        r.atoms.push_back(h);
        ++rep.hydrogens_added;
      }
    }
  }
  if (!rep.missing_atoms.empty())
    diag.warn(DiagCode::MalformedRecord, std::to_string(rep.missing_atoms.size()) +
                                             " heavy atoms of standard residues are missing (e.g. " + rep.missing_atoms[0] +
                                             "); their hydrogens are not placed");
  if (!rep.unknown_residues.empty()) {
    std::string l;
    for (const auto& u : rep.unknown_residues) l += (l.empty() ? "" : ", ") + u;
    diag.warn(DiagCode::FileContentIgnored, "no template for residue(s) " + l + ": left without hydrogens");
  }
  return rep;
}

}  // namespace chemlarp
