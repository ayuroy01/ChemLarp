// Missing heavy atoms: side chains from a rotamer library, missing residues by loop closure.
#include <cstdio>
#include <cstring>
#include <set>
#include <tuple>

#include "protein_internal.hpp"

namespace chemlarp {

using namespace detail;

namespace {

struct SideChainAtom {
  const char *res, *atom, *a, *b, *c;
  double bond, angle, dihedral;
  int chi;
};
struct Rotamer {
  const char* res;
  double freq;
  int nchi;
  double chi[4];
};
struct Chi1Bias {
  const char* res;
  int cls;
  double f[3];
};
#include "protein_geometry.inc"

// Backbone geometry (Engh & Huber, Acta Cryst. A47 (1991) 392).
constexpr double kNCA = 1.458, kCAC = 1.525, kCN = 1.329, kCO = 1.231;
constexpr double kNCAC = 111.0, kCACN = 116.2, kCNCA = 121.7, kCACO = 120.8;

struct Span {
  const SideChainAtom* begin;
  const SideChainAtom* end;
  const Rotamer* rbegin;
  const Rotamer* rend;
};
const std::map<std::string, Span>& spans() {
  static const std::map<std::string, Span> m = [] {
    std::map<std::string, Span> o;
    for (const auto& a : kSideChainAtoms) {
      auto& s = o[a.res];
      if (!s.begin) s.begin = &a;
      s.end = &a + 1;
    }
    for (const auto& r : kRotamers) {
      auto& s = o[r.res];
      if (!s.rbegin) s.rbegin = &r;
      s.rend = &r + 1;
    }
    o["GLY"];
    return o;
  }();
  return m;
}


double angle_diff(double a, double b) {
  double d = std::fmod(a - b, 360.0);
  if (d > 180.0) d -= 360.0;
  if (d < -180.0) d += 360.0;
  return std::fabs(d);
}

struct EnvAtom {
  Vec x;
  std::uint32_t res;
  char el;  // C N O S or X
  bool polar;
};

char el_code(const std::string& e) { return e == "C" || e == "N" || e == "O" || e == "S" || e == "H" ? e[0] : 'X'; }

// Overlap and contact score of two heavy atoms at distance d (weights set on a development set of nine
// structures; see docs/ARCHITECTURE.md).
double pair_score(char e1, char e2, double d) {
  auto rad = [](char e) { return e == 'C' ? 1.7 : e == 'N' ? 1.55 : e == 'O' ? 1.5 : e == 'S' ? 1.8 : e == 'H' ? 1.0 : 1.8; };
  const bool p1 = e1 == 'N' || e1 == 'O', p2 = e2 == 'N' || e2 == 'O';
  double d0 = rad(e1) + rad(e2) - 0.05;      // allowed approach
  if (p1 && p2) d0 = 2.6;                    // hydrogen-bonded pairs
  if (e1 == 'S' && e2 == 'S') d0 = 2.0;      // disulfide
  if (d < d0) {
    const double t = (d0 - d) / 0.4;
    return std::min(30.0 * t * t, 60.0);
  }
  if (p1 && p2 && d < 3.2) return -0.6;      // hydrogen bond
  if (!p1 && !p2 && d < d0 + 1.6) return -0.15;  // packing
  return 0.0;
}

struct Job {
  std::size_t res;                 // residue index
  const Span* span;
  std::vector<int> missing;        // indices into the span of atoms to place
  std::vector<std::array<double, 4>> cand;  // chi sets
  std::vector<double> prior;
  std::vector<std::vector<Vec>> coords;      // per candidate: positions of the missing atoms
  int choice = 0;
  std::vector<std::size_t> near;   // other jobs close enough to interact
  Vec ca;
};

// Positions of the atoms in `missing` for the chi set; existing atoms are read from the residue.
bool build_atoms(const ProteinResidue& r, const Span& sp, const std::vector<int>& missing, const std::array<double, 4>& chi,
                 std::vector<Vec>& out) {
  out.clear();
  std::map<std::string, Vec> placed;
  auto get = [&](const char* n, Vec& x) {
    auto it = placed.find(n);
    if (it != placed.end()) {
      x = it->second;
      return true;
    }
    const ProteinAtom* a = r.atom(n);
    if (!a) return false;
    x = pos(*a);
    return true;
  };
  for (int k : missing) {
    const SideChainAtom& t = sp.begin[k];
    Vec A, B, C;
    if (!get(t.a, A) || !get(t.b, B) || !get(t.c, C)) return false;
    const Vec x = place(A, B, C, t.bond, t.angle, t.chi >= 0 ? chi[t.chi] : t.dihedral);
    placed[t.atom] = x;
    out.push_back(x);
  }
  return true;
}

void recompute_topology(Protein& p) {
  const auto& defs = residue_defs();
  p.chain_breaks.clear();
  std::map<std::string, std::size_t> last;
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    auto& r = p.residues[i];
    if (!defs.count(canonical_name(r.name))) continue;
    r.n_terminal = r.c_terminal = true;
    auto it = last.find(r.chain);
    if (it != last.end()) {
      auto& q = p.residues[it->second];
      const ProteinAtom* N = r.atom("N");
      const ProteinAtom* C = q.atom("C");
      if (N && C && dist(*C, *N) < 2.0) r.n_terminal = q.c_terminal = false;
      else p.chain_breaks.push_back(i);
    }
    last[r.chain] = i;
  }
}

// xorshift generator: the same numbers on every platform.
struct Rng {
  std::uint64_t s;
  double uniform() {
    s ^= s << 13, s ^= s >> 7, s ^= s << 17;
    return double(s >> 11) / 9007199254740992.0;
  }
};

Vec rotate_about(Vec x, Vec origin, Vec axis, double c, double s) {
  const Vec r = sub(x, origin);
  const Vec rot = add(add(mul(r, c), mul(cross(axis, r), s)), mul(axis, dot(axis, r) * (1.0 - c)));
  return add(origin, rot);
}

struct LoopResult {
  std::vector<std::array<Vec, 4>> bb;  // N CA C O per residue
  double closure = 0.0, score = 1e30;
};

}  // namespace

RepairReport repair(Protein& p, const RepairOptions& opt, Diagnostics& diag) {
  RepairReport rep;
  const auto& defs = residue_defs();
  const auto& sp = spans();

  // Residues whose backbone is incomplete cannot anchor a side chain: they are rebuilt whole, like missing
  // residues (a stub at a chain end always).
  std::set<std::tuple<std::string, int, char>> forced;
  if (opt.side_chains || opt.loops) {
    std::vector<std::size_t> remap(p.residues.size(), SIZE_MAX);
    std::vector<ProteinResidue> keep;
    std::map<std::string, int> per_chain;
    std::map<std::string, int> incomplete;
    for (const auto& r : p.residues)
      if (!r.hetero && defs.count(canonical_name(r.name))) {
        ++per_chain[r.chain];
        if (!r.atom("N") || !r.atom("CA") || !r.atom("C")) ++incomplete[r.chain];
      }
    for (std::size_t i = 0; i < p.residues.size(); ++i) {
      auto& r = p.residues[i];
      const bool amino = !r.hetero && defs.count(canonical_name(r.name));
      // Not for CA-only or otherwise skeletal chains (more than a tenth incomplete): those are left as read.
      if (amino && incomplete[r.chain] * 10 <= per_chain[r.chain] && (!r.atom("N") || !r.atom("CA") || !r.atom("C"))) {
        p.missing_residues.push_back({canonical_name(r.name), r.chain, r.number, r.icode});
        forced.insert({r.chain, r.number, r.icode});
        continue;
      }
      remap[i] = keep.size();
      keep.push_back(std::move(r));
    }
    p.residues = std::move(keep);
    if (!forced.empty()) {
      std::vector<std::pair<std::size_t, std::size_t>> ss;
      for (const auto& [a, b] : p.disulfides)
        if (remap[a] != SIZE_MAX && remap[b] != SIZE_MAX) ss.push_back({remap[a], remap[b]});
      p.disulfides = std::move(ss);
      recompute_topology(p);
    }
  }

  // ---- Missing residues: loops and (optionally) chain ends ------------------------------------------------
  if ((opt.loops || opt.termini || !forced.empty()) && !p.missing_residues.empty()) {
    // Insertion point of each missing residue: the index of the first residue of its chain that follows it.
    struct Run {
      std::size_t before = SIZE_MAX, after = SIZE_MAX;  // residue indices with coordinates
      std::vector<MissingResidue> res;
    };
    std::map<std::string, std::vector<std::size_t>> chain_res;
    for (std::size_t i = 0; i < p.residues.size(); ++i)
      if (defs.count(canonical_name(p.residues[i].name))) chain_res[p.residues[i].chain].push_back(i);
    std::map<std::pair<std::string, std::size_t>, Run> runs;
    for (const auto& m : p.missing_residues) {
      auto it = chain_res.find(m.chain);
      if (it == chain_res.end()) {
        continue;  // a chain with no coordinates at all
      }
      const auto& idx = it->second;
      std::size_t k = 0;
      while (k < idx.size()) {
        const auto& r = p.residues[idx[k]];
        if (r.number > m.number || (r.number == m.number && r.icode > m.icode)) break;
        ++k;
      }
      Run& run = runs[{m.chain, k}];
      run.before = k > 0 ? idx[k - 1] : SIZE_MAX;
      run.after = k < idx.size() ? idx[k] : SIZE_MAX;
      run.res.push_back(m);
    }
    // Environment: heavy atoms of everything present.
    std::vector<EnvAtom> env;
    Grid grid(4.0);
    auto add_env = [&](Vec x, std::uint32_t res, char el) {
      grid.add(x, std::uint32_t(env.size()));
      env.push_back({x, res, el, el == 'N' || el == 'O'});
    };
    for (std::size_t i = 0; i < p.residues.size(); ++i)
      for (const auto& a : p.residues[i].atoms)
        if (a.element != "H" && a.element != "D") add_env(pos(a), std::uint32_t(i), el_code(a.element));
    Rng rng{0x9E3779B97F4A7C15ull};
    struct Built {
      std::size_t insert_at;  // index in the original residue list before which to insert
      std::vector<ProteinResidue> res;
    };
    std::vector<Built> built;
    std::set<std::tuple<std::string, int, char>> done;
    auto label = [](const MissingResidue& m) {
      return (m.chain.empty() ? std::string("_") : m.chain) + " " + std::to_string(m.number) + (m.icode != ' ' ? std::string(1, m.icode) : "") + " " + m.name;
    };
    for (auto& [key, run] : runs) {
      std::stable_sort(run.res.begin(), run.res.end(), [](const MissingResidue& a, const MissingResidue& b) {
        return a.number != b.number ? a.number < b.number : a.icode < b.icode;
      });
      const bool internal = run.before != SIZE_MAX && run.after != SIZE_MAX;
      auto is_forced = [&](const MissingResidue& m) { return forced.count({m.chain, m.number, m.icode}) > 0; };
      if (!internal && !opt.termini) {
        // Only the stubs next to the chain: residues that had some atoms.
        if (run.before != SIZE_MAX) {
          std::size_t k = 0;
          while (k < run.res.size() && is_forced(run.res[k])) ++k;
          run.res.resize(k);
        } else {
          std::size_t k = run.res.size();
          while (k > 0 && is_forced(run.res[k - 1])) --k;
          run.res.erase(run.res.begin(), run.res.begin() + long(k));
        }
        if (run.res.empty()) continue;
      }
      const std::size_t n = run.res.size();
      const bool any_forced = std::any_of(run.res.begin(), run.res.end(), is_forced);
      const std::string range = label(run.res.front()) + (n > 1 ? " to " + label(run.res.back()) : "");
      if (internal && !opt.loops && !any_forced) continue;
      bool known = true;
      for (const auto& m : run.res) known = known && defs.count(canonical_name(m.name));
      if (!known) {
        rep.not_built.push_back(range + ": nonstandard residue in the gap");
        continue;
      }
      const ProteinResidue* A = run.before != SIZE_MAX ? &p.residues[run.before] : nullptr;
      const ProteinResidue* B = run.after != SIZE_MAX ? &p.residues[run.after] : nullptr;
      auto has_bb = [](const ProteinResidue* r) { return r && r->atom("N") && r->atom("CA") && r->atom("C"); };
      if ((A && !has_bb(A)) || (B && !has_bb(B))) {
        rep.not_built.push_back(range + ": neighbouring residue has an incomplete backbone");
        continue;
      }
      if (internal) {
        const double gap = norm(sub(pos(*B->atom("N")), pos(*A->atom("C"))));
        if (gap > 3.8 * double(n) + 1.4) {
          char b[96];
          std::snprintf(b, sizeof b, ": the %.1f Å gap is too long for %zu residue(s)", gap, n);
          rep.not_built.push_back(range + b);
          continue;
        }
      }
      // Torsion sampling: residue k has phi[k], psi[k]; psiA is the psi of residue A.
      const int trials = std::max(1, n > 60 ? opt.loop_trials / 4 : n <= 6 ? opt.loop_trials * 3 : opt.loop_trials);
      LoopResult best;
      for (int trial = 0; trial < trials; ++trial) {
        std::vector<double> phi(n), psi(n);
        for (std::size_t k = 0; k < n; ++k) {
          const double u = rng.uniform();
          const bool gly = run.res[k].name == "GLY", pro = run.res[k].name == "PRO";
          double f, s;
          if (u < 0.35) f = -63, s = -43;
          else if (u < 0.75) f = -120, s = 130;
          else f = -75, s = 150;
          if (gly && rng.uniform() < 0.4) f = 80, s = rng.uniform() < 0.5 ? 10 : -170;
          f += 50.0 * (rng.uniform() - 0.5), s += 50.0 * (rng.uniform() - 0.5);
          if (pro) f = -65.0;
          phi[k] = f, psi[k] = s;
        }
        // Atoms: 3 per residue (N, CA, C), then for a closed loop the virtual N and CA of residue B.
        std::vector<Vec> x;
        if (A) {
          const Vec nA = pos(*A->atom("N")), caA = pos(*A->atom("CA")), cA = pos(*A->atom("C"));
          Vec prevN = nA, prevCA = caA, prevC = cA;
          double psi_prev;
          bool psi_fixed = false;
          Vec first_n{};
          if (const ProteinAtom* O = A->atom("O")) {  // N of the next residue is trans to O
            first_n = place(pos(*O), caA, cA, kCN, kCACN, 180.0);
            psi_fixed = true;
          }
          psi_prev = 120.0 + 60.0 * (rng.uniform() - 0.5);
          for (std::size_t k = 0; k <= n; ++k) {
            const bool last = k == n;
            if (last && !B) break;
            Vec N = (k == 0 && psi_fixed) ? first_n : place(prevN, prevCA, prevC, kCN, kCACN, psi_prev);
            Vec CA = place(prevCA, prevC, N, kNCA, kCNCA, 180.0);
            x.push_back(N), x.push_back(CA);
            if (last) break;
            Vec C = place(prevC, N, CA, kCAC, kNCAC, phi[k]);
            x.push_back(C);
            prevN = N, prevCA = CA, prevC = C, psi_prev = psi[k];
          }
          if (B) {
            // Cyclic coordinate descent on phi/psi of the built residues toward N and CA of residue B.
            const Vec tN = pos(*B->atom("N")), tCA = pos(*B->atom("CA"));
            const std::size_t m = x.size();  // 3n + 2
            double err = 1e9;
            for (int sweep = 0; sweep < 1500 && err > 0.01; ++sweep) {
              for (std::size_t k = 0; k < n; ++k) {
                for (int which = 0; which < 2; ++which) {  // phi: axis N-CA; psi: axis CA-C
                  if (which == 0 && run.res[k].name == "PRO") continue;
                  const std::size_t i0 = 3 * k + which, i1 = i0 + 1;
                  const Vec o = x[i1], ax = unit(sub(x[i1], x[i0]));
                  double num = 0.0, den = 0.0;
                  const Vec targets[2] = {tN, tCA};
                  for (int t = 0; t < 2; ++t) {
                    const Vec M = x[m - 2 + t];
                    const Vec rm = sub(M, o), rf = sub(targets[t], o);
                    const Vec r = sub(rm, mul(ax, dot(rm, ax)));
                    const Vec f = sub(rf, mul(ax, dot(rf, ax)));
                    num += dot(f, cross(ax, r));
                    den += dot(f, r);
                  }
                  const double alpha = std::atan2(num, den);
                  if (std::fabs(alpha) < 1e-9) continue;
                  const double c = std::cos(alpha), s = std::sin(alpha);
                  for (std::size_t j = i1 + 1; j < m; ++j) x[j] = rotate_about(x[j], o, ax, c, s);
                }
              }
              err = std::sqrt((dot(sub(x[m - 2], tN), sub(x[m - 2], tN)) + dot(sub(x[m - 1], tCA), sub(x[m - 1], tCA))) / 2.0);
            }
            x.resize(m - 2);
            x.push_back(tN);  // for the last O
            LoopResult cur;
            cur.closure = err;
            cur.bb.resize(n);
            for (std::size_t k = 0; k < n; ++k) {
              const Vec nextN = x[3 * k + 3];
              cur.bb[k] = {x[3 * k], x[3 * k + 1], x[3 * k + 2], place(nextN, x[3 * k + 1], x[3 * k + 2], kCO, kCACO, 180.0)};
            }
            double score = 2000.0 * std::max(0.0, err - 0.02);
            for (std::size_t k = 0; k < n; ++k)
              for (int a = 0; a < 4; ++a) {
                const char el = a == 0 ? 'N' : a == 3 ? 'O' : 'C';
                grid.query(cur.bb[k][a], 4.0, [&](std::uint32_t id, Vec, double d2) {
                  const EnvAtom& e = env[id];
                  if ((e.res == run.before && k == 0) || (e.res == run.after && k + 1 == n)) return;
                  const double v = pair_score(el, e.el, std::sqrt(d2));
                  if (v > 0) score += v;
                });
                for (std::size_t j = 0; j + 1 < k; ++j)
                  for (int b2 = 0; b2 < 4; ++b2) {
                    const double v = pair_score(el, b2 == 0 ? 'N' : b2 == 3 ? 'O' : 'C', norm(sub(cur.bb[k][a], cur.bb[j][b2])));
                    if (v > 0) score += v;
                  }
              }
            cur.score = score;
            if (cur.score < best.score) best = cur;
            if (best.score < 1e-6) break;
            continue;
          }
          // C-terminal tail: no closure.
          LoopResult cur;
          cur.bb.resize(n);
          for (std::size_t k = 0; k < n; ++k) {
            const Vec N = x[3 * k], CA = x[3 * k + 1], C = x[3 * k + 2];
            const Vec nextN = k + 1 < n ? x[3 * k + 3] : place(N, CA, C, kCN, kCACN, psi[k]);
            cur.bb[k] = {N, CA, C, place(nextN, CA, C, kCO, kCACO, 180.0)};
          }
          double score = 0.0;
          for (std::size_t k = 0; k < n; ++k)
            for (int a = 0; a < 4; ++a) {
              const char el = a == 0 ? 'N' : a == 3 ? 'O' : 'C';
              grid.query(cur.bb[k][a], 4.0, [&](std::uint32_t id, Vec, double d2) {
                if (env[id].res == run.before && k == 0) return;
                const double v = pair_score(el, env[id].el, std::sqrt(d2));
                if (v > 0) score += v;
              });
              for (std::size_t j = 0; j + 1 < k; ++j)
                for (int b2 = 0; b2 < 4; ++b2) {
                  const double v = pair_score(el, b2 == 0 ? 'N' : b2 == 3 ? 'O' : 'C', norm(sub(cur.bb[k][a], cur.bb[j][b2])));
                  if (v > 0) score += v;
                }
            }
          cur.score = score;
          if (cur.score < best.score) best = cur;
          if (best.score < 1e-6) break;
        } else {
          // N-terminal tail: built backwards from residue B.
          Vec nextN = pos(*B->atom("N")), nextCA = pos(*B->atom("CA")), nextC = pos(*B->atom("C"));
          LoopResult cur;
          cur.bb.resize(n);
          double phi_next = -90.0 + 60.0 * (rng.uniform() - 0.5);
          for (std::size_t kk = n; kk-- > 0;) {
            const Vec C = place(nextC, nextCA, nextN, kCN, kCNCA, phi_next);
            const Vec CA = place(nextCA, nextN, C, kCAC, kCACN, 180.0);
            const Vec N = place(nextN, C, CA, kNCA, kNCAC, psi[kk]);
            cur.bb[kk] = {N, CA, C, place(nextN, CA, C, kCO, kCACO, 180.0)};
            nextN = N, nextCA = CA, nextC = C, phi_next = phi[kk];
          }
          double score = 0.0;
          for (std::size_t k = 0; k < n; ++k)
            for (int a = 0; a < 4; ++a) {
              const char el = a == 0 ? 'N' : a == 3 ? 'O' : 'C';
              grid.query(cur.bb[k][a], 4.0, [&](std::uint32_t id, Vec, double d2) {
                if (env[id].res == run.after && k + 1 == n) return;
                const double v = pair_score(el, env[id].el, std::sqrt(d2));
                if (v > 0) score += v;
              });
              for (std::size_t j = k + 2; j < n; ++j)
                for (int b2 = 0; b2 < 4; ++b2) {
                  const double v = pair_score(el, b2 == 0 ? 'N' : b2 == 3 ? 'O' : 'C', norm(sub(cur.bb[k][a], cur.bb[j][b2])));
                  if (v > 0) score += v;
                }
            }
          cur.score = score;
          if (cur.score < best.score) best = cur;
          if (best.score < 1e-6) break;
        }
      }
      if (best.bb.empty() || (internal && best.closure > 0.35)) {
        rep.not_built.push_back(range + ": the loop could not be closed");
        continue;
      }
      if (internal && best.closure > 0.1) {
        char b[96];
        std::snprintf(b, sizeof b, ": built, but the closing peptide bond is off by up to %.2f Å", best.closure);
        rep.not_built.push_back(range + b);
      }
      Built b;
      b.insert_at = run.after != SIZE_MAX ? run.after : run.before + 1;
      for (std::size_t k = 0; k < n; ++k) {
        ProteinResidue r;
        r.name = run.res[k].name;
        r.chain = run.res[k].chain;
        r.number = run.res[k].number;
        r.icode = run.res[k].icode;
        r.predicted = true;
        const char* names[4] = {"N", "CA", "C", "O"};
        for (int a = 0; a < 4; ++a) {
          ProteinAtom at;
          at.name = names[a];
          at.element = a == 0 ? "N" : a == 3 ? "O" : "C";
          set_pos(at, best.bb[k][a]);
          at.modelled = true;
          r.atoms.push_back(at);
          add_env(best.bb[k][a], UINT32_MAX, at.element[0]);
          ++rep.atoms_added;
        }
        rep.built_residues.push_back(label(run.res[k]));
        done.insert({r.chain, r.number, r.icode});
        b.res.push_back(std::move(r));
      }
      built.push_back(std::move(b));
    }
    if (!built.empty()) {
      std::stable_sort(built.begin(), built.end(), [](const Built& a, const Built& b) { return a.insert_at < b.insert_at; });
      std::vector<ProteinResidue> out;
      std::vector<std::size_t> remap(p.residues.size());
      std::size_t bi = 0;
      for (std::size_t i = 0; i <= p.residues.size(); ++i) {
        while (bi < built.size() && built[bi].insert_at == i) {
          for (auto& r : built[bi].res) out.push_back(std::move(r));
          ++bi;
        }
        if (i < p.residues.size()) {
          remap[i] = out.size();
          out.push_back(std::move(p.residues[i]));
        }
      }
      p.residues = std::move(out);
      for (auto& [a, b] : p.disulfides) a = remap[a], b = remap[b];
      p.missing_residues.erase(std::remove_if(p.missing_residues.begin(), p.missing_residues.end(),
                                              [&](const MissingResidue& m) { return done.count({m.chain, m.number, m.icode}) > 0; }),
                               p.missing_residues.end());
      recompute_topology(p);
    }
  }

  // ---- Missing side-chain atoms ---------------------------------------------------------------------------
  if (opt.side_chains) {
    std::vector<Job> jobs;
    std::vector<int> job_of(p.residues.size(), -1);
    for (std::size_t i = 0; i < p.residues.size(); ++i) {
      auto& r = p.residues[i];
      const std::string base = canonical_name(r.name);
      auto it = sp.find(base);
      if (it == sp.end() || !defs.count(base)) continue;
      if (!r.atom("N") || !r.atom("CA") || !r.atom("C")) {
        rep.not_built.push_back(residue_label(r) + ": backbone atoms missing");
        continue;
      }
      // Backbone O.
      if (!r.atom("O")) {
        Vec x;
        const ProteinAtom* nn = nullptr;
        for (std::size_t j = i + 1; j < p.residues.size() && j <= i + 1; ++j)
          if (p.residues[j].chain == r.chain && !r.c_terminal) nn = p.residues[j].atom("N");
        if (nn) x = place(pos(*nn), pos(*r.atom("CA")), pos(*r.atom("C")), kCO, kCACO, 180.0);
        else x = place(pos(*r.atom("N")), pos(*r.atom("CA")), pos(*r.atom("C")), kCO, kCACO, -40.0);
        ProteinAtom o;
        o.name = "O", o.element = "O", o.modelled = true;
        set_pos(o, x);
        r.atoms.push_back(o);
        ++rep.atoms_added;
        if (!r.predicted) rep.completed_residues.push_back(residue_label(r));
      }
      const Span& s = it->second;
      if (!s.begin) continue;
      Job j;
      j.res = i;
      j.span = &s;
      for (const SideChainAtom* a = s.begin; a != s.end; ++a)
        if (!r.atom(a->atom)) j.missing.push_back(int(a - s.begin));
      if (j.missing.empty()) continue;
      // Chi angles fixed by atoms that are present.
      std::array<double, 4> have{};
      std::array<bool, 4> known{};
      std::array<bool, 4> needed{};
      for (const SideChainAtom* a = s.begin; a != s.end; ++a) {
        if (a->chi < 0) continue;
        if (r.atom(a->atom) && r.atom(a->a) && r.atom(a->b) && r.atom(a->c)) {
          known[a->chi] = true;
          have[a->chi] = dihedral_deg(pos(*r.atom(a->a)), pos(*r.atom(a->b)), pos(*r.atom(a->c)), pos(*r.atom(a->atom)));
        } else {
          needed[a->chi] = true;
        }
      }
      // Backbone class for the chi1 factor (as in scripts/gen_protein_geometry.py).
      const double* bias = nullptr;
      if (!known[0] && i > 0 && i + 1 < p.residues.size() && !r.n_terminal && !r.c_terminal) {
        const ProteinAtom* pc = p.residues[i - 1].atom("C");
        const ProteinAtom* nn = p.residues[i + 1].atom("N");
        if (pc && nn && dist(*pc, *r.atom("N")) < 1.45 && dist(*r.atom("C"), *nn) < 1.45) {
          const double phi = dihedral_deg(pos(*pc), pos(*r.atom("N")), pos(*r.atom("CA")), pos(*r.atom("C")));
          const double psi = dihedral_deg(pos(*r.atom("N")), pos(*r.atom("CA")), pos(*r.atom("C")), pos(*nn));
          const int cls = phi >= 0 ? 2 : (psi >= -120 && psi < 50 ? 0 : 1);
          for (const auto& b : kChi1Bias)
            if (b.cls == cls && base == b.res) bias = b.f;
        }
      }
      std::vector<std::pair<std::array<double, 4>, double>> cands;
      for (int pass = 0; pass < 2 && cands.empty(); ++pass)
        for (const Rotamer* rot = s.rbegin; rot && rot != s.rend; ++rot) {
          bool ok = true;
          std::array<double, 4> c{};
          for (int k = 0; k < rot->nchi; ++k) {
            if (known[k]) {
              c[k] = have[k];
              if (pass == 0 && angle_diff(rot->chi[k], have[k]) > 45.0) ok = false;
            } else {
              c[k] = rot->chi[k];
            }
          }
          if (!ok) continue;
          double w = rot->freq;
          if (bias) w *= bias[int(std::fmod(c[0] + 360.0, 360.0) / 120.0) % 3];
          bool dup = false;
          for (auto& e : cands) {
            bool same = true;
            for (int k = 0; k < 4; ++k) same = same && (!needed[k] || std::fabs(e.first[k] - c[k]) < 1e-9);
            if (same) e.second += w, dup = true;
          }
          if (!dup) cands.push_back({c, w});
        }
      if (cands.empty()) cands.push_back({std::array<double, 4>{}, 1.0});
      // Sub-rotamers: chi1 and chi2 also 15 degrees to either side of the library value.
      {
        const double SUB = 15.0, SUBW = 0.5;
        const std::size_t n0 = cands.size();
        if (base != "PRO")
          for (std::size_t c = 0; c < n0; ++c)
            for (int d1 = -1; d1 <= 1; ++d1)
              for (int d2 = -1; d2 <= 1; ++d2) {
                if (!d1 && !d2) continue;
                if ((d1 && !needed[0]) || (d2 && !needed[1])) continue;
                auto e = cands[c];
                e.first[0] += SUB * d1, e.first[1] += SUB * d2;
                e.second *= (d1 ? SUBW : 1.0) * (d2 ? SUBW : 1.0);
                cands.push_back(e);
              }
      }
      double total = 0.0;
      for (const auto& c : cands) total += c.second;
      for (const auto& c : cands) {
        std::vector<Vec> xs;
        if (!build_atoms(r, s, j.missing, c.first, xs)) continue;
        j.cand.push_back(c.first);
        j.prior.push_back(-std::log(c.second / total));
        j.coords.push_back(std::move(xs));
      }
      if (j.cand.empty()) {
        rep.not_built.push_back(residue_label(r) + ": side chain cannot be placed");
        continue;
      }
      j.ca = pos(*r.atom("CA"));
      job_of[i] = int(jobs.size());
      jobs.push_back(std::move(j));
    }
    if (!jobs.empty()) {
      std::vector<EnvAtom> env;
      Grid grid(4.0);
      for (std::size_t i = 0; i < p.residues.size(); ++i)
        for (const auto& a : p.residues[i].atoms)
          if (a.element != "H" && a.element != "D") {
            grid.add(pos(a), std::uint32_t(env.size()));
            env.push_back({pos(a), std::uint32_t(i), el_code(a.element), false});
          }
      Grid cagrid(8.0);
      for (std::size_t k = 0; k < jobs.size(); ++k) cagrid.add(jobs[k].ca, std::uint32_t(k));
      for (std::size_t k = 0; k < jobs.size(); ++k)
        cagrid.query(jobs[k].ca, 16.0, [&](std::uint32_t id, Vec, double) {
          if (id != k) jobs[k].near.push_back(id);
        });
      // Score of each candidate against the atoms that are present (computed once).
      auto el_of = [&](const Job& j, std::size_t a) { return el_code(std::string(1, j.span->begin[j.missing[a]].atom[0])); };
      std::vector<std::vector<double>> self(jobs.size());
      for (std::size_t k = 0; k < jobs.size(); ++k) {
        Job& j = jobs[k];
        self[k].resize(j.cand.size());
        for (std::size_t c = 0; c < j.cand.size(); ++c) {
          double e = j.prior[c];
          for (std::size_t a = 0; a < j.missing.size(); ++a) {
            if (std::strcmp(j.span->begin[j.missing[a]].atom, "CB") == 0) continue;  // fixed by the backbone
            const char el = el_of(j, a);
            grid.query(j.coords[c][a], 5.0, [&](std::uint32_t id, Vec, double d2) {
              if (env[id].res == j.res) return;
              e += pair_score(el, env[id].el, std::sqrt(d2));
            });
          }
          self[k][c] = e;
        }
        j.choice = int(std::min_element(self[k].begin(), self[k].end()) - self[k].begin());
      }
      auto pair_energy = [&](const Job& a, int ca, const Job& b, int cb) {
        double e = 0.0;
        for (std::size_t i = 0; i < a.missing.size(); ++i)
          for (std::size_t j = 0; j < b.missing.size(); ++j) {
            const Vec d = sub(a.coords[ca][i], b.coords[cb][j]);
            const double d2 = dot(d, d);
            if (d2 < 25.0) e += pair_score(el_of(a, i), el_of(b, j), std::sqrt(d2));
          }
        return e;
      };
      for (int sweep = 0; sweep < 8; ++sweep) {
        bool changed = false;
        for (std::size_t k = 0; k < jobs.size(); ++k) {
          Job& j = jobs[k];
          if (j.cand.size() < 2 || j.near.empty()) continue;
          int bestc = j.choice;
          double beste = 1e300;
          for (std::size_t c = 0; c < j.cand.size(); ++c) {
            double e = self[k][c];
            for (std::size_t o : j.near) e += pair_energy(j, int(c), jobs[o], jobs[o].choice);
            if (e < beste - 1e-12) beste = e, bestc = int(c);
          }
          if (bestc != j.choice) j.choice = bestc, changed = true;
        }
        if (!changed) break;
      }
      for (auto& j : jobs) {
        auto& r = p.residues[j.res];
        for (std::size_t a = 0; a < j.missing.size(); ++a) {
          const SideChainAtom& t = j.span->begin[j.missing[a]];
          ProteinAtom at;
          at.name = t.atom;
          at.element = std::string(1, t.atom[0]);
          at.modelled = true;
          set_pos(at, j.coords[j.choice][a]);
          r.atoms.push_back(at);
          ++rep.atoms_added;
        }
        if (!r.predicted) rep.completed_residues.push_back(residue_label(r));
      }
    }
  }
  if (!rep.built_residues.empty())
    diag.warn(DiagCode::MalformedRecord, std::to_string(rep.built_residues.size()) +
                                             " residue(s) without coordinates were built from the sequence (" + rep.built_residues.front() +
                                             (rep.built_residues.size() > 1 ? " ...)" : ")") + "; their conformation is a prediction");
  if (!rep.completed_residues.empty())
    diag.warn(DiagCode::MalformedRecord, std::to_string(rep.completed_residues.size()) + " residue(s) had missing atoms rebuilt (e.g. " +
                                             rep.completed_residues.front() + ")");
  return rep;
}

}  // namespace chemlarp
