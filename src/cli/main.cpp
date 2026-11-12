// chemlarp — terminal front end for ChemLarp.
//
//   chemlarp parse <SMILES> [--json]
//   chemlarp canon <SMILES> [--random N] [--kekule]
//   chemlarp kinetics <mechanism> [--thermo therm.dat] [--smiles map.smi] [-m iso-v|iso-p|adiabatic-v|adiabatic-p]
//                     [-T K] [-P bar] [-t s] [--csv file] NAME=moles ...      (mechanism: .ocm, CHEMKIN-II or Cantera YAML)
//   chemlarp embed <SMILES> [--seed N] [--candidates N] [--no-uff] [--sdf file] [--xyz file]
//   chemlarp flash -T K -P bar [--srk] <SMILES> <fraction> ...      (vapour-liquid equilibrium)
//   chemlarp psat <SMILES> -T K [--srk]                              (vapour pressure)
//   chemlarp reaction "<equation>" [-T kelvin] [--json]
//   chemlarp equilibrate [-T kelvin] [-P bar] [--json] <SMILES> <moles> [<SMILES> <moles> ...]

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chemlarp/canonical.hpp"
#include "chemlarp/cip.hpp"
#include "chemlarp/conformer.hpp"
#include "chemlarp/depict.hpp"
#include "chemlarp/uff.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/eos.hpp"
#include "chemlarp/gibbs.hpp"
#include "chemlarp/inchi.hpp"
#include "chemlarp/json.hpp"
#include "chemlarp/fingerprint.hpp"
#include "chemlarp/kinetics.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/reaction.hpp"
#include "chemlarp/smarts.hpp"
#include "chemlarp/smiles.hpp"
#include "chemlarp/thermo.hpp"

using namespace chemlarp;

namespace {

int usage() {
  std::fputs(
      "usage:\n"
      "  chemlarp parse <SMILES> [--json]\n"
      "  chemlarp canon <SMILES> [--random N] [--kekule]\n"
      "  chemlarp embed <SMILES> [--seed N] [--candidates N] [--no-uff] [--sdf file] [--xyz file]\n"
      "  chemlarp flash -T K -P bar [--srk] <SMILES> <fraction> [<SMILES> <fraction> ...]\n"
      "  chemlarp psat <SMILES> -T K [--srk]\n"
      "  chemlarp inchi <SMILES> [--key]\n"
      "  chemlarp read <file.mol|.sdf|.pdb|.cml> [--format mol|pdb|cml] [--keep-h] [--inchi]\n"
      "  chemlarp match <SMILES> <SMARTS> [--chiral]          substructure matches (atom indices)\n"
      "  chemlarp similarity <SMILES> <SMILES> [--radius N] [--bits N]\n"
      "  chemlarp search <file.smi|.sdf|...> [--smarts Q] [--like SMILES] [--min T] [--radius N] [--bits N]\n"
      "  chemlarp depict <SMILES> [-o file.svg] [--mol file.mol] [--cip] [--numbers] [--mono] [--bond px]\n"
      "  chemlarp kinetics <mechanism> [--thermo therm.dat] [--smiles map.smi]\n"
      "                    [-m iso-v|iso-p|adiabatic-v|adiabatic-p] [-T K] [-P bar] [-t s]\n"
      "                    [--rtol x] [--csv file] NAME=moles ...\n"
      "                    mechanism: .ocm, CHEMKIN-II (chem.inp) or Cantera YAML (detected from content)\n"
      "  chemlarp reaction \"<coef> SMILES + ... -> <coef> SMILES + ...\" [-T kelvin] [--json]\n"
      "  chemlarp equilibrate [-T kelvin] [-P bar] [--json] <SMILES> <moles> [<SMILES> <moles> ...]\n",
      stderr);
  return 2;
}

// Prints each diagnostic with a caret under the offending character.
void print_diagnostics(const Diagnostics& d, std::string_view input) {
  for (const auto& x : d.items()) {
    std::fprintf(stderr, "%s[%s %d]: %s\n", to_string(x.severity), to_string(x.code), static_cast<int>(x.code),
                 x.message.c_str());
    if (!input.empty() && x.position != Diagnostic::npos && x.position <= input.size()) {
      std::fprintf(stderr, "    %.*s\n    %*s^\n", static_cast<int>(input.size()), input.data(),
                   static_cast<int>(x.position), "");
    }
  }
}

bool parse_double(const char* s, double& out) {
  char* end = nullptr;
  out = std::strtod(s, &end);
  return end && *end == '\0' && end != s;
}

int cmd_parse(const std::vector<std::string>& args) {
  bool json = false;
  std::string smiles;
  for (const auto& a : args) {
    if (a == "--json") json = true;
    else if (smiles.empty()) smiles = a;
    else return usage();
  }
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  if (json) {
    std::printf("{\"ok\":%s,\"structure\":%s,\"diagnostics\":%s}\n", sp >= 0 ? "true" : "false",
                structure_to_json(st).c_str(), diagnostics_to_json(d).c_str());
    return sp >= 0 ? 0 : 1;
  }
  print_diagnostics(d, smiles);
  if (sp < 0) return 1;
  const Species& s = st.species[static_cast<std::size_t>(sp)];
  std::printf("canonical   %s\nformula     %s\nmolar mass  %.4f g/mol\natoms       %u (heavy+explicit H)\nbonds       %u\n",
              s.canonical_smiles.c_str(), s.formula.c_str(), s.molar_mass, s.atom_count, s.bond_count);
  for (const auto& c : cip_centres(st, static_cast<std::uint32_t>(sp)))
    std::printf("stereo      atom %u: %s\n", c.index, c.label == CipLabel::None ? "undetermined" : to_string(c.label));
  for (const auto& c : cip_double_bonds(st, static_cast<std::uint32_t>(sp)))
    std::printf("stereo      bond %u: %s\n", c.index, c.label == CipLabel::None ? "undetermined" : to_string(c.label));
  std::printf("\n  idx  el  chg  H  arom ring  chiral\n");
  for (std::uint32_t i = s.atom_begin; i < s.atom_begin + s.atom_count; ++i) {
    const Atom a = st.atoms.get(i);
    std::printf("  %3u  %-2.*s  %+3d  %u  %-4s %-4s  %s\n", i, static_cast<int>(element_symbol(a.element).size()),
                element_symbol(a.element).data(), a.formal_charge, a.implicit_h, (a.flags & kAtomAromatic) ? "yes" : "",
                (a.flags & kAtomInRing) ? "yes" : "",
                a.chirality == Chirality::None ? "" : a.chirality == Chirality::Clockwise ? "@@" : "@");
  }
  std::printf("\n  bond  atoms     order\n");
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    std::printf("  %4u  %3u-%-3u   %u%s\n", b, st.bonds.begin[b], st.bonds.end[b], st.bonds.order[b],
                (st.bonds.flags[b] & kBondAromatic) ? " (aromatic)" : "");
  }
  if (auto t = resolve_thermo(st, static_cast<std::uint32_t>(sp), kTref, d)) {
    std::printf("\nthermo (%s) at 298.15 K: H = %.2f kJ/mol, S = %.2f J/mol/K, G = %.2f kJ/mol\n", t->label.c_str(),
                t->h_RT(kTref) * kGasConstant * kTref / 1000, t->s_R(kTref) * kGasConstant,
                t->g_RT(kTref) * kGasConstant * kTref / 1000);
  }
  return 0;
}

int cmd_canon(const std::vector<std::string>& args) {
  std::string smiles;
  unsigned random = 0;
  bool kekule = false;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--kekule") kekule = true;
    else if (args[i] == "--random" && i + 1 < args.size()) random = static_cast<unsigned>(std::strtoul(args[++i].c_str(), nullptr, 10));
    else if (smiles.empty()) smiles = args[i];
    else return usage();
  }
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  print_diagnostics(d, smiles);
  if (sp < 0) return 1;
  const auto s = static_cast<std::uint32_t>(sp);
  SmilesWriteOptions o;
  o.kekule = kekule;
  std::printf("%s\n", write_smiles(st, s, o).c_str());
  for (unsigned i = 0; i < random; ++i) {
    o.canonical = false;
    o.random_seed = i + 1;
    std::printf("  %s\n", write_smiles(st, s, o).c_str());
  }
  return 0;
}

int cmd_kinetics(const std::vector<std::string>& args) {
  std::string path, csv, thermo_path, smiles_path;
  ReactorSetup su;
  su.model = ReactorModel::AdiabaticConstantVolume;
  double P_bar = 1.01325;
  std::vector<std::pair<std::string, double>> feed;
  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string& a = args[i];
    auto next = [&](double& v) { return i + 1 < args.size() && parse_double(args[++i].c_str(), v); };
    if (a == "-m" && i + 1 < args.size()) {
      auto mm = reactor_model_from_string(args[++i]);
      if (!mm) return usage();
      su.model = *mm;
    } else if (a == "-T") { if (!next(su.temperature_K)) return usage(); }
    else if (a == "-P") { if (!next(P_bar)) return usage(); }
    else if (a == "-t") { if (!next(su.t_end)) return usage(); }
    else if (a == "--rtol") { if (!next(su.rtol)) return usage(); }
    else if (a == "--csv" && i + 1 < args.size()) csv = args[++i];
    else if (a == "--thermo" && i + 1 < args.size()) thermo_path = args[++i];
    else if (a == "--smiles" && i + 1 < args.size()) smiles_path = args[++i];
    else if (path.empty()) path = a;
    else {
      const auto eq = a.find('=');
      double v = 0;
      if (eq == std::string::npos || !parse_double(a.c_str() + eq + 1, v)) return usage();
      feed.push_back({a.substr(0, eq), v});
    }
  }
  if (path.empty() || feed.empty()) return usage();
  su.pressure_Pa = P_bar * 1e5;
  Diagnostics d;
  auto mech = load_mechanism(path, d, thermo_path, smiles_path);
  print_diagnostics(d, {});
  if (!mech) return 1;
  d.clear();
  su.moles.assign(mech->species_count(), 0.0);
  for (auto& [name, v] : feed) {
    const std::int32_t sp = mech->find(name);
    if (sp < 0) {
      std::fprintf(stderr, "error: species '%s' is not in the mechanism\n", name.c_str());
      return 1;
    }
    su.moles[static_cast<std::size_t>(sp)] = v;
  }
  for (int k = 1; k <= 20; ++k) su.output_times.push_back(su.t_end * k / 20.0);
  ReactorResult r;
  std::vector<ReactorState> every;
  {
    // Full-resolution run for the CSV / ignition detection, sampled run for the table.
    ReactorSetup all = su;
    all.output_times.clear();
    r = run_reactor(*mech, all, d);
    every = r.samples;
    if (r.ok) {
      Diagnostics d2;
      ReactorResult sampled = run_reactor(*mech, su, d2);
      r.samples = sampled.samples;
    }
  }
  print_diagnostics(d, {});
  if (!r.ok) return 1;
  std::printf("%s: %zu species, %zu reactions, model %s\n", path.c_str(), mech->species_count(),
              mech->reactions.size(), to_string(su.model));
  std::printf("BDF: %zu steps, %zu rejected, %zu Newton failures, %zu f-evals, %zu Jacobians, max order %d\n",
              r.stats.steps, r.stats.rejected, r.stats.newton_failures, r.stats.rhs_evals, r.stats.jacobian_evals,
              r.stats.max_order_used);
  if (r.ignition_time > 0)
    std::printf("ignition (max dT/dt = %.3g K/s) at t = %.6g s\n", r.max_dTdt, r.ignition_time);
  // Show the species that ever exceed 1e-4 mole fraction.
  std::vector<std::size_t> shown;
  for (std::size_t i = 0; i < mech->species_count(); ++i) {
    double mx = 0;
    for (const auto& s : every) {
      double N = 0;
      for (double v : s.n) N += v;
      mx = std::max(mx, s.n[i] / N);
    }
    if (mx > 1e-4) shown.push_back(i);
  }
  std::printf("\n  %-11s %-8s %-9s", "t / s", "T / K", "P / bar");
  for (auto i : shown) std::printf(" %-9s", mech->names[i].c_str());
  std::printf("   (mole fractions)\n");
  for (const auto& s : r.samples) {
    double N = 0;
    for (double v : s.n) N += v;
    std::printf("  %-11.4e %-8.1f %-9.4f", s.t, s.T, s.P / 1e5);
    for (auto i : shown) std::printf(" %-9.3e", s.n[i] / N);
    std::printf("\n");
  }
  if (!csv.empty()) {
    std::ofstream out(csv);
    out << "t_s,T_K,P_Pa";
    for (const auto& nm : mech->names) out << ",n_" << nm;
    out << "\n";
    out.precision(10);
    for (const auto& s : every) {
      out << s.t << ',' << s.T << ',' << s.P;
      for (double v : s.n) out << ',' << v;
      out << "\n";
    }
    std::printf("\nwrote %zu rows to %s\n", every.size(), csv.c_str());
  }
  return 0;
}

int cmd_embed(const std::vector<std::string>& args) {
  std::string smiles, sdf_path, xyz_path;
  EmbedOptions opt;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--no-uff") opt.minimize = false;
    else if (args[i] == "--seed" && i + 1 < args.size()) opt.seed = std::strtoull(args[++i].c_str(), nullptr, 10);
    else if (args[i] == "--candidates" && i + 1 < args.size()) opt.candidates = std::atoi(args[++i].c_str());
    else if (args[i] == "--sdf" && i + 1 < args.size()) sdf_path = args[++i];
    else if (args[i] == "--xyz" && i + 1 < args.size()) xyz_path = args[++i];
    else if (smiles.empty()) smiles = args[i];
    else return usage();
  }
  ChemicalSystemState st;
  Diagnostics d;
  EmbedResult r;
  const std::int32_t sp = embed_smiles(st, smiles, d, &r, opt);
  print_diagnostics(d, smiles);
  if (sp < 0) return 1;
  const auto s = static_cast<std::uint32_t>(sp);
  std::printf("%s: %u atoms (explicit H), %d attempts, distance-geometry residual %.3g\n",
              st.species[s].canonical_smiles.c_str(), st.species[s].atom_count, r.attempts, r.dg_energy);
  if (opt.minimize) {
    Diagnostics q;
    auto ff = UffForceField::build(st, s, q);
    std::printf("UFF: %.3f kcal/mol (%s in %d iterations, max |gradient| %.1e kcal/mol/Å)\n", r.uff.energy_final,
                r.uff.converged ? "converged" : "NOT converged", r.uff.iterations, r.uff.grad_max);
    if (ff) {
      std::vector<double> x;
      for (std::uint32_t i = 0; i < st.species[s].atom_count; ++i) {
        const std::uint32_t g = st.species[s].atom_begin + i;
        x.insert(x.end(), {st.atoms.x[g], st.atoms.y[g], st.atoms.z[g]});
      }
      const auto e = ff->breakdown(x);
      std::printf("     bond %.3f  angle %.3f  torsion %.3f  inversion %.3f  vdW %.3f\n", e.bond, e.angle, e.torsion,
                  e.inversion, e.vdw);
    }
  }
  std::printf("stereo: %zu/%zu centres, %zu/%zu double bonds reproduced in 3D\n", r.stereo_centres_ok, r.stereo_centres,
              r.stereo_bonds_ok, r.stereo_bonds);
  auto save = [](const std::string& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    out << text;
    std::printf("wrote %s\n", path.c_str());
  };
  if (!sdf_path.empty()) save(sdf_path, write_sdf(st, s, smiles));
  if (!xyz_path.empty()) save(xyz_path, write_xyz(st, s, smiles));
  if (sdf_path.empty() && xyz_path.empty()) std::printf("\n%s", write_xyz(st, s, smiles).c_str());
  return 0;
}

int cmd_flash(const std::vector<std::string>& args) {
  double T = 300, P_bar = 1;
  CubicEos eos = CubicEos::PengRobinson;
  std::vector<std::string> smiles;
  std::vector<double> z;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "-T" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), T)) return usage(); }
    else if (args[i] == "-P" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), P_bar)) return usage(); }
    else if (args[i] == "--srk") eos = CubicEos::SoaveRedlichKwong;
    else if (i + 1 < args.size()) {
      double v;
      if (!parse_double(args[i + 1].c_str(), v)) return usage();
      smiles.push_back(args[i]);
      z.push_back(v);
      ++i;
    } else return usage();
  }
  if (smiles.empty()) return usage();
  Diagnostics d;
  auto mix = mixture_from_smiles(smiles, eos, d);
  FlashResult f;
  if (mix) f = flash_tp(*mix, T, P_bar * 1e5, z, d);
  print_diagnostics(d, {});
  if (!f.ok) return 1;
  std::printf("%s flash at %.2f K, %.4g bar: %s", eos == CubicEos::PengRobinson ? "Peng-Robinson" : "SRK", T, P_bar, to_string(f.region));
  if (f.region == PhaseRegion::TwoPhase) std::printf(", vapour fraction %.5f (%d iterations)", f.vapor_fraction, f.iterations);
  std::printf("\n\n  %-18s %-16s %10s %10s %10s %10s\n", "component", "", "feed z", "liquid x", "vapour y", "K = y/x");
  double zs = 0;
  for (double v : z) zs += v;
  for (std::size_t i = 0; i < smiles.size(); ++i)
    std::printf("  %-18s %-16s %10.5f %10.5f %10.5f %10.4g\n", smiles[i].c_str(), mix->components[i].name.c_str(), z[i] / zs, f.x[i], f.y[i], f.K[i]);
  std::printf("\n  Z liquid %.5f   Z vapour %.5f\n", f.liquid.Z, f.vapor.Z);
  return 0;
}

int cmd_inchi(const std::vector<std::string>& args) {
  std::string smiles;
  bool key_only = false;
  for (const auto& a : args) {
    if (a == "--key") key_only = true;
    else if (smiles.empty()) smiles = a;
    else return usage();
  }
  if (smiles.empty()) return usage();
  Diagnostics d;
  const InchiResult r = inchi_from_smiles(smiles, d);
  print_diagnostics(d, smiles);
  if (!r.ok) return 1;
  if (!key_only) std::printf("%s\n", r.inchi.c_str());
  std::printf("%s\n", r.key.c_str());
  return 0;
}

// Prints one line per record: canonical SMILES, a tab and the title (a SMILES file).
int cmd_read(const std::vector<std::string>& args) {
  std::string path;
  ReadOptions opt;
  bool with_inchi = false;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--keep-h") opt.keep_hydrogens = true;
    else if (args[i] == "--inchi") with_inchi = true;
    else if (args[i] == "--format" && i + 1 < args.size()) {
      const std::string f = args[++i];
      if (f == "mol" || f == "sdf") opt.format = MolFormat::Molfile;
      else if (f == "pdb") opt.format = MolFormat::Pdb;
      else if (f == "cml") opt.format = MolFormat::Cml;
      else return usage();
    } else if (path.empty()) path = args[i];
    else return usage();
  }
  if (path.empty()) return usage();
  std::ifstream in(path, std::ios::binary);
  if (!in) { std::fprintf(stderr, "cannot read %s\n", path.c_str()); return 1; }
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  if (opt.format == MolFormat::Auto) opt.format = detect_format(text, path);
  ChemicalSystemState st;
  Diagnostics d;
  const auto records = read_molecules(st, text, d, opt);
  print_diagnostics(d, {});
  int failed = 0;
  for (const MolRecord& r : records) {
    if (r.species < 0) { ++failed; std::printf("*\t%s\n", r.title.c_str()); continue; }
    const Species& s = st.species[static_cast<std::size_t>(r.species)];
    std::printf("%s\t%s", s.canonical_smiles.c_str(), r.title.c_str());
    if (with_inchi) {
      Diagnostics di;
      const InchiResult ir = inchi(st, static_cast<std::uint32_t>(r.species), &di);
      std::printf("\t%s", ir.ok ? ir.key.c_str() : "*");
    }
    std::printf("\n");
  }
  return failed || records.empty() ? 1 : 0;
}

int cmd_match(const std::vector<std::string>& args) {
  std::vector<std::string> pos;
  bool chiral = false;
  for (const auto& a : args) {
    if (a == "--chiral") chiral = true;
    else pos.push_back(a);
  }
  if (pos.size() != 2) return usage();
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, pos[0], d);
  print_diagnostics(d, pos[0]);
  Diagnostics dq;
  const SmartsQuery q = parse_smarts(pos[1], dq);
  print_diagnostics(dq, pos[1]);
  if (sp < 0 || !q.valid()) return 1;
  SmartsMatchOptions o;
  o.use_chirality = chiral;
  o.max_matches = 0;
  const auto matches = match_smarts(st, static_cast<std::uint32_t>(sp), q, o);
  const std::uint32_t base = st.species[static_cast<std::size_t>(sp)].atom_begin;
  std::printf("%zu match%s\n", matches.size(), matches.size() == 1 ? "" : "es");
  for (const auto& m : matches) {
    for (std::size_t k = 0; k < m.size(); ++k) std::printf("%s%u", k ? " " : "", m[k] - base);
    std::printf("\n");
  }
  return 0;
}

bool parse_int_arg(const std::string& s, int& out) {
  char* end = nullptr;
  const long v = std::strtol(s.c_str(), &end, 10);
  if (!end || *end || end == s.c_str() || v <= 0 || v > 1 << 20) return false;
  out = static_cast<int>(v);
  return true;
}

int cmd_similarity(const std::vector<std::string>& args) {
  std::vector<std::string> pos;
  MorganOptions o;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--radius" && i + 1 < args.size()) { if (!parse_int_arg(args[++i], o.radius)) return usage(); }
    else if (args[i] == "--bits" && i + 1 < args.size()) { if (!parse_int_arg(args[++i], o.bits)) return usage(); }
    else pos.push_back(args[i]);
  }
  if (pos.size() != 2) return usage();
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t a = parse_smiles(st, pos[0], d);
  const std::int32_t b = parse_smiles(st, pos[1], d);
  print_diagnostics(d, {});
  if (a < 0 || b < 0) return 1;
  const Fingerprint fa = morgan_fingerprint(st, static_cast<std::uint32_t>(a), o);
  const Fingerprint fb = morgan_fingerprint(st, static_cast<std::uint32_t>(b), o);
  std::printf("tanimoto %.4f\ndice     %.4f\n", tanimoto(fa, fb), dice(fa, fb));
  return 0;
}

// Screens a SMILES file (one per line, optional name after whitespace) or any structure file.
int cmd_search(const std::vector<std::string>& args) {
  std::string path, smarts, like;
  double min_sim = 0.0;
  MorganOptions o;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--smarts" && i + 1 < args.size()) smarts = args[++i];
    else if (args[i] == "--like" && i + 1 < args.size()) like = args[++i];
    else if (args[i] == "--min" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), min_sim)) return usage(); }
    else if (args[i] == "--radius" && i + 1 < args.size()) { if (!parse_int_arg(args[++i], o.radius)) return usage(); }
    else if (args[i] == "--bits" && i + 1 < args.size()) { if (!parse_int_arg(args[++i], o.bits)) return usage(); }
    else if (path.empty()) path = args[i];
    else return usage();
  }
  if (path.empty() || (smarts.empty() && like.empty())) return usage();
  std::ifstream in(path, std::ios::binary);
  if (!in) { std::fprintf(stderr, "cannot read %s\n", path.c_str()); return 1; }
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  Diagnostics d;
  SmartsQuery q;
  if (!smarts.empty()) {
    q = parse_smarts(smarts, d);
    if (!q.valid()) { print_diagnostics(d, smarts); return 1; }
  }
  ChemicalSystemState st;
  Fingerprint ref;
  if (!like.empty()) {
    const std::int32_t r = parse_smiles(st, like, d);
    if (r < 0) { print_diagnostics(d, like); return 1; }
    ref = morgan_fingerprint(st, static_cast<std::uint32_t>(r), o);
  }
  struct Hit { std::string smiles, title; std::size_t count; double sim; };
  std::vector<Hit> hits;
  std::size_t total = 0;
  auto consider = [&](std::uint32_t sp, const std::string& title) {
    ++total;
    std::size_t count = 0;
    if (q.valid()) {
      count = match_smarts(st, sp, q).size();
      if (count == 0) return;
    }
    double sim = 0;
    if (!like.empty()) {
      sim = tanimoto(ref, morgan_fingerprint(st, sp, o));
      if (sim < min_sim) return;
    }
    hits.push_back({st.species[sp].canonical_smiles, title, count, sim});
  };
  const MolFormat fmt = detect_format(text, path);
  if (fmt == MolFormat::Auto) {  // SMILES file
    std::size_t start = 0;
    while (start < text.size()) {
      std::size_t end = text.find('\n', start);
      if (end == std::string::npos) end = text.size();
      std::string line = text.substr(start, end - start);
      start = end + 1;
      while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      const std::size_t cut = line.find_first_of(" \t");
      const std::string smi = line.substr(0, cut);
      const std::string title = cut == std::string::npos ? "" : line.substr(line.find_first_not_of(" \t", cut));
      Diagnostics local;
      const std::int32_t sp = parse_smiles(st, smi, local);
      if (sp >= 0) consider(static_cast<std::uint32_t>(sp), title);
    }
  } else {
    ReadOptions ro;
    ro.format = fmt;
    for (const MolRecord& r : read_molecules(st, text, d, ro))
      if (r.species >= 0) consider(static_cast<std::uint32_t>(r.species), r.title);
  }
  if (!like.empty())
    std::stable_sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.sim > b.sim; });
  for (const Hit& h : hits) {
    std::printf("%s\t%s", h.smiles.c_str(), h.title.c_str());
    if (q.valid()) std::printf("\tmatches=%zu", h.count);
    if (!like.empty()) std::printf("\ttanimoto=%.3f", h.sim);
    std::printf("\n");
  }
  std::fprintf(stderr, "%zu of %zu structures\n", hits.size(), total);
  return 0;
}

int cmd_depict(const std::vector<std::string>& args) {
  std::string smiles, svg_path, mol_path;
  DepictOptions opt;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "-o" && i + 1 < args.size()) svg_path = args[++i];
    else if (args[i] == "--mol" && i + 1 < args.size()) mol_path = args[++i];
    else if (args[i] == "--cip") opt.cip_labels = true;
    else if (args[i] == "--numbers") opt.atom_numbers = true;
    else if (args[i] == "--mono") opt.color = false;
    else if (args[i] == "--bond" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), opt.bond_length)) return usage(); }
    else if (smiles.empty()) smiles = args[i];
    else return usage();
  }
  if (smiles.empty()) return usage();
  Diagnostics d;
  const Depiction dp = depict_smiles(smiles, d, opt);
  print_diagnostics(d, smiles);
  if (!dp.ok) return 1;
  auto save = [](const std::string& path, const std::string& text) {
    std::ofstream f(path, std::ios::binary);
    f << text;
    return static_cast<bool>(f);
  };
  if (svg_path.empty()) std::fputs((dp.svg + "\n").c_str(), stdout);
  else if (!save(svg_path, dp.svg)) { std::fprintf(stderr, "cannot write %s\n", svg_path.c_str()); return 1; }
  if (!mol_path.empty() && !save(mol_path, dp.molfile)) { std::fprintf(stderr, "cannot write %s\n", mol_path.c_str()); return 1; }
  if (!svg_path.empty())
    std::fprintf(stderr, "%s: %.0fx%.0f px, stereo %zu/%zu centres, %zu/%zu double bonds\n", svg_path.c_str(), dp.width, dp.height,
                 dp.stereo_centres_ok, dp.stereo_centres, dp.stereo_bonds_ok, dp.stereo_bonds);
  return 0;
}

int cmd_psat(const std::vector<std::string>& args) {
  double T = 300;
  CubicEos eos = CubicEos::PengRobinson;
  std::string smiles;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "-T" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), T)) return usage(); }
    else if (args[i] == "--srk") eos = CubicEos::SoaveRedlichKwong;
    else if (smiles.empty()) smiles = args[i];
    else return usage();
  }
  Diagnostics d;
  auto mix = mixture_from_smiles({smiles}, eos, d);
  std::optional<double> p;
  if (mix) p = saturation_pressure(mix->components[0], eos, T, d);
  print_diagnostics(d, smiles);
  if (!p) return 1;
  std::printf("%s vapour pressure at %.2f K: %.6g bar (%.6g kPa)\n", mix->components[0].name.c_str(), T, *p / 1e5, *p / 1e3);
  return 0;
}

int cmd_reaction(const std::vector<std::string>& args) {
  bool json = false;
  double T = kTref;
  std::string eq;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--json") json = true;
    else if (args[i] == "-T" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), T)) return usage(); }
    else if (eq.empty()) eq = args[i];
    else return usage();
  }
  ChemicalSystemState st;
  Diagnostics d;
  auto rxn = parse_reaction(st, eq, d);
  bool balanced = rxn && check_balance(st, *rxn, d);
  std::optional<ReactionThermo> th;
  if (balanced) th = reaction_thermo(st, *rxn, T, d);
  if (json) {
    std::printf("{\"ok\":%s,\"balanced\":%s,\"diagnostics\":%s", rxn && !d.has_errors() ? "true" : "false",
                balanced ? "true" : "false", diagnostics_to_json(d).c_str());
    if (th) std::printf(",\"thermo\":{\"T\":%.17g,\"dG_RT\":%.17g,\"dH_J\":%.17g,\"dS_J\":%.17g,\"K\":%.17g}", T, th->dG_RT,
                        th->dH_J, th->dS_J, th->K);
    std::printf("}\n");
    return d.has_errors() ? 1 : 0;
  }
  print_diagnostics(d, eq);
  if (!rxn) return 1;
  std::printf("balanced: %s\n", balanced ? "yes" : "no");
  if (th) {
    std::printf("T = %.2f K   ΔrH° = %.2f kJ/mol   ΔrS° = %.2f J/mol/K   ΔrG° = %.2f kJ/mol   K = %.6g\n", T,
                th->dH_J / 1000, th->dS_J, th->dG_RT * kGasConstant * T / 1000, th->K);
  }
  return d.has_errors() ? 1 : 0;
}

int cmd_equilibrate(const std::vector<std::string>& args) {
  bool json = false;
  double T = 1000.0, P_bar = 1.0;
  std::vector<std::string> rest;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--json") json = true;
    else if (args[i] == "-T" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), T)) return usage(); }
    else if (args[i] == "-P" && i + 1 < args.size()) { if (!parse_double(args[++i].c_str(), P_bar)) return usage(); }
    else rest.push_back(args[i]);
  }
  if (rest.empty() || rest.size() % 2 != 0) return usage();
  ChemicalSystemState st;
  st.conditions = {T, P_bar * 1e5};
  Diagnostics d;
  bool ok = true;
  for (std::size_t i = 0; i < rest.size(); i += 2) {
    double n = 0;
    if (!parse_double(rest[i + 1].c_str(), n)) {
      std::fprintf(stderr, "error: '%s' is not a number of moles\n", rest[i + 1].c_str());
      return 2;
    }
    const std::int32_t sp = parse_smiles(st, rest[i], d);
    if (sp < 0) { ok = false; print_diagnostics(d, rest[i]); d.clear(); continue; }
    st.moles[static_cast<std::size_t>(sp)] = n;
  }
  if (!ok) return 1;
  EquilibriumResult r = equilibrate(st, d);
  if (json) {
    std::printf("{\"ok\":%s,\"equilibrium\":%s,\"diagnostics\":%s}\n", d.has_errors() ? "false" : "true",
                r.n.empty() ? "null" : equilibrium_to_json(st, r).c_str(), diagnostics_to_json(d).c_str());
    return d.has_errors() ? 1 : 0;
  }
  print_diagnostics(d, {});
  if (r.n.empty()) return 1;
  std::printf("T = %.2f K, P = %.4g bar, converged = %s in %d iterations, element error = %.2e\n\n", T, P_bar,
              r.converged ? "yes" : "no", r.iterations, r.max_relative_element_error);
  std::printf("  %-18s %-8s %14s %14s %12s\n", "species", "formula", "n0 / mol", "n_eq / mol", "x_eq");
  for (std::size_t i = 0; i < st.species.size(); ++i) {
    double n0 = 0;
    parse_double(rest[2 * i + 1].c_str(), n0);
    std::printf("  %-18s %-8s %14.6g %14.6g %12.6g\n", st.species[i].smiles.c_str(), st.species[i].formula.c_str(), n0,
                r.n[i], r.mole_fraction[i]);
  }
  std::printf("\n  G/RT = %.6f   element potentials (λ/RT):", r.G_RT);
  for (std::size_t j = 0; j < r.elements.size(); ++j)
    std::printf(" %.*s=%.4f", static_cast<int>(element_symbol(r.elements[j]).size()), element_symbol(r.elements[j]).data(),
                r.element_potentials[j]);
  std::printf("\n");
  return r.converged ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
  SetConsoleOutputCP(CP_UTF8);  // output contains UTF-8 (Δ, °, λ)
#endif
  if (argc < 2) return usage();
  const std::string cmd = argv[1];
  std::vector<std::string> args(argv + 2, argv + argc);
  if (cmd == "parse" && !args.empty()) return cmd_parse(args);
  if (cmd == "reaction" && !args.empty()) return cmd_reaction(args);
  if (cmd == "canon" && !args.empty()) return cmd_canon(args);
  if (cmd == "embed" && !args.empty()) return cmd_embed(args);
  if (cmd == "flash" && !args.empty()) return cmd_flash(args);
  if (cmd == "psat" && !args.empty()) return cmd_psat(args);
  if (cmd == "depict" && !args.empty()) return cmd_depict(args);
  if (cmd == "inchi" && !args.empty()) return cmd_inchi(args);
  if (cmd == "read" && !args.empty()) return cmd_read(args);
  if (cmd == "match" && !args.empty()) return cmd_match(args);
  if (cmd == "similarity" && !args.empty()) return cmd_similarity(args);
  if (cmd == "search" && !args.empty()) return cmd_search(args);
  if (cmd == "kinetics" && !args.empty()) return cmd_kinetics(args);
  if (cmd == "equilibrate") return cmd_equilibrate(args);
  return usage();
}
