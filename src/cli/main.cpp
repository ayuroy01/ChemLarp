// omnichem — terminal front end for OmniChem Core.
//
//   omnichem parse <SMILES> [--json]
//   omnichem canon <SMILES> [--random N] [--kekule]
//   omnichem kinetics <mechanism.ocm> [-m iso-v|iso-p|adiabatic-v|adiabatic-p] [-T K] [-P bar] [-t s] [--csv file] NAME=moles ...
//   omnichem reaction "<equation>" [-T kelvin] [--json]
//   omnichem equilibrate [-T kelvin] [-P bar] [--json] <SMILES> <moles> [<SMILES> <moles> ...]

#include <cstdio>
#include <cstdlib>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include "omnichem/canonical.hpp"
#include "omnichem/cip.hpp"
#include "omnichem/elements.hpp"
#include "omnichem/gibbs.hpp"
#include "omnichem/json.hpp"
#include "omnichem/kinetics.hpp"
#include "omnichem/reaction.hpp"
#include "omnichem/smiles.hpp"
#include "omnichem/thermo.hpp"

using namespace omnichem;

namespace {

int usage() {
  std::fputs(
      "usage:\n"
      "  omnichem parse <SMILES> [--json]\n"
      "  omnichem canon <SMILES> [--random N] [--kekule]\n"
      "  omnichem kinetics <mechanism.ocm> [-m iso-v|iso-p|adiabatic-v|adiabatic-p] [-T K] [-P bar] [-t s]\n"
      "                    [--rtol x] [--csv file] NAME=moles ...\n"
      "  omnichem reaction \"<coef> SMILES + ... -> <coef> SMILES + ...\" [-T kelvin] [--json]\n"
      "  omnichem equilibrate [-T kelvin] [-P bar] [--json] <SMILES> <moles> [<SMILES> <moles> ...]\n",
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
  std::string path, csv;
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
  auto mech = load_mechanism(path, d);
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
  if (cmd == "kinetics" && !args.empty()) return cmd_kinetics(args);
  if (cmd == "equilibrate") return cmd_equilibrate(args);
  return usage();
}
