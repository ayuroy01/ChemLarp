// This translation unit defines the exported C ABI: export unless the static library is being built.
#if !defined(CHEMLARP_STATIC) && !defined(CHEMLARP_BUILD_DLL)
#define CHEMLARP_BUILD_DLL
#endif
#include "chemlarp/chemlarp_c.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <optional>
#include <string>
#include <vector>

#include "chemlarp/canonical.hpp"
#include "chemlarp/cip.hpp"
#include "chemlarp/conformer.hpp"
#include "chemlarp/depict.hpp"
#include "chemlarp/diagnostics.hpp"
#include "chemlarp/eos.hpp"
#include "chemlarp/gibbs.hpp"
#include "chemlarp/inchi.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/json.hpp"
#include "chemlarp/fingerprint.hpp"
#include "chemlarp/kinetics.hpp"
#include "chemlarp/reaction.hpp"
#include "chemlarp/smarts.hpp"
#include "chemlarp/smiles.hpp"
#include "chemlarp/system_state.hpp"

using namespace chemlarp;

struct cl_system {
  ChemicalSystemState state;
  Diagnostics last;
};

namespace {

char* dup(const std::string& s) {
  char* p = static_cast<char*>(std::malloc(s.size() + 1));
  if (!p) return nullptr;
  std::memcpy(p, s.c_str(), s.size() + 1);
  return p;
}

template <class F>
char* guarded(F&& f) noexcept {
  try {
    return dup(f());
  } catch (const std::exception& e) {
    return dup(std::string("{\"ok\":false,\"diagnostics\":[{\"code\":\"InternalError\",\"message\":") +
               json_escape(e.what()) + "}]}");
  } catch (...) {
    return nullptr;
  }
}

}  // namespace

namespace {
// Splits "a b" / "a=b" lines into (name, number) pairs; returns false on a malformed line.
bool parse_pairs(const char* text, std::vector<std::string>& names, std::vector<double>& values, std::string& bad) {
  std::string all = text ? text : "";
  std::size_t start = 0;
  while (start <= all.size()) {
    std::size_t e = all.find('\n', start);
    if (e == std::string::npos) e = all.size();
    std::string line = all.substr(start, e - start);
    start = e + 1;
    for (char& c : line)
      if (c == '\r' || c == '\t') c = ' ';
    const std::size_t a = line.find_first_not_of(' ');
    if (a == std::string::npos) { if (e == all.size()) break; continue; }
    line = line.substr(a, line.find_last_not_of(' ') - a + 1);
    std::size_t sep = line.find_last_of(" =");
    if (sep == std::string::npos) { bad = line; return false; }
    std::string name = line.substr(0, sep), num = line.substr(sep + 1);
    while (!name.empty() && (name.back() == ' ' || name.back() == '=')) name.pop_back();
    char* endp = nullptr;
    const double v = std::strtod(num.c_str(), &endp);
    if (name.empty() || num.empty() || *endp != '\0') { bad = line; return false; }
    names.push_back(name);
    values.push_back(v);
    if (e == all.size()) break;
  }
  return true;
}

char* bad_input(const std::string& line) {
  return dup(std::string("{\"ok\":false,\"diagnostics\":[{\"code\":\"InvalidInput\",\"severity\":\"error\",\"message\":") +
             json_escape("expected 'name amount' but got: " + line) + "}]}");
}
}  // namespace

namespace {
CubicEos eos_from(const char* e) { return (e && std::string(e) == "srk") ? CubicEos::SoaveRedlichKwong : CubicEos::PengRobinson; }
std::string jnum(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.12g", v);
  return std::isfinite(v) ? b : "null";
}
}  // namespace

extern "C" {

const char* cl_version(void) { return "0.1.0"; }

cl_system* cl_system_new(void) { return new (std::nothrow) cl_system(); }
void cl_system_free(cl_system* sys) { delete sys; }

int32_t cl_system_add_smiles(cl_system* sys, const char* smiles) {
  if (!sys || !smiles) return -1;
  sys->last.clear();
  try {
    return parse_smiles(sys->state, smiles, sys->last);
  } catch (...) {
    return -1;
  }
}

int cl_system_set_moles(cl_system* sys, uint32_t species, double moles) {
  if (!sys || species >= sys->state.species.size() || !(moles >= 0.0) || !std::isfinite(moles)) return 0;
  sys->state.moles[species] = moles;
  ++sys->state.revision;
  return 1;
}

int cl_system_set_conditions(cl_system* sys, double temperature_K, double pressure_Pa) {
  if (!sys || !(temperature_K > 0.0) || !(pressure_Pa > 0.0)) return 0;
  sys->state.conditions = {temperature_K, pressure_Pa};
  ++sys->state.revision;
  return 1;
}

int cl_system_equilibrate(cl_system* sys) {
  if (!sys) return 0;
  sys->last.clear();
  try {
    return equilibrate(sys->state, sys->last).converged ? 1 : 0;
  } catch (...) {
    return 0;
  }
}

size_t cl_system_atom_count(const cl_system* s) { return s ? s->state.atoms.size() : 0; }
size_t cl_system_bond_count(const cl_system* s) { return s ? s->state.bonds.size() : 0; }
size_t cl_system_species_count(const cl_system* s) { return s ? s->state.species.size() : 0; }
const uint8_t* cl_system_atom_element(const cl_system* s) { return s ? s->state.atoms.element.data() : nullptr; }
const int8_t* cl_system_atom_charge(const cl_system* s) { return s ? s->state.atoms.formal_charge.data() : nullptr; }
const uint8_t* cl_system_atom_implicit_h(const cl_system* s) { return s ? s->state.atoms.implicit_h.data() : nullptr; }
double* cl_system_atom_x(cl_system* s) { return s ? s->state.atoms.x.data() : nullptr; }
double* cl_system_atom_y(cl_system* s) { return s ? s->state.atoms.y.data() : nullptr; }
double* cl_system_atom_z(cl_system* s) { return s ? s->state.atoms.z.data() : nullptr; }
const uint32_t* cl_system_bond_begin(const cl_system* s) { return s ? s->state.bonds.begin.data() : nullptr; }
const uint32_t* cl_system_bond_end(const cl_system* s) { return s ? s->state.bonds.end.data() : nullptr; }
const uint8_t* cl_system_bond_order(const cl_system* s) { return s ? s->state.bonds.order.data() : nullptr; }
const double* cl_system_moles(const cl_system* s) { return s ? s->state.moles.data() : nullptr; }

char* cl_system_structure_json(const cl_system* s) {
  if (!s) return nullptr;
  return guarded([&] { return structure_to_json(s->state); });
}
char* cl_system_thermo_json(const cl_system* s) {
  if (!s) return nullptr;
  return guarded([&] { return thermo_to_json(s->state); });
}
char* cl_system_diagnostics_json(const cl_system* s) {
  if (!s) return nullptr;
  return guarded([&] { return diagnostics_to_json(s->last); });
}

char* cl_parse_smiles_json(const char* smiles) {
  if (!smiles) return nullptr;
  return guarded([&] {
    ChemicalSystemState st;
    Diagnostics d;
    const bool ok = parse_smiles(st, smiles, d) >= 0;
    return std::string("{\"ok\":") + (ok ? "true" : "false") + ",\"structure\":" + structure_to_json(st) +
           ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_check_reaction_json(const char* equation, double temperature_K) {
  if (!equation) return nullptr;
  return guarded([&] {
    ChemicalSystemState st;
    Diagnostics d;
    auto rxn = parse_reaction(st, equation, d);
    bool balanced = false;
    std::string thermo = "null";
    if (rxn) {
      balanced = check_balance(st, *rxn, d);
      if (balanced && temperature_K > 0.0) {
        if (auto t = reaction_thermo(st, *rxn, temperature_K, d)) {
          char buf[256];
          std::snprintf(buf, sizeof buf, "{\"temperature_K\":%.17g,\"dG_RT\":%.17g,\"dH_J\":%.17g,\"dS_J\":%.17g,\"K\":%.17g}",
                        temperature_K, t->dG_RT, t->dH_J, t->dS_J, t->K);
          thermo = buf;
        }
      }
    }
    return std::string("{\"ok\":") + (rxn && !d.has_errors() ? "true" : "false") +
           ",\"balanced\":" + (balanced ? "true" : "false") + ",\"thermo\":" + thermo +
           ",\"structure\":" + structure_to_json(st) + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_equilibrate_json(const char* const* smiles, const double* moles, size_t count, double temperature_K,
                          double pressure_Pa) {
  if ((!smiles || !moles) && count > 0) return nullptr;
  return guarded([&] {
    ChemicalSystemState st;
    Diagnostics d;
    st.conditions = {temperature_K, pressure_Pa};
    bool parsed = true;
    for (size_t i = 0; i < count; ++i) {
      const std::int32_t sp = parse_smiles(st, smiles[i] ? smiles[i] : "", d);
      if (sp < 0) { parsed = false; continue; }
      st.moles[static_cast<std::size_t>(sp)] = moles[i];
    }
    std::string eq = "null";
    if (parsed) {
      EquilibriumResult r = equilibrate(st, d);
      if (!r.n.empty()) eq = equilibrium_to_json(st, r);
    }
    return std::string("{\"ok\":") + (!d.has_errors() ? "true" : "false") + ",\"equilibrium\":" + eq +
           ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_canonical_smiles_json(const char* smiles, unsigned n_random, uint64_t seed) {
  if (!smiles) return nullptr;
  return guarded([&] {
    ChemicalSystemState st;
    Diagnostics d;
    const std::int32_t sp = parse_smiles(st, smiles, d);
    std::string out = std::string("{\"ok\":") + (sp >= 0 ? "true" : "false");
    if (sp >= 0) {
      const auto s = static_cast<std::uint32_t>(sp);
      out += ",\"canonical\":" + json_escape(st.species[s].canonical_smiles) + ",\"random\":[";
      for (unsigned i = 0; i < n_random; ++i) {
        SmilesWriteOptions o;
        o.canonical = false;
        o.random_seed = seed + i;
        out += (i ? "," : "") + json_escape(write_smiles(st, s, o));
      }
      out += "],\"cip\":{\"centres\":[";
      bool first = true;
      for (const auto& c : cip_centres(st, s)) {
        out += std::string(first ? "" : ",") + "{\"atom\":" + std::to_string(c.index) + ",\"label\":" + json_escape(to_string(c.label)) + "}";
        first = false;
      }
      out += "],\"double_bonds\":[";
      first = true;
      for (const auto& c : cip_double_bonds(st, s)) {
        out += std::string(first ? "" : ",") + "{\"bond\":" + std::to_string(c.index) + ",\"label\":" + json_escape(to_string(c.label)) + "}";
        first = false;
      }
      out += "]}";
    }
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_kinetics_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                       const char* model, double temperature_K, double pressure_Pa,
                       const char* const* names, const double* moles, size_t count, double t_end, size_t n_samples,
                       double rtol) {
  if (!mechanism_text || !model || ((!names || !moles) && count > 0)) return nullptr;
  return guarded([&] {
    Diagnostics d;
    std::string out = "{\"ok\":";
    auto mech = parse_any_mechanism(mechanism_text, thermo_text ? thermo_text : "", smiles_map_text ? smiles_map_text : "",
                                    MechanismFormat::Auto, d);
    auto mm = reactor_model_from_string(model);
    if (!mm) d.error(DiagCode::InvalidReactorSetup, std::string("unknown reactor model '") + model + "'");
    ReactorResult r;
    if (mech && mm) {
      ReactorSetup su;
      su.model = *mm;
      su.temperature_K = temperature_K;
      su.pressure_Pa = pressure_Pa;
      su.t_end = t_end;
      if (rtol > 0) su.rtol = rtol;
      su.moles.assign(mech->species_count(), 0.0);
      bool names_ok = true;
      for (size_t i = 0; i < count; ++i) {
        const std::int32_t sp = mech->find(names[i] ? names[i] : "");
        if (sp < 0) {
          d.error(DiagCode::UnknownSpecies, std::string("species '") + (names[i] ? names[i] : "") + "' is not in the mechanism");
          names_ok = false;
        } else {
          su.moles[static_cast<std::size_t>(sp)] += moles[i];
        }
      }
      for (size_t k = 1; k <= n_samples; ++k) su.output_times.push_back(t_end * static_cast<double>(k) / static_cast<double>(n_samples));
      if (names_ok) r = run_reactor(*mech, su, d);
    }
    out += (r.ok && !d.has_errors()) ? "true" : "false";
    if (mech) {
      out += ",\"species\":[";
      for (std::size_t i = 0; i < mech->names.size(); ++i) out += (i ? "," : "") + json_escape(mech->names[i]);
      out += "]";
    }
    auto num = [](double v) {
      char b[40];
      std::snprintf(b, sizeof b, "%.17g", v);
      return std::string(std::isfinite(v) ? b : "null");
    };
    out += ",\"samples\":[";
    for (std::size_t s = 0; s < r.samples.size(); ++s) {
      const auto& x = r.samples[s];
      out += (s ? "," : "") + std::string("{\"t\":") + num(x.t) + ",\"T\":" + num(x.T) + ",\"P\":" + num(x.P) + ",\"n\":[";
      for (std::size_t i = 0; i < x.n.size(); ++i) out += (i ? "," : "") + num(x.n[i]);
      out += "]}";
    }
    out += "],\"ignition_time\":" + (r.ignition_time > 0 ? num(r.ignition_time) : std::string("null"));
    out += ",\"stats\":{\"steps\":" + std::to_string(r.stats.steps) + ",\"rejected\":" + std::to_string(r.stats.rejected) +
           ",\"newton_failures\":" + std::to_string(r.stats.newton_failures) + ",\"rhs_evals\":" +
           std::to_string(r.stats.rhs_evals) + ",\"max_order\":" + std::to_string(r.stats.max_order_used) + "}";
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_embed_json(const char* smiles, uint64_t seed, int candidates, int minimize) {
  if (!smiles) return nullptr;
  return guarded([&] {
    ChemicalSystemState st;
    Diagnostics d;
    EmbedOptions opt;
    opt.seed = seed;
    if (candidates > 0) opt.candidates = candidates;
    opt.minimize = minimize != 0;
    EmbedResult r;
    const std::int32_t sp = embed_smiles(st, smiles, d, &r, opt);
    std::string out = std::string("{\"ok\":") + (sp >= 0 ? "true" : "false");
    if (sp >= 0) {
      const auto s = static_cast<std::uint32_t>(sp);
      char buf[64];
      std::snprintf(buf, sizeof buf, "%.10g", r.uff.energy_final);
      out += ",\"canonical\":" + json_escape(st.species[s].canonical_smiles) + ",\"sdf\":" + json_escape(write_sdf(st, s, smiles)) +
             ",\"xyz\":" + json_escape(write_xyz(st, s, smiles)) + ",\"energy_kcal_mol\":" + (opt.minimize ? std::string(buf) : "null") +
             ",\"converged\":" + (r.uff.converged ? "true" : "false") + ",\"attempts\":" + std::to_string(r.attempts) +
             ",\"stereo\":{\"centres\":" + std::to_string(r.stereo_centres) + ",\"centres_ok\":" + std::to_string(r.stereo_centres_ok) +
             ",\"bonds\":" + std::to_string(r.stereo_bonds) + ",\"bonds_ok\":" + std::to_string(r.stereo_bonds_ok) + "}";
    }
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}


char* cl_depict_json(const char* smiles, double bond_px, int flags) {
  if (!smiles) return nullptr;
  return guarded([&] {
    Diagnostics d;
    DepictOptions opt;
    if (bond_px > 0) opt.bond_length = bond_px;
    opt.color = (flags & CL_DEPICT_COLOR) != 0;
    opt.cip_labels = (flags & CL_DEPICT_CIP) != 0;
    opt.atom_numbers = (flags & CL_DEPICT_NUMBERS) != 0;
    ChemicalSystemState st;
    const std::int32_t sp = parse_smiles(st, smiles, d);
    Depiction dp;
    if (sp >= 0 && !d.has_errors()) dp = depict(st, static_cast<std::uint32_t>(sp), d, opt);
    std::string out = std::string("{\"ok\":") + (dp.ok ? "true" : "false");
    if (dp.ok) {
      const std::uint32_t begin = st.species[static_cast<std::uint32_t>(sp)].atom_begin;
      out += ",\"svg\":" + json_escape(dp.svg) + ",\"molfile\":" + json_escape(dp.molfile) + ",\"width\":" + jnum(dp.width) +
             ",\"height\":" + jnum(dp.height) + ",\"atoms\":[";
      for (std::size_t i = 0; i < dp.atoms.size(); ++i) {
        const auto& a = dp.atoms[i];
        out += std::string(i ? "," : "") + "{\"atom\":" + std::to_string(a.atom < 0 ? -1 : a.atom - static_cast<std::int32_t>(begin)) +
               ",\"x\":" + jnum(a.x) + ",\"y\":" + jnum(a.y) + ",\"label\":" + json_escape(a.label) + "}";
      }
      out += "],\"bonds\":[";
      for (std::size_t i = 0; i < dp.bonds.size(); ++i) {
        const auto& b = dp.bonds[i];
        out += std::string(i ? "," : "") + "{\"a\":" + std::to_string(b.a) + ",\"b\":" + std::to_string(b.b) +
               ",\"order\":" + std::to_string(b.order) + ",\"stereo\":" + std::to_string(static_cast<int>(b.stereo)) + "}";
      }
      out += "],\"stereo\":{\"centres\":" + std::to_string(dp.stereo_centres) + ",\"centres_ok\":" + std::to_string(dp.stereo_centres_ok) +
             ",\"bonds\":" + std::to_string(dp.stereo_bonds) + ",\"bonds_ok\":" + std::to_string(dp.stereo_bonds_ok) + "}";
    }
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_inchi_json(const char* smiles) {
  if (!smiles) return nullptr;
  return guarded([&] {
    Diagnostics d;
    const InchiResult r = inchi_from_smiles(smiles, d);
    std::string out = std::string("{\"ok\":") + (r.ok ? "true" : "false");
    if (r.ok) out += ",\"inchi\":" + json_escape(r.inchi) + ",\"key\":" + json_escape(r.key);
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_smarts_match_json(const char* smiles, const char* smarts, int use_chirality, int max_matches) {
  if (!smiles || !smarts) return nullptr;
  return guarded([&] {
    Diagnostics d;
    ChemicalSystemState st;
    const std::int32_t sp = parse_smiles(st, smiles, d);
    const SmartsQuery q = parse_smarts(smarts, d);
    if (sp < 0 || !q.valid()) return std::string("{\"ok\":false,\"diagnostics\":") + diagnostics_to_json(d) + "}";
    SmartsMatchOptions o;
    o.use_chirality = use_chirality != 0;
    o.max_matches = max_matches > 0 ? static_cast<std::size_t>(max_matches) : 0;
    const auto matches = match_smarts(st, static_cast<std::uint32_t>(sp), q, o);
    const std::uint32_t base = st.species[static_cast<std::size_t>(sp)].atom_begin;
    std::string m = "[";
    for (std::size_t k = 0; k < matches.size(); ++k) {
      m += k ? ",[" : "[";
      for (std::size_t a = 0; a < matches[k].size(); ++a) m += (a ? "," : "") + std::to_string(matches[k][a] - base);
      m += "]";
    }
    m += "]";
    return "{\"ok\":true,\"count\":" + std::to_string(matches.size()) + ",\"matches\":" + m +
           ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_similarity_json(const char* smiles_a, const char* smiles_b, int radius, int bits) {
  if (!smiles_a || !smiles_b) return nullptr;
  return guarded([&] {
    Diagnostics d;
    ChemicalSystemState st;
    const std::int32_t a = parse_smiles(st, smiles_a, d);
    const std::int32_t b = parse_smiles(st, smiles_b, d);
    if (a < 0 || b < 0) return std::string("{\"ok\":false,\"diagnostics\":") + diagnostics_to_json(d) + "}";
    MorganOptions o;
    if (radius > 0) o.radius = radius;
    if (bits > 0) o.bits = bits;
    const Fingerprint fa = morgan_fingerprint(st, static_cast<std::uint32_t>(a), o);
    const Fingerprint fb = morgan_fingerprint(st, static_cast<std::uint32_t>(b), o);
    return "{\"ok\":true,\"tanimoto\":" + jnum(tanimoto(fa, fb)) + ",\"dice\":" + jnum(dice(fa, fb)) +
           ",\"bits_a\":" + std::to_string(fa.count()) + ",\"bits_b\":" + std::to_string(fb.count()) +
           ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_search_json(const char* smiles_lines, const char* smarts, const char* query_smiles, int radius, int bits) {
  if (!smiles_lines) return nullptr;
  return guarded([&] {
    Diagnostics d;
    const std::string sm = smarts ? smarts : "";
    const std::string qs = query_smiles ? query_smiles : "";
    SmartsQuery q;
    if (!sm.empty()) {
      q = parse_smarts(sm, d);
      if (!q.valid()) return std::string("{\"ok\":false,\"diagnostics\":") + diagnostics_to_json(d) + "}";
    }
    MorganOptions o;
    if (radius > 0) o.radius = radius;
    if (bits > 0) o.bits = bits;
    Fingerprint ref;
    if (!qs.empty()) {
      ChemicalSystemState qst;
      const std::int32_t qsp = parse_smiles(qst, qs, d);
      if (qsp < 0) return std::string("{\"ok\":false,\"diagnostics\":") + diagnostics_to_json(d) + "}";
      ref = morgan_fingerprint(qst, static_cast<std::uint32_t>(qsp), o);
    }
    std::string res = "[";
    const std::string all = smiles_lines;
    std::size_t start = 0, line = 0;
    while (start <= all.size()) {
      std::size_t end = all.find('\n', start);
      if (end == std::string::npos) end = all.size();
      std::string smi = all.substr(start, end - start);
      while (!smi.empty() && (smi.back() == '\r' || smi.back() == ' ')) smi.pop_back();
      ++line;
      if (!smi.empty()) {
        ChemicalSystemState st;
        Diagnostics local;
        const std::int32_t sp = parse_smiles(st, smi, local);
        if (res.size() > 1) res += ',';
        res += "{\"line\":" + std::to_string(line);
        if (sp < 0) {
          res += ",\"ok\":false}";
        } else {
          res += ",\"ok\":true";
          if (q.valid()) res += ",\"count\":" + std::to_string(match_smarts(st, static_cast<std::uint32_t>(sp), q).size());
          if (!qs.empty()) res += ",\"similarity\":" + jnum(tanimoto(ref, morgan_fingerprint(st, static_cast<std::uint32_t>(sp), o)));
          res += "}";
        }
      }
      if (end == all.size()) break;
      start = end + 1;
    }
    res += "]";
    return "{\"ok\":true,\"results\":" + res + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_read_molecules_json(const char* text, const char* format, int keep_hydrogens, int max_records) {
  if (!text) return nullptr;
  return guarded([&] {
    const std::string f = format ? format : "";
    ReadOptions opt;
    opt.keep_hydrogens = keep_hydrogens != 0;
    opt.max_records = max_records > 0 ? static_cast<std::size_t>(max_records) : 0;
    if (f == "mol" || f == "sdf") opt.format = MolFormat::Molfile;
    else if (f == "pdb") opt.format = MolFormat::Pdb;
    else if (f == "cml") opt.format = MolFormat::Cml;
    else opt.format = detect_format(text, f);
    ChemicalSystemState st;
    Diagnostics d;
    const auto records = read_molecules(st, text, d, opt);
    bool any = false;
    std::string recs = "[";
    for (const MolRecord& r : records) {
      if (recs.size() > 1) recs += ',';
      recs += std::string("{\"ok\":") + (r.species >= 0 ? "true" : "false") + ",\"line\":" + std::to_string(r.line) +
              ",\"title\":" + json_escape(r.title) + ",\"dimensions\":" + std::to_string(r.dimensions);
      if (r.species >= 0) {
        any = true;
        const Species& s = st.species[static_cast<std::size_t>(r.species)];
        recs += ",\"smiles\":" + json_escape(s.canonical_smiles) + ",\"formula\":" + json_escape(s.formula) +
                ",\"molar_mass\":" + jnum(s.molar_mass) + ",\"atoms\":" + std::to_string(s.atom_count);
      }
      recs += ",\"properties\":[";
      for (std::size_t k = 0; k < r.properties.size(); ++k)
        recs += (k ? "," : "") + std::string("[") + json_escape(r.properties[k].first) + "," + json_escape(r.properties[k].second) + "]";
      recs += "]}";
    }
    recs += "]";
    return std::string("{\"ok\":") + (any ? "true" : "false") + ",\"format\":" + json_escape(to_string(opt.format)) +
           ",\"records\":" + recs + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_canonical_smiles32_json(const char* smiles, unsigned n_random, unsigned seed) {
  return cl_canonical_smiles_json(smiles, n_random, seed);
}

char* cl_embed32_json(const char* smiles, unsigned seed, int candidates, int minimize) {
  return cl_embed_json(smiles, seed, candidates, minimize);
}

char* cl_equilibrate_text_json(const char* species_lines, double temperature_K, double pressure_Pa) {
  std::vector<std::string> names;
  std::vector<double> moles;
  std::string bad;
  if (!parse_pairs(species_lines, names, moles, bad)) return bad_input(bad);
  std::vector<const char*> ptrs;
  for (const auto& n : names) ptrs.push_back(n.c_str());
  return cl_equilibrate_json(ptrs.data(), moles.data(), ptrs.size(), temperature_K, pressure_Pa);
}

char* cl_kinetics_text_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                            const char* model, double temperature_K, double pressure_Pa, const char* feed_lines,
                            double t_end, size_t n_samples, double rtol) {
  std::vector<std::string> names;
  std::vector<double> moles;
  std::string bad;
  if (!parse_pairs(feed_lines, names, moles, bad)) return bad_input(bad);
  std::vector<const char*> ptrs;
  for (const auto& n : names) ptrs.push_back(n.c_str());
  return cl_kinetics_json(mechanism_text, thermo_text, smiles_map_text, model, temperature_K, pressure_Pa, ptrs.data(),
                          moles.data(), ptrs.size(), t_end, n_samples, rtol);
}


char* cl_flash_text_json(const char* component_lines, double temperature_K, double pressure_Pa, const char* eos) {
  std::vector<std::string> smiles;
  std::vector<double> z;
  std::string bad;
  if (!parse_pairs(component_lines, smiles, z, bad)) return bad_input(bad);
  return guarded([&] {
    Diagnostics d;
    std::string out = "{\"ok\":";
    auto mix = mixture_from_smiles(smiles, eos_from(eos), d);
    FlashResult f;
    if (mix) f = flash_tp(*mix, temperature_K, pressure_Pa, z, d);
    out += (f.ok && !d.has_errors()) ? "true" : "false";
    if (f.ok) {
      double zs = 0;
      for (double v : z) zs += v;
      out += ",\"region\":" + json_escape(to_string(f.region)) + ",\"vapor_fraction\":" + jnum(f.vapor_fraction) +
             ",\"converged\":" + (f.converged ? "true" : "false") + ",\"iterations\":" + std::to_string(f.iterations) +
             ",\"components\":[";
      for (std::size_t i = 0; i < smiles.size(); ++i) {
        out += (i ? "," : "") + std::string("{\"smiles\":") + json_escape(smiles[i]) + ",\"name\":" +
               json_escape(mix->components[i].name) + ",\"z\":" + jnum(z[i] / zs) + ",\"x\":" + jnum(f.x[i]) +
               ",\"y\":" + jnum(f.y[i]) + ",\"K\":" + jnum(f.K[i]) + "}";
      }
      out += "],\"liquid\":{\"Z\":" + jnum(f.liquid.Z) + ",\"molar_volume\":" + jnum(f.liquid.molar_volume) +
             "},\"vapor\":{\"Z\":" + jnum(f.vapor.Z) + ",\"molar_volume\":" + jnum(f.vapor.molar_volume) + "}";
    }
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* cl_saturation_pressure_json(const char* smiles, double temperature_K, const char* eos) {
  if (!smiles) return nullptr;
  return guarded([&] {
    Diagnostics d;
    auto mix = mixture_from_smiles({smiles}, eos_from(eos), d);
    std::optional<double> p;
    if (mix) p = saturation_pressure(mix->components[0], eos_from(eos), temperature_K, d);
    std::string out = std::string("{\"ok\":") + (p ? "true" : "false");
    if (mix) out += ",\"name\":" + json_escape(mix->components[0].name);
    if (p) out += ",\"pressure_Pa\":" + jnum(*p);
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

void cl_free_string(char* s) { std::free(s); }

}  // extern "C"
