// This translation unit defines the exported C ABI: export unless the static library is being built.
#if !defined(OMNICHEM_STATIC) && !defined(OMNICHEM_BUILD_DLL)
#define OMNICHEM_BUILD_DLL
#endif
#include "omnichem/omnichem_c.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <optional>
#include <string>
#include <vector>

#include "omnichem/canonical.hpp"
#include "omnichem/cip.hpp"
#include "omnichem/conformer.hpp"
#include "omnichem/depict.hpp"
#include "omnichem/diagnostics.hpp"
#include "omnichem/eos.hpp"
#include "omnichem/gibbs.hpp"
#include "omnichem/inchi.hpp"
#include "omnichem/molio.hpp"
#include "omnichem/json.hpp"
#include "omnichem/kinetics.hpp"
#include "omnichem/reaction.hpp"
#include "omnichem/smiles.hpp"
#include "omnichem/system_state.hpp"

using namespace omnichem;

struct oc_system {
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

const char* oc_version(void) { return "0.1.0"; }

oc_system* oc_system_new(void) { return new (std::nothrow) oc_system(); }
void oc_system_free(oc_system* sys) { delete sys; }

int32_t oc_system_add_smiles(oc_system* sys, const char* smiles) {
  if (!sys || !smiles) return -1;
  sys->last.clear();
  try {
    return parse_smiles(sys->state, smiles, sys->last);
  } catch (...) {
    return -1;
  }
}

int oc_system_set_moles(oc_system* sys, uint32_t species, double moles) {
  if (!sys || species >= sys->state.species.size() || !(moles >= 0.0) || !std::isfinite(moles)) return 0;
  sys->state.moles[species] = moles;
  ++sys->state.revision;
  return 1;
}

int oc_system_set_conditions(oc_system* sys, double temperature_K, double pressure_Pa) {
  if (!sys || !(temperature_K > 0.0) || !(pressure_Pa > 0.0)) return 0;
  sys->state.conditions = {temperature_K, pressure_Pa};
  ++sys->state.revision;
  return 1;
}

int oc_system_equilibrate(oc_system* sys) {
  if (!sys) return 0;
  sys->last.clear();
  try {
    return equilibrate(sys->state, sys->last).converged ? 1 : 0;
  } catch (...) {
    return 0;
  }
}

size_t oc_system_atom_count(const oc_system* s) { return s ? s->state.atoms.size() : 0; }
size_t oc_system_bond_count(const oc_system* s) { return s ? s->state.bonds.size() : 0; }
size_t oc_system_species_count(const oc_system* s) { return s ? s->state.species.size() : 0; }
const uint8_t* oc_system_atom_element(const oc_system* s) { return s ? s->state.atoms.element.data() : nullptr; }
const int8_t* oc_system_atom_charge(const oc_system* s) { return s ? s->state.atoms.formal_charge.data() : nullptr; }
const uint8_t* oc_system_atom_implicit_h(const oc_system* s) { return s ? s->state.atoms.implicit_h.data() : nullptr; }
double* oc_system_atom_x(oc_system* s) { return s ? s->state.atoms.x.data() : nullptr; }
double* oc_system_atom_y(oc_system* s) { return s ? s->state.atoms.y.data() : nullptr; }
double* oc_system_atom_z(oc_system* s) { return s ? s->state.atoms.z.data() : nullptr; }
const uint32_t* oc_system_bond_begin(const oc_system* s) { return s ? s->state.bonds.begin.data() : nullptr; }
const uint32_t* oc_system_bond_end(const oc_system* s) { return s ? s->state.bonds.end.data() : nullptr; }
const uint8_t* oc_system_bond_order(const oc_system* s) { return s ? s->state.bonds.order.data() : nullptr; }
const double* oc_system_moles(const oc_system* s) { return s ? s->state.moles.data() : nullptr; }

char* oc_system_structure_json(const oc_system* s) {
  if (!s) return nullptr;
  return guarded([&] { return structure_to_json(s->state); });
}
char* oc_system_thermo_json(const oc_system* s) {
  if (!s) return nullptr;
  return guarded([&] { return thermo_to_json(s->state); });
}
char* oc_system_diagnostics_json(const oc_system* s) {
  if (!s) return nullptr;
  return guarded([&] { return diagnostics_to_json(s->last); });
}

char* oc_parse_smiles_json(const char* smiles) {
  if (!smiles) return nullptr;
  return guarded([&] {
    ChemicalSystemState st;
    Diagnostics d;
    const bool ok = parse_smiles(st, smiles, d) >= 0;
    return std::string("{\"ok\":") + (ok ? "true" : "false") + ",\"structure\":" + structure_to_json(st) +
           ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* oc_check_reaction_json(const char* equation, double temperature_K) {
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

char* oc_equilibrate_json(const char* const* smiles, const double* moles, size_t count, double temperature_K,
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

char* oc_canonical_smiles_json(const char* smiles, unsigned n_random, uint64_t seed) {
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

char* oc_kinetics_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
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

char* oc_embed_json(const char* smiles, uint64_t seed, int candidates, int minimize) {
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


char* oc_depict_json(const char* smiles, double bond_px, int flags) {
  if (!smiles) return nullptr;
  return guarded([&] {
    Diagnostics d;
    DepictOptions opt;
    if (bond_px > 0) opt.bond_length = bond_px;
    opt.color = (flags & OC_DEPICT_COLOR) != 0;
    opt.cip_labels = (flags & OC_DEPICT_CIP) != 0;
    opt.atom_numbers = (flags & OC_DEPICT_NUMBERS) != 0;
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

char* oc_inchi_json(const char* smiles) {
  if (!smiles) return nullptr;
  return guarded([&] {
    Diagnostics d;
    const InchiResult r = inchi_from_smiles(smiles, d);
    std::string out = std::string("{\"ok\":") + (r.ok ? "true" : "false");
    if (r.ok) out += ",\"inchi\":" + json_escape(r.inchi) + ",\"key\":" + json_escape(r.key);
    return out + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
  });
}

char* oc_read_molecules_json(const char* text, const char* format, int keep_hydrogens, int max_records) {
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

char* oc_canonical_smiles32_json(const char* smiles, unsigned n_random, unsigned seed) {
  return oc_canonical_smiles_json(smiles, n_random, seed);
}

char* oc_embed32_json(const char* smiles, unsigned seed, int candidates, int minimize) {
  return oc_embed_json(smiles, seed, candidates, minimize);
}

char* oc_equilibrate_text_json(const char* species_lines, double temperature_K, double pressure_Pa) {
  std::vector<std::string> names;
  std::vector<double> moles;
  std::string bad;
  if (!parse_pairs(species_lines, names, moles, bad)) return bad_input(bad);
  std::vector<const char*> ptrs;
  for (const auto& n : names) ptrs.push_back(n.c_str());
  return oc_equilibrate_json(ptrs.data(), moles.data(), ptrs.size(), temperature_K, pressure_Pa);
}

char* oc_kinetics_text_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                            const char* model, double temperature_K, double pressure_Pa, const char* feed_lines,
                            double t_end, size_t n_samples, double rtol) {
  std::vector<std::string> names;
  std::vector<double> moles;
  std::string bad;
  if (!parse_pairs(feed_lines, names, moles, bad)) return bad_input(bad);
  std::vector<const char*> ptrs;
  for (const auto& n : names) ptrs.push_back(n.c_str());
  return oc_kinetics_json(mechanism_text, thermo_text, smiles_map_text, model, temperature_K, pressure_Pa, ptrs.data(),
                          moles.data(), ptrs.size(), t_end, n_samples, rtol);
}


char* oc_flash_text_json(const char* component_lines, double temperature_K, double pressure_Pa, const char* eos) {
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

char* oc_saturation_pressure_json(const char* smiles, double temperature_K, const char* eos) {
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

void oc_free_string(char* s) { std::free(s); }

}  // extern "C"
