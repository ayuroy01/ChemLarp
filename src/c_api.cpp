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
#include <string>

#include "omnichem/canonical.hpp"
#include "omnichem/cip.hpp"
#include "omnichem/diagnostics.hpp"
#include "omnichem/gibbs.hpp"
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

char* oc_kinetics_json(const char* mechanism_text, const char* model, double temperature_K, double pressure_Pa,
                       const char* const* names, const double* moles, size_t count, double t_end, size_t n_samples,
                       double rtol) {
  if (!mechanism_text || !model || ((!names || !moles) && count > 0)) return nullptr;
  return guarded([&] {
    Diagnostics d;
    std::string out = "{\"ok\":";
    auto mech = parse_mechanism(mechanism_text, d);
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

void oc_free_string(char* s) { std::free(s); }

}  // extern "C"
