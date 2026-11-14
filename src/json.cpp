#include "chemlarp/json.hpp"

#include <cmath>
#include <cstdio>

#include "chemlarp/cip.hpp"
#include "chemlarp/elements.hpp"

namespace chemlarp {
namespace {

std::string num(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

template <class Vec, class F>
std::string json_array(const Vec& v, F&& fmt) {
  std::string s = "[";
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (i) s += ',';
    s += fmt(v[i]);
  }
  return s + "]";
}

auto as_int = [](auto v) { return std::to_string(static_cast<long long>(v)); };

const char* phase_name(Phase p) {
  switch (p) {
    case Phase::Gas: return "gas";
    case Phase::Liquid: return "liquid";
    case Phase::Solid: return "solid";
    case Phase::Aqueous: return "aqueous";
  }
  return "gas";
}

}  // namespace

std::string json_escape(std::string_view s) {
  std::string out = "\"";
  for (char ch : s) {
    const unsigned char c = static_cast<unsigned char>(ch);
    switch (ch) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (c < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out += ch;
        }
    }
  }
  return out + "\"";
}

std::string diagnostics_to_json(const Diagnostics& d) {
  std::string s = "[";
  bool first = true;
  for (const auto& x : d.items()) {
    if (!first) s += ',';
    first = false;
    s += "{\"code\":" + json_escape(to_string(x.code));
    s += ",\"code_id\":" + std::to_string(static_cast<int>(x.code));
    s += ",\"severity\":" + json_escape(to_string(x.severity));
    s += ",\"position\":" + (x.position == Diagnostic::npos ? std::string("null") : std::to_string(x.position));
    s += ",\"atom\":" + (x.atom < 0 ? std::string("null") : std::to_string(x.atom));
    s += ",\"species\":" + (x.species < 0 ? std::string("null") : std::to_string(x.species));
    if (x.code == DiagCode::ElementImbalance) {
      s += ",\"element\":" + json_escape(element_symbol(x.element)) + ",\"delta\":" + num(x.delta);
    } else if (x.code == DiagCode::ChargeImbalance) {
      s += ",\"delta\":" + num(x.delta);
    }
    s += ",\"message\":" + json_escape(x.message) + "}";
  }
  return s + "]";
}

std::string structure_to_json(const ChemicalSystemState& st) {
  const auto& a = st.atoms;
  const auto& b = st.bonds;
  std::string s = "{\"atoms\":{";
  s += "\"count\":" + std::to_string(a.size());
  s += ",\"element\":" + json_array(a.element, as_int);
  s += ",\"formal_charge\":" + json_array(a.formal_charge, as_int);
  s += ",\"implicit_h\":" + json_array(a.implicit_h, as_int);
  s += ",\"flags\":" + json_array(a.flags, as_int);
  s += ",\"isotope\":" + json_array(a.isotope, as_int);
  s += ",\"chirality\":" + json_array(a.chirality, as_int);
  s += ",\"atom_map\":" + json_array(a.atom_map, as_int);
  s += ",\"species\":" + json_array(a.species, as_int);
  s += ",\"x\":" + json_array(a.x, num) + ",\"y\":" + json_array(a.y, num) + ",\"z\":" + json_array(a.z, num);
  s += "},\"bonds\":{";
  s += "\"count\":" + std::to_string(b.size());
  s += ",\"begin\":" + json_array(b.begin, as_int);
  s += ",\"end\":" + json_array(b.end, as_int);
  s += ",\"order\":" + json_array(b.order, as_int);
  s += ",\"flags\":" + json_array(b.flags, as_int);
  s += ",\"species\":" + json_array(b.species, as_int);
  s += "},\"species\":[";
  for (std::size_t i = 0; i < st.species.size(); ++i) {
    const Species& sp = st.species[i];
    if (i) s += ',';
    s += "{\"smiles\":" + json_escape(sp.smiles) + ",\"canonical_smiles\":" + json_escape(sp.canonical_smiles) +
         ",\"formula\":" + json_escape(sp.formula);
    s += ",\"atom_begin\":" + std::to_string(sp.atom_begin) + ",\"atom_count\":" + std::to_string(sp.atom_count);
    s += ",\"bond_begin\":" + std::to_string(sp.bond_begin) + ",\"bond_count\":" + std::to_string(sp.bond_count);
    s += ",\"net_charge\":" + std::to_string(sp.net_charge) + ",\"molar_mass\":" + num(sp.molar_mass);
    s += ",\"phase\":" + json_escape(phase_name(sp.phase)) + ",\"composition\":{";
    bool first = true;
    for (std::size_t z = 1; z < kElementSlots; ++z) {
      if (sp.composition[z] == 0) continue;
      if (!first) s += ',';
      first = false;
      s += json_escape(element_symbol(static_cast<std::uint8_t>(z))) + ":" + std::to_string(sp.composition[z]);
    }
    s += "}}";
  }
  auto ref = [](std::uint32_t r) {
    if (r == kImplicitHRef) return std::string("\"H\"");
    if (r == kLonePairRef) return std::string("\"LP\"");
    return std::to_string(r);
  };
  s += "],\"stereo\":{\"tetrahedral\":[";
  for (std::size_t i = 0; i < st.tetrahedral.size(); ++i) {
    const auto& t = st.tetrahedral[i];
    if (i) s += ',';
    s += "{\"atom\":" + std::to_string(t.atom) + ",\"refs\":[" + ref(t.refs[0]) + "," + ref(t.refs[1]) + "," +
         ref(t.refs[2]) + "," + ref(t.refs[3]) + "],\"parity\":" +
         json_escape(t.parity == Chirality::Clockwise ? "@@" : "@") + ",\"cip\":" + json_escape(to_string(cip_label(st, t))) + "}";
  }
  s += "],\"double_bonds\":[";
  for (std::size_t i = 0; i < st.double_bonds.size(); ++i) {
    const auto& d = st.double_bonds[i];
    if (i) s += ',';
    s += "{\"bond\":" + std::to_string(d.bond) + ",\"ref_begin\":" + std::to_string(d.ref_begin) +
         ",\"ref_end\":" + std::to_string(d.ref_end) + ",\"cis\":" + (d.cis ? "true" : "false") +
         ",\"cip\":" + json_escape(to_string(cip_label(st, d))) + "}";
  }
  s += "],\"axial\":[";
  for (std::size_t i = 0; i < st.axial.size(); ++i) {
    const auto& ax = st.axial[i];
    if (i) s += ',';
    s += "{\"begin\":" + std::to_string(ax.begin) + ",\"end\":" + std::to_string(ax.end) + ",\"centre\":" +
         (ax.centre == kNoAtom ? std::string("null") : std::to_string(ax.centre)) + ",\"refs\":[" + ref(ax.refs[0]) + "," +
         ref(ax.refs[1]) + "," + ref(ax.refs[2]) + "," + ref(ax.refs[3]) + "],\"parity\":" +
         json_escape(ax.parity == Chirality::Clockwise ? "@@" : "@") + ",\"cip\":" + json_escape(to_string(cip_label(st, ax))) + "}";
  }
  return s + "]}}";
}

std::string thermo_to_json(const ChemicalSystemState& st) {
  std::string s = "{\"temperature_K\":" + num(st.conditions.temperature_K);
  s += ",\"pressure_Pa\":" + num(st.conditions.pressure_Pa);
  s += ",\"moles\":" + json_array(st.moles, num);
  s += ",\"g_RT\":" + json_array(st.g_RT, num) + "}";
  return s;
}

std::string equilibrium_to_json(const ChemicalSystemState& st, const EquilibriumResult& r) {
  std::string s = "{\"converged\":" + std::string(r.converged ? "true" : "false");
  s += ",\"iterations\":" + std::to_string(r.iterations);
  s += ",\"total_moles\":" + num(r.total_moles);
  s += ",\"G_RT\":" + num(r.G_RT);
  s += ",\"max_relative_element_error\":" + num(r.max_relative_element_error);
  s += ",\"species\":[";
  for (std::size_t i = 0; i < r.n.size() && i < st.species.size(); ++i) {
    if (i) s += ',';
    s += "{\"smiles\":" + json_escape(st.species[i].smiles) + ",\"formula\":" + json_escape(st.species[i].formula);
    s += ",\"moles\":" + num(r.n[i]) + ",\"mole_fraction\":" + num(r.mole_fraction[i]) + "}";
  }
  s += "],\"element_potentials\":{";
  for (std::size_t j = 0; j < r.elements.size(); ++j) {
    if (j) s += ',';
    s += json_escape(element_symbol(r.elements[j])) + ":" + num(r.element_potentials[j]);
  }
  return s + "}}";
}

}  // namespace chemlarp
