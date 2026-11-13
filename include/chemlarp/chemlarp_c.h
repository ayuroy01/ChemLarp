/* ChemLarp — stable C ABI for Python (ctypes/cffi), Node (N-API/ffi) and others.
 *
 * Two layers:
 *  1. Handle API: an opaque cl_system owning a ChemicalSystemState. The column
 *     accessors return pointers *into* the SoA storage — wrap them with
 *     numpy.ctypeslib.as_array / a JS ArrayBuffer for zero-copy views. Pointers
 *     are invalidated by the next call that adds species to the same system.
 *  2. One-shot JSON API for convenience. Returned strings are malloc'd; release
 *     them with cl_free_string.
 */
#ifndef CHEMLARP_C_H
#define CHEMLARP_C_H

#include <stddef.h>
#include <stdint.h>

/* Symbol visibility. Building the shared library defines CHEMLARP_BUILD_DLL;
 * linking the static library defines CHEMLARP_STATIC. */
#if defined(_WIN32) && !defined(CHEMLARP_STATIC)
#  if defined(CHEMLARP_BUILD_DLL)
#    define CHEMLARP_API __declspec(dllexport)
#  else
#    define CHEMLARP_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define CHEMLARP_API __attribute__((visibility("default")))
#else
#  define CHEMLARP_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cl_system cl_system;

CHEMLARP_API const char* cl_version(void);

/* ---- handle API ---- */
CHEMLARP_API cl_system* cl_system_new(void);
CHEMLARP_API void cl_system_free(cl_system* sys);

/* Parse and append one species. Returns its index, or -1 (see cl_system_diagnostics_json). */
CHEMLARP_API int32_t cl_system_add_smiles(cl_system* sys, const char* smiles);
CHEMLARP_API int cl_system_set_moles(cl_system* sys, uint32_t species, double moles);
CHEMLARP_API int cl_system_set_conditions(cl_system* sys, double temperature_K, double pressure_Pa);
/* Gibbs-minimise in place; returns 1 if converged, 0 otherwise. */
CHEMLARP_API int cl_system_equilibrate(cl_system* sys);

CHEMLARP_API size_t cl_system_atom_count(const cl_system* sys);
CHEMLARP_API size_t cl_system_bond_count(const cl_system* sys);
CHEMLARP_API size_t cl_system_species_count(const cl_system* sys);
CHEMLARP_API const uint8_t* cl_system_atom_element(const cl_system* sys);
CHEMLARP_API const int8_t* cl_system_atom_charge(const cl_system* sys);
CHEMLARP_API const uint8_t* cl_system_atom_implicit_h(const cl_system* sys);
CHEMLARP_API double* cl_system_atom_x(cl_system* sys);
CHEMLARP_API double* cl_system_atom_y(cl_system* sys);
CHEMLARP_API double* cl_system_atom_z(cl_system* sys);
CHEMLARP_API const uint32_t* cl_system_bond_begin(const cl_system* sys);
CHEMLARP_API const uint32_t* cl_system_bond_end(const cl_system* sys);
CHEMLARP_API const uint8_t* cl_system_bond_order(const cl_system* sys);
CHEMLARP_API const double* cl_system_moles(const cl_system* sys);

/* JSON snapshots (malloc'd; free with cl_free_string). */
CHEMLARP_API char* cl_system_structure_json(const cl_system* sys);
CHEMLARP_API char* cl_system_thermo_json(const cl_system* sys);
CHEMLARP_API char* cl_system_diagnostics_json(const cl_system* sys); /* diagnostics of the most recent call */

/* ---- one-shot JSON API ---- */
CHEMLARP_API char* cl_parse_smiles_json(const char* smiles);
/* {"ok", "canonical", "random":[...n_random forms], "cip":{"centres":[{atom,label}],"double_bonds":[{bond,label}]}, "diagnostics"} */
CHEMLARP_API char* cl_canonical_smiles_json(const char* smiles, unsigned n_random, uint64_t seed);
CHEMLARP_API char* cl_check_reaction_json(const char* equation, double temperature_K);
CHEMLARP_API char* cl_equilibrate_json(const char* const* smiles, const double* moles, size_t count, double temperature_K,
                          double pressure_Pa);

/* Batch-reactor kinetics. `mechanism_text` is .ocm, CHEMKIN-II or Cantera YAML source (detected
 * from content); `thermo_text` (CHEMKIN therm.dat) and `smiles_map_text` ("NAME SMILES" lines) may
 * be NULL. `model` is "iso-v", "iso-p",
 * "adiabatic-v" or "adiabatic-p"; species are given by mechanism name. Returns
 * {"ok", "species":[names], "samples":[{t,T,P,n:[...]}], "ignition_time", "stats", "diagnostics"}
 * with `n_samples` evenly spaced samples in (0, t_end]. */
CHEMLARP_API char* cl_kinetics_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                                    const char* model, double temperature_K, double pressure_Pa,
                       const char* const* names, const double* moles, size_t count, double t_end, size_t n_samples,
                       double rtol);

/* 3D structure: parse, add hydrogens, distance-geometry embed, UFF-minimise.
 * {"ok", "canonical", "sdf", "xyz", "energy_kcal_mol", "converged", "attempts",
 *  "stereo": {"centres", "centres_ok", "bonds", "bonds_ok"}, "diagnostics"} */
CHEMLARP_API char* cl_embed_json(const char* smiles, uint64_t seed, int candidates, int minimize);

/* String-only variants for JavaScript / WebAssembly callers (no pointer arrays).
 * `species_lines`: one "SMILES moles" pair per line.  `feed_lines`: "NAME=moles" or "NAME moles" per line. */
CHEMLARP_API char* cl_canonical_smiles32_json(const char* smiles, unsigned n_random, unsigned seed);
CHEMLARP_API char* cl_embed32_json(const char* smiles, unsigned seed, int candidates, int minimize);
CHEMLARP_API char* cl_equilibrate_text_json(const char* species_lines, double temperature_K, double pressure_Pa);
CHEMLARP_API char* cl_kinetics_text_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                                         const char* model, double temperature_K, double pressure_Pa,
                                         const char* feed_lines, double t_end, size_t n_samples, double rtol);

/* Vapour-liquid equilibrium with a cubic equation of state. `eos` is "pr" (Peng-Robinson) or "srk".
 * Flash: `component_lines` has one "SMILES mole_fraction" per line. Returns
 * {"ok","region","vapor_fraction","converged","iterations","components":[{smiles,name,z,x,y,K}],
 *  "liquid":{"Z","molar_volume"},"vapor":{...},"diagnostics"}.
 * Saturation: {"ok","name","pressure_Pa","diagnostics"}. */
CHEMLARP_API char* cl_flash_text_json(const char* component_lines, double temperature_K, double pressure_Pa, const char* eos);
CHEMLARP_API char* cl_saturation_pressure_json(const char* smiles, double temperature_K, const char* eos);

/* Vapour-liquid equilibrium with UNIFAC activity coefficients and Antoine vapour pressures (modified Raoult).
 * `mode` is "bubble-T" or "dew-T" (value = pressure in Pa) or "bubble-P" or "dew-P" (value = temperature in K);
 * `component_lines` has one "SMILES mole_fraction" per line (liquid x for bubble, vapour y for dew). Returns
 * {"ok","T_K","P_Pa","extrapolated","iterations","components":[{smiles,name,x,y,gamma}],"diagnostics"}. */
/* Pure-component properties (constants from published compilations; see scripts/gen_property_data.py).
 * With temperature_K > 0 also "at":{T_K,P_Pa,cp_ideal_gas_JmolK,H_ideal_gas_Jmol (H(T) - H(298.15 K)),
 * S_ideal_gas_JmolK,Hvap_Jmol,V_liquid_m3mol,rho_liquid_kgm3,viscosity_liquid_Pas,viscosity_gas_Pas,
 * conductivity_liquid_WmK,conductivity_gas_WmK,surface_tension_Nm,cp_liquid_JmolK} (gas values at low pressure); pressure_Pa <= 0 means 1 bar. Missing values
 * are null; a correlation that cannot be evaluated adds a warning. */
CHEMLARP_API char* cl_pure_component_json(const char* smiles, double temperature_K, double pressure_Pa);

CHEMLARP_API char* cl_vle_text_json(const char* component_lines, const char* mode, double value);

/* Steady-state flowsheet: units connected by named streams, recycles converged by tear-stream iteration
 * (Wegstein or direct substitution). `flowsheet_text` uses the format documented in chemlarp/flowsheet.hpp
 * (components, model, feed, guess, tear, method, tolerance, max-iterations, and one line per unit:
 * mixer, splitter, heater, valve, flash, separator, pump, compressor, exchanger, reactor). Returns
 * {"ok","converged","iterations","method","tears","order","residuals","components","streams","units",
 * "balance","diagnostics"}; see flowsheet.hpp for the fields. */
CHEMLARP_API char* cl_flowsheet_json(const char* flowsheet_text);

/* 2D structure drawing. bond_px: bond length in pixels (<= 0 for 30). flags: CL_DEPICT_* bits.
 * {"ok","svg","molfile","width","height","atoms":[{atom,x,y,label}],"bonds":[{a,b,order,stereo}],
 *  "stereo":{"centres","centres_ok","bonds","bonds_ok"},"diagnostics"}.
 * atom is the 0-based atom index in the input, or -1 for a hydrogen drawn to carry a wedge;
 * stereo is 0, 1 (wedge, narrow end at a) or 6 (hash). */
#define CL_DEPICT_COLOR 1
#define CL_DEPICT_CIP 2
#define CL_DEPICT_NUMBERS 4
CHEMLARP_API char* cl_depict_json(const char* smiles, double bond_px, int flags);

/* Standard InChI and InChIKey of a SMILES string (clean-room generator; see inchi.hpp for limits).
 * JSON: {"ok","inchi","key","diagnostics"}. */
CHEMLARP_API char* cl_inchi_json(const char* smiles);

/* Systematic IUPAC name of a SMILES string (2013 Recommendations; see iupac.hpp for scope).
 * JSON: {"ok","name","diagnostics"}; ok false with a NamingUnsupported diagnostic outside the scope. */
CHEMLARP_API char* cl_iupac_name_json(const char* smiles);

/* Reads a MOL/SDF, PDB or CML file (see molio.hpp). `format` is "mol", "sdf", "pdb", "cml", a file name
 * whose extension decides, or "" to detect from the content. max_records 0 = all.
 * JSON: {"ok","format","records":[{"ok","line","title","smiles","formula","molar_mass","atoms","dimensions",
 *        "properties":[[name,value],...]}],"diagnostics"}. "ok" is true when at least one record was read. */
CHEMLARP_API char* cl_read_molecules_json(const char* text, const char* format, int keep_hydrogens, int max_records);

/* SMARTS matches in one SMILES (see smarts.hpp). Atoms are 0-based indices in SMILES order.
 * JSON: {"ok","count","matches":[[atom,...],...],"diagnostics"}. max_matches 0 = no limit. */
CHEMLARP_API char* cl_smarts_match_json(const char* smiles, const char* smarts, int use_chirality, int max_matches);

/* Morgan fingerprint similarity of two SMILES (radius 2 and 2048 bits when given 0).
 * JSON: {"ok","tanimoto","dice","bits_a","bits_b","diagnostics"}. */
CHEMLARP_API char* cl_similarity_json(const char* smiles_a, const char* smiles_b, int radius, int bits);

/* Screens one SMILES per line: SMARTS match count (when smarts is non-empty) and Tanimoto to
 * query_smiles (when non-empty). JSON: {"ok","results":[{"line","ok","count","similarity"}],"diagnostics"};
 * a line that does not parse has "ok": false. */
/* Quantum chemistry single point (see qc.hpp). `input` is a SMILES (embedded to 3D with explicit
 * hydrogens first) or a 3D molfile / SDF record. method: "hf", "svwn", "pbe", "b3lyp"; basis:
 * "STO-3G", "3-21G", "6-31G", "6-31G*"; multiplicity 0 = lowest; grid_level 1-3; memory_mb <= 0
 * for the default limit on stored two-electron integrals.
 * JSON: {"ok","converged","method","basis","energy","nuclear_repulsion","one_electron","two_electron",
 * "xc","iterations","n_basis","n_alpha","n_beta","multiplicity","charge","homo","lumo",
 * "orbitals_alpha":[...],"orbitals_beta":[...],"mulliken":[...],"elements":[...],"xyz":[[x,y,z],...],
 * "dipole":[x,y,z],"dipole_total","s2","electrons_on_grid","diagnostics"} (hartree, angstrom, debye). */
CHEMLARP_API char* cl_qc_json(const char* input, const char* method, const char* basis, int multiplicity, int grid_level,
                              double memory_mb);

CHEMLARP_API char* cl_search_json(const char* smiles_lines, const char* smarts, const char* query_smiles, int radius, int bits);

CHEMLARP_API void cl_free_string(char* s);

#ifdef __cplusplus
}
#endif

#endif /* CHEMLARP_C_H */
