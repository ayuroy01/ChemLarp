/* OmniChem Core — stable C ABI for Python (ctypes/cffi), Node (N-API/ffi) and others.
 *
 * Two layers:
 *  1. Handle API: an opaque oc_system owning a ChemicalSystemState. The column
 *     accessors return pointers *into* the SoA storage — wrap them with
 *     numpy.ctypeslib.as_array / a JS ArrayBuffer for zero-copy views. Pointers
 *     are invalidated by the next call that adds species to the same system.
 *  2. One-shot JSON API for convenience. Returned strings are malloc'd; release
 *     them with oc_free_string.
 */
#ifndef OMNICHEM_C_H
#define OMNICHEM_C_H

#include <stddef.h>
#include <stdint.h>

/* Symbol visibility. Building the shared library defines OMNICHEM_BUILD_DLL;
 * linking the static library defines OMNICHEM_STATIC. */
#if defined(_WIN32) && !defined(OMNICHEM_STATIC)
#  if defined(OMNICHEM_BUILD_DLL)
#    define OMNICHEM_API __declspec(dllexport)
#  else
#    define OMNICHEM_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define OMNICHEM_API __attribute__((visibility("default")))
#else
#  define OMNICHEM_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct oc_system oc_system;

OMNICHEM_API const char* oc_version(void);

/* ---- handle API ---- */
OMNICHEM_API oc_system* oc_system_new(void);
OMNICHEM_API void oc_system_free(oc_system* sys);

/* Parse and append one species. Returns its index, or -1 (see oc_system_diagnostics_json). */
OMNICHEM_API int32_t oc_system_add_smiles(oc_system* sys, const char* smiles);
OMNICHEM_API int oc_system_set_moles(oc_system* sys, uint32_t species, double moles);
OMNICHEM_API int oc_system_set_conditions(oc_system* sys, double temperature_K, double pressure_Pa);
/* Gibbs-minimise in place; returns 1 if converged, 0 otherwise. */
OMNICHEM_API int oc_system_equilibrate(oc_system* sys);

OMNICHEM_API size_t oc_system_atom_count(const oc_system* sys);
OMNICHEM_API size_t oc_system_bond_count(const oc_system* sys);
OMNICHEM_API size_t oc_system_species_count(const oc_system* sys);
OMNICHEM_API const uint8_t* oc_system_atom_element(const oc_system* sys);
OMNICHEM_API const int8_t* oc_system_atom_charge(const oc_system* sys);
OMNICHEM_API const uint8_t* oc_system_atom_implicit_h(const oc_system* sys);
OMNICHEM_API double* oc_system_atom_x(oc_system* sys);
OMNICHEM_API double* oc_system_atom_y(oc_system* sys);
OMNICHEM_API double* oc_system_atom_z(oc_system* sys);
OMNICHEM_API const uint32_t* oc_system_bond_begin(const oc_system* sys);
OMNICHEM_API const uint32_t* oc_system_bond_end(const oc_system* sys);
OMNICHEM_API const uint8_t* oc_system_bond_order(const oc_system* sys);
OMNICHEM_API const double* oc_system_moles(const oc_system* sys);

/* JSON snapshots (malloc'd; free with oc_free_string). */
OMNICHEM_API char* oc_system_structure_json(const oc_system* sys);
OMNICHEM_API char* oc_system_thermo_json(const oc_system* sys);
OMNICHEM_API char* oc_system_diagnostics_json(const oc_system* sys); /* diagnostics of the most recent call */

/* ---- one-shot JSON API ---- */
OMNICHEM_API char* oc_parse_smiles_json(const char* smiles);
/* {"ok", "canonical", "random":[...n_random forms], "cip":{"centres":[{atom,label}],"double_bonds":[{bond,label}]}, "diagnostics"} */
OMNICHEM_API char* oc_canonical_smiles_json(const char* smiles, unsigned n_random, uint64_t seed);
OMNICHEM_API char* oc_check_reaction_json(const char* equation, double temperature_K);
OMNICHEM_API char* oc_equilibrate_json(const char* const* smiles, const double* moles, size_t count, double temperature_K,
                          double pressure_Pa);

/* Batch-reactor kinetics. `mechanism_text` is .ocm, CHEMKIN-II or Cantera YAML source (detected
 * from content); `thermo_text` (CHEMKIN therm.dat) and `smiles_map_text` ("NAME SMILES" lines) may
 * be NULL. `model` is "iso-v", "iso-p",
 * "adiabatic-v" or "adiabatic-p"; species are given by mechanism name. Returns
 * {"ok", "species":[names], "samples":[{t,T,P,n:[...]}], "ignition_time", "stats", "diagnostics"}
 * with `n_samples` evenly spaced samples in (0, t_end]. */
OMNICHEM_API char* oc_kinetics_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                                    const char* model, double temperature_K, double pressure_Pa,
                       const char* const* names, const double* moles, size_t count, double t_end, size_t n_samples,
                       double rtol);

/* 3D structure: parse, add hydrogens, distance-geometry embed, UFF-minimise.
 * {"ok", "canonical", "sdf", "xyz", "energy_kcal_mol", "converged", "attempts",
 *  "stereo": {"centres", "centres_ok", "bonds", "bonds_ok"}, "diagnostics"} */
OMNICHEM_API char* oc_embed_json(const char* smiles, uint64_t seed, int candidates, int minimize);

/* String-only variants for JavaScript / WebAssembly callers (no pointer arrays).
 * `species_lines`: one "SMILES moles" pair per line.  `feed_lines`: "NAME=moles" or "NAME moles" per line. */
OMNICHEM_API char* oc_canonical_smiles32_json(const char* smiles, unsigned n_random, unsigned seed);
OMNICHEM_API char* oc_embed32_json(const char* smiles, unsigned seed, int candidates, int minimize);
OMNICHEM_API char* oc_equilibrate_text_json(const char* species_lines, double temperature_K, double pressure_Pa);
OMNICHEM_API char* oc_kinetics_text_json(const char* mechanism_text, const char* thermo_text, const char* smiles_map_text,
                                         const char* model, double temperature_K, double pressure_Pa,
                                         const char* feed_lines, double t_end, size_t n_samples, double rtol);

/* Vapour-liquid equilibrium with a cubic equation of state. `eos` is "pr" (Peng-Robinson) or "srk".
 * Flash: `component_lines` has one "SMILES mole_fraction" per line. Returns
 * {"ok","region","vapor_fraction","converged","iterations","components":[{smiles,name,z,x,y,K}],
 *  "liquid":{"Z","molar_volume"},"vapor":{...},"diagnostics"}.
 * Saturation: {"ok","name","pressure_Pa","diagnostics"}. */
OMNICHEM_API char* oc_flash_text_json(const char* component_lines, double temperature_K, double pressure_Pa, const char* eos);
OMNICHEM_API char* oc_saturation_pressure_json(const char* smiles, double temperature_K, const char* eos);

/* 2D structure drawing. bond_px: bond length in pixels (<= 0 for 30). flags: OC_DEPICT_* bits.
 * {"ok","svg","molfile","width","height","atoms":[{atom,x,y,label}],"bonds":[{a,b,order,stereo}],
 *  "stereo":{"centres","centres_ok","bonds","bonds_ok"},"diagnostics"}.
 * atom is the 0-based atom index in the input, or -1 for a hydrogen drawn to carry a wedge;
 * stereo is 0, 1 (wedge, narrow end at a) or 6 (hash). */
#define OC_DEPICT_COLOR 1
#define OC_DEPICT_CIP 2
#define OC_DEPICT_NUMBERS 4
OMNICHEM_API char* oc_depict_json(const char* smiles, double bond_px, int flags);

/* Standard InChI and InChIKey of a SMILES string (clean-room generator; see inchi.hpp for limits).
 * JSON: {"ok","inchi","key","diagnostics"}. */
OMNICHEM_API char* oc_inchi_json(const char* smiles);

/* Reads a MOL/SDF, PDB or CML file (see molio.hpp). `format` is "mol", "sdf", "pdb", "cml", a file name
 * whose extension decides, or "" to detect from the content. max_records 0 = all.
 * JSON: {"ok","format","records":[{"ok","line","title","smiles","formula","molar_mass","atoms","dimensions",
 *        "properties":[[name,value],...]}],"diagnostics"}. "ok" is true when at least one record was read. */
OMNICHEM_API char* oc_read_molecules_json(const char* text, const char* format, int keep_hydrogens, int max_records);

/* SMARTS matches in one SMILES (see smarts.hpp). Atoms are 0-based indices in SMILES order.
 * JSON: {"ok","count","matches":[[atom,...],...],"diagnostics"}. max_matches 0 = no limit. */
OMNICHEM_API char* oc_smarts_match_json(const char* smiles, const char* smarts, int use_chirality, int max_matches);

/* Morgan fingerprint similarity of two SMILES (radius 2 and 2048 bits when given 0).
 * JSON: {"ok","tanimoto","dice","bits_a","bits_b","diagnostics"}. */
OMNICHEM_API char* oc_similarity_json(const char* smiles_a, const char* smiles_b, int radius, int bits);

/* Screens one SMILES per line: SMARTS match count (when smarts is non-empty) and Tanimoto to
 * query_smiles (when non-empty). JSON: {"ok","results":[{"line","ok","count","similarity"}],"diagnostics"};
 * a line that does not parse has "ok": false. */
OMNICHEM_API char* oc_search_json(const char* smiles_lines, const char* smarts, const char* query_smiles, int radius, int bits);

OMNICHEM_API void oc_free_string(char* s);

#ifdef __cplusplus
}
#endif

#endif /* OMNICHEM_C_H */
