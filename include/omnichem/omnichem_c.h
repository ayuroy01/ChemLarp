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

OMNICHEM_API void oc_free_string(char* s);

#ifdef __cplusplus
}
#endif

#endif /* OMNICHEM_C_H */
