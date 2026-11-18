// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
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
 * be NULL; `mechanism_text` and `model` may not, and NULL is returned if they are. `model` is "iso-v", "iso-p",
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

/* Aqueous solution chemistry (pH, speciation, saturation indices, equilibrium with minerals and gases) with
 * data from PHREEQC's phreeqc.dat. `solution_text` uses the format documented in chemlarp/aqueous.hpp
 * (temperature, units, pH, component totals, add <formula>, equilibrium <phase>). Returns
 * {"ok","T_K","pH","ionic_strength","water_activity","A_dh","B_dh","charge_balance_eq","alkalinity_eq",
 * "iterations","totals","species","saturation","phases","diagnostics"}. */
CHEMLARP_API char* cl_aqueous_json(const char* solution_text);

/* Molecular dynamics of one molecule with UFF (SMILES embedded in 3D with explicit hydrogens and minimized).
 * thermostat: "none" (NVE), "bussi" or "langevin"; tau_fs: coupling time (Bussi) or 1/friction (Langevin);
 * compute_modes != 0 adds the harmonic normal-mode frequencies. Returns {"ok","elements","bonds",
 * "minimized_energy","dof","mean_temperature_K","mean_potential","drift_kcal_mol_ps","fluctuation",
 * "frames":[{t,K,U,E,conserved,T}],"positions":[[3N per stored frame]],"position_stride","modes":[cm-1],
 * "spectrum":{"max_cm1","intensity"},"xyz","diagnostics"}; units Å, fs, kcal/mol, K. force_field: "uff"
 * (default), "mmff94" or "mmff94s". */
CHEMLARP_API char* cl_md_json(const char* smiles, double temperature_K, unsigned steps, double dt_fs,
                              const char* thermostat, double tau_fs, unsigned sample_every, unsigned seed,
                              int compute_modes, const char* force_field);

/* MMFF94 / MMFF94s for a SMILES (embedded in 3D with explicit hydrogens; minimized when minimize != 0).
 * variant: "MMFF94" (default) or "MMFF94s". Returns {"ok","variant","elements","bonds","types","symbols",
 * "charges","formal_charges","energy":{total,bond,angle,stretch_bend,out_of_plane,torsion,vdw,electrostatic},
 * "minimized","converged","iterations","empirical_parameters","sdf","xyz","diagnostics"}. */
CHEMLARP_API char* cl_mmff_json(const char* smiles, const char* variant, int minimize);

/* Protein from PDB or mmCIF text: chains, sequences, disulfides; missing side-chain atoms and missing loop
 * residues built; pKa of each titratable group from the structure; protonation states at the pH with the
 * hydrogen-bond network optimised (Asn/Gln/His flips, His tautomers, OH/SH/NH3 rotors); hydrogens rebuilt.
 * Returns {"ok","title","chains":[{id,sequence,residues,pI,molecular_weight}],"residues":[titratable residues
 * and termini {chain,number,icode,name,state,charge,n_terminal,c_terminal,predicted}],"pka":[{chain,number,
 * name,group,reference,pKa,desolvation,hydrogen_bonds,charges,buried,features}],"disulfides","chain_breaks",
 * "missing_atoms","unknown_residues","flipped","undetermined","repair":{atoms_added,completed_residues,
 * built_residues,not_built,missing_residues_listed,missing_residues_left},"atoms","hydrogens_added",
 * "hydrogens_removed","net_charge","pI","titration":[[pH,charge]],"titration_structure":[[pH,charge]],
 * "format":"pdb"|"cif","pdb","cif","diagnostics"}. The structure is in "pdb", or in "cif" when the PDB
 * format cannot hold it (more than 99,999 atoms or chain names of several characters).
 * cl_protein_json_opts flags (sum): 1 leave missing atoms, 2 do not build loops, 4 also build missing chain
 * ends, 8 no flips, 16 reference pKa values instead of structure-based ones. */
CHEMLARP_API char* cl_protein_json(const char* pdb_text, double pH);
CHEMLARP_API char* cl_protein_json_opts(const char* text, double pH, int flags);

/* Lennard-Jones fluid in reduced units (4 cells^3 particles, truncated and shifted at `cutoff`), Bussi
 * thermostat unless thermostat is "none" or "langevin". Returns {"ok","particles","box","energy_per_particle",
 * "energy_error","pressure","pressure_error","temperature","drift","diagnostics"}. */
CHEMLARP_API char* cl_lj_json(int cells, double density, double temperature, double cutoff, double dt,
                              unsigned equilibration, unsigned production, const char* thermostat, unsigned seed);

/* MBAR on a matrix of reduced potentials (units of kT). Text: the sample count of each of the K states on the
 * first line, then K lines of N reduced potentials u_k(x_n), the samples grouped by the state they came from;
 * '#' lines are skipped. Returns {"ok","converged","iterations","f":[K],"df":[[K]],"ddf":[[K]],"overlap":[[K]],
 * "min_neighbour_overlap","diagnostics"}; df[i][j] = f_j - f_i. */
CHEMLARP_API char* cl_mbar_json(const char* text);

/* Virtual lab (include/chemlarp/lab.hpp). Mixing: "portion <label> <mL> [<T C>]" lines, each followed by
 * "<reagent id or formula> <mol/L>" lines. Returns {"ok","volume_mL","T_initial_C","T_final_C","heat_released_kJ",
 * "pH","ionic_strength","precipitates":[{phase,formula,mol,grams}],"gases":[{phase,mol,mL}],"major_species",
 * "not_modelled","hazards":[{kind,severity,message,source,substances}],"diagnostics"}. */
CHEMLARP_API char* cl_lab_mix_json(const char* text);
/* Titration: "analyte <mL> [<T C>]", solute lines, "titrant <reagent> <mol/L>", "max <mL>", "points <n>".
 * Returns {"ok","curve":[[mL,pH]],"equivalence":[{volume_mL,pH}],"not_modelled","hazards","diagnostics"}. */
CHEMLARP_API char* cl_lab_titration_json(const char* text);
/* The reagent catalogue with reactive groups and GHS data: {"ok","reagents":[...],"ghs_compounds"}. */
CHEMLARP_API char* cl_lab_reagents_json(void);
/* GHS classification of a structure from the bundled PubChem table: {"ok","found","name","cid","signal",
 * "pictograms","hazards":[{code,percent,text}],"diagnostics"?}. */
CHEMLARP_API char* cl_ghs_json(const char* smiles);

/* Excess chemical potential of a Lennard-Jones particle in the LJ fluid (reduced units) by alchemical coupling
 * over lambda windows (MBAR, BAR, TI, EXP) and by Widom insertion (insertions > 0). Zero or negative arguments
 * take the defaults. keep_matrix adds "u_kn" and "n_k". Returns {"ok","particles","box","temperature",
 * "density","mu_ex":{"mbar"|"bar"|"ti"|"exp"|"widom":{"kT":{df,error},"epsilon":{df,error}}},
 * "mbar_converged","min_neighbour_overlap","windows":[{lambda,samples,inefficiency,dudl,dudl_error}],
 * "overlap","diagnostics"}. */
CHEMLARP_API char* cl_lj_solvation_json(int cells, double density, double temperature, double cutoff, double dt,
                                        unsigned equilibration, unsigned production, unsigned insertions, unsigned seed,
                                        int keep_matrix);

/* TIP3P water box (solvation.hpp): NPT Langevin dynamics with PME. Zero or negative arguments take the defaults
 * (400 waters, 298.15 K, 1 atm, 20 ps equilibration, 100 ps production; pressure < 0 runs at constant volume).
 * A negative or non-finite argument, or one past what the engine can finish (20000 waters, 10000 K, 1 us), is
 * rejected with an InvalidInput diagnostic rather than silently falling back to the defaults, which would start
 * a simulation lasting minutes.
 * Returns {"ok","waters","box","density","density_error","potential_per_molecule","potential_error",
 * "temperature","barostat_acceptance","g_oo":[[r,g]...],"diagnostics"}. */
CHEMLARP_API char* cl_water_box_json(int waters, double temperature, double pressure_atm, double equilibration_ps,
                                     double production_ps, unsigned seed);

/* Hydration free energy of `solute` (SMILES of a neutral molecule, or "water") in TIP3P water by alchemical
 * decoupling (solvation.hpp), kcal/mol. Defaults as for cl_water_box_json; production_ps is per window.
 * Returns {"ok","solute","solute_atoms","waters","box","dg_hyd":{"mbar"|"bar"|"ti":{dg,error}},"dg_elec",
 * "dg_vdw","min_neighbour_overlap","windows":[{lambda_elec,lambda_vdw,samples,inefficiency,dudl_elec,
 * dudl_vdw}],"diagnostics"}. Minutes to hours of computation. */
CHEMLARP_API char* cl_hydration_json(const char* solute, int waters, double temperature, double equilibration_ps,
                                     double production_ps, unsigned seed);

/* Protein-ligand docking (docking.hpp): a flexible ligand into a rigid receptor, with an empirical score in
 * the style of AutoDock Vina, receptor terms on grid maps, and Monte Carlo with a local L-BFGS.
 * `receptor` is a PDB or mmCIF file. `ligand` is a SMILES string or MOL/SDF text (detected from the
 * content); from SMILES the conformation is generated, from MOL/SDF it is kept. `reference` is optional and
 * is either the name of a hetero residue of the receptor file ("MK1") or MOL/SDF text of a known pose: it
 * places the search box when none is given, clears the hetero groups sharing that site, and (as MOL/SDF) is
 * the pose the RMSD of each result is measured against.
 * The box is the centre and the edge lengths in angstrom; all six zero means "around the reference ligand".
 * Zero means the default for runs (8 independent Monte Carlo runs), steps (200 per run plus 100 per
 * rotatable bond), modes (9 poses), pH (7.4) and seed (42). A negative or non-finite argument, or one past
 * what the engine can finish (a box edge over 200 A, over 256 runs, over 10^6 steps, over 100 poses, a pH
 * outside -2 to 16), is rejected with an InvalidInput diagnostic rather than silently falling back to the
 * default, which could start a search lasting hours.
 * `flags` (sum): 1 drop the waters the structure resolves (they are kept by default, since a resolved water
 * is part of the site), 2 complete missing side chains first, 4 take the donors from the hydrogens already
 * in the file instead of predicting the protonation at the pH.
 * Returns {"ok","smiles","heavy_atoms","rotatable_bonds","receptor_atoms","box":{"centre","size"},"runs",
 * "steps","evaluations","seconds","poses":[{"affinity","intermolecular","internal","rmsd","terms":{...},
 * "coordinates":[[x,y,z]...]}],"reference":{"atoms","smiles","affinity","rmsd"},"receptor":{...},"sdf",
 * "diagnostics"}. Affinities are kcal/mol, best pose first. Seconds to minutes of computation. */
#define CL_DOCK_REMOVE_WATERS 1
#define CL_DOCK_REPAIR 2
#define CL_DOCK_FILE_HYDROGENS 4
CHEMLARP_API char* cl_dock_json(const char* receptor, const char* ligand, const char* reference, double centre_x,
                                double centre_y, double centre_z, double size_x, double size_y, double size_z,
                                int runs, int steps, int modes, double pH, unsigned seed, int flags);

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

/* Threads used by the parallel kernels (quantum chemistry integrals, J/K build and DFT grid). Results do
 * not depend on the number. Default: CHEMLARP_THREADS from the environment, else the hardware
 * concurrency; WebAssembly builds always use 1. cl_set_threads(n <= 0) restores the default. */
CHEMLARP_API void cl_set_threads(int n);
CHEMLARP_API int cl_threads(void);

CHEMLARP_API void cl_free_string(char* s);

#ifdef __cplusplus
}
#endif

#endif /* CHEMLARP_C_H */
