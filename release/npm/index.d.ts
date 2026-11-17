// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Every call returns the JSON object of the C API (include/chemlarp/chemlarp_c.h): {ok, ..., diagnostics}.
export interface Diagnostic { code: string; severity: 'error' | 'warning'; message: string; position?: number | null }
export interface Result { ok: boolean; diagnostics: Diagnostic[]; [key: string]: unknown }
export interface ChemLarp {
  version(): string;
  parse(smiles: string): Result;
  canonical(smiles: string, random?: number, seed?: number): Result;
  inchi(smiles: string): Result;
  iupacName(smiles: string): Result;
  depict(smiles: string, bondPx?: number, flags?: number): Result;
  embed(smiles: string, options?: { seed?: number; candidates?: number; minimize?: boolean }): Result;
  readMolecules(text: string, options?: { format?: string; keepHydrogens?: boolean; max?: number }): Result;
  smartsMatch(smiles: string, smarts: string, options?: { chiral?: boolean; max?: number }): Result;
  similarity(a: string, b: string, options?: { radius?: number; bits?: number }): Result;
  checkReaction(equation: string, temperature?: number): Result;
  equilibrate(speciesLines: string, temperature: number, pressure: number): Result;
  qc(input: string, options?: { method?: string; basis?: string; multiplicity?: number; grid?: number; memoryMb?: number }): Result;
  mmff(smiles: string, options?: { variant?: string; minimize?: boolean }): Result;
  pureComponent(smiles: string, temperature?: number, pressure?: number): Result;
  mbar(uKn: number[][], nK: number[]): Result;
}
export function load(): Promise<ChemLarp>;
