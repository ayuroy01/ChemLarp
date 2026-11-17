// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
'use strict';
// ChemLarp for Node.js: the C++ engine compiled to WebAssembly (no native build, any OS).
//
//   const chemlarp = await require('chemlarp').load();
//   chemlarp.canonical('C1=CC=CC=C1').canonical      // 'c1ccccc1'
//
// Every function returns the parsed JSON the C API returns ({ok, ..., diagnostics}); see
// include/chemlarp/chemlarp_c.h for the fields. The signatures are those of web/engine-worker.js.
const path = require('path');
const createChemLarp = require('./chemlarp.js');

function wrap(M) {
  function call(name, types, args) {
    const ptr = M.ccall(name, 'number', types, args);
    if (!ptr) throw new Error(name + ' returned nothing');
    try {
      return JSON.parse(M.UTF8ToString(ptr));
    } finally {
      M.ccall('cl_free_string', null, ['number'], [ptr]);
    }
  }
  return {
    version: () => M.ccall('cl_version', 'string', [], []),
    parse: (smiles) => call('cl_parse_smiles_json', ['string'], [smiles]),
    canonical: (smiles, random = 0, seed = 1) => call('cl_canonical_smiles32_json', ['string', 'number', 'number'], [smiles, random, seed]),
    inchi: (smiles) => call('cl_inchi_json', ['string'], [smiles]),
    iupacName: (smiles) => call('cl_iupac_name_json', ['string'], [smiles]),
    depict: (smiles, bondPx = 30, flags = 1) => call('cl_depict_json', ['string', 'number', 'number'], [smiles, bondPx, flags]),
    // Adds `elements` and `coordinates` ([x, y, z] in Å per atom) from the XYZ text, as the Python binding does.
    embed: (smiles, { seed = 42, candidates = 5, minimize = true } = {}) => {
      const r = call('cl_embed32_json', ['string', 'number', 'number', 'number'], [smiles, seed, candidates, minimize ? 1 : 0]);
      if (r.ok && typeof r.xyz === 'string') {
        const rows = r.xyz.split('\n').slice(2).filter((l) => l.trim());
        r.elements = rows.map((l) => l.trim().split(/\s+/)[0]);
        r.coordinates = rows.map((l) => l.trim().split(/\s+/).slice(1, 4).map(Number));
      }
      return r;
    },
    readMolecules: (text, { format = '', keepHydrogens = false, max = 0 } = {}) =>
      call('cl_read_molecules_json', ['string', 'string', 'number', 'number'], [text, format, keepHydrogens ? 1 : 0, max]),
    smartsMatch: (smiles, smarts, { chiral = false, max = 1000 } = {}) =>
      call('cl_smarts_match_json', ['string', 'string', 'number', 'number'], [smiles, smarts, chiral ? 1 : 0, max]),
    similarity: (a, b, { radius = 0, bits = 0 } = {}) => call('cl_similarity_json', ['string', 'string', 'number', 'number'], [a, b, radius, bits]),
    checkReaction: (equation, temperature = 298.15) => call('cl_check_reaction_json', ['string', 'number'], [equation, temperature]),
    equilibrate: (speciesLines, temperature, pressure) => call('cl_equilibrate_text_json', ['string', 'number', 'number'], [speciesLines, temperature, pressure]),
    qc: (input, { method = 'hf', basis = 'STO-3G', multiplicity = 0, grid = 2, memoryMb = 384 } = {}) =>
      call('cl_qc_json', ['string', 'string', 'string', 'number', 'number', 'number'], [input, method, basis, multiplicity, grid, memoryMb]),
    mmff: (smiles, { variant = 'MMFF94', minimize = true } = {}) => call('cl_mmff_json', ['string', 'string', 'number'], [smiles, variant, minimize ? 1 : 0]),
    pureComponent: (smiles, temperature = 0, pressure = 0) => call('cl_pure_component_json', ['string', 'number', 'number'], [smiles, temperature, pressure]),
    mbar: (uKn, nK) => {
      const lines = [nK.join(' ')].concat(uKn.map((row) => row.map((x) => String(x)).join(' ')));
      return call('cl_mbar_json', ['string'], [lines.join('\n')]);
    },
  };
}

let loading = null;
// Loads the WebAssembly module once; resolves to the API object.
function load() {
  if (!loading) loading = createChemLarp({ locateFile: (f) => path.join(__dirname, f) }).then(wrap);
  return loading;
}

module.exports = { load };
