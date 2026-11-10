// Runs the WebAssembly engine off the main thread so long calculations
// (kinetics, 3D embedding) never freeze the interface.
importScripts('omnichem.js');

const ready = createOmniChem({ locateFile: (path) => path });

// Every exported engine call returns a malloc'd JSON string that we must free.
function callJson(M, name, types, args) {
  const ptr = M.ccall(name, 'number', types, args);
  if (!ptr) throw new Error(name + ' returned nothing');
  try {
    return JSON.parse(M.UTF8ToString(ptr));
  } finally {
    M.ccall('oc_free_string', null, ['number'], [ptr]);
  }
}

const calls = {
  version: (M) => M.ccall('oc_version', 'string', [], []),
  parse: (M, a) => callJson(M, 'oc_parse_smiles_json', ['string'], [a.smiles]),
  canonical: (M, a) => callJson(M, 'oc_canonical_smiles32_json', ['string', 'number', 'number'], [a.smiles, a.random || 0, a.seed || 1]),
  reaction: (M, a) => callJson(M, 'oc_check_reaction_json', ['string', 'number'], [a.equation, a.temperature]),
  equilibrate: (M, a) => callJson(M, 'oc_equilibrate_text_json', ['string', 'number', 'number'], [a.species, a.temperature, a.pressure]),
  kinetics: (M, a) => callJson(M, 'oc_kinetics_text_json',
    ['string', 'string', 'string', 'string', 'number', 'number', 'string', 'number', 'number', 'number'],
    [a.mechanism, a.thermo || '', a.smiles || '', a.model, a.temperature, a.pressure, a.feed, a.tEnd, a.samples, a.rtol || 0]),
  flash: (M, a) => callJson(M, 'oc_flash_text_json', ['string', 'number', 'number', 'string'], [a.mixture, a.temperature, a.pressure, a.eos]),
  psat: (M, a) => callJson(M, 'oc_saturation_pressure_json', ['string', 'number', 'string'], [a.smiles, a.temperature, a.eos]),
  depict: (M, a) => callJson(M, 'oc_depict_json', ['string', 'number', 'number'], [a.smiles, a.bond || 30, a.flags == null ? 1 : a.flags]),
  embed: (M, a) => callJson(M, 'oc_embed32_json', ['string', 'number', 'number', 'number'], [a.smiles, a.seed || 42, a.candidates || 5, a.minimize ? 1 : 0]),
};

self.onmessage = async (event) => {
  const { id, op, args } = event.data;
  try {
    const M = await ready;
    if (!calls[op]) throw new Error('unknown operation ' + op);
    const t0 = performance.now();
    const result = calls[op](M, args || {});
    self.postMessage({ id, ok: true, result, ms: performance.now() - t0 });
  } catch (err) {
    self.postMessage({ id, ok: false, error: String(err && err.message ? err.message : err) });
  }
};
