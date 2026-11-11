// OmniChem web app: workspace around one open compound, backed by the
// WebAssembly engine running in engine-worker.js.
(function () {
  'use strict';
  const BUILD = '__BUILD_ID__';
  const $ = (id) => document.getElementById(id);
  const app = document.querySelector('.app');

  // ---- engine ------------------------------------------------------------------
  const worker = new Worker('engine-worker.js');
  const pending = new Map();
  let nextId = 1;
  worker.onmessage = (e) => {
    const p = pending.get(e.data.id);
    if (!p) return;
    pending.delete(e.data.id);
    e.data.ok ? p.resolve(e.data) : p.reject(new Error(e.data.error));
  };
  worker.onerror = (e) => { $('engine-status').textContent = 'Engine failed to load: ' + (e.message || 'unknown error'); };
  const engine = (op, args) => new Promise((resolve, reject) => {
    const id = nextId++;
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, op, args });
  });
  engine('version').then((r) => {
    const s = $('engine-status');
    s.textContent = 'Engine ' + r.result + ' · offline ready';
    s.classList.add('ready');
  }).catch((err) => { $('engine-status').textContent = 'Engine failed to load: ' + err.message; });

  // ---- element data ----------------------------------------------------------------
  // Standard atomic weights (IUPAC) and monoisotopic masses of the most abundant isotope.
  const AVG = { H: 1.008, B: 10.81, C: 12.011, N: 14.007, O: 15.999, F: 18.998403163, Na: 22.98976928, Mg: 24.305, Al: 26.9815385, Si: 28.085, P: 30.973761998, S: 32.06, Cl: 35.45, K: 39.0983, Ca: 40.078, Fe: 55.845, Cu: 63.546, Zn: 65.38, Se: 78.971, Br: 79.904, I: 126.90447, Ar: 39.948, He: 4.002602, Li: 6.94 };
  const MONO = { H: 1.00782503207, B: 11.0093054, C: 12, N: 14.0030740048, O: 15.99491461956, F: 18.99840322, Na: 22.9897692809, Mg: 23.9850417, Al: 26.98153863, Si: 27.9769265325, P: 30.97376163, S: 31.972071, Cl: 34.96885268, K: 38.96370668, Ca: 39.96259098, Fe: 55.9349375, Cu: 62.9295975, Zn: 63.9291422, Se: 79.9165213, Br: 78.9183371, I: 126.904473, Ar: 39.9623831225, He: 4.00260325, Li: 7.01600455 };
  const ELEMENT_COLOR = (sym) => `var(--el-${['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I'].includes(sym) ? sym : 'other'})`;
  const ELECTRON = 0.00054857990946;

  // ---- DOM helpers -------------------------------------------------------------------
  function h(tag, attrs, ...kids) {
    const n = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v == null || v === false) continue;
      if (k === 'class') n.className = v;
      else if (k === 'text') n.textContent = v;
      else if (k === 'html') n.innerHTML = v; // only for markup this file builds from trusted pieces
      else if (k.startsWith('on')) n.addEventListener(k.slice(2), v);
      else n.setAttribute(k, v);
    }
    for (const c of kids.flat()) if (c != null && c !== false) n.append(c);
    return n;
  }
  const icon = (name) => {
    const s = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    const u = document.createElementNS('http://www.w3.org/2000/svg', 'use');
    u.setAttribute('href', '#i-' + name);
    s.append(u);
    return s;
  };
  const fmt = (v, digits = 4) => {
    if (v == null || !isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(digits - 1).replace('e+', 'e');
    return (+v.toPrecision(digits)).toLocaleString('en-US', { maximumFractionDigits: 10 });
  };
  const fixed = (v, d) => (v == null || !isFinite(v) ? '—' : v.toFixed(d));
  function formulaNode(formula, cls = 'f') {
    // Hill formula "C9H8O4", "H4N+" or "C6H8N+" -> subscripts and a superscript charge.
    const n = h('span', { class: cls });
    const m = /^(.*?)(\d*[+-])?$/.exec(formula || '');
    for (const [, sym, count] of (m[1] || '').matchAll(/([A-Z][a-z]?)(\d*)/g)) {
      n.append(sym);
      if (count) n.append(h('sub', { text: count }));
    }
    if (m[2]) n.append(h('sup', { text: m[2] }));
    return n;
  }
  const metric = (label, value, unit) => h('div', { class: 'metric' }, h('div', { class: 'label', text: label }),
    h('div', { class: 'value' }, value, unit ? h('span', { class: 'unit', text: unit }) : null));
  function busy(out, text) { out.replaceChildren(h('div', { class: 'busy' }, h('div', { class: 'spinner' }), h('span', { text }))); }
  function fail(out, err) { out.replaceChildren(h('div', { class: 'diag' }, h('div', { class: 'code', text: 'ERROR' }), h('div', { text: err.message }))); }
  function diagnostics(list, input) {
    return (list || []).map((d) => {
      const box = h('div', { class: 'diag' + (d.severity === 'warning' ? ' warning' : '') },
        h('div', { class: 'code', text: (d.severity === 'warning' ? 'WARNING · ' : 'ERROR · ') + d.code }),
        h('div', { text: d.message }));
      if (input && d.position != null && d.position <= input.length) {
        box.append(h('pre', {}, input + '\n' + ' '.repeat(d.position), h('span', { class: 'caret', text: '▲' })));
      }
      return box;
    });
  }
  function status(kind, title, sub) {
    return h('div', { class: 'status ' + kind }, icon(kind === 'ok' ? 'check' : 'alert'), h('span', { text: title }), sub ? h('span', { class: 'sub', text: sub }) : null);
  }
  let toastTimer = 0;
  function toast(text) {
    const t = $('toast');
    t.textContent = text;
    t.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => { t.hidden = true; }, 1600);
  }
  async function copy(text) {
    try { await navigator.clipboard.writeText(text); toast('Copied'); }
    catch { toast('Copy not available here'); }
  }
  function download(name, text, type) {
    const url = URL.createObjectURL(new Blob([text], { type }));
    const a = h('a', { href: url, download: name });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  const store = {
    get(k, d) { try { const v = localStorage.getItem('omnichem.' + k); return v ? JSON.parse(v) : d; } catch { return d; } },
    set(k, v) { try { localStorage.setItem('omnichem.' + k, JSON.stringify(v)); } catch { /* storage unavailable */ } },
  };
  function segmented(group, onChange) {
    group.addEventListener('click', (e) => {
      const b = e.target.closest('button[role="radio"]');
      if (!b) return;
      group.querySelectorAll('button[role="radio"]').forEach((x) => x.setAttribute('aria-checked', String(x === b)));
      onChange(b);
    });
  }
  const checked = (group) => group.querySelector('[aria-checked="true"]');

  // ---- the open compound ----------------------------------------------------------------
  let compound = null;  // {input, canonical, formula, composition, mass, charge, atoms, bonds, cip, random}
  let model3d = { for: null, sdf: null, xyz: null };

  function properties(c) {
    const comp = c.composition;
    const syms = Object.keys(comp);
    const known = syms.every((s) => AVG[s] != null && MONO[s] != null);
    const avg = syms.reduce((s, e) => s + (AVG[e] || 0) * comp[e], 0);
    const mono = known ? syms.reduce((s, e) => s + MONO[e] * comp[e], 0) - c.charge * ELECTRON : null;
    const order = ['C', 'H', ...syms.filter((s) => s !== 'C' && s !== 'H').sort()].filter((s) => comp[s]);
    const ea = known ? order.map((s) => ({ sym: s, pct: (100 * AVG[s] * comp[s]) / avg })) : [];
    const n = (s) => comp[s] || 0;
    const halo = n('F') + n('Cl') + n('Br') + n('I');
    const dbe = n('C') + n('Si') - (n('H') + halo) / 2 + (n('N') + n('P')) / 2 + 1;
    return { mono, ea, dbe, mz: mono != null && c.charge ? mono / Math.abs(c.charge) : null };
  }

  function renderInspector() {
    const body = $('inspector-body');
    if (!compound) { body.replaceChildren(h('p', { class: 'muted', text: 'Open a compound to see its properties.' })); return; }
    const c = compound;
    const p = properties(c);
    const stereoCount = c.cip.centres.length + c.cip.double_bonds.length;
    body.replaceChildren(
      h('div', { class: 'thumb', id: 'insp-drawing', role: 'img', 'aria-label': 'Structure drawing' }),
      h('div', {}, formulaNode(c.formula, 'formula')),
      h('dl', { class: 'props' },
        h('dt', { text: 'Molecular Weight' }), h('dd', {}, fixed(c.mass, 4), h('span', { class: 'unit', text: 'g/mol' })),
        h('dt', { text: 'Exact Mass' }), h('dd', {}, p.mono != null ? fixed(p.mono, 6) : '—', h('span', { class: 'unit', text: 'Da' })),
        p.mz != null ? [h('dt', { text: 'm/z' }), h('dd', { text: fixed(p.mz, 6) })] : null,
        h('dt', { text: 'Degree of Unsaturation' }), h('dd', { text: Number.isInteger(p.dbe) ? String(p.dbe) : p.dbe.toFixed(1) }),
        h('dt', { text: 'Heavy Atoms' }), h('dd', { text: String(c.heavy) }),
        h('dt', { text: 'Net Charge' }), h('dd', { text: c.charge > 0 ? '+' + c.charge : String(c.charge) }),
        h('dt', { text: 'Stereo Elements' }), h('dd', { text: String(stereoCount) })),
      h('div', { class: 'inspector-section' },
        h('h3', { text: 'Elemental Analysis' }),
        h('div', { class: 'mono', style: 'font-size:12.5px;line-height:1.6' }, p.ea.length ? p.ea.map((e) => `${e.sym}, ${e.pct.toFixed(2)}`).join('; ') : '—'),
        h('div', { class: 'composition', role: 'img', 'aria-label': 'Mass composition' },
          p.ea.map((e) => h('span', { style: `width:${e.pct}%;background:${ELEMENT_COLOR(e.sym)}`, title: `${e.sym} ${e.pct.toFixed(2)}%` }))),
        h('div', { class: 'legend' }, p.ea.map((e) => h('span', {}, h('i', { style: `background:${ELEMENT_COLOR(e.sym)}` }), `${e.sym} ${e.pct.toFixed(1)}%`)))),
      h('div', { class: 'inspector-section' },
        h('h3', { text: 'Canonical SMILES' }),
        h('div', { class: 'identifier' }, h('span', { class: 'text', text: c.canonical }),
          h('button', { class: 'icon-btn', type: 'button', title: 'Copy', 'aria-label': 'Copy canonical SMILES', onclick: () => copy(c.canonical) }, icon('copy')))),
      vapourPressureSection(c),
      h('div', { class: 'actions' },
        h('button', { class: 'btn', type: 'button', onclick: () => { location.hash = 'model3d'; } , text: '3D Model' }),
        h('button', { class: 'btn', type: 'button', onclick: () => { location.hash = 'structure'; }, text: 'Structure' })),
    );
    paintThumb(c);
  }

  function vapourPressureSection(c) {
    const T = h('input', { type: 'number', value: store.get('psatT', 298.15), step: 'any', min: '1', 'aria-label': 'Temperature in kelvin' });
    const result = h('div', { class: 'mono', style: 'font-size:12.5px;min-height:1.6em' });
    async function run(e) {
      if (e) e.preventDefault();
      store.set('psatT', +T.value);
      result.textContent = '…';
      try {
        const { result: r } = await engine('psat', { smiles: c.canonical, temperature: +T.value, eos: 'pr' });
        if (compound !== c) return;
        const d0 = r.diagnostics[0];
        result.textContent = r.ok ? `${fmt(r.pressure_Pa / 1e5, 5)} bar  (${fmt(r.pressure_Pa / 1e3, 5)} kPa)`
          : d0 && d0.code === 'NoCriticalData' ? 'Not in the critical-property database yet.'
          : d0 && d0.code === 'SupercriticalTemperature' ? 'Above the critical temperature: no vapour pressure.'
          : d0 ? d0.message : 'Not available';
      } catch (err) { result.textContent = err.message; }
    }
    queueMicrotask(run);
    return h('form', { class: 'inspector-section', onsubmit: run },
      h('h3', { text: 'Vapour Pressure (Peng–Robinson)' }),
      h('div', { style: 'display:flex;gap:8px;align-items:center' },
        h('span', { class: 'unit-input', style: 'flex:1' }, T, h('span', { class: 'unit', text: 'K' })),
        h('button', { class: 'btn', type: 'submit', text: 'Compute' })),
      result);
  }

  // ---- 2D drawing ----------------------------------------------------------------------
  const drawOpts = Object.assign({ color: true, cip: false, numbers: false }, store.get('draw', {}));
  const drawFlags = () => (drawOpts.color ? 1 : 0) | (drawOpts.cip ? 2 : 0) | (drawOpts.numbers ? 4 : 0);
  const depictCache = new Map();
  async function depiction(c, flags = drawFlags(), bond = 30) {
    const key = `${c.canonical}|${flags}|${bond}`;
    if (!depictCache.has(key)) {
      if (depictCache.size > 40) depictCache.clear();
      depictCache.set(key, engine('depict', { smiles: c.input, bond, flags }).then((x) => x.result));
    }
    return depictCache.get(key);
  }
  // Standalone SVG for files: black ink instead of currentColor.
  const standaloneSvg = (svg) => '<?xml version="1.0" encoding="UTF-8"?>\n' + svg.replace(/currentColor/g, '#000000');
  function savePng(r, name) {
    const scale = 4;
    const img = new Image();
    const url = URL.createObjectURL(new Blob([standaloneSvg(r.svg)], { type: 'image/svg+xml' }));
    img.onload = () => {
      const cv = document.createElement('canvas');
      cv.width = Math.ceil(r.width * scale);
      cv.height = Math.ceil(r.height * scale);
      const g = cv.getContext('2d');
      g.fillStyle = '#ffffff';
      g.fillRect(0, 0, cv.width, cv.height);
      g.drawImage(img, 0, 0, cv.width, cv.height);
      URL.revokeObjectURL(url);
      cv.toBlob((blob) => {
        const a = h('a', { href: URL.createObjectURL(blob), download: name });
        document.body.append(a);
        a.click();
        a.remove();
      }, 'image/png');
    };
    img.src = url;
  }
  async function paintThumb(c) {
    const host = $('insp-drawing');
    if (!host) return;
    try {
      const r = await depiction(c, drawOpts.color ? 1 : 0);
      if (compound !== c || !r.ok) return;
      host.innerHTML = r.svg;  // engine-generated SVG; labels are element symbols and the title is escaped
    } catch { /* the thumbnail is optional */ }
  }
  function drawingCard(c) {
    const host = h('div', { class: 'drawing', role: 'img', 'aria-label': `Structure of ${c.formula}` }, h('div', { class: 'spinner' }));
    const note = h('span', { class: 'meta' });
    let current = null;
    async function paint() {
      try {
        const r = await depiction(c, drawFlags(), 36);
        if (compound !== c) return;
        current = r;
        if (!r.ok) { host.replaceChildren(...diagnostics(r.diagnostics)); return; }
        host.innerHTML = r.svg;
        const st = r.stereo;
        note.textContent = st.centres + st.bonds ? `Stereo shown: ${st.centres_ok}/${st.centres} centres · ${st.bonds_ok}/${st.bonds} double bonds` : 'No stereo elements';
      } catch (err) { host.replaceChildren(h('p', { class: 'muted', text: err.message })); }
    }
    const toggle = (label, key) => h('button', {
      class: 'chip', type: 'button', 'aria-pressed': String(drawOpts[key]), text: label,
      onclick: (e) => {
        drawOpts[key] = !drawOpts[key];
        e.currentTarget.setAttribute('aria-pressed', String(drawOpts[key]));
        store.set('draw', drawOpts);
        paint();
        if (key === 'color') paintThumb(c);
      },
    });
    const base = (c.formula || 'structure').replace(/[^A-Za-z0-9+-]/g, '');
    const exp = (label, fn) => h('button', { class: 'chip', type: 'button', text: label, onclick: () => current && current.ok && fn(current) });
    paint();
    return h('div', { class: 'glass-card card' },
      h('div', { class: 'card-head' }, h('h2', { text: 'Structure' }), note),
      host,
      h('div', { class: 'drawing-bar' },
        h('div', { class: 'samples', role: 'group', 'aria-label': 'Drawing options' },
          toggle('Colour', 'color'), toggle('Stereo labels', 'cip'), toggle('Atom numbers', 'numbers')),
        h('div', { class: 'samples', role: 'group', 'aria-label': 'Export' },
          h('span', { class: 'samples-label', text: 'Export' }),
          exp('SVG', (r) => download(base + '.svg', standaloneSvg(r.svg), 'image/svg+xml')),
          exp('PNG', (r) => savePng(r, base + '.png')),
          exp('MOL', (r) => download(base + '.mol', r.molfile, 'chemical/x-mdl-molfile')),
          exp('Copy SVG', (r) => copy(standaloneSvg(r.svg))))));
  }

  function renderStructure(r, input) {
    const out = $('structure-out');
    if (!r.ok) {
      out.replaceChildren(status('bad', 'This structure cannot be read', `${r.diagnostics.filter((d) => d.severity === 'error').length} problem(s)`), ...diagnostics(r.diagnostics, input));
      return;
    }
    const c = compound;
    const cip = c.cip;
    const stereo = [
      ...cip.centres.map((x) => h('span', { class: 'stereo' }, `Atom ${x.atom + 1}`, h('b', { text: x.label || '?' }))),
      ...cip.double_bonds.map((x) => h('span', { class: 'stereo' }, `Bond ${x.bond + 1}`, h('b', { text: x.label || '?' }))),
    ];
    out.replaceChildren(
      status('ok', 'Valid structure', `${c.heavy} heavy atoms · ${c.bonds} bonds`),
      drawingCard(c),
      h('div', { class: 'glass-card card' },
        h('div', { class: 'card-head' }, h('h2', { text: 'Identifiers' })),
        h('dl', { class: 'props' },
          h('dt', { text: 'Formula' }), h('dd', {}, formulaNode(c.formula)),
          h('dt', { text: 'Canonical SMILES' }), h('dd', {}, h('span', { class: 'identifier' }, h('span', { class: 'text', text: c.canonical }),
            h('button', { class: 'icon-btn', type: 'button', 'aria-label': 'Copy canonical SMILES', onclick: () => copy(c.canonical) }, icon('copy')))),
          c.inchi ? [
            h('dt', { text: 'InChI' }), h('dd', {}, h('span', { class: 'identifier' }, h('span', { class: 'text', text: c.inchi }),
              h('button', { class: 'icon-btn', type: 'button', 'aria-label': 'Copy InChI', onclick: () => copy(c.inchi) }, icon('copy')))),
            h('dt', { text: 'InChIKey' }), h('dd', {}, h('span', { class: 'identifier' }, h('span', { class: 'text', text: c.inchiKey }),
              h('button', { class: 'icon-btn', type: 'button', 'aria-label': 'Copy InChIKey', onclick: () => copy(c.inchiKey) }, icon('copy')))),
          ] : null,
          h('dt', { text: 'As entered' }), h('dd', { class: 'mono', text: input }))),
      h('div', { class: 'glass-card card' },
        h('div', { class: 'card-head' }, h('h2', { text: 'Stereochemistry' }), h('span', { class: 'meta', text: 'Cahn–Ingold–Prelog' })),
        stereo.length ? h('div', { class: 'stereo-list' }, stereo) : h('p', { class: 'muted', style: 'margin:0', text: 'No stereocentres or stereo double bonds specified.' })),
      c.random.length ? h('div', { class: 'glass-card card' },
        h('div', { class: 'card-head' }, h('h2', { text: 'Equivalent SMILES' }), h('span', { class: 'meta', text: 'All give the same canonical form' })),
        h('div', { class: 'alt-list' }, c.random.map((s) => h('div', { text: s })))) : null,
      ...diagnostics(r.diagnostics, input),
    );
  }

  async function openCompound(input, { remember = true } = {}) {
    input = (input || '').trim();
    if (!input) return;
    $('compound').value = input;
    busy($('structure-out'), 'Reading structure…');
    try {
      const [p, cn, ik] = await Promise.all([engine('parse', { smiles: input }), engine('canonical', { smiles: input, random: 3, seed: 7 }),
        engine('inchi', { smiles: input }).catch(() => ({ result: { ok: false } }))]);
      const r = p.result;
      if (!r.ok) {
        compound = null;
        renderInspector();
        renderStructure(r, input);
        if (current() !== 'structure') toast('Structure has errors: see Structure');
        return;
      }
      const s = r.structure.species[0];
      const atoms = r.structure.atoms;
      compound = {
        input, canonical: s.canonical_smiles, formula: s.formula, composition: s.composition, mass: s.molar_mass,
        charge: s.net_charge, heavy: atoms.element.filter((z) => z !== 1).length, bonds: s.bond_count,
        cip: cn.result.ok ? cn.result.cip : { centres: [], double_bonds: [] }, random: cn.result.ok ? cn.result.random : [],
        inchi: ik.result.ok ? ik.result.inchi : null, inchiKey: ik.result.ok ? ik.result.key : null,
      };
      renderInspector();
      renderStructure(r, input);
      if (remember) remember_(compound);
      store.set('last', input);
      if (current() === 'model3d') build3d();
    } catch (err) { fail($('structure-out'), err); }
  }

  // ---- structure files (MOL/SDF, PDB, CML) ---------------------------------------------
  async function openFile(file) {
    if (!file) return;
    if (file.size > 25e6) { toast('File is larger than 25 MB'); return; }
    location.hash = 'structure';
    const out = $('structure-out');
    busy(out, `Reading ${file.name}…`);
    try {
      const text = await file.text();
      const r = (await engine('read', { text, format: file.name, max: 500 })).result;
      const good = r.ok ? r.records.filter((x) => x.ok) : [];
      const warnings = r.diagnostics.filter((d) => d.severity === 'warning').length;
      if (good.length === 1 && r.records.length === 1) {
        await openCompound(good[0].smiles);
        if (warnings) toast(`${file.name}: ${warnings} warning(s) — see the file's notes below`);
        if (warnings) out.append(...diagnostics(r.diagnostics));
        return;
      }
      if (!good.length) {
        out.replaceChildren(status('bad', `Nothing could be read from ${file.name}`, r.format === 'auto' ? 'Not a MOL, SDF, PDB or CML file' : `${r.format} file`),
          ...diagnostics(r.diagnostics));
        return;
      }
      const failed = r.records.length - good.length;
      out.replaceChildren(
        status(failed ? 'warn' : 'ok', `${r.records.length} records in ${file.name}`, failed ? `${failed} could not be read` : 'Choose one to open'),
        h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'Records' }), h('span', { class: 'meta', text: r.records.length >= 500 ? 'First 500 shown' : r.format.toUpperCase() })),
          h('div', { class: 'record-list' }, r.records.map((x, i) => x.ok
            ? h('button', { type: 'button', class: 'record', onclick: () => openCompound(x.smiles), title: x.smiles },
              h('span', { class: 'r-title', text: x.title || `Record ${i + 1}` }), h('span', { class: 'r-formula' }, formulaNode(x.formula)),
              h('span', { class: 'r-smiles', text: x.smiles }))
            : h('div', { class: 'record bad' }, h('span', { class: 'r-title', text: x.title || `Record ${i + 1}` }), h('span', { class: 'muted', text: `line ${x.line}: could not be read` }))))),
        ...diagnostics(r.diagnostics));
    } catch (err) { fail(out, err); }
  }

  // ---- recents --------------------------------------------------------------------------
  function remember_(c) {
    const list = store.get('recents', []).filter((x) => x.canonical !== c.canonical);
    list.unshift({ input: c.input, canonical: c.canonical, formula: c.formula });
    store.set('recents', list.slice(0, 8));
    renderRecents();
  }
  function renderRecents() {
    const list = store.get('recents', []);
    const ul = $('recent-list');
    if (!list.length) { ul.replaceChildren(h('li', { class: 'empty', text: 'Compounds you open appear here.' })); return; }
    ul.replaceChildren(...list.map((x) => h('li', {}, h('button', { type: 'button', onclick: () => openCompound(x.input), title: x.canonical },
      h('span', { class: 'r-formula' }, formulaNode(x.formula)), h('span', { class: 'r-smiles', text: x.canonical })))));
  }

  // ---- views ------------------------------------------------------------------------------
  const views = ['structure', 'model3d', 'reaction', 'equilibrium', 'phases', 'kinetics', 'roadmap'];
  const autorun = { reaction: 'rx-form', equilibrium: 'eq-form', phases: 'ph-form', kinetics: 'kin-form' };
  const ran = new Set();
  const current = () => (views.includes(location.hash.slice(1)) ? location.hash.slice(1) : 'structure');
  function show() {
    const name = current();
    for (const v of views) $(v).hidden = v !== name;
    document.querySelectorAll('[data-tool]').forEach((a) => (a.dataset.tool === name ? a.setAttribute('aria-current', 'page') : a.removeAttribute('aria-current')));
    $('content').scrollTop = 0;
    app.classList.remove('show-inspector');
    if (name === 'model3d') build3d();
    if (name === 'roadmap') loadRoadmap();
    if (autorun[name] && !ran.has(name)) { ran.add(name); $(autorun[name]).requestSubmit(); }
  }
  window.addEventListener('hashchange', show);

  // ---- 3D ---------------------------------------------------------------------------------
  let viewer = null;
  async function build3d() {
    if (!viewer) viewer = new MoleculeViewer($('viewer'));
    requestAnimationFrame(() => viewer.draw());
    if (!compound) { $('d3-out').replaceChildren(h('p', { class: 'muted', text: 'Open a compound to build its 3D model.' })); return; }
    if (model3d.for === compound.canonical) return;
    const target = compound.canonical;
    model3d = { for: target, sdf: null, xyz: null };
    $('viewer-empty').hidden = false;
    $('d3-sdf').disabled = true;
    busy($('d3-out'), 'Embedding and minimising…');
    try {
      const { result: r, ms } = await engine('embed', { smiles: compound.input, seed: 42, candidates: 5, minimize: true });
      if (!compound || compound.canonical !== target) return;
      $('viewer-empty').hidden = true;
      if (!r.ok) { $('d3-out').replaceChildren(status('bad', 'No 3D model for this compound'), ...diagnostics(r.diagnostics)); return; }
      model3d = { for: target, sdf: r.sdf, xyz: r.xyz };
      viewer.labels = checked($('model3d').querySelector('.segmented')).dataset.style === 'labels';
      viewer.load(r.sdf);
      $('d3-sdf').disabled = false;
      const st = r.stereo;
      const stereoOk = st.centres_ok === st.centres && st.bonds_ok === st.bonds;
      $('d3-out').replaceChildren(
        st.centres + st.bonds === 0 ? status('ok', '3D model ready', 'No stereo elements to reproduce')
          : status(stereoOk ? 'ok' : 'warn', stereoOk ? 'Stereochemistry reproduced in 3D' : 'Stereochemistry only partly reproduced',
            `${st.centres_ok}/${st.centres} centres · ${st.bonds_ok}/${st.bonds} double bonds`),
        h('div', { class: 'metrics' },
          metric('UFF energy', fmt(r.energy_kcal_mol, 5), 'kcal/mol'),
          metric('Minimisation', r.converged ? 'Converged' : 'Stopped'),
          metric('Conformer attempts', String(r.attempts)),
          metric('Compute time', String(Math.round(ms)), 'ms')),
        h('div', { class: 'samples' }, h('span', { class: 'samples-label', text: 'Export' }),
          h('button', { class: 'chip', type: 'button', onclick: () => download(fileBase() + '.sdf', r.sdf, 'chemical/x-mdl-sdfile'), text: 'SDF (MDL molfile)' }),
          h('button', { class: 'chip', type: 'button', onclick: () => download(fileBase() + '.xyz', r.xyz, 'chemical/x-xyz'), text: 'XYZ' })),
        ...diagnostics(r.diagnostics));
    } catch (err) { fail($('d3-out'), err); }
  }
  const fileBase = () => (compound ? compound.formula : 'structure');
  segmented($('model3d').querySelector('.segmented'), (b) => { if (viewer) { viewer.labels = b.dataset.style === 'labels'; viewer.draw(); } });
  $('d3-reset').addEventListener('click', () => { if (viewer && model3d.sdf) viewer.load(model3d.sdf); });
  $('d3-sdf').addEventListener('click', () => model3d.sdf && download(fileBase() + '.sdf', model3d.sdf, 'chemical/x-mdl-sdfile'));

  // ---- Reaction ---------------------------------------------------------------------------
  $('rx-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const equation = $('rx-eq').value.trim();
    const T = +$('rx-t').value;
    const out = $('rx-out');
    busy(out, 'Checking stoichiometry…');
    try {
      const { result: r } = await engine('reaction', { equation, temperature: T });
      const kids = [r.balanced ? status('ok', 'Balanced', 'Every element and charge is conserved') : status('bad', 'Not balanced')];
      if (r.thermo) {
        const t = r.thermo;
        const dG = (t.dG_RT * 8.314462618 * t.temperature_K) / 1000;
        kids.push(h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'Standard reaction thermodynamics' }),
            h('span', { class: 'meta', text: `${fmt(t.temperature_K, 6)} K (${fixed(t.temperature_K - 273.15, 2)} °C) · 1 bar` })),
          h('div', { class: 'metrics' },
            metric('ΔrH°', fmt(t.dH_J / 1000, 5), 'kJ/mol'),
            metric('ΔrS°', fmt(t.dS_J, 5), 'J/(mol·K)'),
            metric('ΔrG°', fmt(dG, 5), 'kJ/mol'),
            metric('K', fmt(t.K, 4))),
          h('p', { class: 'muted', style: 'margin:0', text: dG < 0 ? 'Spontaneous under standard conditions (ΔrG° < 0).' : 'Not spontaneous under standard conditions (ΔrG° > 0).' })));
      }
      out.replaceChildren(...kids, ...diagnostics(r.diagnostics, equation));
    } catch (err) { fail(out, err); }
  });

  // ---- Equilibrium -------------------------------------------------------------------------
  function barRow(label, sub, values, frac) {
    return h('tr', {}, h('td', {}, label), h('td', { class: 'mono' }, sub), ...values.map((v) => h('td', { class: 'num', text: v })),
      h('td', { class: 'num' }, h('div', { class: 'bar-cell' }, h('i', { style: `width:${Math.max(1, frac * 120)}px` }))));
  }
  $('eq-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('eq-out');
    busy(out, 'Minimising Gibbs energy…');
    try {
      const { result: r } = await engine('equilibrate', { species: $('eq-species').value, temperature: +$('eq-t').value, pressure: +$('eq-p').value * 1e5 });
      if (!r.equilibrium) { out.replaceChildren(status('bad', 'Could not solve'), ...diagnostics(r.diagnostics)); return; }
      const q = r.equilibrium;
      const max = Math.max(...q.species.map((s) => s.mole_fraction));
      out.replaceChildren(
        q.converged ? status('ok', 'Equilibrium found', `${q.iterations} Newton iterations`) : status('bad', 'Did not converge'),
        h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'Equilibrium composition' }), h('span', { class: 'meta', text: `${fmt(q.total_moles, 5)} mol total` })),
          h('div', { class: 'table-wrap' }, h('table', {},
            h('thead', {}, h('tr', {}, h('th', { text: 'Formula' }), h('th', { text: 'SMILES' }), h('th', { class: 'num', text: 'Amount (mol)' }), h('th', { class: 'num', text: 'Mole fraction' }), h('th', {}))),
            h('tbody', {}, q.species.map((s) => barRow(formulaNode(s.formula), s.smiles, [fmt(s.moles), fmt(s.mole_fraction)], s.mole_fraction / max))))))
        , ...diagnostics(r.diagnostics));
    } catch (err) { fail(out, err); }
  });

  // ---- Phase behaviour ---------------------------------------------------------------------
  const regionText = { 'two-phase': 'Two phases: vapour + liquid', liquid: 'Single phase: liquid', vapor: 'Single phase: vapour', supercritical: 'Single phase: supercritical fluid' };
  segmented($('ph-eos'), () => $('ph-form').requestSubmit());
  $('ph-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('ph-out');
    busy(out, 'Running stability test and flash…');
    try {
      const eos = checked($('ph-eos')).dataset.eos;
      const { result: r } = await engine('flash', { mixture: $('ph-mix').value, temperature: +$('ph-t').value, pressure: +$('ph-p').value * 1e5, eos });
      if (!r.ok) { out.replaceChildren(status('bad', 'Could not flash this mixture'), ...diagnostics(r.diagnostics)); return; }
      const two = r.region === 'two-phase';
      out.replaceChildren(
        status(two ? 'warn' : 'ok', regionText[r.region] || r.region, `${eos === 'pr' ? 'Peng–Robinson' : 'SRK'}${two ? ` · ${r.iterations} iterations` : ''}`),
        h('div', { class: 'metrics' },
          metric('Vapour fraction', two ? fixed(r.vapor_fraction, 4) : r.region === 'vapor' ? '1' : '0'),
          metric('Z (liquid)', two || r.region === 'liquid' ? fixed(r.liquid.Z, 4) : '—'),
          metric('Z (vapour)', two || r.region !== 'liquid' ? fixed(r.vapor.Z, 4) : '—')),
        h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'Phase compositions' }), h('span', { class: 'meta', text: 'mole fractions' })),
          h('div', { class: 'table-wrap' }, h('table', {},
            h('thead', {}, h('tr', {}, h('th', { text: 'Component' }), h('th', { text: 'SMILES' }), h('th', { class: 'num', text: 'Feed z' }), h('th', { class: 'num', text: 'Liquid x' }), h('th', { class: 'num', text: 'Vapour y' }), h('th', { class: 'num', text: 'K = y/x' }))),
            h('tbody', {}, r.components.map((c) => h('tr', {}, h('td', { text: c.name }), h('td', { class: 'mono', text: c.smiles }),
              h('td', { class: 'num', text: fixed(c.z, 4) }), h('td', { class: 'num', text: two ? fixed(c.x, 4) : '—' }),
              h('td', { class: 'num', text: two ? fixed(c.y, 4) : '—' }), h('td', { class: 'num', text: two ? fmt(c.K, 4) : '—' }))))))),
        ...diagnostics(r.diagnostics));
    } catch (err) { fail(out, err); }
  });

  // ---- Batch reactor --------------------------------------------------------------------------
  $('kin-mech').addEventListener('change', () => { $('kin-custom').hidden = $('kin-mech').value !== 'custom'; });
  const readFile = (input) => (input.files[0] ? input.files[0].text() : Promise.resolve(''));
  let lastKin = null;
  async function mechanismText() {
    const v = $('kin-mech').value;
    if (v !== 'custom') return { mechanism: await (await fetch(v)).text(), thermo: '' };
    if (!$('kin-file').files[0]) throw new Error('Choose a mechanism file first.');
    return { mechanism: await readFile($('kin-file')), thermo: await readFile($('kin-thermo')) };
  }
  function drawKinetics() {
    if (!lastKin) return;
    const { samples, species } = lastKin;
    const tms = samples.map((s) => s.t * 1000);
    Charts.lineChart($('kin-chart-t'), [{ name: 'Temperature (K)', points: samples.map((s, i) => [tms[i], s.T]) }], { xLabel: 'time (ms)', yLabel: 'K' });
    const totals = samples.map((s) => s.n.reduce((a, b) => a + b, 0));
    const shown = species.map((name, i) => ({ name, i, max: Math.max(...samples.map((s, k) => s.n[i] / totals[k])) })).filter((x) => x.max > 1e-6);
    Charts.lineChart($('kin-chart-x'), shown.map((x) => ({ name: x.name, points: samples.map((s, k) => [tms[k], s.n[x.i] / totals[k]]) })),
      { xLabel: 'time (ms)', yLabel: 'mole fraction', logY: checked($('kin-scale')).dataset.log === '1' });
  }
  segmented($('kin-scale'), drawKinetics);
  $('kin-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('kin-out');
    busy(out, 'Integrating stiff kinetics…');
    try {
      const { mechanism, thermo } = await mechanismText();
      const { result: r, ms } = await engine('kinetics', {
        mechanism, thermo, model: $('kin-model').value, temperature: +$('kin-t').value, pressure: +$('kin-p').value * 1e5,
        feed: $('kin-feed').value, tEnd: +$('kin-time').value / 1000, samples: 220,
      });
      if (!r.ok) { $('kin-charts').hidden = true; out.replaceChildren(status('bad', 'The run did not complete'), ...diagnostics(r.diagnostics)); return; }
      const end = r.samples[r.samples.length - 1];
      out.replaceChildren(
        r.ignition_time ? status('warn', 'Ignition', `at ${fmt(r.ignition_time * 1000, 4)} ms`) : status('ok', 'Run complete', 'no ignition in this window'),
        h('div', { class: 'metrics' },
          metric('Final temperature', fmt(end.T, 5), 'K'),
          metric('Final pressure', fmt(end.P / 1e5, 4), 'bar'),
          metric('Solver steps', String(r.stats.steps)),
          metric('Compute time', String(Math.round(ms)), 'ms')),
        ...diagnostics(r.diagnostics));
      lastKin = { samples: r.samples, species: r.species };
      $('kin-charts').hidden = false;
      drawKinetics();
    } catch (err) { fail(out, err); }
  });
  window.addEventListener('resize', () => { if (lastKin && current() === 'kinetics') drawKinetics(); });

  // ---- Roadmap ---------------------------------------------------------------------------------
  let roadmap = null;
  async function fetchRoadmap() {
    if (!roadmap) roadmap = await (await fetch('data/roadmap.json')).json();
    return roadmap;
  }
  function roadmapStats(road) {
    const w = road.weights;
    const all = road.areas.flatMap((a) => a.modules.map((m) => ({ ...m, area: a.name, replaces: a.replaces, weight: w[m.size] })));
    const sum = (ms, st) => ms.filter((m) => !st || m.status === st).reduce((s, m) => s + m.weight, 0);
    return { all, sum, total: sum(all), done: sum(all, 'done'), prog: sum(all, 'in_progress') };
  }
  async function loadRoadmap() {
    const out = $('roadmap-out');
    busy(out, 'Loading roadmap…');
    try {
      const road = await fetchRoadmap();
      const { all, sum, total, done, prog } = roadmapStats(road);
      const bar = (d, p, t) => h('div', { class: 'composition' },
        h('span', { style: `width:${(100 * d) / t}%;background:var(--ok)` }), h('span', { style: `width:${(100 * p) / t}%;background:var(--warn)` }));
      const order = ['in_progress', 'done', 'planned'];
      const label = { done: 'Done', in_progress: 'In progress', planned: 'Planned' };
      out.replaceChildren(
        h('div', { class: 'metrics' },
          metric('Complete', ((100 * done) / total).toFixed(1), '%'),
          metric('In progress', ((100 * prog) / total).toFixed(1), '%'),
          metric('Modules done', `${all.filter((m) => m.status === 'done').length} / ${all.length}`)),
        h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'By product area' }), h('span', { class: 'meta', text: `Updated ${road.updated}` })),
          road.areas.map((a) => {
            const ms = all.filter((m) => m.area === a.name);
            const t = sum(ms), d = sum(ms, 'done'), p = sum(ms, 'in_progress');
            return h('div', { style: 'display:grid;gap:6px' },
              h('div', { class: 'card-head' }, h('span', { style: 'font-weight:600', text: a.name }), h('span', { class: 'meta', text: `${a.replaces} · ${Math.round((100 * d) / t)}%` })),
              bar(d, p, t));
          }),
          h('div', { class: 'legend' }, h('span', {}, h('i', { style: 'background:var(--ok)' }), 'Done'), h('span', {}, h('i', { style: 'background:var(--warn)' }), 'In progress'), h('span', {}, h('i', { style: 'background:var(--hairline)' }), 'Planned'))),
        h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'Modules' })),
          h('div', { class: 'table-wrap' }, h('table', {},
            h('thead', {}, h('tr', {}, h('th', { text: 'Module' }), h('th', { text: 'Area' }), h('th', { text: 'Size' }), h('th', { text: 'Status' }))),
            h('tbody', {}, [...all].sort((a, b) => order.indexOf(a.status) - order.indexOf(b.status)).map((m) =>
              h('tr', {}, h('td', { text: m.name }), h('td', { text: m.area }), h('td', { text: m.size }),
                h('td', {}, h('span', { class: 'stereo', style: `color:var(--${m.status === 'done' ? 'ok' : m.status === 'in_progress' ? 'warn' : 'ink-3'})`, text: label[m.status] }))))))))
      );
    } catch (err) { fail(out, err); }
  }
  fetchRoadmap().then((road) => { const s = roadmapStats(road); $('roadmap-pct').textContent = ((100 * s.done) / s.total).toFixed(0) + '%'; }).catch(() => {});

  // ---- wiring -------------------------------------------------------------------------------
  $('compound-form').addEventListener('submit', (e) => { e.preventDefault(); openCompound($('compound').value); });
  $('open-file').addEventListener('click', () => $('file-input').click());
  $('file-input').addEventListener('change', (e) => { openFile(e.target.files[0]); e.target.value = ''; });
  document.addEventListener('dragover', (e) => {
    if ([...(e.dataTransfer?.types || [])].includes('Files')) { e.preventDefault(); document.body.classList.add('dropping'); }
  });
  document.addEventListener('dragleave', (e) => { if (!e.relatedTarget) document.body.classList.remove('dropping'); });
  document.addEventListener('drop', (e) => {
    if (!e.dataTransfer?.files?.length) return;
    e.preventDefault();
    document.body.classList.remove('dropping');
    openFile(e.dataTransfer.files[0]);
  });
  document.addEventListener('click', (e) => {
    const chip = e.target.closest('[data-compound]');
    if (chip) { openCompound(chip.dataset.compound); return; }
    const fill = e.target.closest('[data-fill]');
    if (fill) { const f = $(fill.dataset.fill); f.value = fill.dataset.value; f.form.requestSubmit(); return; }
    const add = e.target.closest('[data-add-compound]');
    if (add) {
      if (!compound) { toast('Open a compound first'); return; }
      const ta = $(add.dataset.addCompound);
      ta.value = ta.value.replace(/\s*$/, '') + `\n${compound.canonical} ${add.dataset.amount}`;
      toast(`Added ${compound.formula}`);
    }
  });
  $('copy-smiles').addEventListener('click', () => (compound ? copy(compound.canonical) : toast('Open a compound first')));
  function toggleInspector() {
    const wide = window.matchMedia('(min-width: 1181px)').matches;
    if (wide) app.classList.toggle('no-inspector');
    else app.classList.toggle('show-inspector');
    const visible = wide ? !app.classList.contains('no-inspector') : app.classList.contains('show-inspector');
    $('toggle-inspector').setAttribute('aria-pressed', String(visible));
  }
  $('toggle-inspector').addEventListener('click', toggleInspector);
  document.addEventListener('keydown', (e) => {
    const mod = e.metaKey || e.ctrlKey;
    if (!mod) return;
    if (e.key === 'k') { e.preventDefault(); $('compound').focus(); $('compound').select(); }
    else if (e.key === 'i') { e.preventDefault(); toggleInspector(); }
    else if (/^[1-6]$/.test(e.key)) {
      const a = document.querySelector(`.tool-list a[data-key="${e.key}"]`);
      if (a) { e.preventDefault(); location.hash = a.dataset.tool; }
    }
  });

  // ---- install + offline ----------------------------------------------------------------------
  let installEvent = null;
  window.addEventListener('beforeinstallprompt', (e) => { e.preventDefault(); installEvent = e; $('install-btn').hidden = false; });
  $('install-btn').addEventListener('click', async () => {
    if (!installEvent) return;
    installEvent.prompt();
    await installEvent.userChoice;
    installEvent = null;
    $('install-btn').hidden = true;
  });
  const local = ['localhost', '127.0.0.1'].includes(location.hostname);
  if ('serviceWorker' in navigator && location.protocol !== 'file:' && !local) navigator.serviceWorker.register('sw.js').catch(() => {});
  console.info('OmniChem web build', BUILD);

  renderRecents();
  show();
  openCompound(store.get('last', 'CC(=O)Oc1ccccc1C(=O)O'), { remember: true });
})();
