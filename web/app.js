// ChemLarp web app: workspace around one open compound, backed by the
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
    // Settings saved before the rename used the 'omnichem.' prefix; read them as a fallback.
    get(k, d) { try { const v = localStorage.getItem('chemlarp.' + k) ?? localStorage.getItem('omnichem.' + k); return v ? JSON.parse(v) : d; } catch { return d; } },
    set(k, v) { try { localStorage.setItem('chemlarp.' + k, JSON.stringify(v)); } catch { /* storage unavailable */ } },
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
          h('dt', { text: 'IUPAC name' }), c.iupac
            ? h('dd', {}, h('span', { class: 'identifier' }, h('span', { class: 'text', text: c.iupac }),
              h('button', { class: 'icon-btn', type: 'button', 'aria-label': 'Copy IUPAC name', onclick: () => copy(c.iupac) }, icon('copy'))))
            : h('dd', { class: 'muted', text: c.iupacWhy || 'Not available' }),
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

  // Short reason shown when no IUPAC name is produced (the namer reports, it never guesses).
  function namingReason(diags) {
    const d = (diags || []).find((x) => x.code === 'NamingUnsupported');
    if (!d) return 'Not available';
    const m = d.message.replace(/^no IUPAC name: /, '');
    return `Not named: ${m.length > 90 ? m.slice(0, 87) + '…' : m}`;
  }

  async function openCompound(input, { remember = true } = {}) {
    input = (input || '').trim();
    if (!input) return;
    $('compound').value = input;
    busy($('structure-out'), 'Reading structure…');
    try {
      const [p, cn, ik, nm] = await Promise.all([engine('parse', { smiles: input }), engine('canonical', { smiles: input, random: 3, seed: 7 }),
        engine('inchi', { smiles: input }).catch(() => ({ result: { ok: false } })),
        engine('iupac', { smiles: input }).catch(() => ({ result: { ok: false, diagnostics: [] } }))]);
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
        iupac: nm.result.ok ? nm.result.name : null,
        iupacWhy: nm.result.ok ? null : namingReason(nm.result.diagnostics),
      };
      renderInspector();
      renderStructure(r, input);
      if (remember) remember_(compound);
      store.set('last', input);
      if (current() === 'model3d') build3d();
      if (current() === 'draw' && editor && !editor.atoms.length) edLoadCompound();
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
      lastFile = good.length ? { name: file.name, records: good } : lastFile;
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

  // ---- search: SMARTS, similarity, screening ----------------------------------------------
  let lastFile = null;  // records of the last structure file opened

  // Draws the open compound with the atoms (and bonds between them) of `atoms` highlighted.
  function highlightDrawing(c, matches) {
    const host = h('div', { class: 'drawing', role: 'img', 'aria-label': `Structure of ${c.formula} with matches highlighted` }, h('div', { class: 'spinner' }));
    let current = null;
    const paint = (sel) => {
      if (!current) return;
      const svg = host.querySelector('svg');
      if (!svg) return;
      svg.querySelector('g.hl')?.remove();
      const ns = 'http://www.w3.org/2000/svg';
      const g = document.createElementNS(ns, 'g');
      g.setAttribute('class', 'hl');
      const on = new Set(sel.flat());
      const pos = new Map(current.atoms.map((a, i) => [a.atom, { x: a.x, y: a.y, i }]));
      const idx = new Map(current.atoms.map((a, i) => [i, a.atom]));
      for (const b of current.bonds) {
        const A = idx.get(b.a), B = idx.get(b.b);
        if (!on.has(A) || !on.has(B)) continue;
        const pa = pos.get(A), pb = pos.get(B);
        const l = document.createElementNS(ns, 'line');
        for (const [k, v] of Object.entries({ x1: pa.x, y1: pa.y, x2: pb.x, y2: pb.y })) l.setAttribute(k, v);
        g.append(l);
      }
      for (const a of on) {
        const p = pos.get(a);
        if (!p) continue;
        const ci = document.createElementNS(ns, 'circle');
        ci.setAttribute('cx', p.x); ci.setAttribute('cy', p.y); ci.setAttribute('r', 11);
        g.append(ci);
      }
      svg.insertBefore(g, svg.firstChild);
    };
    depiction(c, drawFlags() & ~2).then((r) => {
      if (!r.ok) { host.replaceChildren(h('p', { class: 'muted', text: 'No drawing for this structure.' })); return; }
      current = r;
      host.innerHTML = r.svg;  // engine-generated SVG
      paint(matches);
    }).catch(() => host.replaceChildren());
    return { host, paint };
  }

  async function runSearch() {
    const out = $('sr-out');
    if (!compound) { out.replaceChildren(status('warn', 'Open a compound first', 'Search works on the open compound')); return; }
    const smarts = $('sr-smarts').value.trim();
    const like = $('sr-like').value.trim();
    const chiral = checked($('sr-chiral')).dataset.chiral === '1';
    busy(out, 'Searching…');
    try {
      const cards = [];
      if (smarts) {
        const r = (await engine('smarts', { smiles: compound.input, smarts, chiral, max: 500 })).result;
        if (!r.ok) {
          cards.push(status('bad', 'This SMARTS cannot be read'), ...diagnostics(r.diagnostics, smarts));
        } else {
          const d = highlightDrawing(compound, r.matches);
          const list = r.matches.slice(0, 60).map((m, i) => h('button', {
            type: 'button', class: 'chip', title: 'Highlight this match only',
            onclick: (e) => {
              const btn = e.currentTarget;
              const was = btn.getAttribute('aria-pressed') === 'true';
              btn.parentNode.querySelectorAll('.chip').forEach((x) => x.setAttribute('aria-pressed', 'false'));
              btn.setAttribute('aria-pressed', String(!was));
              d.paint(was ? r.matches : [m]);
            },
            text: `${i + 1}: ${m.map((a) => a + 1).join(' ')}`,
          }));
          cards.push(
            status(r.count ? 'ok' : 'warn', r.count ? `${r.count} match${r.count === 1 ? '' : 'es'} in ${compound.formula}` : `No match in ${compound.formula}`,
              r.count ? 'Matched atoms are highlighted; numbers are atom positions' : smarts),
            h('div', { class: 'glass-card card' },
              h('div', { class: 'card-head' }, h('h2', { text: 'Substructure' }), h('span', { class: 'meta mono', text: smarts })),
              d.host,
              r.count ? h('div', { class: 'samples match-list', role: 'group', 'aria-label': 'Matches' }, list) : null,
              r.count > 60 ? h('p', { class: 'muted', text: `First 60 of ${r.count} shown.` }) : null));
        }
      }
      if (like) {
        const r = (await engine('similarity', { a: compound.input, b: like })).result;
        if (!r.ok) cards.push(status('bad', 'The comparison SMILES cannot be read'), ...diagnostics(r.diagnostics, like));
        else cards.push(h('div', { class: 'glass-card card' },
          h('div', { class: 'card-head' }, h('h2', { text: 'Similarity' }), h('span', { class: 'meta', text: 'Morgan radius 2 (ECFP4), 2048 bits' })),
          h('div', { class: 'metrics' }, metric('Tanimoto', fmt(r.tanimoto, 3)), metric('Dice', fmt(r.dice, 3)),
            metric('Bits set', `${r.bits_a} / ${r.bits_b}`))));
      }
      if (lastFile && lastFile.records.length > 1) {
        const query = like || compound.input;
        const r = (await engine('search', { lines: lastFile.records.map((x) => x.smiles).join('\n'), smarts, query })).result;
        if (r.ok) {
          const rows = r.results.map((x) => ({ ...x, rec: lastFile.records[x.line - 1] }))
            .filter((x) => x.ok && (!smarts || x.count > 0))
            .sort((a, b) => b.similarity - a.similarity);
          cards.push(h('div', { class: 'glass-card card' },
            h('div', { class: 'card-head' }, h('h2', { text: `Screen ${lastFile.name}` }),
              h('span', { class: 'meta', text: `${rows.length} of ${lastFile.records.length} records${smarts ? ' match' : ''} · ranked by similarity to ${like ? 'the comparison SMILES' : 'the open compound'}` })),
            rows.length ? h('div', { class: 'record-list' }, rows.slice(0, 100).map((x) =>
              h('button', { type: 'button', class: 'record', onclick: () => { openCompound(x.rec.smiles); location.hash = 'structure'; }, title: x.rec.smiles },
                h('span', { class: 'r-title', text: x.rec.title || `Record ${x.line}` }),
                h('span', { class: 'r-formula', text: `${fmt(x.similarity, 3)}${smarts ? ` · ${x.count}×` : ''}` }),
                h('span', { class: 'r-smiles', text: x.rec.smiles }))))
              : h('p', { class: 'muted', text: 'No record contains this substructure.' })));
        }
      } else {
        cards.push(h('p', { class: 'muted', text: 'Open an SD file (folder button or drag and drop) to screen its records here.' }));
      }
      out.replaceChildren(...cards);
    } catch (err) { fail(out, err); }
  }

  // ---- draw: structure editor -------------------------------------------------------------
  let editor = null;
  let edSmiles = null;  // canonical SMILES of the drawing, when valid
  let edTimer = null;
  function edSide(nodes) { $('ed-side').replaceChildren(...nodes); }
  function edHelp() {
    return h('details', { class: 'ed-help' }, h('summary', { text: 'Shortcuts' }),
      h('ul', {}, [
        'Click empty space or an atom to add a bond; drag to aim it.',
        'Click a bond to raise its order; with a wedge tool, click again to flip it.',
        'Hover an atom and type C N O S P F I H, L for Cl, B for Br, or + and −.',
        'Hover a bond and press 1, 2 or 3; Delete removes what is under the pointer.',
        'Scroll to pan, pinch or ⌘-scroll to zoom; ⌘Z undoes.',
      ].map((t) => h('li', { text: t }))));
  }
  async function edCheck() {
    if (!editor) return;
    if (!editor.atoms.length) {
      edSmiles = null;
      edSide([status('warn', 'Empty canvas', 'Draw a structure, or load the open compound'), edActions(), edHelp()]);
      return;
    }
    try {
      const r = (await engine('read', { text: editor.toMolfile(), format: 'mol' })).result;
      const rec = r.records && r.records[0];
      const nums = [];
      for (const d of r.diagnostics) for (const m of d.message.matchAll(/atom (\d+)/g)) nums.push(Number(m[1]));
      editor.setErrors(nums);
      if (rec && rec.ok) {
        edSmiles = rec.smiles;
        edSide([
          status('ok', 'Valid structure', `${rec.formula} · ${fmt(rec.molar_mass, 2)} g/mol`),
          h('div', { class: 'field' }, h('label', { text: 'Canonical SMILES' }), h('div', { class: 'ed-smiles mono', text: rec.smiles })),
          edActions(),
          ...diagnostics(r.diagnostics.filter((d) => d.severity === 'warning')),
          edHelp(),
        ]);
      } else {
        edSmiles = null;
        edSide([status('bad', 'Not a valid structure yet', 'Problem atoms are marked in red'), edActions(), ...diagnostics(r.diagnostics), edHelp()]);
      }
    } catch (err) { edSide([status('bad', 'Could not check the drawing', err.message), edActions()]); }
  }
  function edActions() {
    return h('div', { class: 'ed-side-actions' },
      h('button', { type: 'button', class: 'btn btn-primary', disabled: !edSmiles, onclick: () => { openCompound(edSmiles); location.hash = 'structure'; }, text: 'Open in ChemLarp' }),
      h('button', { type: 'button', class: 'btn', disabled: !compound, onclick: edLoadCompound, text: 'Load open compound' }),
      h('button', { type: 'button', class: 'btn', disabled: !edSmiles, onclick: () => copy(edSmiles), text: 'Copy SMILES' }),
      h('button', { type: 'button', class: 'btn', disabled: !editor || !editor.atoms.length, onclick: () => download('structure.mol', editor.toMolfile(), 'chemical/x-mdl-molfile'), text: 'Molfile' }));
  }
  async function edLoadCompound() {
    if (!compound) { toast('Open a compound first'); return; }
    try {
      const r = await depiction(compound, 1);
      if (r.ok) editor.loadMolfile(r.molfile);
    } catch (err) { toast(err.message); }
  }
  async function edClean() {
    if (!edSmiles) { toast('Clean needs a valid structure'); return; }
    const r = (await engine('depict', { smiles: edSmiles, bond: 30, flags: 1 })).result;
    if (r.ok) editor.loadMolfile(r.molfile);
  }
  function setupEditor() {
    if (editor) { requestAnimationFrame(() => editor.render()); return; }
    editor = new StructureEditor($('ed-host'), { onChange: () => { clearTimeout(edTimer); edTimer = setTimeout(edCheck, 220); } });
    const bar = document.querySelector('.editor-toolbar');
    // ring button icons
    bar.querySelectorAll('polygon.ed-ring-icon').forEach((pg) => {
      const n = Number(pg.dataset.n), r = 8.5;
      pg.setAttribute('points', Array.from({ length: n }, (_, k) => {
        const a = -Math.PI / 2 + (2 * Math.PI * k) / n;
        return `${12 + r * Math.cos(a)},${12.5 + r * Math.sin(a)}`;
      }).join(' '));
    });
    bar.querySelectorAll('[data-aromatic] svg').forEach((sv) => {
      const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
      c.setAttribute('cx', 12); c.setAttribute('cy', 12.5); c.setAttribute('r', 4.6);
      sv.append(c);
    });
    const pick = (b) => {
      bar.querySelectorAll('[aria-pressed]').forEach((x) => x.setAttribute('aria-pressed', 'false'));
      b.setAttribute('aria-pressed', 'true');
    };
    bar.querySelectorAll('.ed-btn').forEach((b) => { if (!['undo', 'redo', 'clean', 'clear', 'other'].includes(b.dataset.ed)) b.setAttribute('aria-pressed', b.getAttribute('aria-pressed') || 'false'); });
    bar.addEventListener('click', (e) => {
      const b = e.target.closest('.ed-btn');
      if (!b) return;
      const k = b.dataset.ed;
      if (k === 'undo') return editor.undo();
      if (k === 'redo') return editor.redo();
      if (k === 'clean') return edClean();
      if (k === 'clear') return editor.clear();
      if (k === 'other') {
        const el = (prompt('Element symbol (e.g. Si, Se, Na, B):') || '').trim();
        if (!/^[A-Z][a-z]?$/.test(el)) return;
        editor.element = el;
        editor.tool = { kind: 'element' };
        pick(b);
        return;
      }
      pick(b);
      if (k === 'bond') editor.tool = { kind: 'bond', order: Number(b.dataset.order || 1), stereo: Number(b.dataset.stereo || 0) };
      else if (k === 'chain') editor.tool = { kind: 'chain' };
      else if (k === 'ring') editor.tool = { kind: 'ring', size: Number(b.dataset.size), aromatic: b.hasAttribute('data-aromatic') };
      else if (k === 'element') { editor.element = b.dataset.el; editor.tool = { kind: 'element' }; }
      else if (k === 'charge') editor.tool = { kind: 'charge', delta: Number(b.dataset.delta) };
      else if (k === 'erase') editor.tool = { kind: 'erase' };
      else if (k === 'select') editor.tool = { kind: 'select' };
      editor.svg.focus({ preventScroll: true });
    });
    if (compound) edLoadCompound(); else edCheck();
  }

  // ---- quantum chemistry ---------------------------------------------------------------------
  const HARTREE_EV = 27.211386245988;
  function orbitalDiagram(r) {
    const spins = r.orbitals_beta && r.orbitals_beta.length ? [['α', r.orbitals_alpha, r.n_alpha], ['β', r.orbitals_beta, r.n_beta]] : [['', r.orbitals_alpha, r.n_alpha]];
    const pick = (orb, nocc) => orb.map((e, i) => ({ e: e * HARTREE_EV, occ: i < nocc, i })).filter((o) => o.i >= nocc - 6 && o.i < nocc + 4);
    const cols = spins.map(([label, orb, nocc]) => ({ label, levels: pick(orb, nocc) }));
    const all = cols.flatMap((c) => c.levels.map((l) => l.e));
    const lo = Math.min(...all), hi = Math.max(...all), span = Math.max(1e-6, hi - lo);
    const W = 300, H = 260, top = 14, bottom = 16, left = 46;
    const y = (e) => top + (H - top - bottom) * (1 - (e - lo) / span);
    const colW = (W - left - 10) / cols.length;
    const parts = [`<svg viewBox="0 0 ${W} ${H}" class="orbital-diagram" role="img" aria-label="Orbital energy levels">`];
    for (let k = 0; k <= 4; k++) {
      const e = lo + (span * k) / 4;
      parts.push(`<text x="${left - 6}" y="${y(e) + 4}" text-anchor="end" class="od-axis">${e.toFixed(1)}</text>`);
    }
    parts.push(`<text x="4" y="${top - 2}" class="od-axis">eV</text>`);
    cols.forEach((c, ci) => {
      const x0 = left + ci * colW + 10, x1 = left + (ci + 1) * colW - 10;
      if (c.label) parts.push(`<text x="${(x0 + x1) / 2}" y="${H - 2}" text-anchor="middle" class="od-axis">${c.label}</text>`);
      for (const l of c.levels) parts.push(`<line x1="${x0}" x2="${x1}" y1="${y(l.e)}" y2="${y(l.e)}" class="${l.occ ? 'od-occ' : 'od-virt'}"><title>${l.occ ? 'occupied' : 'virtual'} ${l.e.toFixed(3)} eV</title></line>`);
    });
    parts.push('</svg>');
    return h('div', { class: 'orbital-wrap', html: parts.join('') });  // built only from numbers above
  }
  function chargeDrawing(c, r) {
    // heavy atoms come first in the embedded structure, in SMILES order; add each H to its nearest heavy atom
    const heavy = r.elements.map((e, i) => (e !== 'H' ? i : -1)).filter((i) => i >= 0);
    const q = new Map(heavy.map((i) => [i, r.mulliken[i]]));
    r.elements.forEach((e, i) => {
      if (e !== 'H') return;
      let best = -1, bd = Infinity;
      for (const k of heavy) {
        const d = Math.hypot(...[0, 1, 2].map((t) => r.xyz[i][t] - r.xyz[k][t]));
        if (d < bd) { bd = d; best = k; }
      }
      if (best >= 0) q.set(best, q.get(best) + r.mulliken[i]);
    });
    const host = h('div', { class: 'drawing', role: 'img', 'aria-label': 'Structure with Mulliken charges' }, h('div', { class: 'spinner' }));
    depiction(c, drawFlags() & ~2, 46).then((d) => {
      if (!d.ok) { host.replaceChildren(); return; }
      host.innerHTML = d.svg;  // engine-generated SVG
      const svg = host.querySelector('svg');
      const ns = 'http://www.w3.org/2000/svg';
      const g = document.createElementNS(ns, 'g');
      const cx = d.atoms.reduce((s, a) => s + a.x, 0) / d.atoms.length, cy = d.atoms.reduce((s, a) => s + a.y, 0) / d.atoms.length;
      for (const a of d.atoms) {
        if (!q.has(a.atom)) continue;
        const v = q.get(a.atom);
        const halo = document.createElementNS(ns, 'circle');
        halo.setAttribute('cx', a.x); halo.setAttribute('cy', a.y);
        halo.setAttribute('r', 6 + Math.min(10, Math.abs(v) * 18));
        halo.setAttribute('class', v < 0 ? 'q-neg' : 'q-pos');
        g.append(halo);
        // label outward from the molecule's centre, clear of the atom symbol
        let ux = a.x - cx, uy = a.y - cy;
        const len = Math.hypot(ux, uy);
        if (len < 1e-6) { ux = 0; uy = -1; } else { ux /= len; uy /= len; }
        const t = document.createElementNS(ns, 'text');
        t.setAttribute('x', a.x + ux * 21); t.setAttribute('y', a.y + uy * 21 + 3.5);
        t.setAttribute('text-anchor', 'middle');
        t.setAttribute('class', 'q-label');
        t.textContent = (v > 0 ? '+' : '') + v.toFixed(2);
        g.append(t);
      }
      svg.insertBefore(g, svg.firstChild);
      svg.append(...[...g.querySelectorAll('text')]);
      svg.style.overflow = 'visible';  // labels may sit just outside the drawing's box
      host.style.padding = '18px 12px';
    }).catch(() => host.replaceChildren());
    return host;
  }
  async function runQc() {
    const out = $('qc-out');
    if (!compound) { out.replaceChildren(status('warn', 'Open a compound first', 'Quantum chemistry works on the open compound')); return; }
    const method = checked($('qc-method')).dataset.method, basis = $('qc-basis').value, mult = Number($('qc-mult').value);
    busy(out, `Running ${method.toUpperCase()}/${basis}…`);
    const t0 = performance.now();
    try {
      const r = (await engine('qc', { input: compound.input, method, basis, mult, grid: 2, memory: 384 })).result;
      if (!r.ok) { out.replaceChildren(status('bad', 'The calculation could not run'), ...diagnostics(r.diagnostics)); return; }
      const secs = (performance.now() - t0) / 1000;
      const gap = (r.lumo - r.homo) * HARTREE_EV;
      const metrics = [
        metric('Total energy', fmt(r.energy, 10), 'Eₕ'),
        metric('HOMO', fmt(r.homo * HARTREE_EV, 4), 'eV'),
        metric('LUMO', fmt(r.lumo * HARTREE_EV, 4), 'eV'),
        metric('Gap', fmt(gap, 4), 'eV'),
        metric('Dipole', fmt(r.dipole_total, 4), 'D'),
      ];
      if (r.multiplicity > 1) metrics.push(metric('⟨S²⟩', fmt(r.s2, 4)));
      const comps = [['Nuclear repulsion', r.nuclear_repulsion], ['One-electron', r.one_electron], ['Two-electron', r.two_electron]];
      if (method !== 'hf') comps.push(['Exchange–correlation', r.xc]);
      out.replaceChildren(
        status(r.converged ? 'ok' : 'warn', r.converged ? `${r.method}/${r.basis} converged` : `${r.method}/${r.basis} did not converge`,
          `${r.iterations} iterations · ${r.n_basis} basis functions · ${secs.toFixed(1)} s`),
        h('div', { class: 'glass-card card' }, h('div', { class: 'card-head' }, h('h2', { text: 'Energy' }),
          h('span', { class: 'meta', text: `charge ${r.charge} · multiplicity ${r.multiplicity}` })), h('div', { class: 'metrics' }, metrics)),
        h('div', { class: 'qc-grid' },
          h('div', { class: 'glass-card card' }, h('div', { class: 'card-head' }, h('h2', { text: 'Mulliken charges' }),
            h('span', { class: 'meta', text: 'hydrogens summed into their atoms' })), chargeDrawing(compound, r)),
          h('div', { class: 'glass-card card' }, h('div', { class: 'card-head' }, h('h2', { text: 'Orbital energies' }),
            h('span', { class: 'meta', text: 'filled: occupied' })), orbitalDiagram(r))),
        h('div', { class: 'glass-card card' }, h('div', { class: 'card-head' }, h('h2', { text: 'Energy components' }), h('span', { class: 'meta', text: 'hartree' })),
          h('table', { class: 'data' }, h('tbody', {}, comps.map(([k, v]) => h('tr', {}, h('td', { text: k }), h('td', { class: 'num', text: fmt(v, 8) })))))),
        ...diagnostics(r.diagnostics.filter((d) => d.severity === 'warning')));
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
  const views = ['structure', 'model3d', 'draw', 'search', 'reaction', 'equilibrium', 'phases', 'quantum', 'kinetics', 'roadmap'];
  const autorun = { search: 'sr-form', reaction: 'rx-form', equilibrium: 'eq-form', phases: 'ph-form', kinetics: 'kin-form' };
  const ran = new Set();
  const current = () => (views.includes(location.hash.slice(1)) ? location.hash.slice(1) : 'structure');
  function show() {
    const name = current();
    for (const v of views) $(v).hidden = v !== name;
    document.querySelectorAll('[data-tool]').forEach((a) => (a.dataset.tool === name ? a.setAttribute('aria-current', 'page') : a.removeAttribute('aria-current')));
    $('content').scrollTop = 0;
    app.classList.remove('show-inspector');
    if (name === 'model3d') build3d();
    if (name === 'draw') setupEditor();
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
  const phModel = () => checked($('ph-model')).dataset.model;
  const phKind = () => checked($('ph-kind')).dataset.kind;
  const phGiven = () => checked($('ph-given')).dataset.given;
  // Each model keeps its own inputs, so switching back and forth loses nothing.
  const phSaved = { unifac: { mix: 'CCO 0.3\nO 0.7', T: '350', P: '1.01325' } };
  let phShown = 'eos', lastVle = null;
  function phLayout() {
    const m = phModel(), act = m === 'unifac', bubble = phKind() === 'bubble';
    if (m !== phShown) {
      phSaved[phShown] = { mix: $('ph-mix').value, T: $('ph-t').value, P: $('ph-p').value };
      const v = phSaved[m];
      if (v) { $('ph-mix').value = v.mix; $('ph-t').value = v.T; $('ph-p').value = v.P; }
      phShown = m;
      $('ph-out').replaceChildren();
      $('ph-charts').hidden = true;
      lastVle = null;
    }
    $('ph-act-fields').hidden = !act; $('ph-act-samples').hidden = !act;
    $('ph-eos-field').hidden = act; $('ph-eos-samples').hidden = act;
    // The variable being solved for is an output, not an input.
    $('ph-t').closest('.field').hidden = act && phGiven() === 'P';
    $('ph-p').closest('.field').hidden = act && phGiven() === 'T';
    $('ph-submit').textContent = act ? (bubble ? 'Find bubble point' : 'Find dew point') : 'Flash';
    $('ph-mix-label').replaceChildren(act ? (bubble ? 'Liquid composition' : 'Vapour composition') : 'Mixture',
      ', one per line as ', h('code', { text: 'SMILES mole_fraction' }));
  }
  const phRerun = () => { phLayout(); $('ph-form').requestSubmit(); };
  segmented($('ph-model'), phRerun);
  segmented($('ph-kind'), phRerun);
  segmented($('ph-given'), phRerun);
  phLayout();

  function drawVle() {
    if (!lastVle) return;
    const { bubble, dew, tie, az, given, a } = lastVle;
    const markers = az ? [{ x: az.x, y: az.v, label: `azeotrope ${fixed(az.x, 3)}` }] : [];
    Charts.lineChart($('ph-chart'), [
      { name: given === 'P' ? 'Bubble T (liquid)' : 'Bubble P (liquid)', points: bubble },
      { name: given === 'P' ? 'Dew T (vapour)' : 'Dew P (vapour)', points: dew },
      { name: 'This tie line', points: tie, color: getComputedStyle(document.documentElement).getPropertyValue('--ink-3').trim() || '#888' },
    ], { xLabel: `mole fraction ${a.name || a.smiles}`, yLabel: given === 'P' ? 'K' : 'bar', xRange: [0, 1], markers });
  }

  // Binary T-x-y (given P) or P-x-y (given T) diagram from bubble points across the composition range.
  async function phaseDiagram(r, given, value) {
    const [a, b] = r.components;
    const mode = given === 'P' ? 'bubble-T' : 'bubble-P';
    const val = (o) => (given === 'P' ? o.T_K : o.P_Pa / 1e5);
    const xs = Array.from({ length: 51 }, (_, i) => i / 50);
    const pts = (await Promise.all(xs.map((x) =>
      engine('vle', { mixture: `${a.smiles} ${x}\n${b.smiles} ${1 - x}`, mode, value }).then((o) => o.result, () => null))))
      .filter((o) => o && o.ok);
    if (pts.length < 3) { $('ph-charts').hidden = true; return; }
    const bubble = pts.map((o) => [o.components[0].x, val(o)]);
    const dew = pts.map((o) => [o.components[0].y, val(o)]).sort((p, q) => p[0] - q[0]);
    // An azeotrope is where y − x changes sign between two interior points.
    let az = null;
    const inner = pts.filter((o) => o.components[0].x > 0 && o.components[0].x < 1);
    for (let i = 1; i < inner.length && !az; ++i) {
      const p = inner[i - 1], q = inner[i];
      const d1 = p.components[0].y - p.components[0].x, d2 = q.components[0].y - q.components[0].x;
      if (d1 * d2 < 0) {
        const t = d1 / (d1 - d2);
        az = { x: p.components[0].x + t * (q.components[0].x - p.components[0].x), v: val(p) + t * (val(q) - val(p)) };
      }
    }
    const v = val(r);
    lastVle = { bubble, dew, tie: [[r.components[0].x, v], [r.components[0].y, v]], az, given, a };
    $('ph-chart-title').textContent = given === 'P' ? `T–x–y at ${fmt(value / 1e5, 5)} bar` : `P–x–y at ${fmt(value, 5)} K`;
    $('ph-chart-meta').textContent = az
      ? `azeotrope: x = ${fixed(az.x, 3)} ${a.name || a.smiles}, ${given === 'P' ? fixed(az.v, 2) + ' K' : fmt(az.v, 4) + ' bar'}`
      : 'no azeotrope';
    $('ph-charts').hidden = false;
    drawVle();
  }

  async function runVle(out) {
    const kind = phKind(), given = phGiven();
    const value = given === 'P' ? +$('ph-p').value * 1e5 : +$('ph-t').value;
    busy(out, `Finding the ${kind} point…`);
    const { result: r } = await engine('vle', { mixture: $('ph-mix').value, mode: `${kind}-${given === 'P' ? 'T' : 'P'}`, value });
    if (!r.ok) {
      $('ph-charts').hidden = true;
      lastVle = null;
      out.replaceChildren(status('bad', `Could not find the ${kind} point`), ...diagnostics(r.diagnostics));
      return;
    }
    const found = given === 'P' ? `${fixed(r.T_K, 2)} K (${fixed(r.T_K - 273.15, 2)} °C)` : `${fmt(r.P_Pa / 1e5, 5)} bar`;
    out.replaceChildren(
      status(r.extrapolated ? 'warn' : 'ok', `${kind === 'bubble' ? 'Bubble' : 'Dew'} point: ${found}`,
        `UNIFAC + Antoine, ideal vapour${r.extrapolated ? ' · vapour pressure extrapolated' : ''}`),
      h('div', { class: 'metrics' },
        metric('Temperature', fixed(r.T_K, 2), 'K'),
        metric('Pressure', fmt(r.P_Pa / 1e5, 5), 'bar'),
        metric('Iterations', String(r.iterations))),
      h('div', { class: 'glass-card card' },
        h('div', { class: 'card-head' }, h('h2', { text: 'Phase compositions' }), h('span', { class: 'meta', text: 'mole fractions' })),
        h('div', { class: 'table-wrap' }, h('table', {},
          h('thead', {}, h('tr', {}, h('th', { text: 'Component' }), h('th', { text: 'SMILES' }), h('th', { class: 'num', text: 'Liquid x' }),
            h('th', { class: 'num', text: 'Vapour y' }), h('th', { class: 'num', text: 'γ' }), h('th', { class: 'num', text: 'K = y/x' }))),
          h('tbody', {}, r.components.map((c) => h('tr', {}, h('td', { text: c.name }), h('td', { class: 'mono', text: c.smiles }),
            h('td', { class: 'num', text: fixed(c.x, 4) }), h('td', { class: 'num', text: fixed(c.y, 4) }),
            h('td', { class: 'num', text: fixed(c.gamma, 4) }), h('td', { class: 'num', text: c.x > 0 ? fmt(c.y / c.x, 4) : '—' }))))))),
      ...diagnostics(r.diagnostics));
    if (r.components.length === 2) await phaseDiagram(r, given, value);
    else { $('ph-charts').hidden = true; lastVle = null; }
  }

  segmented($('ph-eos'), () => $('ph-form').requestSubmit());
  $('ph-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('ph-out');
    if (phModel() === 'unifac') { try { await runVle(out); } catch (err) { fail(out, err); } return; }
    $('ph-charts').hidden = true;
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
  window.addEventListener('resize', () => { if (lastVle && current() === 'phases') drawVle(); });

  // ---- Roadmap ---------------------------------------------------------------------------------
  // The roadmap changes with every deploy: fetch it fresh (no HTTP cache) whenever it is shown, and when
  // the app comes back into view, so a long-open tab or an installed app never shows an old figure.
  async function fetchRoadmap() {
    return (await fetch('data/roadmap.json', { cache: 'no-store' })).json();
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
  function refreshRoadmapPct() {
    fetchRoadmap().then((road) => {
      const s = roadmapStats(road);
      $('roadmap-pct').textContent = ((100 * s.done) / s.total).toFixed(1) + '%';
    }).catch(() => {});
  }
  refreshRoadmapPct();
  document.addEventListener('visibilitychange', () => {
    if (document.visibilityState !== 'visible') return;
    refreshRoadmapPct();
    if (current() === 'roadmap') loadRoadmap();
  });

  // ---- wiring -------------------------------------------------------------------------------
  $('compound-form').addEventListener('submit', (e) => { e.preventDefault(); openCompound($('compound').value); });
  $('open-file').addEventListener('click', () => $('file-input').click());
  $('sr-form').addEventListener('submit', (e) => { e.preventDefault(); runSearch(); });
  $('qc-form').addEventListener('submit', (e) => { e.preventDefault(); runQc(); });
  segmented($('qc-method'), () => {});
  segmented($('sr-chiral'), () => $('sr-form').requestSubmit());
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
    else if (/^[1-9]$/.test(e.key)) {
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
  if ('serviceWorker' in navigator && location.protocol !== 'file:' && !local) {
    // A new deploy installs a new worker, which takes control at once (skipWaiting + claim); reload then,
    // so the page runs the new build instead of waiting for the next visit. Check for one whenever the
    // app comes back into view.
    const hadController = !!navigator.serviceWorker.controller;
    let reloading = false;
    navigator.serviceWorker.addEventListener('controllerchange', () => {
      if (!hadController || reloading) return;
      reloading = true;
      location.reload();
    });
    navigator.serviceWorker.register('sw.js').then((reg) => {
      document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'visible') reg.update().catch(() => {}); });
    }).catch(() => {});
  }
  console.info('ChemLarp web build', BUILD);

  renderRecents();
  show();
  openCompound(store.get('last', 'CC(=O)Oc1ccccc1C(=O)O'), { remember: true });
})();
