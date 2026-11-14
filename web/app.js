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
    const stereoCount = c.cip.centres.length + c.cip.double_bonds.length + (c.cip.axes || []).length;
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
        const axes = (c.cip && c.cip.axes || []).length;
        const axisNote = axes ? `${axes} chirality axis${axes === 1 ? '' : 'es'} (labelled, not drawn with wedges)` : '';
        note.textContent = st.centres + st.bonds
          ? `Stereo shown: ${st.centres_ok}/${st.centres} centres · ${st.bonds_ok}/${st.bonds} double bonds${axisNote ? ' · ' + axisNote : ''}`
          : axisNote || 'No stereo elements';
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
      ...(cip.axes || []).map((x) => h('span', { class: 'stereo' }, x.centre !== null ? `Axis at atom ${x.centre + 1}` : `Axis ${x.begin + 1}–${x.end + 1}`, h('b', { text: x.label || '?' }))),
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
        cip: cn.result.ok ? cn.result.cip : { centres: [], double_bonds: [], axes: [] }, random: cn.result.ok ? cn.result.random : [],
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
  const views = ['structure', 'model3d', 'dynamics', 'draw', 'search', 'reaction', 'equilibrium', 'phases', 'aqueous', 'protein', 'quantum', 'kinetics', 'flowsheet', 'roadmap'];
  const autorun = { search: 'sr-form', reaction: 'rx-form', equilibrium: 'eq-form', phases: 'ph-form', aqueous: 'aq-form', protein: 'pr-form', kinetics: 'kin-form', flowsheet: 'fs-form' };
  const ran = new Set();
  const current = () => (views.includes(location.hash.slice(1)) ? location.hash.slice(1) : 'structure');
  function show() {
    const name = current();
    for (const v of views) $(v).hidden = v !== name;
    document.querySelectorAll('[data-tool]').forEach((a) => (a.dataset.tool === name ? a.setAttribute('aria-current', 'page') : a.removeAttribute('aria-current')));
    // The phone tab bar scrolls sideways: keep the current tool in view.
    const tab = document.querySelector('.tabbar a[aria-current]');
    if (tab && tab.offsetParent) tab.parentElement.scrollTo({ left: tab.offsetLeft - (tab.parentElement.clientWidth - tab.offsetWidth) / 2 });
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
    const ffSel = checked($('d3-ff')).dataset.ff;
    if (model3d.for === compound.canonical + '|' + ffSel) return;
    const target = compound.canonical + '|' + ffSel;
    model3d = { for: target, sdf: null, xyz: null };
    $('viewer-empty').hidden = false;
    $('d3-sdf').disabled = true;
    busy($('d3-out'), 'Embedding and minimising…');
    if (ffSel !== 'uff') {
      try {
        const { result: r, ms } = await engine('mmff', { smiles: compound.input, variant: ffSel, minimize: true });
        if (!compound || compound.canonical + '|' + checked($('d3-ff')).dataset.ff !== target) return;
        $('viewer-empty').hidden = true;
        if (!r.ok) { $('d3-out').replaceChildren(status('bad', `No ${ffSel} model for this compound`), ...diagnostics(r.diagnostics)); return; }
        model3d = { for: target, sdf: r.sdf, xyz: r.xyz };
        viewer.labels = checked($('model3d').querySelector('.segmented')).dataset.style === 'labels';
        viewer.load(r.sdf);
        $('d3-sdf').disabled = false;
        const E = r.energy;
        const terms = [['Bond stretching', E.bond], ['Angle bending', E.angle], ['Stretch–bend', E.stretch_bend],
          ['Out-of-plane', E.out_of_plane], ['Torsion', E.torsion], ['van der Waals', E.vdw], ['Electrostatic', E.electrostatic]];
        $('d3-out').replaceChildren(
          status(r.converged ? 'ok' : 'warn', `${r.variant} minimum`, `${r.iterations} iterations${r.empirical_parameters ? ` · ${r.empirical_parameters} empirical-rule parameters` : ''}`),
          h('div', { class: 'metrics' },
            metric(`${r.variant} energy`, fmt(E.total, 6), 'kcal/mol'),
            metric('Net charge', fmt(r.charges.reduce((a, b) => a + b, 0), 3), 'e'),
            metric('Compute time', String(Math.round(ms)), 'ms')),
          h('div', { class: 'glass-card card' },
            h('div', { class: 'card-head' }, h('h2', { text: 'Energy terms' }), h('span', { class: 'meta', text: 'kcal/mol' })),
            h('div', { class: 'table-wrap' }, h('table', {}, h('tbody', {}, terms.map(([k, v]) => h('tr', {}, h('td', { text: k }), h('td', { class: 'num', text: fixed(v, 4) }))))))),
          h('div', { class: 'glass-card card' },
            h('div', { class: 'card-head' }, h('h2', { text: 'Atom types and partial charges' }), h('span', { class: 'meta', text: 'MMFF bond charge increments' })),
            h('div', { class: 'table-wrap' }, h('table', {},
              h('thead', {}, h('tr', {}, h('th', { text: 'Atom' }), h('th', { text: 'Element' }), h('th', { class: 'num', text: 'Type' }), h('th', { text: 'Symbol' }), h('th', { class: 'num', text: 'Charge' }))),
              h('tbody', {}, r.elements.map((el, i) => h('tr', {}, h('td', { text: String(i + 1) }), h('td', { text: el }),
                h('td', { class: 'num', text: String(r.types[i]) }), h('td', { class: 'mono', text: r.symbols[i] }),
                h('td', { class: 'num ' + (r.charges[i] > 0.05 ? 'q-pos' : r.charges[i] < -0.05 ? 'q-neg' : ''), text: fixed(r.charges[i], 4) }))))))),
          h('div', { class: 'samples' }, h('span', { class: 'samples-label', text: 'Export' }),
            h('button', { class: 'chip', type: 'button', onclick: () => download(fileBase() + '.sdf', r.sdf, 'chemical/x-mdl-sdfile'), text: 'SDF (MDL molfile)' }),
            h('button', { class: 'chip', type: 'button', onclick: () => download(fileBase() + '.xyz', r.xyz, 'chemical/x-xyz'), text: 'XYZ' })),
          ...diagnostics(r.diagnostics));
      } catch (err) { fail($('d3-out'), err); }
      return;
    }
    try {
      const { result: r, ms } = await engine('embed', { smiles: compound.input, seed: 42, candidates: 5, minimize: true });
      if (!compound || compound.canonical + '|uff' !== target) return;
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
  segmented($('d3-ff'), () => build3d());
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

  // ---- Molecular dynamics -----------------------------------------------------------------------
  let mdViewer = null, mdRun = null, mdTimer = 0;
  function mdShowFrame(i) {
    if (!mdRun) return;
    const k = Math.max(0, Math.min(mdRun.positions.length - 1, i));
    $('md-frame').value = String(k);
    mdViewer.setPositions(mdRun.positions[k]);
    const t = mdRun.frames[Math.min(mdRun.frames.length - 1, k * mdRun.position_stride)].t;
    $('md-tlabel').textContent = `${fmt(t, 4)} fs`;
  }
  function mdStop() { clearInterval(mdTimer); mdTimer = 0; $('md-play').textContent = '▶'; }
  $('md-play').addEventListener('click', () => {
    if (!mdRun) return;
    if (mdTimer) { mdStop(); return; }
    $('md-play').textContent = '❚❚';
    mdTimer = setInterval(() => {
      if (current() !== 'dynamics') { mdStop(); return; }
      const next = (+$('md-frame').value + 1) % mdRun.positions.length;
      mdShowFrame(next);
    }, 40);
  });
  $('md-frame').addEventListener('input', () => { mdStop(); mdShowFrame(+$('md-frame').value); });
  $('md-xyz').addEventListener('click', () => { if (mdRun) download(fileBase() + '-md.xyz', mdRun.xyz, 'chemical/x-xyz'); });
  function drawMd() {
    if (!mdRun) return;
    const ps = (f) => f.t / 1000;
    const fr = mdRun.frames;
    const thermo = mdRun.thermostat !== 'none';
    Charts.lineChart($('md-chart-e'), [
      { name: 'potential', points: fr.map((f) => [ps(f), f.U]) },
      { name: 'kinetic', points: fr.map((f) => [ps(f), f.K]) },
      { name: 'total', points: fr.map((f) => [ps(f), f.E]) },
      ...(thermo ? [{ name: 'conserved', points: fr.map((f) => [ps(f), f.conserved]) }] : []),
    ], { xLabel: 'time (ps)', yLabel: 'kcal/mol' });
    Charts.lineChart($('md-chart-t'), [{ name: 'T (K)', points: fr.map((f) => [ps(f), f.T]) }], { xLabel: 'time (ps)', yLabel: 'K' });
    const sp = mdRun.spectrum.intensity, mx = mdRun.spectrum.max_cm1;
    const pts = sp.map((v, i) => [(mx * i) / (sp.length - 1), v]);
    Charts.lineChart($('md-chart-s'), [{ name: 'spectrum', points: pts }],
      { xLabel: 'wavenumber (cm⁻¹)', yLabel: 'intensity', xRange: [0, mx], markers: (mdRun.modes || []).filter((f) => f > 0 && f < mx).map((f) => ({ x: f, y: 0 })) });
  }
  window.addEventListener('resize', () => { if (mdRun && current() === 'dynamics') drawMd(); });
  $('md-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('md-out');
    mdStop();
    if (!compound) { out.replaceChildren(h('p', { class: 'muted', text: 'Open a compound to run its dynamics.' })); return; }
    const dt = +$('md-dt').value, ps = +$('md-time').value;
    const steps = Math.round((ps * 1000) / dt);
    if (!(steps > 0) || steps > 400000) { fail(out, new Error('Choose a duration and time step giving 1 to 400,000 steps.')); return; }
    const every = Math.max(1, Math.round(steps / 1500));
    const thermostat = $('md-thermo').value;
    busy(out, `Running ${steps.toLocaleString('en-US')} steps…`);
    try {
      const { result: r, ms } = await engine('md', {
        smiles: compound.input, temperature: +$('md-t').value, steps, dt, thermostat, tau: +$('md-tau').value, every, modes: $('md-modes').checked,
        forceField: $('md-ff').value,
      });
      if (!r.ok) { $('md-results').hidden = true; mdRun = null; out.replaceChildren(status('bad', 'The dynamics could not run'), ...diagnostics(r.diagnostics)); return; }
      mdRun = { ...r, thermostat };
      const nve = thermostat === 'none';
      out.replaceChildren(
        status('ok', `${(steps * dt / 1000).toLocaleString('en-US', { maximumFractionDigits: 3 })} ps of dynamics`,
          `${r.elements.length} atoms · ${r.dof} kinetic degrees of freedom · ${nve ? 'NVE' : thermostat === 'bussi' ? 'Bussi NVT' : 'Langevin NVT'} · ${$('md-ff').selectedOptions[0].text.split(' ')[0]}`),
        h('div', { class: 'metrics' },
          metric('Mean temperature', fixed(r.mean_temperature_K, 1), 'K'),
          metric('Mean potential', fmt(r.mean_potential - r.minimized_energy, 4), 'kcal/mol above min.'),
          metric(nve ? 'Energy drift' : 'Conserved-energy drift', fmt(r.drift_kcal_mol_ps, 2), 'kcal/mol/ps'),
          metric('Compute time', String(Math.round(ms)), 'ms')),
        ...diagnostics(r.diagnostics));
      $('md-results').hidden = false;
      if (!mdViewer) mdViewer = new MoleculeViewer($('md-viewer'));
      const p0 = r.positions[0];
      mdViewer.loadAtoms({
        atoms: r.elements.map((el, i) => ({ el, x: p0[3 * i], y: p0[3 * i + 1], z: p0[3 * i + 2] })),
        bonds: r.bonds.map(([a, b, order]) => ({ a, b, order })),
      });
      $('md-frame').max = String(r.positions.length - 1);
      mdShowFrame(0);
      $('md-modes-card').hidden = !r.modes.length;
      $('md-modes-list').replaceChildren(...r.modes.map((f) => h('span', { class: 'chip' + (f < 0 ? ' md-imag' : ''), text: f < 0 ? `${fixed(-f, 1)}i` : fixed(f, 1) })));
      drawMd();
    } catch (err) { fail(out, err); }
  });

  // ---- Protein ----------------------------------------------------------------------------------
  let prViewer = null, prLast = null, prText = null, prName = '1CRN';
  // Covalent bonds by distance on a 2.2 Å grid (PDB files carry no bond list for standard residues).
  function proteinBonds(atoms) {
    const cell = 2.2, grid = new Map(), bonds = [];
    const key = (x, y, z) => `${x},${y},${z}`;
    atoms.forEach((a, i) => {
      const k = key(Math.floor(a.x / cell), Math.floor(a.y / cell), Math.floor(a.z / cell));
      (grid.get(k) || grid.set(k, []).get(k)).push(i);
    });
    atoms.forEach((a, i) => {
      const cx = Math.floor(a.x / cell), cy = Math.floor(a.y / cell), cz = Math.floor(a.z / cell);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        for (const j of grid.get(key(cx + dx, cy + dy, cz + dz)) || []) {
          if (j <= i) continue;
          const b = atoms[j], hyd = a.el === 'H' || b.el === 'H';
          if (hyd && a.el === b.el) continue;
          const max = hyd ? 1.3 : (a.el === 'S' || b.el === 'S') ? 2.2 : 1.9;
          const d = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
          if (d > 0.4 && d < max) bonds.push({ a: i, b: j, order: 1 });
        }
      }
    });
    return bonds;
  }
  // Atoms of the result: PDB columns, or the mmCIF atom_site rows the engine writes (fixed column order).
  function prAtoms() {
    const atoms = [];
    const fix = (el) => (el.length === 2 ? el[0] + el[1].toLowerCase() : el);
    if (prLast.format === 'pdb') {
      for (const l of prLast.pdb.split('\n')) {
        if (!l.startsWith('ATOM') && !l.startsWith('HETATM')) continue;
        atoms.push({ el: fix(l.slice(76, 78).trim()), name: l.slice(12, 16).trim(), res: l.slice(17, 20).trim(), chain: l[21], num: +l.slice(22, 26),
          x: +l.slice(30, 38), y: +l.slice(38, 46), z: +l.slice(46, 54), predicted: +l.slice(54, 60) === 0 });
      }
    } else {
      for (const l of prLast.cif.split('\n')) {
        if (!l.startsWith('ATOM') && !l.startsWith('HETATM')) continue;
        const t = l.match(/"[^"]*"|'[^']*'|\S+/g);
        atoms.push({ el: fix(t[2]), name: t[3].replace(/^["']|["']$/g, ''), res: t[5], chain: t[6], num: +t[8], x: +t[10], y: +t[11], z: +t[12], predicted: false });
      }
    }
    return atoms;
  }
  const PR_TRACE_ABOVE = 12000;  // heavy atoms; larger structures are drawn as a chain trace
  function prShowStructure() {
    if (!prLast) return;
    const withH = checked($('pr-show')).dataset.h === '1';
    const all = prAtoms().filter((a) => a.res !== 'HOH' && a.res !== 'WAT' && a.res !== 'DOD');
    const heavy = all.filter((a) => a.el !== 'H' && a.el !== 'D');
    if (!prViewer) prViewer = new MoleculeViewer($('pr-viewer'));
    const trace = heavy.length > PR_TRACE_ABOVE;
    $('pr-show').hidden = trace;
    if (trace) {
      // One point per residue (CA of proteins, P of nucleic acids), joined along each chain.
      const pts = heavy.filter((a) => a.name === 'CA' || a.name === 'P');
      const bonds = [];
      for (let i = 1; i < pts.length; i++) {
        const a = pts[i - 1], b = pts[i];
        if (a.chain === b.chain && a.name === b.name && Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z) < (a.name === 'P' ? 8 : 4.3)) bonds.push({ a: i - 1, b: i, order: 1 });
      }
      prViewer.loadAtoms({ atoms: pts.map((a) => ({ el: a.name === 'P' ? 'P' : 'C', x: a.x, y: a.y, z: a.z })), bonds });
      return;
    }
    const atoms = (withH ? all : heavy).map((a) => ({ el: a.el, x: a.x, y: a.y, z: a.z }));
    prViewer.loadAtoms({ atoms, bonds: proteinBonds(atoms) });
  }
  segmented($('pr-show'), prShowStructure);
  function drawTitration() {
    if (!prLast) return;
    const pH = +$('pr-ph').value;
    // Charge at the chosen pH, interpolated on the 0.25-step curve from the structure-based pKa values.
    const t = prLast.titration_structure, k = Math.min(t.length - 2, Math.max(0, Math.floor(pH / 0.25)));
    const q = isFinite(pH) ? [pH, t[k][1] + (t[k + 1][1] - t[k][1]) * (pH - t[k][0]) / (t[k + 1][0] - t[k][0])] : null;
    Charts.lineChart($('pr-titration'), [{ name: 'from the structure', points: t }, { name: 'from the sequence', points: prLast.titration }], {
      xLabel: 'pH', yLabel: 'net charge', xRange: [0, 14],
      markers: q ? [{ x: q[0], y: q[1], label: `pH ${fixed(q[0], 2)}` }] : [],
    });
  }
  function seqBlock(seq) {
    // Groups of ten with the residue number of each group's first residue.
    const rows = [];
    for (let i = 0; i < seq.length; i += 50) {
      const groups = [];
      for (let k = i; k < Math.min(seq.length, i + 50); k += 10) groups.push(seq.slice(k, k + 10));
      rows.push(h('div', { class: 'pr-seq-row' }, h('span', { class: 'pr-seq-n', text: String(i + 1) }), h('span', { text: groups.join(' ') })));
    }
    return h('div', { class: 'pr-seq mono' }, ...rows);
  }
  // PDB titles are upper case: sentence case, cut at a word near 70 characters.
  function prTitle(t) {
    let s = t.toLowerCase().replace(/^./, (c) => c.toUpperCase());
    if (s.length > 72) s = s.slice(0, s.lastIndexOf(' ', 70)) + '…';
    return s;
  }
  async function prLoadDefault() {
    if (prText) return;
    try {
      const r = await fetch('data/1crn.pdb');
      if (r.ok) { prText = await r.text(); prName = '1CRN'; $('pr-id').value = '1CRN'; }
    } catch { /* offline without cache */ }
  }
  $('pr-file').addEventListener('change', async () => {
    const f = $('pr-file').files[0];
    if (!f) return;
    prText = await f.text();
    prName = f.name.replace(/\.(pdb|ent|cif|mmcif|txt)$/i, '');
    $('pr-id').value = '';
    $('pr-form').requestSubmit();
  });
  async function prFetch() {
    const id = $('pr-id').value.trim().toUpperCase();
    const out = $('pr-out');
    if (!/^[0-9][A-Z0-9]{3}$/.test(id)) { out.replaceChildren(status('bad', 'A PDB ID has four characters, e.g. 1CRN')); return; }
    busy(out, `Downloading ${id} from RCSB…`);
    try {
      // Very large entries have no PDB-format file: fall back to mmCIF.
      let r = await fetch(`https://files.rcsb.org/download/${id}.pdb`);
      if (r.status === 404) r = await fetch(`https://files.rcsb.org/download/${id}.cif`);
      if (!r.ok) throw new Error(r.status === 404 ? `RCSB has no entry ${id}` : `RCSB answered ${r.status}`);
      prText = await r.text();
      prName = id;
      $('pr-form').requestSubmit();
    } catch (err) { fail(out, err); }
  }
  $('pr-fetch').addEventListener('click', prFetch);
  $('pr-id').addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); prFetch(); } });
  $('pr-save').addEventListener('click', () => prLast && (prLast.format === 'pdb'
    ? download(`${prName}_H.pdb`, prLast.pdb, 'chemical/x-pdb') : download(`${prName}_H.cif`, prLast.cif, 'chemical/x-mmcif')));
  $('pr-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('pr-out');
    await prLoadDefault();
    if (!prText) { out.replaceChildren(status('bad', 'Fetch a PDB ID or open a file first')); return; }
    const pH = +$('pr-ph').value;
    if (!(pH >= 0 && pH <= 14)) { out.replaceChildren(status('bad', 'pH must be between 0 and 14')); return; }
    busy(out, 'Protonating…');
    try {
      const flags = ($('pr-repair').checked ? 0 : 1) | ($('pr-loops').checked ? 0 : 2) | ($('pr-termini').checked ? 4 : 0) | ($('pr-flips').checked ? 0 : 8);
      const { result: r, ms } = await engine('protein', { pdb: prText, pH, flags });
      if (!r.ok) {
        $('pr-results').hidden = true;
        out.replaceChildren(status('bad', 'The structure could not be read'), ...diagnostics(r.diagnostics));
        return;
      }
      prLast = r;
      const nres = r.chains.reduce((s, c) => s + c.sequence.length, 0);
      out.replaceChildren(
        status('ok', prName + (r.title ? ' · ' + prTitle(r.title) : ''),
          `${r.chains.length} chain${r.chains.length === 1 ? '' : 's'} · ${nres} residues · ${r.disulfides.length} disulfide${r.disulfides.length === 1 ? '' : 's'}`),
        h('div', { class: 'metrics' },
          metric(`Net charge at pH ${fixed(pH, 1)}`, (r.net_charge > 0 ? '+' : '') + r.net_charge),
          metric('Isoelectric point', fixed(r.pI, 2), 'from the sequence'),
          metric('Hydrogens added', r.hydrogens_added.toLocaleString('en-US'), r.hydrogens_removed ? `${r.hydrogens_removed.toLocaleString('en-US')} replaced` : ''),
          ...(r.repair.atoms_added ? [metric('Atoms built', r.repair.atoms_added.toLocaleString('en-US'),
            [r.repair.completed_residues.length ? `${r.repair.completed_residues.length} side chains` : '', r.repair.built_residues.length ? `${r.repair.built_residues.length} residues predicted` : ''].filter(Boolean).join(' · '))] : []),
          metric('Compute time', String(Math.round(ms)), 'ms')),
        ...diagnostics(r.diagnostics));
      $('pr-chains').replaceChildren(...r.chains.map((c) => h('div', { class: 'pr-chain' },
        h('div', { class: 'pr-chain-head' }, h('strong', { text: `Chain ${c.id || '–'}` }),
          h('span', { class: 'meta', text: `${c.sequence.length} residues · pI ${fixed(c.pI, 2)} · ${fmt(c.molecular_weight / 1000, 5)} kDa` }),
          h('button', { type: 'button', class: 'icon-btn', title: 'Copy sequence', 'aria-label': `Copy chain ${c.id} sequence`, onclick: () => copy(c.sequence) }, icon('copy'))),
        seqBlock(c.sequence))),
        ...[r.disulfides.length ? h('p', { class: 'meta', text: 'Disulfides: ' + r.disulfides.map(([a, n, b, m]) => `Cys ${a}${n}–Cys ${b}${m}`).join(', ') }) : null,
        r.chain_breaks.length ? h('p', { class: 'meta', text: 'Chain breaks before: ' + r.chain_breaks.join(', ') }) : null].filter(Boolean));
      const charged = (q) => h('td', { class: 'num ' + (q > 0 ? 'pr-pos' : q < 0 ? 'pr-neg' : ''), text: q > 0 ? `+${q}` : String(q) });
      // pKa of each residue's groups (side chain, and the terminus it carries).
      const pk = new Map();
      for (const s of r.pka) {
        const key = `${s.chain}|${s.number}|${s.name}`;
        (pk.get(key) || pk.set(key, []).get(key)).push(s);
      }
      const pkCell = (x) => {
        const list = pk.get(`${x.chain}|${x.number}|${x.name}`) || [];
        const txt = list.map((s) => (list.length > 1 ? (s.group === 'N+' ? 'N-term ' : s.group === 'C-' ? 'C-term ' : 'side chain ') : '') + fixed(s.pKa, 2)).join(' · ');
        const shift = list.reduce((m, s) => (Math.abs(s.pKa - s.reference) > Math.abs(m) ? s.pKa - s.reference : m), 0);
        return [h('td', { class: 'num', text: txt || '–' }), h('td', { class: 'num meta', text: Math.abs(shift) >= 0.05 ? (shift > 0 ? '+' : '−') + fixed(Math.abs(shift), 2) : '' })];
      };
      const listNote = (label, items) => (items.length ? h('p', { class: 'meta', text: `${label}: ${items.slice(0, 40).join(', ')}${items.length > 40 ? ` and ${items.length - 40} more` : ''}` }) : null);
      $('pr-notes').replaceChildren(...[
        listNote('Flipped by 180° for better hydrogen bonds', r.flipped),
        listNote('Not decided by the structure (kept as read; histidine as HIE)', r.undetermined.map((s) => s.split(':')[0])),
        listNote('Built from the sequence, conformation predicted', r.repair.built_residues),
        listNote('Not built', r.repair.not_built),
        r.repair.missing_residues_left ? h('p', { class: 'meta', text: `${r.repair.missing_residues_left} residues listed by the file have no coordinates (chain ends are built only on request).` }) : null,
      ].filter(Boolean));
      const PR_ROWS = 1500;
      $('pr-residues').replaceChildren(h('table', {},
        h('thead', {}, h('tr', {}, ...['Residue', 'Chain', 'pKa', 'Shift', 'State', 'Charge'].map((t, i) => h('th', { class: i === 2 || i === 3 || i === 5 ? 'num' : '', text: t })))),
        h('tbody', {}, ...r.residues.slice(0, PR_ROWS).map((x) => h('tr', {},
          h('td', { text: `${x.name[0]}${x.name.slice(1).toLowerCase()} ${x.number}${x.icode}${x.n_terminal ? ' · N-term' : ''}${x.c_terminal ? ' · C-term' : ''}${x.predicted ? ' · built' : ''}` }),
          h('td', { text: x.chain || '–' }), ...pkCell(x), h('td', { class: 'mono', text: x.state }), charged(x.charge))))),
        ...(r.residues.length > PR_ROWS ? [h('p', { class: 'meta', text: `First ${PR_ROWS.toLocaleString('en-US')} of ${r.residues.length.toLocaleString('en-US')} residues shown; the saved file has all.` })] : []));
      $('pr-results').hidden = false;
      prShowStructure();
      drawTitration();
    } catch (err) { fail(out, err); }
  });
  window.addEventListener('resize', () => { if (prLast && current() === 'protein') drawTitration(); });

  // ---- Aqueous chemistry ------------------------------------------------------------------------
  const AQ_EXAMPLES = [
    { name: 'Hard water · CaCl₂ + NaHCO₃', text: `# 1 mmol/kg CaCl2 and 2 mmol/kg NaHCO3 at 25 °C
add CaCl2 1
add NaHCO3 2
` },
    { name: 'Calcite with atmospheric CO₂', text: `# Limestone water open to air (10^-3.5 atm CO2)
temperature 25
equilibrium Calcite 0
equilibrium CO2(g) -3.5
` },
    { name: 'Rain water', text: `# Pure water in equilibrium with atmospheric CO2
equilibrium CO2(g) -3.5
` },
    { name: 'Gypsum solubility', text: `equilibrium Gypsum 0
` },
    { name: 'Phosphate buffer', text: `add NaH2PO4 10
add Na2HPO4 10
` },
    { name: 'Ammonia solution', text: `# 10 mmol/kg NH3: the charge balance sets the pH
N(-3) 10
` },
    { name: 'Seawater', text: `units mmol/kgw
pH 8.22 charge
Ca 10.6
Mg 54.7
Na 485.8
K 10.6
Cl 565.8
S(6) 29.3
C(4) 2.4
` },
    { name: 'Iron(III) hydrolysis', text: `add FeCl3 0.1
` },
  ];
  $('aq-example').replaceChildren(...AQ_EXAMPLES.map((x, i) => h('option', { value: String(i), text: x.name })));
  const aqSaved = store.get('aqueous', null);
  $('aq-text').value = aqSaved ?? AQ_EXAMPLES[0].text;
  if (aqSaved != null) $('aq-example').append(h('option', { value: 'saved', text: 'Your last solution' })), ($('aq-example').value = 'saved');
  $('aq-example').addEventListener('change', () => {
    const x = AQ_EXAMPLES[+$('aq-example').value];
    if (x) { $('aq-text').value = x.text; $('aq-form').requestSubmit(); }
  });
  $('aq-text').addEventListener('input', () => store.set('aqueous', $('aq-text').value));
  let lastAq = null, aqDist = null;
  function aqSpeciesTable() {
    const all = checked($('aq-all')).dataset.all === '1';
    const rows = lastAq.species.filter((s) => all || s.molality >= 1e-9);
    $('aq-species').replaceChildren(h('table', {},
      h('thead', {}, h('tr', {}, h('th', { text: 'Species' }), h('th', { class: 'num', text: 'Molality' }),
        h('th', { class: 'num', text: 'Activity' }), h('th', { class: 'num', text: 'log γ' }), h('th', { class: 'num', text: 'log a' }))),
      h('tbody', {}, rows.map((s) => h('tr', {}, h('td', { class: 'mono', text: s.name }), h('td', { class: 'num', text: fmt(s.molality, 4) }),
        h('td', { class: 'num', text: fmt(s.activity, 4) }), h('td', { class: 'num', text: fixed(s.log_gamma, 4) }),
        h('td', { class: 'num', text: s.activity > 0 ? fixed(Math.log10(s.activity), 3) : '—' }))))));
  }
  segmented($('aq-all'), () => lastAq && aqSpeciesTable());
  // Species distribution against pH for one component: the solved composition (component totals after any
  // phase reactions) at fixed pH from 0 to 14.
  async function aqDistribution() {
    const comp = $('aq-dist-comp').value;
    if (!comp || !lastAq) return;
    const base = `temperature ${(lastAq.T_K - 273.15).toFixed(4)}\nunits mol/kgw\n` +
      lastAq.totals.map((t) => `${t.component} ${t.molality.toPrecision(12)}`).join('\n');
    const pts = [];
    for (let p = 0; p <= 14.0001; p += 0.25) {
      const { result: r } = await engine('solution', { text: `${base}\npH ${p.toFixed(2)}\n` });
      if (r.ok) pts.push([p, r]);
    }
    const series = new Map();
    for (const [p, r] of pts) {
      const tot = (r.totals.find((t) => t.component === comp) || {}).molality;
      if (!(tot > 0)) continue;
      for (const s of r.species) {
        const c = s.components[comp];
        if (!c) continue;
        if (!series.has(s.name)) series.set(s.name, []);
        series.get(s.name).push([p, (c * s.molality) / tot]);
      }
    }
    aqDist = [...series.entries()].filter(([, v]) => Math.max(...v.map((q) => q[1])) > 0.01).map(([name, points]) => ({ name, points }));
    drawAqDist();
  }
  function drawAqDist() {
    if (!aqDist) return;
    const markers = lastAq && lastAq.pH >= 0 && lastAq.pH <= 14 ? [{ x: lastAq.pH, y: 0, label: `pH ${fixed(lastAq.pH, 2)}` }] : [];
    Charts.lineChart($('aq-dist'), aqDist, { xLabel: 'pH', yLabel: 'fraction of ' + $('aq-dist-comp').value, xRange: [0, 14], markers });
  }
  $('aq-dist-comp').addEventListener('change', () => aqDistribution().catch(() => {}));
  $('aq-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('aq-out');
    busy(out, 'Solving the speciation…');
    try {
      const text = $('aq-text').value;
      const { result: r, ms } = await engine('solution', { text });
      if (!r.ok) {
        $('aq-results').hidden = true;
        lastAq = null;
        out.replaceChildren(status('bad', 'The solution could not be solved'), ...diagnostics(r.diagnostics));
        return;
      }
      lastAq = r;
      const sat = r.saturation.filter((s) => s.si > 0).map((s) => s.phase);
      out.replaceChildren(
        status(r.ionic_strength > 0.7 ? 'warn' : 'ok', `pH ${fixed(r.pH, 3)}`,
          `${fixed(r.T_K - 273.15, 1)} °C · ${sat.length ? 'supersaturated: ' + sat.slice(0, 4).join(', ') + (sat.length > 4 ? '…' : '') : 'no supersaturated phase'}`),
        h('div', { class: 'metrics' },
          metric('pH', fixed(r.pH, 3)),
          metric('Ionic strength', fmt(r.ionic_strength, 4), 'mol/kg'),
          metric('Alkalinity', fmt(r.alkalinity_eq * 1e3, 4), 'meq/kg'),
          metric('Compute time', String(Math.round(ms)), 'ms')),
        ...diagnostics(r.diagnostics));
      $('aq-results').hidden = false;
      aqSpeciesTable();
      $('aq-phases-card').hidden = !r.phases.length;
      $('aq-phases').replaceChildren(h('table', {},
        h('thead', {}, h('tr', {}, h('th', { text: 'Phase' }), h('th', { class: 'num', text: 'Target SI' }),
          h('th', { class: 'num', text: 'Dissolved' }), h('th', { class: 'num', text: 'Remaining' }))),
        h('tbody', {}, r.phases.map((p) => h('tr', {}, h('td', { text: p.phase }), h('td', { class: 'num', text: fixed(p.target_si, 2) }),
          h('td', { class: 'num', text: fmt(p.dissolved_mol, 4) }), h('td', { class: 'num', text: fmt(p.remaining_mol, 4) }))))));
      const si = [...r.saturation].sort((a, b) => b.si - a.si);
      $('aq-si').replaceChildren(h('table', {},
        h('thead', {}, h('tr', {}, h('th', { text: 'Phase' }), h('th', { text: 'Formula' }), h('th', { class: 'num', text: 'SI' }),
          h('th', { class: 'num', text: 'log IAP' }), h('th', { class: 'num', text: 'log K' }))),
        h('tbody', {}, si.map((s) => h('tr', {}, h('td', { text: s.phase }), h('td', { class: 'mono', text: s.formula }),
          h('td', { class: 'num ' + (s.si > 0.005 ? 'aq-over' : s.si < -0.005 ? 'aq-under' : ''), text: fixed(s.si, 3) }),
          h('td', { class: 'num', text: fixed(s.log_iap, 3) }), h('td', { class: 'num', text: fixed(s.log_k, 3) }))))));
      const comps = r.totals.map((t) => t.component);
      const prev = $('aq-dist-comp').value;
      $('aq-dist-comp').replaceChildren(...comps.map((c) => h('option', { value: c, text: c })));
      const pick = comps.includes(prev) ? prev : comps.find((c) => /^(C\(4\)|P|N\(-3\)|S\(-2\)|Fe\(3\)|Al)$/.test(c)) || comps[0];
      if (pick) { $('aq-dist-comp').value = pick; aqDistribution().catch(() => {}); }
    } catch (err) { fail(out, err); }
  });
  window.addEventListener('resize', () => { if (aqDist && current() === 'aqueous') drawAqDist(); });

  // ---- Flowsheet ------------------------------------------------------------------------------
  const FS_EXAMPLES = [
    { name: 'Flash with liquid recycle · benzene/toluene', text: `# Part of the flash liquid is recycled to the feed.
components c1ccccc1 Cc1ccccc1
feed F T=300 P=1atm flows c1ccccc1=50 Cc1ccccc1=50
mixer M F R -> S1
heater H S1 -> S2 T=368
flash D S2 -> V L T=368 P=1atm
splitter SP L -> R B fractions=0.6,0.4
` },
    { name: 'Propane dehydrogenation · Felder & Rousseau 4.7-2', text: `# 95 % overall conversion; the separator keeps 0.555 % of the
# propane and 5 % of the propylene for recycle.
components CCC C=CC [HH]
feed F T=300 P=1bar flows CCC=100
mixer M F R -> Q1
reactor RX Q1 -> Q2 key=CCC X=0.095391 T=800 nu CCC=-1 C=CC=1 [HH]=1
separator S Q2 -> P R T=300 split CCC=0.00555 C=CC=0.95 [HH]=1
` },
    { name: 'Methanol loop with purge · CO₂ + 3 H₂', text: `# Reactor feed 28 % CO2, 70 % H2, 2 % argon; 60 % single-pass H2
# conversion; the condenser takes out methanol and water.
components O=C=O [HH] [Ar] CO O
feed F T=300 P=50bar flows O=C=O=15.7193 [HH]=45.4386 [Ar]=0.245614
mixer M F R -> Q1
reactor RX Q1 -> Q2 key=[HH] X=0.6 T=500 nu O=C=O=-1 [HH]=-3 CO=1 O=1
separator C Q2 -> G L T=300 split O=C=O=1 [HH]=1 [Ar]=1
splitter SP G -> R PURGE fractions=0.877193,0.122807
` },
    { name: 'Hydrodealkylation · feed–effluent exchanger', text: `# Toluene + H2 -> benzene + CH4. The reactor effluent preheats the feed.
components Cc1ccccc1 [HH] c1ccccc1 C
feed F T=300 P=25bar flows Cc1ccccc1=10 [HH]=50
exchanger E Q2 F -> PROD Q0 Thot=500
heater H Q0 -> Q1 T=890
reactor RX Q1 -> Q2 key=Cc1ccccc1 X=0.75 nu Cc1ccccc1=-1 [HH]=-1 c1ccccc1=1 C=1
` },
    { name: 'Distillation · benzene/toluene', text: `# 15 equilibrium stages (1 = total condenser, 15 = reboiler).
components c1ccccc1 Cc1ccccc1
feed F T=360 P=1atm flows c1ccccc1=50 Cc1ccccc1=50
column T1 F -> D B stages=15 feed=8 R=2 D=50 P=1atm
` },
    { name: 'Distillation · ethanol/water (UNIFAC)', text: `# The distillate cannot pass the azeotrope (about 89 % ethanol).
components CCO O
model unifac
feed F T=350 P=1atm flows CCO=10 O=90
column T1 F -> D B stages=20 feed=12 R=3 D=11 P=1atm
` },
    { name: 'Water pump, boiler and throttle', text: `# No recycle: one pass through the units.
components O
feed W T=25C P=1bar flows O=1
pump P1 W -> W2 P=10bar eff=0.75
heater B W2 -> S T=500
valve V S -> S2 P=2bar
cooler C S2 -> L T=330
` },
  ];
  $('fs-example').replaceChildren(...FS_EXAMPLES.map((x, i) => h('option', { value: String(i), text: x.name })));
  const fsSaved = store.get('flowsheet', null);
  $('fs-text').value = fsSaved ?? FS_EXAMPLES[0].text;
  if (fsSaved != null) $('fs-example').append(h('option', { value: 'saved', text: 'Your last flowsheet' })), ($('fs-example').value = 'saved');
  $('fs-example').addEventListener('change', () => {
    const x = FS_EXAMPLES[+$('fs-example').value];
    if (x) { $('fs-text').value = x.text; $('fs-form').requestSubmit(); }
  });

  const SVG_NS = 'http://www.w3.org/2000/svg';
  function svg(tag, attrs, ...kids) {
    const n = document.createElementNS(SVG_NS, tag);
    for (const [k, v] of Object.entries(attrs || {})) if (v != null) n.setAttribute(k, v);
    for (const c of kids.flat()) if (c != null) n.append(c);
    return n;
  }
  const kW = (w) => fmt(w / 1e3, 4);
  const streamTip = (s) => `${s.name}: ${fixed(s.T_K, 2)} K, ${fmt(s.P_Pa / 1e5, 4)} bar, ${fmt(s.flows.reduce((a, b) => a + b, 0), 4)} mol/s, vapour fraction ${fixed(s.vapor_fraction, 3)}`;

  // Layered drawing: units in columns by longest path along non-tear streams, feeds just before the unit they
  // enter, products just after the unit they leave; rows ordered by a barycentre sweep. Tear streams run back
  // underneath, dashed.
  function flowsheetDiagram(r) {
    const tears = new Set(r.tears);
    const edges = r.streams.map((s) => ({ s, from: s.from === 'feed' ? 'in:' + s.name : s.from, to: s.to || 'out:' + s.name, tear: tears.has(s.name) }));
    const layer = new Map();
    for (const u of r.order) {
      let L = 1;
      for (const e of edges) if (e.to === u && !e.tear && !e.from.startsWith('in:')) L = Math.max(L, layer.get(e.from) + 1);
      layer.set(u, L);
    }
    for (const e of edges) {
      if (e.from.startsWith('in:')) layer.set(e.from, (layer.get(e.to) ?? 1) - 1);
      if (e.to.startsWith('out:')) layer.set(e.to, layer.get(e.from) + 1);
    }
    const nodes = [...layer.keys()];
    const maxL = Math.max(...layer.values());
    const cols = Array.from({ length: maxL + 1 }, () => []);
    const rank = (n) => (n.startsWith('in:') || n.startsWith('out:') ? r.order.length : r.order.indexOf(n));
    nodes.sort((a, b) => rank(a) - rank(b)).forEach((n) => cols[layer.get(n)].push(n));
    const pos = new Map();
    const place = () => cols.forEach((c) => c.forEach((n, i) => pos.set(n, i - (c.length - 1) / 2)));
    place();
    for (let sweep = 0; sweep < 3; sweep++) {
      for (const c of cols) {
        const bary = new Map(c.map((n) => {
          const nb = edges.filter((e) => !e.tear && (e.to === n || e.from === n)).map((e) => pos.get(e.to === n ? e.from : e.to));
          return [n, nb.length ? nb.reduce((a, b) => a + b, 0) / nb.length : pos.get(n)];
        }));
        c.sort((a, b) => bary.get(a) - bary.get(b));
      }
      place();
    }
    const colW = 168, rowH = 92, UW = 112, UH = 50, TW = 64, TH = 26, padX = 18;
    const rows = Math.max(...cols.map((c) => c.length));
    const top = 24 + (rows - 1) / 2 * rowH;
    const box = (n) => {
      const term = n.includes(':');
      const w = term ? TW : UW, hh = term ? TH : UH;
      const cx = padX + layer.get(n) * colW + UW / 2, cy = top + pos.get(n) * rowH;
      return { x: cx - w / 2, y: cy - hh / 2, w, h: hh, cx, cy };
    };
    // Ports: inlets on the left, outlets on the right, in stream order.
    const ins = new Map(), outs = new Map();
    for (const e of edges) {
      (ins.get(e.to) || ins.set(e.to, []).get(e.to)).push(e);
      (outs.get(e.from) || outs.set(e.from, []).get(e.from)).push(e);
    }
    const port = (n, list, e, right) => {
      const b = box(n), k = list.indexOf(e), m = list.length;
      return [right ? b.x + b.w : b.x, b.y + b.h * (k + 1) / (m + 1)];
    };
    const width = padX * 2 + maxL * colW + UW;
    let bottom = top + (rows - 1) / 2 * rowH + UH / 2 + 22;
    const parts = [];
    const labels = [];
    let back = 0;
    for (const e of edges) {
      const [x1, y1] = port(e.from, outs.get(e.from), e, true);
      const [x2, y2] = port(e.to, ins.get(e.to), e, false);
      let d, lx, ly;
      if (layer.get(e.to) > layer.get(e.from)) {
        const dx = Math.max(28, (x2 - x1) / 2);
        d = `M${x1} ${y1} C${x1 + dx} ${y1} ${x2 - dx} ${y2} ${x2} ${y2}`;
        lx = (x1 + x2) / 2; ly = (y1 + y2) / 2;
      } else {
        const yb = bottom + back++ * 22;
        d = `M${x1} ${y1} H${x1 + 16} V${yb} H${x2 - 16} V${y2} H${x2}`;
        lx = (x1 + x2) / 2; ly = yb;
      }
      parts.push(svg('path', { d, class: 'fs-edge' + (e.tear ? ' tear' : ''), 'marker-end': `url(#fs-arrow${e.tear ? '-tear' : ''})` },
        svg('title', {}, streamTip(e.s))));
      if (!e.from.startsWith('in:') && !e.to.startsWith('out:')) labels.push([lx, ly, e.s.name, e.tear]);
    }
    const height = bottom + Math.max(0, back - 1) * 22 + (back ? 22 : 0);
    for (const [x, y, t, tear] of labels) {
      const w = 8 + t.length * 7.2;
      parts.push(svg('g', { class: 'fs-label' + (tear ? ' tear' : '') }, svg('rect', { x: x - w / 2, y: y - 9, width: w, height: 18, rx: 9 }),
        svg('text', { x, y: y + 4, 'text-anchor': 'middle' }, t)));
    }
    const unitOf = new Map(r.units.map((u) => [u.name, u]));
    for (const n of nodes) {
      const b = box(n);
      if (n.includes(':')) {
        const s = edges.find((e) => e.from === n || e.to === n).s;
        parts.push(svg('g', { class: 'fs-term' }, svg('title', {}, streamTip(s)),
          svg('rect', { x: b.x, y: b.y, width: b.w, height: b.h, rx: b.h / 2 }),
          svg('text', { x: b.cx, y: b.cy + 4, 'text-anchor': 'middle' }, s.name)));
      } else {
        const u = unitOf.get(n);
        const tip = `${u.name} (${u.type})` + (u.duty_W ? `, duty ${kW(u.duty_W)} kW` : '') + (u.work_W ? `, work ${kW(u.work_W)} kW` : '');
        parts.push(svg('g', { class: 'fs-unit' }, svg('title', {}, tip),
          svg('rect', { x: b.x, y: b.y, width: b.w, height: b.h, rx: 12 }),
          svg('text', { x: b.cx, y: b.cy - 2, 'text-anchor': 'middle', class: 'name' }, u.name),
          svg('text', { x: b.cx, y: b.cy + 14, 'text-anchor': 'middle', class: 'type' }, u.type)));
      }
    }
    const arrow = (id) => svg('marker', { id, viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: 'auto-start-reverse' },
      svg('path', { d: 'M0 1 L9 5 L0 9 Z' }));
    // Scales down to 70 % to fit narrow windows, then scrolls.
    return svg('svg', { class: 'fs-svg', viewBox: `0 0 ${width} ${height}`, role: 'img', 'aria-label': 'Process flow diagram',
      style: `width: 100%; max-width: ${width}px; min-width: ${Math.round(width * 0.7)}px; height: auto` },
      svg('defs', {}, arrow('fs-arrow'), arrow('fs-arrow-tear')), ...parts);
  }

  const UNIT_VALUES = {
    outlet_T_K: ['outlet', (v) => fixed(v, 2) + ' K'], T_K: ['T', (v) => fixed(v, 2) + ' K'],
    vapor_fraction: ['vapour fraction', (v) => fixed(v, 4)], exchanger_duty_W: ['duty', (v) => kW(v) + ' kW'],
    min_approach_K: ['min. approach', (v) => fmt(v, 4) + ' K'], UA_WK: ['UA', (v) => fmt(v / 1e3, 4) + ' kW/K'],
    isentropic_work_W: ['isentropic', (v) => kW(v) + ' kW'], extent_mols: ['extent', (v) => fmt(v, 5) + ' mol/s'],
    heat_of_reaction_Jmol: ['ΔrH°', (v) => fmt(v / 1e3, 4) + ' kJ/mol'],
    condenser_duty_W: ['condenser', (v) => kW(v) + ' kW'], reboiler_duty_W: ['reboiler', (v) => kW(v) + ' kW'],
    top_T_K: ['top', (v) => fixed(v, 2) + ' K'], bottom_T_K: ['bottom', (v) => fixed(v, 2) + ' K'],
    newton_iterations: ['Newton iterations', (v) => String(v)],
  };
  let lastFs = null;
  function fsStreamTable() {
    const r = lastFs;
    const frac = checked($('fs-basis')).dataset.basis === 'frac';
    const head = h('tr', {}, h('th', { text: 'Stream' }), h('th', { text: 'From → to' }), h('th', { class: 'num', text: 'T (K)' }),
      h('th', { class: 'num', text: 'P (bar)' }), h('th', { class: 'num', text: 'Vapour' }), h('th', { class: 'num', text: 'mol/s' }),
      h('th', { class: 'num', text: 'kg/s' }), h('th', { class: 'num', text: 'H (kW)' }),
      ...r.components.map((c) => h('th', { class: 'num', title: c.smiles, text: c.name })));
    const rows = r.streams.map((s) => {
      const F = s.flows.reduce((a, b) => a + b, 0);
      return h('tr', { class: r.tears.includes(s.name) ? 'fs-tear-row' : null }, h('td', { class: 'mono', text: s.name }),
        h('td', { text: `${s.from === 'feed' ? 'feed' : s.from} → ${s.to || 'product'}` }),
        h('td', { class: 'num', text: fixed(s.T_K, 2) }), h('td', { class: 'num', text: fmt(s.P_Pa / 1e5, 4) }),
        h('td', { class: 'num', text: F > 0 ? fixed(s.vapor_fraction, 3) : '—' }), h('td', { class: 'num', text: fmt(F, 5) }),
        h('td', { class: 'num', text: fmt(s.mass_kgs, 4) }), h('td', { class: 'num', text: fmt(s.H_W / 1e3, 5) }),
        ...s.flows.map((n) => h('td', { class: 'num', text: frac ? (F > 0 ? fixed(n / F, 4) : '—') : fmt(n, 5) })));
    });
    $('fs-streams').replaceChildren(h('table', { class: 'fs-table' }, h('thead', {}, head), h('tbody', {}, rows)));
  }
  segmented($('fs-basis'), () => lastFs && fsStreamTable());
  function drawColumnProfile() {
    const col = lastFs && lastFs.units.find((u) => u.profile);
    $('fs-col-card').hidden = !col;
    if (!col) return;
    $('fs-col-title').textContent = `Column ${col.name}: stage profile`;
    const st = col.profile;
    Charts.lineChart($('fs-col-t'), [{ name: 'Temperature (K)', points: st.map((s, j) => [j + 1, s.T_K]) }],
      { xLabel: 'stage', yLabel: 'K' });
    Charts.lineChart($('fs-col-x'), lastFs.components.map((c, i) => ({ name: c.name, points: st.map((s, j) => [j + 1, s.x[i]]) })),
      { xLabel: 'stage', yLabel: 'liquid mole fraction' });
  }
  function drawFsChart() {
    drawColumnProfile();
    if (!lastFs || !lastFs.residuals.length) return;
    Charts.lineChart($('fs-chart'), [{ name: 'tear change', points: lastFs.residuals.map((v, i) => [i + 1, Math.max(v, 1e-16)]) }],
      { xLabel: 'iteration', yLabel: 'relative change', logY: true });
  }
  $('fs-text').addEventListener('input', () => store.set('flowsheet', $('fs-text').value));
  $('fs-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('fs-out');
    busy(out, 'Solving the flowsheet…');
    try {
      const { result: r, ms } = await engine('flowsheet', { text: $('fs-text').value });
      const errors = r.diagnostics.filter((d) => d.severity !== 'warning');
      const warnings = r.diagnostics.filter((d) => d.severity === 'warning');
      const warnBox = warnings.length ? h('details', { class: 'fs-warnings' }, h('summary', { text: `${warnings.length} note${warnings.length > 1 ? 's' : ''} on the property models` }), ...diagnostics(warnings)) : null;
      if (!r.streams) {
        $('fs-results').hidden = true;
        lastFs = null;
        out.replaceChildren(...[status('bad', 'The flowsheet could not be solved'), ...diagnostics(errors), warnBox].filter(Boolean));
        return;
      }
      const recycle = r.tears.length ? `${r.method === 'direct' ? 'direct substitution' : 'Wegstein'} · tear ${r.tears.join(', ')}` : 'no recycle';
      const massErr = r.balance.mass_in_kgs > 0 ? Math.abs(r.balance.mass_out_kgs - r.balance.mass_in_kgs) / r.balance.mass_in_kgs : 0;
      out.replaceChildren(
        r.converged ? status('ok', `Converged in ${r.iterations} iteration${r.iterations > 1 ? 's' : ''}`, recycle)
          : status('bad', `Not converged after ${r.iterations} iterations`, recycle),
        h('div', { class: 'metrics' },
          metric('Mass in', fmt(r.balance.mass_in_kgs, 5), 'kg/s'),
          metric('Mass balance', massErr < 1e-12 ? '0' : fmt(massErr, 2), 'relative'),
          metric('Energy balance', fmt(r.balance.energy_imbalance_W, 2), 'W'),
          metric('Compute time', String(Math.round(ms)), 'ms')),
        ...diagnostics(errors), ...(warnBox ? [warnBox] : []));
      lastFs = r;
      $('fs-results').hidden = false;
      $('fs-diagram').replaceChildren(flowsheetDiagram(r));
      fsStreamTable();
      $('fs-units').replaceChildren(h('table', {},
        h('thead', {}, h('tr', {}, h('th', { text: 'Unit' }), h('th', { text: 'Type' }), h('th', { class: 'num', text: 'Duty (kW)' }),
          h('th', { class: 'num', text: 'Work (kW)' }), h('th', { text: 'Results' }))),
        h('tbody', {}, r.units.map((u) => h('tr', {}, h('td', { class: 'mono', text: u.name }), h('td', { text: u.type }),
          h('td', { class: 'num', text: u.duty_W ? kW(u.duty_W) : '—' }), h('td', { class: 'num', text: u.work_W ? kW(u.work_W) : '—' }),
          h('td', { class: 'fs-values', text: Object.entries(u.values).map(([k, v]) => (UNIT_VALUES[k] ? `${UNIT_VALUES[k][0]} ${UNIT_VALUES[k][1](v)}` : `${k} ${fmt(v, 4)}`)).join(' · ') || '—' }))))));
      $('fs-conv-card').hidden = !r.residuals.length;
      drawFsChart();
    } catch (err) { fail(out, err); }
  });
  window.addEventListener('resize', () => { if (lastFs && current() === 'flowsheet') drawFsChart(); });

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
