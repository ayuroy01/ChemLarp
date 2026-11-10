// OmniChem web app: UI glue around the WebAssembly engine (engine-worker.js).
(function () {
  'use strict';
  const BUILD = '__BUILD_ID__';
  const $ = (id) => document.getElementById(id);

  // ---- engine (Web Worker) --------------------------------------------------
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
  function engine(op, args) {
    return new Promise((resolve, reject) => {
      const id = nextId++;
      pending.set(id, { resolve, reject });
      worker.postMessage({ id, op, args });
    });
  }
  engine('version').then((r) => {
    $('engine-status').textContent = 'Engine ready · works offline';
    $('version').textContent = r.result;
  }).catch((err) => { $('engine-status').textContent = 'Engine failed to load: ' + err.message; });

  // ---- small DOM helpers ----------------------------------------------------
  function el(tag, attrs, ...children) {
    const n = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === 'class') n.className = v;
      else if (k === 'text') n.textContent = v;
      else n.setAttribute(k, v);
    }
    for (const c of children) if (c != null) n.append(c);
    return n;
  }
  const num = (v, digits = 4) => (v == null || !isFinite(v) ? '—' : Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(digits - 1) : (+v.toPrecision(digits)).toString());
  function facts(pairs) {
    const dl = el('dl', { class: 'facts' });
    for (const [k, v] of pairs) dl.append(el('div', {}, el('dt', { text: k }), el('dd', { text: v })));
    return dl;
  }
  function busy(out, text) { out.replaceChildren(el('p', { class: 'busy', text })); }
  function fail(out, err) { out.replaceChildren(el('div', { class: 'diag' }, el('div', { text: 'Something went wrong: ' + err.message }))); }

  // Engine diagnostics, with a caret under the character they point at.
  function diagnostics(list, input) {
    const box = el('div', { class: 'result' });
    for (const d of list || []) {
      const item = el('div', { class: 'diag' + (d.severity === 'warning' ? ' warning' : '') },
        el('div', { class: 'code', text: `${d.severity} · ${d.code}` }),
        el('div', { text: d.message }));
      if (input && d.position != null && d.position <= input.length) {
        item.append(el('pre', { text: input + '\n' + ' '.repeat(d.position) + '^' }));
      }
      box.append(item);
    }
    return box;
  }

  function download(name, text, type) {
    const url = URL.createObjectURL(new Blob([text], { type }));
    const a = el('a', { href: url, download: name });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  // ---- tabs (hash routing) --------------------------------------------------
  const tabs = ['molecule', 'model3d', 'reaction', 'equilibrium', 'kinetics', 'phases', 'progress'];
  // Each tool runs its example the first time it is opened, so it never starts empty.
  const autorun = { molecule: 'mol-form', model3d: 'd3-form', reaction: 'rx-form', equilibrium: 'eq-form', kinetics: 'kin-form', phases: 'ph-form' };
  const ran = new Set();
  function show() {
    const name = tabs.includes(location.hash.slice(1)) ? location.hash.slice(1) : 'molecule';
    for (const t of tabs) $(t).hidden = t !== name;
    document.querySelectorAll('.tabs a').forEach((a) => (a.dataset.tab === name ? a.setAttribute('aria-current', 'page') : a.removeAttribute('aria-current')));
    if (name === 'model3d' && viewer) viewer.draw();
    if (name === 'progress') loadProgress();
    if (autorun[name] && !ran.has(name)) {
      ran.add(name);
      $(autorun[name]).requestSubmit();
    }
  }
  window.addEventListener('hashchange', show);

  // Example buttons fill their target input.
  document.querySelectorAll('.examples').forEach((box) => {
    box.addEventListener('click', (e) => {
      const b = e.target.closest('button[data-smiles]');
      if (!b) return;
      $(box.dataset.target).value = b.dataset.smiles;
      $(box.dataset.target).form.requestSubmit();
    });
  });

  // ---- Molecule -------------------------------------------------------------
  $('mol-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const smiles = $('mol-smiles').value.trim();
    const out = $('mol-out');
    busy(out, 'Analysing…');
    try {
      const [p, c] = await Promise.all([engine('parse', { smiles }), engine('canonical', { smiles, random: 3, seed: 7 })]);
      const r = p.result;
      out.replaceChildren();
      if (!r.ok) {
        out.append(el('p', { class: 'bad-line', text: 'This SMILES has problems:' }), diagnostics(r.diagnostics, smiles));
        return;
      }
      const s = r.structure.species[0];
      const cip = c.result.ok ? c.result.cip : { centres: [], double_bonds: [] };
      const stereo = [...cip.centres.map((x) => `atom ${x.atom}: ${x.label || 'undetermined'}`), ...cip.double_bonds.map((x) => `bond ${x.bond}: ${x.label || 'undetermined'}`)];
      const card = el('div', { class: 'card' },
        el('div', { class: 'big', text: s.formula }),
        facts([
          ['Canonical SMILES', s.canonical_smiles],
          ['Molar mass', num(s.molar_mass, 6) + ' g/mol'],
          ['Atoms (heavy)', String(s.atom_count)],
          ['Bonds', String(s.bond_count)],
          ['Net charge', String(s.net_charge)],
          ['Stereo (CIP)', stereo.length ? stereo.join(', ') : 'none'],
        ]));
      const actions = el('div', { class: 'row tools' });
      const to3d = el('button', { type: 'button', class: 'secondary', text: 'Open in 3D' });
      to3d.addEventListener('click', () => { $('d3-smiles').value = smiles; location.hash = 'model3d'; $('d3-form').requestSubmit(); });
      actions.append(to3d);
      card.append(actions);
      out.append(card);
      if (c.result.ok && c.result.random.length) {
        out.append(el('div', { class: 'card' }, el('div', { class: 'code', text: 'Other valid ways to write it (all give the same canonical SMILES)' }),
          ...c.result.random.map((x) => el('div', { class: 'mono', text: x }))));
      }
      if (r.diagnostics.length) out.append(diagnostics(r.diagnostics, smiles));
    } catch (err) { fail(out, err); }
  });

  // ---- 3D -------------------------------------------------------------------
  let viewer = null;
  let last3d = null;
  $('d3-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const smiles = $('d3-smiles').value.trim();
    const out = $('d3-out');
    busy(out, 'Building 3D structure…');
    try {
      const { result: r, ms } = await engine('embed', { smiles, seed: 42, candidates: 5, minimize: true });
      out.replaceChildren();
      if (!r.ok) { out.append(diagnostics(r.diagnostics, smiles)); return; }
      if (!viewer) viewer = new MoleculeViewer($('viewer'));
      viewer.labels = $('d3-labels').checked;
      viewer.load(r.sdf);
      $('viewer-empty').hidden = true;
      last3d = { smiles, sdf: r.sdf, xyz: r.xyz };
      $('d3-sdf').disabled = $('d3-xyz').disabled = false;
      const st = r.stereo;
      out.append(el('div', { class: 'card' }, facts([
        ['UFF energy', num(r.energy_kcal_mol, 5) + ' kcal/mol'],
        ['Minimisation', r.converged ? 'converged' : 'not converged'],
        ['Stereocentres reproduced', `${st.centres_ok} of ${st.centres}`],
        ['Stereo double bonds', `${st.bonds_ok} of ${st.bonds}`],
        ['Attempts', String(r.attempts)],
        ['Time', `${Math.round(ms)} ms`],
      ])));
      if (r.diagnostics.length) out.append(diagnostics(r.diagnostics, smiles));
    } catch (err) { fail(out, err); }
  });
  $('d3-labels').addEventListener('change', () => { if (viewer) { viewer.labels = $('d3-labels').checked; viewer.draw(); } });
  $('d3-sdf').addEventListener('click', () => last3d && download('structure.sdf', last3d.sdf, 'chemical/x-mdl-sdfile'));
  $('d3-xyz').addEventListener('click', () => last3d && download('structure.xyz', last3d.xyz, 'chemical/x-xyz'));

  // ---- Reaction -------------------------------------------------------------
  $('rx-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const equation = $('rx-eq').value.trim();
    const T = +$('rx-t').value;
    const out = $('rx-out');
    busy(out, 'Checking…');
    try {
      const { result: r } = await engine('reaction', { equation, temperature: T });
      out.replaceChildren();
      out.append(el('p', { class: r.balanced ? 'ok-line' : 'bad-line', text: r.balanced ? 'Balanced: every element and charge is conserved.' : 'Not balanced.' }));
      if (r.thermo) {
        const t = r.thermo;
        out.append(el('div', { class: 'card' }, facts([
          ['ΔrH°', num(t.dH_J / 1000, 5) + ' kJ/mol'],
          ['ΔrS°', num(t.dS_J, 5) + ' J/(mol·K)'],
          ['ΔrG°', num(t.dG_RT * 8.314462618 * t.temperature_K / 1000, 5) + ' kJ/mol'],
          ['K (1 bar)', num(t.K, 4)],
          ['Temperature', num(t.temperature_K, 6) + ' K'],
        ])));
      }
      if (r.diagnostics.length) out.append(diagnostics(r.diagnostics, equation));
    } catch (err) { fail(out, err); }
  });

  // ---- Equilibrium ----------------------------------------------------------
  $('eq-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('eq-out');
    busy(out, 'Minimising Gibbs energy…');
    try {
      const { result: r } = await engine('equilibrate', { species: $('eq-species').value, temperature: +$('eq-t').value, pressure: +$('eq-p').value * 1e5 });
      out.replaceChildren();
      $('eq-chart').hidden = true;
      if (!r.equilibrium) { out.append(diagnostics(r.diagnostics)); return; }
      const q = r.equilibrium;
      out.append(el('p', { class: q.converged ? 'ok-line' : 'bad-line', text: q.converged ? `Converged in ${q.iterations} iterations.` : 'Did not converge.' }));
      const tbody = el('tbody');
      for (const s of q.species) tbody.append(el('tr', {}, el('td', { class: 'mono', text: s.smiles }), el('td', { text: s.formula }), el('td', { class: 'num', text: num(s.moles) }), el('td', { class: 'num', text: num(s.mole_fraction) })));
      out.append(el('div', { class: 'card table-wrap' }, el('table', {},
        el('thead', {}, el('tr', {}, el('th', { text: 'Species' }), el('th', { text: 'Formula' }), el('th', { class: 'num', text: 'Moles' }), el('th', { class: 'num', text: 'Mole fraction' }))),
        tbody)));
      Charts.barChart($('eq-chart'), q.species.map((s) => ({ label: s.formula, value: s.mole_fraction })), { caption: 'Mole fraction at equilibrium', format: (v) => num(v, 3) });
      if (r.diagnostics.length) out.append(diagnostics(r.diagnostics));
    } catch (err) { fail(out, err); }
  });

  // ---- Kinetics -------------------------------------------------------------
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
      { xLabel: 'time (ms)', yLabel: 'mole fraction', logY: $('kin-log').checked });
  }
  $('kin-log').addEventListener('change', drawKinetics);
  $('kin-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('kin-out');
    busy(out, 'Integrating…');
    try {
      const { mechanism, thermo } = await mechanismText();
      const { result: r, ms } = await engine('kinetics', {
        mechanism, thermo, model: $('kin-model').value, temperature: +$('kin-t').value, pressure: +$('kin-p').value * 1e5,
        feed: $('kin-feed').value, tEnd: +$('kin-time').value / 1000, samples: 200,
      });
      out.replaceChildren();
      if (!r.ok) {
        ['kin-chart-t', 'kin-chart-x', 'kin-scale'].forEach((id) => ($(id).hidden = true));
        out.append(el('p', { class: 'bad-line', text: 'The run did not complete:' }), diagnostics(r.diagnostics));
        return;
      }
      const end = r.samples[r.samples.length - 1];
      out.append(el('div', { class: 'card' }, facts([
        ['Ignition time', r.ignition_time ? num(r.ignition_time * 1000, 4) + ' ms' : 'no ignition'],
        ['Final temperature', num(end.T, 5) + ' K'],
        ['Final pressure', num(end.P / 1e5, 4) + ' bar'],
        ['Solver steps', String(r.stats.steps)],
        ['Time', `${Math.round(ms)} ms`],
      ])));
      lastKin = { samples: r.samples, species: r.species };
      $('kin-scale').hidden = false;
      drawKinetics();
      if (r.diagnostics.length) out.append(diagnostics(r.diagnostics));
    } catch (err) { fail(out, err); }
  });

  // ---- Phases ---------------------------------------------------------------
  const regionText = { 'two-phase': 'Two phases: liquid + vapour', liquid: 'Single phase: liquid', vapor: 'Single phase: vapour', supercritical: 'Single phase: supercritical fluid' };
  $('ph-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('ph-out');
    busy(out, 'Flashing…');
    try {
      const { result: r } = await engine('flash', { mixture: $('ph-mix').value, temperature: +$('ph-t').value, pressure: +$('ph-p').value * 1e5, eos: $('ph-eos').value });
      out.replaceChildren();
      if (!r.ok) { out.append(el('p', { class: 'bad-line', text: 'Could not flash this mixture:' }), diagnostics(r.diagnostics)); return; }
      const two = r.region === 'two-phase';
      out.append(el('div', { class: 'card' },
        el('p', {}, el('span', { class: 'region', text: regionText[r.region] || r.region })),
        facts([
          ['Vapour fraction', two ? num(r.vapor_fraction, 5) : r.region === 'vapor' ? '1' : '0'],
          ['Z liquid', two || r.region === 'liquid' ? num(r.liquid.Z, 5) : '—'],
          ['Z vapour', two || r.region !== 'liquid' ? num(r.vapor.Z, 5) : '—'],
          ['Iterations', String(r.iterations)],
        ])));
      const tbody = el('tbody');
      for (const c of r.components) tbody.append(el('tr', {}, el('td', { text: c.name }), el('td', { class: 'mono', text: c.smiles }),
        el('td', { class: 'num', text: num(c.z, 4) }), el('td', { class: 'num', text: two ? num(c.x, 4) : '—' }),
        el('td', { class: 'num', text: two ? num(c.y, 4) : '—' }), el('td', { class: 'num', text: two ? num(c.K, 4) : '—' })));
      out.append(el('div', { class: 'card table-wrap' }, el('table', {},
        el('thead', {}, el('tr', {}, el('th', { text: 'Component' }), el('th', { text: 'SMILES' }), el('th', { class: 'num', text: 'Feed' }),
          el('th', { class: 'num', text: 'Liquid x' }), el('th', { class: 'num', text: 'Vapour y' }), el('th', { class: 'num', text: 'K = y/x' }))),
        tbody)));
      if (r.diagnostics.length) out.append(diagnostics(r.diagnostics));
    } catch (err) { fail(out, err); }
  });
  $('vp-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const out = $('vp-out');
    busy(out, 'Solving…');
    try {
      const smiles = $('vp-smiles').value.trim();
      const { result: r } = await engine('psat', { smiles, temperature: +$('vp-t').value, eos: $('ph-eos').value });
      out.replaceChildren();
      if (!r.ok) { out.append(diagnostics(r.diagnostics, smiles)); return; }
      out.append(el('div', { class: 'card' }, facts([
        ['Compound', r.name],
        ['Vapour pressure', num(r.pressure_Pa / 1e5, 5) + ' bar'],
        ['', num(r.pressure_Pa / 1e3, 5) + ' kPa'],
      ])));
    } catch (err) { fail(out, err); }
  });

  // ---- Progress -------------------------------------------------------------
  let progressLoaded = false;
  async function loadProgress() {
    if (progressLoaded) return;
    const out = $('prog-out');
    busy(out, 'Loading roadmap…');
    try {
      const road = await (await fetch('data/roadmap.json')).json();
      progressLoaded = true;
      const w = road.weights;
      const all = road.areas.flatMap((a) => a.modules.map((m) => ({ ...m, area: a.name, replaces: a.replaces, weight: w[m.size] })));
      const sum = (ms, st) => ms.filter((m) => !st || m.status === st).reduce((s, m) => s + m.weight, 0);
      const total = sum(all), done = sum(all, 'done'), prog = sum(all, 'in_progress');
      const bar = (d, p, t) => el('div', { class: 'meter' },
        el('div', { style: `width:${(100 * d) / t}%;background:var(--ok)` }),
        el('div', { style: `width:${(100 * p) / t}%;background:var(--warn)` }));
      out.replaceChildren(
        el('div', { class: 'card' },
          el('div', { class: 'big', text: `${((100 * done) / total).toFixed(1)}% built` }),
          el('p', { class: 'hint', text: `${all.filter((m) => m.status === 'done').length} of ${all.length} modules done, ${((100 * prog) / total).toFixed(1)}% more in progress. Updated ${road.updated}.` }),
          bar(done, prog, total),
          el('div', { class: 'legend' }, el('span', {}, el('i', { style: 'background:var(--ok)' }), 'Done'), el('span', {}, el('i', { style: 'background:var(--warn)' }), 'In progress'), el('span', {}, el('i', { style: 'background:var(--line)' }), 'Planned'))),
        el('div', { class: 'card' }, ...road.areas.map((a) => {
          const t = sum(a.modules.map((m) => ({ ...m, weight: w[m.size] })));
          const d = sum(a.modules.map((m) => ({ ...m, weight: w[m.size] })), 'done');
          const p = sum(a.modules.map((m) => ({ ...m, weight: w[m.size] })), 'in_progress');
          return el('div', { class: 'area-row' },
            el('div', { class: 'area-head' }, el('strong', { text: a.name }), el('span', { text: `${a.replaces} · ${Math.round((100 * d) / t)}% done` })),
            bar(d, p, t));
        })),
        el('div', { class: 'card table-wrap' }, el('table', {},
          el('thead', {}, el('tr', {}, el('th', { text: 'Module' }), el('th', { text: 'Area' }), el('th', { text: 'Size' }), el('th', { text: 'Status' }))),
          el('tbody', {}, ...all.sort((a, b) => ['in_progress', 'done', 'planned'].indexOf(a.status) - ['in_progress', 'done', 'planned'].indexOf(b.status)).map((m) =>
            el('tr', {}, el('td', { text: m.name }), el('td', { text: m.area }), el('td', { text: m.size }),
              el('td', {}, el('span', { class: 'pill ' + m.status, text: m.status.replace('_', ' ') }))))))));
    } catch (err) { fail(out, err); }
  }

  // ---- install + offline ----------------------------------------------------
  let installEvent = null;
  window.addEventListener('beforeinstallprompt', (e) => { e.preventDefault(); installEvent = e; $('install-btn').hidden = false; });
  $('install-btn').addEventListener('click', async () => {
    if (!installEvent) return;
    installEvent.prompt();
    await installEvent.userChoice;
    installEvent = null;
    $('install-btn').hidden = true;
  });
  if ('serviceWorker' in navigator && location.protocol !== 'file:') {
    navigator.serviceWorker.register('sw.js').catch(() => {});
  }
  console.info('OmniChem web build', BUILD);

  show();
})();
