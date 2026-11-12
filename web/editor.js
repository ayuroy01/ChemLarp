// ChemLarp structure editor: draw and edit a molecule on an SVG canvas.
//
// Model: atoms {id, el, x, y, charge, iso} and bonds {id, a, b, order 1-3,
// stereo 0 | 1 (wedge) | 6 (hash)}, with a wedge's narrow end at bond.a.
// Coordinates are canvas units (one bond = BOND). The editor exchanges
// structures as MDL V2000 molfiles, so the engine's reader decides validity,
// stereo from wedges, and the canonical SMILES.
'use strict';

(function () {
  const BOND = 40;
  const SVGNS = 'http://www.w3.org/2000/svg';
  const VALENCE = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1] };
  const Z = { H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, Si: 14, P: 15, S: 16, Cl: 17, Br: 35, I: 53 };
  const BY_Z = Object.fromEntries(Object.entries(Z).map(([k, v]) => [v, k]));
  const COLOR = { O: 'var(--el-O)', N: 'var(--el-N)', S: 'var(--el-S)', P: 'var(--el-P)', F: 'var(--el-F)', Cl: 'var(--el-Cl)', Br: 'var(--el-Br)', I: 'var(--el-I)' };
  const KEY_ELEMENT = { c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', f: 'F', i: 'I', h: 'H', l: 'Cl', b: 'Br' };

  const dist = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
  const angle = (p, q) => Math.atan2(q.y - p.y, q.x - p.x);
  const snap = (a) => Math.round(a / (Math.PI / 6)) * (Math.PI / 6);

  // Implicit hydrogens by the MDL rule the engine uses: B C N O P S F Cl Br I fill
  // their lowest default valence; a charge shifts the element to its isoelectronic
  // neighbour (N+ like C, O- like F). Other elements get none.
  function implicitH(el, charge, bondSum) {
    if (!VALENCE[el]) return 0;
    let table = VALENCE[el];
    if (charge) {
      const shifted = BY_Z[Z[el] - charge];
      if (!shifted || !VALENCE[shifted]) return 0;
      table = VALENCE[shifted];
    }
    for (const v of table) if (v >= bondSum) return v - bondSum;
    return 0;
  }

  class StructureEditor {
    constructor(host, { onChange } = {}) {
      this.host = host;
      this.onChange = onChange || (() => {});
      this.atoms = [];
      this.bonds = [];
      this.nextId = 1;
      this.tool = { kind: 'bond', order: 1, stereo: 0 };
      this.element = 'C';
      this.selection = new Set();
      this.errors = new Set();
      this.view = { x: -200, y: -150, scale: 1 };
      this.undoStack = [];
      this.redoStack = [];
      this.hover = null;  // {type: 'atom'|'bond', id}
      this.drag = null;
      this.svg = document.createElementNS(SVGNS, 'svg');
      this.svg.setAttribute('class', 'editor-canvas');
      this.svg.setAttribute('tabindex', '0');
      this.svg.setAttribute('role', 'application');
      this.svg.setAttribute('aria-label', 'Structure drawing canvas');
      host.append(this.svg);
      this.bind();
      new ResizeObserver(() => this.render()).observe(host);
      this.render();
    }

    // ---- model helpers ---------------------------------------------------------------------
    atom(id) { return this.atoms.find((a) => a.id === id); }
    bondBetween(a, b) { return this.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a)); }
    neighbours(id) { return this.bonds.filter((b) => b.a === id || b.b === id).map((b) => (b.a === id ? b.b : b.a)); }
    bondSum(id) { return this.bonds.filter((b) => b.a === id || b.b === id).reduce((s, b) => s + b.order, 0); }
    addAtom(x, y, el = 'C') { const a = { id: this.nextId++, el, x, y, charge: 0, iso: 0 }; this.atoms.push(a); return a; }
    addBond(a, b, order = 1, stereo = 0) {
      if (a === b) return null;
      const old = this.bondBetween(a, b);
      if (old) return old;
      const bond = { id: this.nextId++, a, b, order, stereo };
      this.bonds.push(bond);
      return bond;
    }
    atomNear(p, r = BOND * 0.3, except = null) {
      let best = null, bd = r;
      for (const a of this.atoms) {
        if (a.id === except) continue;
        const d = dist(a, p);
        if (d < bd) { bd = d; best = a; }
      }
      return best;
    }
    bondNear(p, r = BOND * 0.22) {
      let best = null, bd = r;
      for (const b of this.bonds) {
        const A = this.atom(b.a), B = this.atom(b.b);
        const L2 = (B.x - A.x) ** 2 + (B.y - A.y) ** 2;
        const t = Math.max(0.15, Math.min(0.85, ((p.x - A.x) * (B.x - A.x) + (p.y - A.y) * (B.y - A.y)) / L2));
        const d = Math.hypot(p.x - (A.x + t * (B.x - A.x)), p.y - (A.y + t * (B.y - A.y)));
        if (d < bd) { bd = d; best = b; }
      }
      return best;
    }
    // Free direction for a new bond from atom `a`: away from neighbours, zig-zag on chains.
    freeAngle(a) {
      const nb = this.neighbours(a.id).map((id) => angle(a, this.atom(id)));
      if (nb.length === 0) return -Math.PI / 6;
      if (nb.length === 1) {
        const cands = [nb[0] + (2 * Math.PI) / 3, nb[0] - (2 * Math.PI) / 3];
        // continue a zig-zag: the new bond runs parallel to the bond two steps back (Q->P)
        const prev = this.neighbours(a.id)[0];
        const pp = this.neighbours(prev).filter((id) => id !== a.id)[0];
        if (pp !== undefined) {
          const want = angle(this.atom(pp), this.atom(prev));
          const diff = (c) => Math.abs(Math.atan2(Math.sin(c - want), Math.cos(c - want)));
          const pick = cands.slice().sort((x, y) => diff(x) - diff(y))[0];
          if (this.crowding(a, pick) < 0.5) return pick;
        }
        // least crowded; on a tie the more horizontal direction, which draws the usual zig-zag
        return cands.map((c) => [c, this.crowding(a, c), Math.abs(Math.sin(c))])
          .sort((x, y) => (Math.abs(x[1] - y[1]) > 1e-6 ? x[1] - y[1] : x[2] - y[2]))[0][0];
      }
      const s = nb.map((x) => (x + 2 * Math.PI) % (2 * Math.PI)).sort((x, y) => x - y);
      let best = 0, gap = -1;
      for (let i = 0; i < s.length; i++) {
        const g = (i + 1 < s.length ? s[i + 1] : s[0] + 2 * Math.PI) - s[i];
        if (g > gap) { gap = g; best = s[i] + g / 2; }
      }
      return best;
    }
    crowding(a, ang) {
      const p = { x: a.x + BOND * Math.cos(ang), y: a.y + BOND * Math.sin(ang) };
      return this.atoms.reduce((s, q) => s + (q.id === a.id ? 0 : 1 / Math.max(1, dist(p, q) ** 2)), 0);
    }

    // ---- history ------------------------------------------------------------------------------
    snapshot() { return JSON.stringify({ atoms: this.atoms, bonds: this.bonds, nextId: this.nextId }); }
    restore(s) { const o = JSON.parse(s); this.atoms = o.atoms; this.bonds = o.bonds; this.nextId = o.nextId; }
    commit(before) {
      if (before === this.snapshot()) return;
      this.undoStack.push(before);
      if (this.undoStack.length > 200) this.undoStack.shift();
      this.redoStack = [];
      this.changed();
    }
    undo() { if (!this.undoStack.length) return; this.redoStack.push(this.snapshot()); this.restore(this.undoStack.pop()); this.selection.clear(); this.changed(); }
    redo() { if (!this.redoStack.length) return; this.undoStack.push(this.snapshot()); this.restore(this.redoStack.pop()); this.selection.clear(); this.changed(); }
    changed() { this.errors.clear(); this.render(); this.onChange(this); }
    clear() { const before = this.snapshot(); this.atoms = []; this.bonds = []; this.selection.clear(); this.commit(before); }

    // ---- edits ---------------------------------------------------------------------------------
    deleteAtom(id) {
      this.atoms = this.atoms.filter((a) => a.id !== id);
      this.bonds = this.bonds.filter((b) => b.a !== id && b.b !== id);
      this.selection.delete(id);
    }
    deleteSelection() {
      if (!this.selection.size) return;
      const before = this.snapshot();
      for (const id of [...this.selection]) this.deleteAtom(id);
      this.commit(before);
    }
    setElement(id, el) {
      const before = this.snapshot();
      const a = this.atom(id);
      a.el = el;
      this.commit(before);
    }
    adjustCharge(id, delta) {
      const before = this.snapshot();
      const a = this.atom(id);
      a.charge = Math.max(-4, Math.min(4, a.charge + delta));
      this.commit(before);
    }
    // Clicking a bond with a bond tool: cycle order, or set / flip a wedge.
    applyBondTool(b) {
      const before = this.snapshot();
      const t = this.tool;
      if (t.stereo) {
        if (b.order !== 1) b.order = 1;
        if (b.stereo === t.stereo) [b.a, b.b] = [b.b, b.a];  // second click flips the narrow end
        b.stereo = t.stereo;
      } else if (t.order === 1) {
        if (b.stereo) b.stereo = 0;
        else b.order = b.order === 3 ? 1 : b.order + 1;
      } else {
        b.stereo = 0;
        b.order = b.order === t.order ? 1 : t.order;
      }
      this.commit(before);
    }
    // New atom bonded to `a` in direction `ang`, merging with an atom already there.
    grow(a, ang, el = 'C') {
      const p = { x: a.x + BOND * Math.cos(ang), y: a.y + BOND * Math.sin(ang) };
      const target = this.atomNear(p, BOND * 0.3, a.id) || this.addAtom(p.x, p.y, el);
      const bond = this.addBond(a.id, target.id, this.tool.stereo ? 1 : this.tool.order, this.tool.stereo || 0);
      return { target, bond };
    }
    // Ring of n atoms (benzene when aromatic) sharing the edge p->q, built on the side away from crowding.
    ringOnEdge(pa, pb, n, aromatic, side) {
      const mx = (pa.x + pb.x) / 2, my = (pa.y + pb.y) / 2;
      const ux = (pb.x - pa.x) / BOND, uy = (pb.y - pa.y) / BOND;
      const nx = -uy * side, ny = ux * side;
      const apothem = BOND / (2 * Math.tan(Math.PI / n));
      const c = { x: mx + nx * apothem, y: my + ny * apothem };
      const R = BOND / (2 * Math.sin(Math.PI / n));
      const start = angle(c, pa), step = (2 * Math.PI) / n;
      const dir = Math.sign((pb.x - c.x) * Math.sin(start) - (pb.y - c.y) * Math.cos(start)) >= 0 ? -1 : 1;
      const ids = [];
      for (let k = 0; k < n; k++) {
        const p = { x: c.x + R * Math.cos(start + dir * k * step), y: c.y + R * Math.sin(start + dir * k * step) };
        const near = this.atomNear(p, BOND * 0.3);
        ids.push((near || this.addAtom(p.x, p.y)).id);
      }
      const ring = [];
      for (let k = 0; k < n; k++) ring.push(this.addBond(ids[k], ids[(k + 1) % n], 1, 0));
      if (aromatic) this.kekulize(ring);
      return ids;
    }
    // Alternating double bonds around a new benzene ring: keep double bonds the ring
    // already has, and never give an atom two double bonds.
    kekulize(ring) {
      const inRing = new Set(ring);
      const outsideDouble = (id) => this.bonds.some((b) => !inRing.has(b) && (b.a === id || b.b === id) && b.order === 2);
      for (const offset of [0, 1]) {
        const picks = ring.filter((_, k) => k % 2 === offset);
        const keepsExisting = ring.every((b) => b.order !== 2 || picks.includes(b));
        const fits = picks.every((b) => !outsideDouble(b.a) && !outsideDouble(b.b));
        if (keepsExisting && fits) { for (const b of picks) b.order = 2; return; }
      }
      for (const b of ring) {  // fused to a system that fixes the pattern: fill what fits
        const busy = (id) => this.bonds.some((x) => x !== b && (x.a === id || x.b === id) && x.order === 2);
        if (b.order === 1 && !busy(b.a) && !busy(b.b)) b.order = 2;
      }
    }
    // Ring of n atoms with one vertex at atom `a`, its centre in direction `dir` from `a`.
    ringAtAtom(a, dir, n, aromatic) {
      const R = BOND / (2 * Math.sin(Math.PI / n));
      const c = { x: a.x + R * Math.cos(dir), y: a.y + R * Math.sin(dir) };
      const ids = [a.id];
      for (let k = 1; k < n; k++) {
        const t = dir + Math.PI + (2 * Math.PI * k) / n;
        const p = { x: c.x + R * Math.cos(t), y: c.y + R * Math.sin(t) };
        const near = this.atomNear(p, BOND * 0.3);
        ids.push((near || this.addAtom(p.x, p.y)).id);
      }
      const ring = [];
      for (let k = 0; k < n; k++) ring.push(this.addBond(ids[k], ids[(k + 1) % n], 1, 0));
      if (aromatic) this.kekulize(ring);
      return ids;
    }
    sideAwayFromCrowd(pa, pb) {
      const mx = (pa.x + pb.x) / 2, my = (pa.y + pb.y) / 2;
      const ux = (pb.x - pa.x), uy = (pb.y - pa.y);
      let s = 0;
      for (const q of this.atoms) s += Math.sign(ux * (q.y - my) - uy * (q.x - mx));
      return s > 0 ? -1 : 1;
    }

    // ---- molfile I/O ------------------------------------------------------------------------
    toMolfile(title = 'ChemLarp editor') {
      const ids = new Map(this.atoms.map((a, i) => [a.id, i + 1]));
      const f = (v) => v.toFixed(4).padStart(10);
      const lines = [title, '  ChemLarp          2D', '',
        `${String(this.atoms.length).padStart(3)}${String(this.bonds.length).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`];
      for (const a of this.atoms)
        lines.push(`${f((a.x / BOND) * 1.5)}${f((-a.y / BOND) * 1.5)}${f(0)} ${a.el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
      for (const b of this.bonds)
        lines.push(`${String(ids.get(b.a)).padStart(3)}${String(ids.get(b.b)).padStart(3)}${String(b.order).padStart(3)}${String(b.stereo).padStart(3)}`);
      const chg = this.atoms.filter((a) => a.charge).map((a) => [ids.get(a.id), a.charge]);
      for (let i = 0; i < chg.length; i += 8) {
        const part = chg.slice(i, i + 8);
        lines.push(`M  CHG${String(part.length).padStart(3)}` + part.map(([k, v]) => ` ${String(k).padStart(3)} ${String(v).padStart(3)}`).join(''));
      }
      const iso = this.atoms.filter((a) => a.iso).map((a) => [ids.get(a.id), a.iso]);
      for (let i = 0; i < iso.length; i += 8) {
        const part = iso.slice(i, i + 8);
        lines.push(`M  ISO${String(part.length).padStart(3)}` + part.map(([k, v]) => ` ${String(k).padStart(3)} ${String(v).padStart(3)}`).join(''));
      }
      lines.push('M  END');
      return lines.join('\n') + '\n';
    }

    // Loads a V2000 molfile (as written by the engine's depiction), scaled to the canvas bond length.
    loadMolfile(text, { record = true } = {}) {
      const L = text.split(/\r?\n/);
      const na = parseInt(L[3].slice(0, 3), 10), nb = parseInt(L[3].slice(3, 6), 10);
      if (!(na >= 0) || !(nb >= 0)) throw new Error('not a V2000 molfile');
      const before = this.snapshot();
      const raw = [];
      for (let i = 0; i < na; i++) {
        const l = L[4 + i];
        raw.push({ x: parseFloat(l.slice(0, 10)), y: parseFloat(l.slice(10, 20)), el: l.slice(31, 34).trim() || 'C' });
      }
      const bonds = [];
      for (let i = 0; i < nb; i++) {
        const l = L[4 + na + i];
        bonds.push({ a: parseInt(l.slice(0, 3), 10) - 1, b: parseInt(l.slice(3, 6), 10) - 1, order: parseInt(l.slice(6, 9), 10), stereo: parseInt(l.slice(9, 12), 10) || 0 });
      }
      const charge = new Array(na).fill(0), iso = new Array(na).fill(0);
      for (let i = 4 + na + nb; i < L.length && !L[i].startsWith('M  END'); i++) {
        const tag = L[i].slice(0, 6);
        if (tag !== 'M  CHG' && tag !== 'M  ISO') continue;
        const w = L[i].slice(6).trim().split(/\s+/).map(Number);
        for (let k = 0; k < w[0]; k++) (tag === 'M  CHG' ? charge : iso)[w[1 + 2 * k] - 1] = w[2 + 2 * k];
      }
      let sum = 0;
      for (const b of bonds) sum += Math.hypot(raw[b.a].x - raw[b.b].x, raw[b.a].y - raw[b.b].y);
      const scale = bonds.length ? BOND / (sum / bonds.length) : BOND / 1.5;
      this.atoms = [];
      this.bonds = [];
      this.nextId = 1;
      const ids = raw.map((r, i) => {
        const a = this.addAtom(r.x * scale, -r.y * scale, r.el);
        a.charge = charge[i];
        a.iso = iso[i];
        return a.id;
      });
      for (const b of bonds) {
        if (b.order < 1 || b.order > 3) b.order = 1;  // aromatic (4) is not written by the engine
        this.addBond(ids[b.a], ids[b.b], b.order, b.stereo === 1 || b.stereo === 6 ? b.stereo : 0);
      }
      this.selection.clear();
      this.fit();
      if (record) this.commit(before); else this.changed();
    }

    // ---- view ----------------------------------------------------------------------------------
    fit() {
      const r = this.host.getBoundingClientRect();
      const w = Math.max(200, r.width), h = Math.max(200, r.height);
      if (!this.atoms.length) { this.view = { x: -w / 2, y: -h / 2, scale: 1 }; return; }
      const xs = this.atoms.map((a) => a.x), ys = this.atoms.map((a) => a.y);
      const bw = Math.max(...xs) - Math.min(...xs) + 3 * BOND, bh = Math.max(...ys) - Math.min(...ys) + 3 * BOND;
      const scale = Math.min(1.6, Math.max(0.3, Math.min(w / bw, h / bh)));
      const cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
      this.view = { x: cx - w / (2 * scale), y: cy - h / (2 * scale), scale };
    }
    toWorld(ev) {
      const r = this.svg.getBoundingClientRect();
      return { x: this.view.x + (ev.clientX - r.left) / this.view.scale, y: this.view.y + (ev.clientY - r.top) / this.view.scale };
    }

    // ---- input ---------------------------------------------------------------------------------
    bind() {
      const s = this.svg;
      s.addEventListener('pointerdown', (ev) => this.down(ev));
      s.addEventListener('pointermove', (ev) => this.move(ev));
      s.addEventListener('pointerup', (ev) => this.up(ev));
      s.addEventListener('pointerleave', () => { this.hover = null; this.render(); });
      s.addEventListener('wheel', (ev) => {
        ev.preventDefault();
        if (ev.ctrlKey || ev.metaKey) {
          const p = this.toWorld(ev);
          const k = Math.exp(-ev.deltaY * 0.01);
          const scale = Math.min(4, Math.max(0.2, this.view.scale * k));
          this.view.x = p.x - (p.x - this.view.x) * (this.view.scale / scale);
          this.view.y = p.y - (p.y - this.view.y) * (this.view.scale / scale);
          this.view.scale = scale;
        } else {
          this.view.x += ev.deltaX / this.view.scale;
          this.view.y += ev.deltaY / this.view.scale;
        }
        this.render();
      }, { passive: false });
      s.addEventListener('keydown', (ev) => this.key(ev));
    }

    hit(p) {
      const a = this.atomNear(p);
      if (a) return { type: 'atom', id: a.id };
      const b = this.bondNear(p);
      if (b) return { type: 'bond', id: b.id };
      return null;
    }

    down(ev) {
      if (ev.button !== 0 && ev.pointerType === 'mouse') {
        if (ev.button === 1) this.drag = { pan: true, start: { x: ev.clientX, y: ev.clientY }, view: { ...this.view } };
        return;
      }
      this.svg.focus({ preventScroll: true });
      this.svg.setPointerCapture(ev.pointerId);
      const p = this.toWorld(ev);
      const h = this.hit(p);
      this.drag = { start: p, hit: h, moved: false, before: this.snapshot(), screen: { x: ev.clientX, y: ev.clientY } };
      if (this.tool.kind === 'select') {
        if (h && h.type === 'atom') {
          if (!this.selection.has(h.id)) {
            if (!ev.shiftKey) this.selection.clear();
            this.selection.add(h.id);
          }
          this.drag.moveAtoms = [...this.selection].map((id) => ({ a: this.atom(id), x: this.atom(id).x, y: this.atom(id).y }));
        } else if (h && h.type === 'bond') {
          const b = this.bonds.find((x) => x.id === h.id);
          if (!ev.shiftKey) this.selection.clear();
          this.selection.add(b.a).add(b.b);
          this.drag.moveAtoms = [...this.selection].map((id) => ({ a: this.atom(id), x: this.atom(id).x, y: this.atom(id).y }));
        } else {
          if (!ev.shiftKey) this.selection.clear();
          this.drag.box = { x0: p.x, y0: p.y, x1: p.x, y1: p.y };
        }
        this.render();
      }
    }

    move(ev) {
      const p = this.toWorld(ev);
      if (this.drag && this.drag.pan) {
        this.view.x = this.drag.view.x - (ev.clientX - this.drag.start.x) / this.view.scale;
        this.view.y = this.drag.view.y - (ev.clientY - this.drag.start.y) / this.view.scale;
        this.render();
        return;
      }
      if (!this.drag) {
        const h = this.hit(p);
        if (JSON.stringify(h) !== JSON.stringify(this.hover)) { this.hover = h; this.render(); }
        return;
      }
      const d = this.drag;
      if (Math.hypot(ev.clientX - d.screen.x, ev.clientY - d.screen.y) > 4) d.moved = true;
      if (!d.moved) return;
      if (this.tool.kind === 'select') {
        if (d.moveAtoms) {
          for (const m of d.moveAtoms) { m.a.x = m.x + p.x - d.start.x; m.a.y = m.y + p.y - d.start.y; }
        } else if (d.box) {
          d.box.x1 = p.x; d.box.y1 = p.y;
          const [x0, x1] = [Math.min(d.box.x0, d.box.x1), Math.max(d.box.x0, d.box.x1)];
          const [y0, y1] = [Math.min(d.box.y0, d.box.y1), Math.max(d.box.y0, d.box.y1)];
          this.selection = new Set(this.atoms.filter((a) => a.x >= x0 && a.x <= x1 && a.y >= y0 && a.y <= y1).map((a) => a.id));
        }
        this.render();
        return;
      }
      if (this.tool.kind === 'bond' || this.tool.kind === 'chain') {
        d.preview = this.previewDrag(p);
        this.render();
      }
    }

    // Where a drag from the pressed atom (or empty point) would draw to.
    previewDrag(p) {
      const d = this.drag;
      const from = d.hit && d.hit.type === 'atom' ? this.atom(d.hit.id) : d.start;
      const ang = snap(angle(from, p));
      if (this.tool.kind === 'chain') {
        const n = Math.max(1, Math.round(dist(from, p) / (BOND * 0.87)));
        const pts = [];
        let cur = { x: from.x, y: from.y };
        for (let k = 0; k < n; k++) {
          const turn = (k % 2 === 0 ? 1 : -1) * (Math.PI / 6);
          cur = { x: cur.x + BOND * Math.cos(ang + turn), y: cur.y + BOND * Math.sin(ang + turn) };
          pts.push(cur);
        }
        return { from, pts };
      }
      const over = this.atomNear(p, BOND * 0.35, d.hit && d.hit.type === 'atom' ? d.hit.id : null);
      const to = over || { x: from.x + BOND * Math.cos(ang), y: from.y + BOND * Math.sin(ang) };
      return { from, pts: [to], target: over };
    }

    up(ev) {
      const d = this.drag;
      this.drag = null;
      if (!d) return;
      if (d.pan) return;
      const p = this.toWorld(ev);
      const t = this.tool;
      if (t.kind === 'select') {
        if (d.moveAtoms && d.moved) {
          for (const m of d.moveAtoms) {  // dropping an atom onto another merges them
            const other = this.atomNear(m.a, BOND * 0.25, m.a.id);
            if (other && d.moveAtoms.length === 1) this.mergeAtoms(m.a.id, other.id);
          }
          this.commit(d.before);
        } else this.render();
        return;
      }
      const h = d.hit;
      if (t.kind === 'bond' || t.kind === 'chain') {
        if (d.moved && d.preview) {
          const from = h && h.type === 'atom' ? this.atom(h.id) : this.addAtom(d.start.x, d.start.y);
          let prev = from;
          for (const q of d.preview.pts) {
            const target = (d.preview.target && q === d.preview.pts[0] && t.kind === 'bond') ? d.preview.target
              : (this.atomNear(q, BOND * 0.3, prev.id) || this.addAtom(q.x, q.y));
            this.addBond(prev.id, target.id, t.kind === 'chain' ? 1 : (t.stereo ? 1 : t.order), t.kind === 'chain' ? 0 : (t.stereo || 0));
            prev = target;
          }
        } else if (h && h.type === 'bond') {
          this.applyBondTool(this.bonds.find((b) => b.id === h.id));
          return;
        } else if (h && h.type === 'atom') {
          const a = this.atom(h.id);
          this.grow(a, this.freeAngle(a));
        } else {
          const a = this.addAtom(p.x, p.y);
          this.grow(a, -Math.PI / 6);
        }
        this.commit(d.before);
        return;
      }
      if (t.kind === 'element') {
        if (h && h.type === 'atom') { this.setElement(h.id, this.element); return; }
        if (!h) { const before = this.snapshot(); this.addAtom(p.x, p.y, this.element); this.commit(before); }
        return;
      }
      if (t.kind === 'charge') {
        if (h && h.type === 'atom') this.adjustCharge(h.id, t.delta);
        return;
      }
      if (t.kind === 'erase') {
        if (h && h.type === 'atom') { const before = this.snapshot(); this.deleteAtom(h.id); this.commit(before); }
        if (h && h.type === 'bond') { const before = this.snapshot(); this.bonds = this.bonds.filter((b) => b.id !== h.id); this.commit(before); }
        return;
      }
      if (t.kind === 'ring') {
        const before = this.snapshot();
        if (h && h.type === 'bond') {
          const b = this.bonds.find((x) => x.id === h.id);
          const A = this.atom(b.a), B = this.atom(b.b);
          this.ringOnEdge(A, B, t.size, t.aromatic, this.sideAwayFromCrowd(A, B));
        } else if (h && h.type === 'atom') {
          let a = this.atom(h.id);
          let dir = -Math.PI / 2;
          if (this.neighbours(a.id).length) {  // substituent ring on a new bond, pointing away from the molecule
            dir = this.freeAngle(a);
            const saved = this.tool;
            this.tool = { kind: 'bond', order: 1, stereo: 0 };
            a = this.grow(a, dir).target;
            this.tool = saved;
          }
          this.ringAtAtom(a, dir, t.size, t.aromatic);
        } else {
          const pa = { x: p.x - BOND / 2, y: p.y + BOND / (2 * Math.tan(Math.PI / t.size)) };
          const pb = { x: p.x + BOND / 2, y: pa.y };
          this.ringOnEdge(pa, pb, t.size, t.aromatic, -1);
        }
        this.commit(before);
      }
    }

    mergeAtoms(keep, drop) {
      for (const b of this.bonds) {
        if (b.a === keep && b.b === drop) continue;
        if (b.a === drop) b.a = keep;
        if (b.b === drop) b.b = keep;
      }
      this.bonds = this.bonds.filter((b, i, all) => b.a !== b.b && all.findIndex((x) => (x.a === b.a && x.b === b.b) || (x.a === b.b && x.b === b.a)) === i);
      this.atoms = this.atoms.filter((a) => a.id !== drop);
      this.selection.delete(drop);
    }

    key(ev) {
      const mod = ev.metaKey || ev.ctrlKey;
      if (mod && ev.key.toLowerCase() === 'z') { ev.preventDefault(); ev.shiftKey ? this.redo() : this.undo(); return; }
      if (mod && ev.key.toLowerCase() === 'y') { ev.preventDefault(); this.redo(); return; }
      if (mod && ev.key.toLowerCase() === 'a') { ev.preventDefault(); this.selection = new Set(this.atoms.map((a) => a.id)); this.render(); return; }
      if (mod) return;
      if (ev.key === 'Delete' || ev.key === 'Backspace') {
        ev.preventDefault();
        if (this.selection.size) { this.deleteSelection(); return; }
        if (this.hover && this.hover.type === 'atom') { const b = this.snapshot(); this.deleteAtom(this.hover.id); this.hover = null; this.commit(b); }
        else if (this.hover && this.hover.type === 'bond') { const b = this.snapshot(); this.bonds = this.bonds.filter((x) => x.id !== this.hover.id); this.hover = null; this.commit(b); }
        return;
      }
      if (ev.key === 'Escape') { this.selection.clear(); this.render(); return; }
      const h = this.hover;
      if (!h) return;
      if (h.type === 'atom') {
        const el = KEY_ELEMENT[ev.key.toLowerCase()];
        if (el) { ev.preventDefault(); this.setElement(h.id, el); return; }
        if (ev.key === '+' || ev.key === '=') { ev.preventDefault(); this.adjustCharge(h.id, 1); return; }
        if (ev.key === '-') { ev.preventDefault(); this.adjustCharge(h.id, -1); return; }
      }
      if (h.type === 'bond' && ['1', '2', '3'].includes(ev.key)) {
        ev.preventDefault();
        const before = this.snapshot();
        const b = this.bonds.find((x) => x.id === h.id);
        b.order = Number(ev.key);
        b.stereo = 0;
        this.commit(before);
      }
    }

    // ---- rendering ---------------------------------------------------------------------------
    label(a) {
      const deg = this.neighbours(a.id).length;
      if (a.el === 'C' && deg > 0 && !a.charge && !a.iso) return null;
      const h = implicitH(a.el, a.charge, this.bondSum(a.id));
      let s = (a.iso ? String(a.iso) : '') + a.el;
      if (h) s += 'H' + (h > 1 ? h : '');
      return s;
    }
    render() {
      const r = this.host.getBoundingClientRect();
      const w = Math.max(1, r.width), hgt = Math.max(1, r.height);
      const v = this.view;
      this.svg.setAttribute('viewBox', `${v.x} ${v.y} ${w / v.scale} ${hgt / v.scale}`);
      const parts = [];
      const el = (tag, attrs, text) => {
        const n = document.createElementNS(SVGNS, tag);
        for (const [k, val] of Object.entries(attrs)) n.setAttribute(k, val);
        if (text != null) n.textContent = text;
        parts.push(n);
        return n;
      };
      const labels = new Map(this.atoms.map((a) => [a.id, this.label(a)]));
      // selection and error halos
      for (const a of this.atoms) {
        if (this.errors.has(a.id)) el('circle', { cx: a.x, cy: a.y, r: 13, class: 'ed-error' });
        if (this.selection.has(a.id)) el('circle', { cx: a.x, cy: a.y, r: 11, class: 'ed-selected' });
      }
      // ring centres for inner double-bond lines
      const ringSide = (b) => {
        const A = this.atom(b.a), B = this.atom(b.b);
        const common = this.neighbours(b.a).filter((x) => x !== b.b).flatMap((x) => this.neighbours(x))
          .filter((x) => this.neighbours(b.b).includes(x));
        let s = 0;
        for (const id of [...this.neighbours(b.a).filter((x) => x !== b.b), ...this.neighbours(b.b).filter((x) => x !== b.a)]) {
          const q = this.atom(id);
          s += Math.sign((B.x - A.x) * (q.y - A.y) - (B.y - A.y) * (q.x - A.x)) * (common.length ? 2 : 1);
        }
        return s;
      };
      for (const b of this.bonds) {
        const A = this.atom(b.a), B = this.atom(b.b);
        if (!A || !B) continue;
        const len = dist(A, B) || 1;
        const ux = (B.x - A.x) / len, uy = (B.y - A.y) / len;
        const ga = labels.get(A.id) ? 9 : 0, gb = labels.get(B.id) ? 9 : 0;
        const x1 = A.x + ux * ga, y1 = A.y + uy * ga, x2 = B.x - ux * gb, y2 = B.y - uy * gb;
        const cls = 'ed-bond' + (this.hover && this.hover.type === 'bond' && this.hover.id === b.id ? ' hover' : '');
        if (b.stereo === 1) {
          el('polygon', { points: `${x1},${y1} ${x2 - uy * 4.5},${y2 + ux * 4.5} ${x2 + uy * 4.5},${y2 - ux * 4.5}`, class: cls + ' wedge' });
        } else if (b.stereo === 6) {
          const n = 7;
          for (let k = 1; k <= n; k++) {
            const t = k / n, half = 4.5 * t;
            const cx = x1 + (x2 - x1) * t, cy = y1 + (y2 - y1) * t;
            el('line', { x1: cx - uy * half, y1: cy + ux * half, x2: cx + uy * half, y2: cy - ux * half, class: cls });
          }
        } else if (b.order === 1) {
          el('line', { x1, y1, x2, y2, class: cls });
        } else if (b.order === 2) {
          const side = ringSide(b);
          if (side === 0 || labels.get(A.id) || labels.get(B.id)) {
            for (const o of [-2.6, 2.6]) el('line', { x1: x1 - uy * o, y1: y1 + ux * o, x2: x2 - uy * o, y2: y2 + ux * o, class: cls });
          } else {
            const o = 5.5 * Math.sign(side);
            el('line', { x1, y1, x2, y2, class: cls });
            el('line', { x1: x1 + ux * 5 - uy * o, y1: y1 + uy * 5 + ux * o, x2: x2 - ux * 5 - uy * o, y2: y2 - uy * 5 + ux * o, class: cls });
          }
        } else {
          for (const o of [-4.5, 0, 4.5]) el('line', { x1: x1 - uy * o, y1: y1 + ux * o, x2: x2 - uy * o, y2: y2 + ux * o, class: cls });
        }
      }
      for (const a of this.atoms) {
        const lab = labels.get(a.id);
        const hov = this.hover && this.hover.type === 'atom' && this.hover.id === a.id;
        if (hov) el('circle', { cx: a.x, cy: a.y, r: 10, class: 'ed-hover' });
        if (lab) {
          el('circle', { cx: a.x, cy: a.y, r: 8.5, class: 'ed-label-bg' });
          // hydrogens on the side away from the bonds
          const nb = this.neighbours(a.id).map((id) => this.atom(id));
          const left = nb.length && nb.reduce((s, q) => s + (q.x - a.x), 0) > 0;
          const t = el('text', { x: a.x, y: a.y + 5, class: 'ed-label', 'text-anchor': left ? 'end' : 'middle', fill: COLOR[a.el] || 'currentColor' }, lab);
          if (left) t.setAttribute('x', a.x + 5.5);
        }
        if (a.charge) {
          const sign = a.charge > 0 ? '+' : '−';
          el('text', { x: a.x + 8, y: a.y - 6, class: 'ed-charge' }, (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + sign);
        }
      }
      // drag previews
      const d = this.drag;
      if (d && d.preview) {
        let prev = d.preview.from;
        for (const q of d.preview.pts) { el('line', { x1: prev.x, y1: prev.y, x2: q.x, y2: q.y, class: 'ed-preview' }); prev = q; }
      }
      if (d && d.box) {
        el('rect', { x: Math.min(d.box.x0, d.box.x1), y: Math.min(d.box.y0, d.box.y1), width: Math.abs(d.box.x1 - d.box.x0), height: Math.abs(d.box.y1 - d.box.y0), class: 'ed-box' });
      }
      if (!this.atoms.length) {
        el('text', { x: v.x + w / (2 * v.scale), y: v.y + hgt / (2 * v.scale), class: 'ed-hint', 'text-anchor': 'middle' }, 'Click or drag to draw a bond');
      }
      this.svg.replaceChildren(...parts);
    }

    setErrors(atomNumbers) {  // 1-based atom numbers from engine diagnostics
      this.errors = new Set(atomNumbers.map((k) => this.atoms[k - 1] && this.atoms[k - 1].id).filter(Boolean));
      this.render();
    }
  }

  window.StructureEditor = StructureEditor;
  window.StructureEditor.implicitH = implicitH;
})();
