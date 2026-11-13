// Minimal ball-and-stick viewer on a 2D canvas: orthographic projection,
// painter's-algorithm depth sort, drag to rotate, wheel/pinch to zoom.
(function () {
  const CPK = {
    H: '#e8e8e8', C: '#5a5f63', N: '#3060e0', O: '#e0302a', F: '#70d040', Cl: '#30c030', Br: '#a52a2a', I: '#8a3fb0',
    S: '#e8c020', P: '#ff8000', B: '#e0a080', Si: '#c09060',
  };
  const RADIUS = { H: 0.25, C: 0.38, N: 0.36, O: 0.35, F: 0.33, Cl: 0.45, Br: 0.5, I: 0.56, S: 0.48, P: 0.48, B: 0.38, Si: 0.5 };

  // Parse a V2000 molfile into {atoms:[{el,x,y,z}], bonds:[{a,b,order}]}.
  function parseSdf(text) {
    const lines = text.split(/\r?\n/);
    const counts = lines[3] || '';
    const na = parseInt(counts.slice(0, 3), 10), nb = parseInt(counts.slice(3, 6), 10);
    const atoms = [], bonds = [];
    for (let i = 0; i < na; i++) {
      const l = lines[4 + i];
      atoms.push({ x: +l.slice(0, 10), y: +l.slice(10, 20), z: +l.slice(20, 30), el: l.slice(31, 34).trim() });
    }
    for (let i = 0; i < nb; i++) {
      const l = lines[4 + na + i];
      bonds.push({ a: parseInt(l.slice(0, 3), 10) - 1, b: parseInt(l.slice(3, 6), 10) - 1, order: parseInt(l.slice(6, 9), 10) });
    }
    return { atoms, bonds };
  }

  const mul = (A, B) => {
    const C = new Array(9).fill(0);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) C[3 * i + j] += A[3 * i + k] * B[3 * k + j];
    return C;
  };
  const rotX = (t) => [1, 0, 0, 0, Math.cos(t), -Math.sin(t), 0, Math.sin(t), Math.cos(t)];
  const rotY = (t) => [Math.cos(t), 0, Math.sin(t), 0, 1, 0, -Math.sin(t), 0, Math.cos(t)];

  class MoleculeViewer {
    constructor(canvas) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.mol = null;
      this.R = [1, 0, 0, 0, 1, 0, 0, 0, 1];
      this.zoom = 1;
      this.labels = false;
      this.pointers = new Map();
      this.bind();
      new ResizeObserver(() => this.draw()).observe(canvas);
    }

    load(sdf) { this.loadAtoms(parseSdf(sdf)); }

    // m: {atoms:[{el,x,y,z}], bonds:[{a,b,order}]}; centred on its geometric centre.
    loadAtoms(m) {
      const c = m.atoms.reduce((s, a) => [s[0] + a.x, s[1] + a.y, s[2] + a.z], [0, 0, 0]).map((v) => v / m.atoms.length);
      m.atoms.forEach((a) => { a.x -= c[0]; a.y -= c[1]; a.z -= c[2]; });
      this.centre = c;
      this.extent = Math.max(1.5, ...m.atoms.map((a) => Math.hypot(a.x, a.y, a.z) + 0.6));
      this.mol = m;
      this.R = mul(rotX(-0.35), rotY(0.5));
      this.zoom = 1;
      this.draw();
    }

    // New coordinates (flat x0 y0 z0 x1 ...) for the loaded atoms, keeping the view: trajectory frames.
    setPositions(xyz) {
      if (!this.mol || xyz.length !== 3 * this.mol.atoms.length) return;
      const c = this.centre || [0, 0, 0];
      this.mol.atoms.forEach((a, i) => { a.x = xyz[3 * i] - c[0]; a.y = xyz[3 * i + 1] - c[1]; a.z = xyz[3 * i + 2] - c[2]; });
      this.draw();
    }

    bind() {
      const c = this.canvas;
      c.addEventListener('pointerdown', (e) => { c.setPointerCapture(e.pointerId); this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY }); });
      const end = (e) => { this.pointers.delete(e.pointerId); this.pinch = null; };
      c.addEventListener('pointerup', end);
      c.addEventListener('pointercancel', end);
      c.addEventListener('pointermove', (e) => {
        if (!this.pointers.has(e.pointerId)) return;
        const prev = this.pointers.get(e.pointerId);
        this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
        if (this.pointers.size === 2) {
          const [p, q] = [...this.pointers.values()];
          const d = Math.hypot(p.x - q.x, p.y - q.y);
          if (this.pinch) this.zoom = Math.min(6, Math.max(0.3, this.zoom * d / this.pinch));
          this.pinch = d;
        } else {
          const dx = e.clientX - prev.x, dy = e.clientY - prev.y;
          this.R = mul(mul(rotX(dy * 0.01), rotY(dx * 0.01)), this.R);
        }
        this.draw();
      });
      c.addEventListener('wheel', (e) => {
        e.preventDefault();
        this.zoom = Math.min(6, Math.max(0.3, this.zoom * Math.exp(-e.deltaY * 0.0015)));
        this.draw();
      }, { passive: false });
    }

    draw() {
      const c = this.canvas, ctx = this.ctx;
      const dpr = window.devicePixelRatio || 1;
      const w = c.clientWidth, h = c.clientHeight;
      if (!w || !h) return;
      c.width = Math.round(w * dpr);
      c.height = Math.round(h * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, w, h);
      if (!this.mol) return;
      const style = getComputedStyle(document.documentElement);
      const fg = style.getPropertyValue('--fg').trim() || '#000';
      const scale = (Math.min(w, h) / 2 / this.extent) * this.zoom;
      const R = this.R;
      const P = this.mol.atoms.map((a) => ({
        sx: w / 2 + scale * (R[0] * a.x + R[1] * a.y + R[2] * a.z),
        sy: h / 2 - scale * (R[3] * a.x + R[4] * a.y + R[5] * a.z),
        z: R[6] * a.x + R[7] * a.y + R[8] * a.z,
        el: a.el,
      }));
      const items = [];
      this.mol.bonds.forEach((b) => items.push({ kind: 'bond', z: (P[b.a].z + P[b.b].z) / 2 - 0.01, b }));
      P.forEach((p, i) => items.push({ kind: 'atom', z: p.z, i }));
      items.sort((u, v) => u.z - v.z);
      ctx.lineCap = 'round';
      for (const it of items) {
        if (it.kind === 'bond') {
          const p = P[it.b.a], q = P[it.b.b];
          const dx = q.sx - p.sx, dy = q.sy - p.sy, len = Math.hypot(dx, dy) || 1;
          const nx = -dy / len, ny = dx / len;
          const n = Math.min(3, Math.max(1, it.b.order));
          const gap = 0.09 * scale;
          ctx.lineWidth = Math.max(1.5, 0.11 * scale / (n > 1 ? 1.6 : 1));
          for (let k = 0; k < n; k++) {
            const off = (k - (n - 1) / 2) * gap;
            const mx = (p.sx + q.sx) / 2 + nx * off, my = (p.sy + q.sy) / 2 + ny * off;
            ctx.strokeStyle = CPK[p.el] || '#999';
            ctx.beginPath(); ctx.moveTo(p.sx + nx * off, p.sy + ny * off); ctx.lineTo(mx, my); ctx.stroke();
            ctx.strokeStyle = CPK[q.el] || '#999';
            ctx.beginPath(); ctx.moveTo(mx, my); ctx.lineTo(q.sx + nx * off, q.sy + ny * off); ctx.stroke();
          }
        } else {
          const p = P[it.i];
          const r = (RADIUS[p.el] || 0.45) * scale * 0.75;
          const g = ctx.createRadialGradient(p.sx - r * 0.35, p.sy - r * 0.35, r * 0.1, p.sx, p.sy, r);
          g.addColorStop(0, '#ffffff');
          g.addColorStop(0.35, CPK[p.el] || '#c060c0');
          g.addColorStop(1, shade(CPK[p.el] || '#c060c0'));
          ctx.fillStyle = g;
          ctx.beginPath(); ctx.arc(p.sx, p.sy, r, 0, 2 * Math.PI); ctx.fill();
          if (this.labels && p.el !== 'H') {
            ctx.fillStyle = fg;
            ctx.font = `600 ${Math.max(10, r * 0.9)}px system-ui, sans-serif`;
            ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            ctx.fillText(p.el, p.sx, p.sy - r - 8);
          }
        }
      }
    }
  }

  function shade(hex) {
    const v = parseInt(hex.slice(1), 16);
    const f = (s) => Math.round(((v >> s) & 255) * 0.45);
    return `rgb(${f(16)},${f(8)},${f(0)})`;
  }

  window.MoleculeViewer = MoleculeViewer;
  window.parseSdf = parseSdf;
})();
