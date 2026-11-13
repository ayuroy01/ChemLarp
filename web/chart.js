// Small canvas line and bar charts that follow the page's colour tokens.
(function () {
  const PALETTE = ['#14837a', '#d1495b', '#3d6fd8', '#e0a100', '#7b52c4', '#2a9d3f', '#c45a1a', '#4b8fb3', '#a03f8f'];

  function token(name, fallback) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
  }

  function niceTicks(lo, hi, n) {
    if (!(hi > lo)) { hi = lo + 1; }
    const step0 = (hi - lo) / Math.max(1, n);
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= step0) || 10 * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
    return out;
  }

  function fmt(v) {
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e4 || a < 1e-3) return v.toExponential(0).replace('e+', 'e');
    return +v.toPrecision(3) + '';
  }

  // The drawn height comes from data-height (CSS pixels); the canvas's own width/height
  // attributes are the device-pixel buffer and are rewritten on every draw.
  function setup(canvas) {
    const dpr = window.devicePixelRatio || 1;
    const h = parseFloat(canvas.dataset.height) || 240;
    canvas.style.height = h + 'px';
    const w = canvas.clientWidth;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx, w, h };
  }

  // series: [{name, points:[[x,y],...]}]; opts: {xLabel, yLabel, logY}
  function lineChart(canvas, series, opts = {}) {
    canvas.hidden = false;
    const { ctx, w, h } = setup(canvas);
    const fg = token('--ink', '#111'), muted = token('--ink-3', '#666'), line = token('--hairline', '#ddd');
    const tf = (y) => (opts.logY ? Math.log10(Math.max(y, 1e-12)) : y);
    const xs = series.flatMap((s) => s.points.map((p) => p[0]));
    const ys = series.flatMap((s) => s.points.map((p) => tf(p[1])));
    let x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
    if (opts.xRange) [x0, x1] = opts.xRange;
    if (opts.logY) { y0 = Math.max(Math.floor(y0), Math.ceil(y1) - 10); y1 = Math.ceil(y1); }
    else { const pad = (y1 - y0) * 0.06 || 1; y0 -= pad; y1 += pad; }
    const yt = opts.logY ? niceTicks(y0, y1, 5).filter((v) => Number.isInteger(v)) : niceTicks(y0, y1, 5);
    ctx.font = '12px system-ui, sans-serif';
    const yLab = (v) => (opts.logY ? '1e' + v : fmt(v));
    const left = 14 + Math.max(...yt.map((v) => ctx.measureText(yLab(v)).width)) + (opts.yLabel ? 16 : 0);
    // Legend entries are sized by their text; count the rows they wrap into.
    const legendW = series.map((s) => 16 + ctx.measureText(s.name).width + 18);
    let legendRows = 1;
    for (let i = 0, lx = left; i < series.length; ++i) {
      if (lx + legendW[i] > w && lx > left) { legendRows++; lx = left; }
      lx += legendW[i];
    }
    const top = 10 + legendRows * 18, right = 14, bottom = 38;
    const X = (x) => left + ((x - x0) / (x1 - x0 || 1)) * (w - left - right);
    const Y = (y) => top + (1 - (tf(y) - y0) / (y1 - y0 || 1)) * (h - top - bottom);
    const Yraw = (v) => top + (1 - (v - y0) / (y1 - y0 || 1)) * (h - top - bottom);
    ctx.strokeStyle = line; ctx.lineWidth = 1; ctx.fillStyle = muted;
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const v of yt) {
      if (v < y0 || v > y1) continue;
      const py = Yraw(v);
      ctx.beginPath(); ctx.moveTo(left, py); ctx.lineTo(w - right, py); ctx.stroke();
      ctx.fillText(yLab(v), left - 6, py);
    }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    const xt = niceTicks(x0, x1, Math.max(2, Math.floor((w - left) / 90)));
    for (const v of xt) {
      if (v < x0 || v > x1) continue;
      ctx.fillText(fmt(v), X(v), h - bottom + 6);
    }
    if (opts.xLabel) { ctx.fillText(opts.xLabel, left + (w - left - right) / 2, h - 16); }
    if (opts.yLabel) {
      ctx.save(); ctx.translate(12, top + (h - top - bottom) / 2); ctx.rotate(-Math.PI / 2);
      ctx.textBaseline = 'middle'; ctx.fillText(opts.yLabel, 0, 0); ctx.restore();
    }
    series.forEach((s, i) => {
      ctx.strokeStyle = s.color || PALETTE[i % PALETTE.length];
      ctx.lineWidth = 2;
      ctx.beginPath();
      s.points.forEach(([x, y], k) => (k ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y))));
      ctx.stroke();
    });
    // Markers: [{x, y, label}] as dots with a label.
    for (const m of opts.markers || []) {
      ctx.fillStyle = fg;
      ctx.beginPath(); ctx.arc(X(m.x), Y(m.y), 4, 0, 2 * Math.PI); ctx.fill();
      if (m.label) {
        const right = X(m.x) > (left + w) / 2;
        ctx.textAlign = right ? 'right' : 'left'; ctx.textBaseline = 'bottom';
        ctx.fillText(m.label, X(m.x) + (right ? -8 : 8), Y(m.y) - 6);
      }
    }
    // Legend.
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillStyle = fg;
    let lx = left, ly = 8;
    series.forEach((s, i) => {
      if (lx + legendW[i] > w && lx > left) { lx = left; ly += 18; }
      ctx.fillStyle = s.color || PALETTE[i % PALETTE.length];
      ctx.fillRect(lx, ly - 4, 12, 8);
      ctx.fillStyle = fg;
      ctx.fillText(s.name, lx + 16, ly);
      lx += legendW[i];
    });
  }

  // bars: [{label, value}] drawn horizontally.
  function barChart(canvas, bars, opts = {}) {
    canvas.hidden = false;
    const rowH = 24;
    canvas.dataset.height = String(bars.length * rowH + 36);
    const { ctx, w, h } = setup(canvas);
    const fg = token('--ink', '#111'), muted = token('--ink-3', '#666'), accent = token('--accent', '#14534d');
    ctx.font = '13px ui-monospace, Menlo, monospace';
    const left = 16 + Math.max(...bars.map((b) => ctx.measureText(b.label).width));
    const right = 70, max = Math.max(...bars.map((b) => b.value), 1e-12);
    bars.forEach((b, i) => {
      const y = 10 + i * rowH;
      ctx.fillStyle = fg; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      ctx.fillText(b.label, left - 8, y + rowH / 2 - 2);
      const bw = (b.value / max) * (w - left - right);
      ctx.fillStyle = accent; ctx.fillRect(left, y + 3, Math.max(1, bw), rowH - 10);
      ctx.fillStyle = muted; ctx.textAlign = 'left'; ctx.font = '12px system-ui, sans-serif';
      ctx.fillText(opts.format ? opts.format(b.value) : fmt(b.value), left + bw + 6, y + rowH / 2 - 2);
      ctx.font = '13px ui-monospace, Menlo, monospace';
    });
    ctx.fillStyle = muted; ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'left';
    if (opts.caption) ctx.fillText(opts.caption, left, h - 12);
  }

  window.Charts = { lineChart, barChart, palette: PALETTE };
})();
