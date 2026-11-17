// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
'use strict';
// Smoke test of the packed npm module (run from an installation of the tarball, or from build/npm).
const assert = require('assert');
require('./index.js').load().then((c) => {
  assert.strictEqual(c.version(), require('./package.json').version);
  assert.strictEqual(c.canonical('C1=CC=CC=C1').canonical, 'c1ccccc1');
  assert.strictEqual(c.inchi('CCO').key, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N');
  assert.strictEqual(c.iupacName('CC(C)Cc1ccc(cc1)C(C)C(=O)O').name, '2-[4-(2-methylpropyl)phenyl]propanoic acid');
  const e = c.embed('C[C@@H](C(=O)O)N');
  assert.ok(e.ok && e.coordinates.length === 13);
  const m = c.mbar([[0, 0.5, 2, 0.1], [0.2, 1.0, 0.1, 0]], [2, 2]);
  assert.ok(m.ok && Math.abs(m.f[1] - -0.2639) < 1e-3);
  assert.ok(!c.parse('C1CC').ok);
  console.log('chemlarp', c.version(), 'npm package: all checks passed');
}).catch((err) => { console.error(err); process.exit(1); });
