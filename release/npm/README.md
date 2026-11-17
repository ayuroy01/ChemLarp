# chemlarp

The ChemLarp chemistry engine for Node.js, compiled to WebAssembly: nothing to compile, the same on every
operating system.

```js
const chemlarp = await require('chemlarp').load();

chemlarp.canonical('C1=CC=CC=C1').canonical;               // 'c1ccccc1'
chemlarp.inchi('CCO').key;                                  // 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
chemlarp.iupacName('CC(C)Cc1ccc(cc1)C(C)C(=O)O').name;      // '2-[4-(2-methylpropyl)phenyl]propanoic acid'
chemlarp.embed('C[C@@H](C(=O)O)N').coordinates;             // 3D coordinates, stereo kept
chemlarp.qc('O', { method: 'b3lyp', basis: '6-31G*' }).energy;
```

Every call returns the JSON object of the C API, `{ ok, ..., diagnostics }`; when `ok` is false, the
diagnostics say what was wrong and where. The fields of each result are documented in
[include/chemlarp/chemlarp_c.h](https://github.com/ayuroy01/ChemLarp/blob/main/include/chemlarp/chemlarp_c.h).

The browser app with the same engine: https://ayuroy01.github.io/ChemLarp/

Licence: AGPL-3.0-only (commercial licences available; see LICENSING.md in the repository).
