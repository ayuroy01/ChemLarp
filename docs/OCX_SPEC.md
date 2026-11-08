# OmniChem Exchange (`.ocx`) — v0.1

An `.ocx` file is a ZIP container. Its first entry is `manifest.json`, followed by one or more payloads. Payloads use JSON in v0.1. A Protocol Buffers encoding of the same schemas is planned as `encoding: "protobuf"`; readers must check the `encoding` field rather than infer it from the file extension.

## manifest.json

```json
{
  "format": "ocx",
  "version": "0.1",
  "created": "2026-11-08T12:00:00+00:00",
  "generator": "omnichem-core 0.1.0",
  "payloads": [
    {"path": "structure.json", "kind": "structure", "encoding": "json", "sha256": "…"},
    {"path": "thermo.json",    "kind": "thermo",    "encoding": "json", "sha256": "…"}
  ]
}
```

Readers must verify each payload's `sha256` and must ignore unknown payload kinds.

## structure.json

This payload mirrors `ChemicalSystemState` one-to-one. Columns are stored structure-of-arrays, so a reader can map them straight into typed arrays.

| field | type | notes |
|---|---|---|
| `atoms.count` | int | N |
| `atoms.element` | u8[N] | atomic number |
| `atoms.formal_charge` | i8[N] | |
| `atoms.implicit_h` | u8[N] | |
| `atoms.flags` | u8[N] | bit0 aromatic, bit1 bracket, bit2 in-ring |
| `atoms.isotope` | u16[N] | 0 = natural abundance |
| `atoms.chirality` | i8[N] | 0 none, 1 `@`, 2 `@@` (relative to neighbour order in `bonds`) |
| `atoms.atom_map` | u16[N] | SMILES atom class |
| `atoms.species` | u32[N] | owning species index |
| `atoms.x/y/z` | f64[N] | Å; all zeros until 3-D embedding exists |
| `bonds.count` | int | M |
| `bonds.begin/end` | u32[M] | atom indices |
| `bonds.order` | u8[M] | Kekulé order 1–4 |
| `bonds.flags` | u8[M] | bit0 aromatic, bit1 ring-closure, bit2 `/`, bit3 `\`, bit4 in-ring |
| `bonds.species` | u32[M] | |
| `species[]` | object | `smiles, canonical_smiles, formula, atom_begin, atom_count, bond_begin, bond_count, net_charge, molar_mass, phase, composition{symbol: count}` |
| `stereo.tetrahedral[]` | object | `atom, refs[4]` (atom ids, `"H"` implicit H, `"LP"` lone pair), `parity` `"@"`/`"@@"` (viewed from `refs[0]`), `cip` `"R"`/`"S"`/`""` |
| `stereo.double_bonds[]` | object | `bond, ref_begin, ref_end, cis` (refs are neighbours of `bonds.begin/end`), `cip` `"E"`/`"Z"`/`""` |

Invariant: each species owns the contiguous ranges `[atom_begin, atom_begin+atom_count)` and `[bond_begin, bond_begin+bond_count)`. `scripts/ocx_pack.py` validates the column lengths and bond index ranges.

## thermo.json

```json
{"temperature_K": 1000, "pressure_Pa": 100000, "moles": [..S..], "g_RT": [..S..]}
```

Indices in `moles` and `g_RT` match `structure.species`. `g_RT` is null for a species with no thermo bound.

## process.json (reserved)

This payload will hold the flowsheet topology (unit operations, streams referencing species indices, specs). It is not defined yet.
