# Docking reference complex

One protein–ligand complex, so that the redocking test can run anywhere the repository is checked out —
including in CI, under the sanitizers — without downloading anything.

| File | What it is |
|---|---|
| `3PTB.pdb` | PDB entry 3PTB, bovine β-trypsin with benzamidine bound (Bode & Schwager, *J. Mol. Biol.* 98 (1975) 693), as deposited: `files.rcsb.org/download/3PTB.pdb` |
| `3PTB_BEN.sdf` | the benzamidine (chemical component BEN, residue A 1) of that entry at its crystal coordinates, from the RCSB model server with `{"label_comp_id":"BEN"}`: `models.rcsb.org/v1/3ptb/atoms?encoding=sdf&copy_all_categories=false` — the SDF keeps the query it was made with in its data items |

Protein Data Bank data carry no restrictions (CC0); see docs/AUDITS.md, "Bundled data provenance".

`tests/test_docking.cpp` (`docking_redocks_benzamidine_into_trypsin`) reads both: it takes only the ligand's
canonical SMILES from the SDF, generates a conformation from scratch, docks it into the structure with the
ligand's own site cleared, and requires the top-ranked pose within 2 Å of the deposited one. The crystal
coordinates are the reference for the RMSD and never an input to the search.

The other eleven complexes of the redocking set are not in the repository. `scripts/validate_docking.py`
names them with the URLs to fetch them from, and `scripts/validate.py` looks for them under its `dock_dir`
data path. To run the whole set against this directory plus downloads:

```sh
python3 scripts/validate_docking.py <directory with all twelve> --runs 16
```
