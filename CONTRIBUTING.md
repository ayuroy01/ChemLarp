# Contributing

Thank you for helping. Bug reports, chemistry corrections and test cases with references are the most useful
contributions: open an [issue](https://github.com/ayuroy01/ChemLarp/issues) with the input (a SMILES string or
file), what you got, what you expected, and a source if you have one.

## Code contributions

Before writing code, please open an issue to discuss the change.

1. Follow the existing style and architecture (see [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)).
2. Add tests (`tests/`), and validation against an independent reference where the change affects results.
3. `make -j8 && build/chemlarp_tests` must pass.
4. Every new source file starts with the licence header used in the other files.

## Licence of contributions (required)

ChemLarp is licensed under the AGPL-3.0 and also offered under commercial licences (see
[LICENSING.md](LICENSING.md)). For that to remain possible, each contribution needs these terms.

By submitting a contribution (a pull request, patch or code in an issue), you agree that:

1. you wrote it, or have the right to submit it, and it does not include code copied from elsewhere without its
   licence being compatible and stated;
2. it is licensed to everyone under the AGPL-3.0-only; and
3. you grant the ChemLarp copyright holder a perpetual, worldwide, non-exclusive, royalty-free, irrevocable licence to
   use, modify, sublicense and distribute your contribution under any licence, including commercial licences.
   You keep the copyright of your contribution.

Sign off each commit to confirm this (`git commit -s`), which adds a line like:

```
Signed-off-by: Your Name <you@example.com>
```

Pull requests whose commits are not signed off cannot be merged.
