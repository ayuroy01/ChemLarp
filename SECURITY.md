# Security policy

## Supported versions

The latest release and the `main` branch receive security fixes. Older versions do not.

| Version | Supported |
|---|---|
| 0.1.x | yes |
| earlier | no |

## Reporting a vulnerability

Please **do not open a public issue** for a security problem.

Report it privately through GitHub: go to the
[Security tab](https://github.com/ayuroy01/ChemLarp/security/advisories/new) and open a draft advisory
("Report a vulnerability"). That page is private between you and the maintainer.

Please include:

- what the problem is and how serious you think it is;
- the exact input that triggers it (a SMILES string, a file, a URL), small enough to reproduce quickly;
- the version (`chemlarp --version`, or the commit) and the platform;
- a stack trace or sanitizer output if you have one.

**What to expect.** This is a one-person project, so replies are best effort: an acknowledgement within about
a week, and an assessment within about a month. Fixed problems are credited to the reporter in the release
notes unless you ask otherwise. There is no bug-bounty programme.

## What counts as a vulnerability

ChemLarp is a local program and a web app that runs entirely in your browser. There is no ChemLarp server, no
account and no telemetry, so the interesting attack surface is **untrusted input**:

- a crash, hang, out-of-bounds access or excessive memory use when reading a structure file (MOL, SDF, PDB,
  mmCIF, CML), a mechanism (`.ocm`, CHEMKIN, Cantera YAML), a flowsheet, a lab or aqueous input, a SMILES or
  SMARTS string, or JSON passed to the C API;
- script injection in the web app from content in a file the user opens (names, titles, headers);
- the local servers in the desktop and Android apps serving files outside the bundle, or binding to anything
  other than `127.0.0.1`;
- a packaging or build problem that lets a release ship something other than the published source.

**Not vulnerabilities:** wrong chemistry or numerical results (report those as ordinary issues, they are very
welcome), running out of memory on a genuinely huge input, missing code signing on the apps (known, see
docs/AUDITS.md), and anything requiring the attacker to already control the machine.

## Verifying a release

Every release lists `SHA256SUMS`. Check a download against it:

```bash
shasum -a 256 -c SHA256SUMS --ignore-missing
```

Releases are not yet signed with a code-signing certificate or Sigstore; this is tracked in
[docs/AUDITS.md](docs/AUDITS.md) (audit 2.8).
