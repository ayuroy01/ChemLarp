#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
"""GHS hazard classifications of common laboratory chemicals from PubChem -> src/ghs_data.inc.

    python3 scripts/gen_ghs_data.py [--cache DIR] [--names FILE]

For each name in NAMES: the PubChem CID (PUG REST, name lookup), its SMILES, and the "GHS Classification"
section of the compound record (PUG-View), which PubChem aggregates from the notifications to the ECHA C&L
Inventory and other sources. Kept per compound: the signal word, the pictogram codes (GHS01-GHS09), and the
hazard statements (H-codes) with the share of notifications that give each one, from the first source
listed (usually the ECHA C&L Inventory summary). The H-statement texts are the standard UN GHS phrases as
PubChem prints them. Entries are keyed by ChemLarp's canonical SMILES of the PubChem structure.

Network use: about two requests per compound, at most four per second; responses are cached in --cache so a
re-run does not fetch again. Run from the repository after building build/chemlarp.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = "chemlarp-data-builder/0.1 (offline table generation)"

NAMES = """
water; hydrochloric acid; sulfuric acid; nitric acid; phosphoric acid; acetic acid; formic acid; citric acid;
oxalic acid; hydrofluoric acid; perchloric acid; trifluoroacetic acid; benzoic acid; boric acid; lactic acid;
sodium hydroxide; potassium hydroxide; ammonia; calcium hydroxide; sodium carbonate; sodium bicarbonate;
potassium carbonate; triethylamine; pyridine; sodium chloride; potassium chloride; calcium chloride;
magnesium sulfate; sodium sulfate; copper(II) sulfate; silver nitrate; potassium permanganate;
potassium dichromate; sodium hypochlorite; calcium hypochlorite; potassium iodide; sodium iodide;
potassium bromide; barium chloride; lead(II) nitrate; iron(III) chloride; iron(II) sulfate; zinc chloride;
ammonium chloride; ammonium nitrate; ammonium sulfate; sodium nitrite; sodium nitrate; potassium nitrate;
sodium azide; potassium cyanide; sodium cyanide; sodium sulfide; sodium thiosulfate; potassium chlorate;
sodium chlorate; sodium perchlorate; potassium hydrogen phthalate; sodium acetate; sodium phosphate;
disodium hydrogen phosphate; sodium dihydrogen phosphate; magnesium chloride; aluminium chloride;
cobalt(II) chloride; nickel(II) chloride; mercury(II) chloride; cadmium chloride; chromium(III) chloride;
manganese(II) sulfate; zinc sulfate; copper(II) chloride; lithium chloride; cesium chloride;
hydrogen peroxide; sodium peroxide; potassium persulfate; ammonium persulfate; iodine; bromine;
methanol; ethanol; 2-propanol; 1-butanol; acetone; acetonitrile; dichloromethane; chloroform;
carbon tetrachloride; diethyl ether; tetrahydrofuran; toluene; benzene; hexane; heptane; ethyl acetate;
N,N-dimethylformamide; dimethyl sulfoxide; 1,4-dioxane; p-xylene; cyclohexane; pentane; ethylene glycol;
glycerol; acetic anhydride; acetyl chloride; thionyl chloride; phosphorus pentachloride;
phosphoryl chloride; oxalyl chloride; benzoyl chloride; p-toluenesulfonyl chloride;
3-chloroperoxybenzoic acid; dicyclohexylcarbodiimide; sodium borohydride; lithium aluminium hydride;
n-butyllithium; sodium hydride; sodium; potassium; lithium; magnesium; zinc; iron; copper; aluminium;
phenol; aniline; formaldehyde; acetaldehyde; benzaldehyde; nitrobenzene; nitromethane; hydrazine;
hydroxylamine; urea; glycine; glucose; sucrose; caffeine; aspirin; acetaminophen; salicylic acid;
ethylenediaminetetraacetic acid; disodium EDTA; phenolphthalein; methyl orange; bromothymol blue;
carbon dioxide; chlorine; hydrogen sulfide; sulfur dioxide; carbon monoxide; hydrogen; oxygen; nitrogen;
helium; argon; methane; propane; butane; acetylene; ethylene; nitrogen dioxide; hydrogen cyanide;
nitric oxide; ozone; picric acid; ethidium bromide; acrylamide; sodium dodecyl sulfate; tris(hydroxymethyl)aminomethane;
2-mercaptoethanol; dithiothreitol; imidazole; guanidine hydrochloride; potassium ferricyanide;
potassium ferrocyanide; calcium carbonate; magnesium oxide; calcium oxide; silica gel; activated charcoal;
starch; sodium silicate; phosphorus pentoxide; sulfur; red phosphorus; white phosphorus;
benzyl alcohol; cyclohexanone; styrene; acrolein; acrylonitrile; vinyl chloride; ethylene oxide;
propylene oxide; epichlorohydrin; diethylamine; piperidine; morpholine; ethylenediamine; ethanolamine;
triethanolamine; dimethylamine; methylamine; trimethylamine; pyrrolidine; naphthalene; anthracene;
""".replace("\n", " ")


def get(url, cache, key):
    path = os.path.join(cache, re.sub(r"[^A-Za-z0-9_.-]+", "_", key))
    if os.path.exists(path):
        return open(path, "rb").read()
    # curl rather than urllib: it verifies TLS against the system certificate store, which some Python builds
    # (python.org's on macOS) do not use.
    for attempt in range(4):
        p = subprocess.run(["curl", "-sS", "-m", "60", "-A", UA, "-w", "\n%{http_code}", url], capture_output=True)
        body, _, code = p.stdout.rpartition(b"\n")
        if p.returncode == 0 and code == b"200":
            open(path, "wb").write(body)
            time.sleep(0.25)
            return body
        if p.returncode == 0 and code == b"404":
            open(path, "wb").write(b"")
            return b""
        time.sleep(2 + 3 * attempt)
    return b""


def strings(node):
    """All StringWithMarkup strings under a PUG-View node."""
    out = []
    if isinstance(node, dict):
        v = node.get("Value", {})
        for s in v.get("StringWithMarkup", []):
            out.append((s.get("String", ""), [m.get("URL", "") for m in s.get("Markup", [])]))
        for k in ("Information", "Section"):
            for x in node.get(k, []):
                out += strings(x)
    return out


def ghs(record):
    """Signal word, pictograms, hazard statements of the first source of the GHS Classification section."""
    found = None

    def walk(sec):
        nonlocal found
        if found is not None:
            return
        if sec.get("TOCHeading") == "GHS Classification":
            found = sec
            return
        for s in sec.get("Section", []):
            walk(s)

    for s in record.get("Record", {}).get("Section", []):
        walk(s)
    if not found:
        return None
    infos = found.get("Information", [])
    first_ref = infos[0].get("ReferenceNumber") if infos else None
    signal, pictos, hazards = "", [], []
    for info in infos:
        if info.get("ReferenceNumber") != first_ref:
            continue
        name = info.get("Name", "")
        for text, urls in strings(info):
            if name == "Signal":
                signal = signal or text.strip()
            elif name == "Pictogram(s)":
                for u in urls:
                    m = re.search(r"(GHS0\d)", u)
                    if m and m.group(1) not in pictos:
                        pictos.append(m.group(1))
            elif name == "GHS Hazard Statements":
                m = re.match(r"(H\d{3}[A-Za-z+H\d]*)\s*(?:\(([\d.]+)%\))?:?\s*(.*)", text.strip())
                if m:
                    pct = float(m.group(2)) if m.group(2) else 100.0
                    phrase = re.sub(r"\s*\[.*\]$", "", m.group(3)).strip()
                    hazards.append((m.group(1), pct, phrase))
    source = ""
    for ref in record.get("Record", {}).get("Reference", []):
        if ref.get("ReferenceNumber") == first_ref:
            source = ref.get("SourceName", "")
    return {"signal": signal, "pictograms": pictos, "hazards": hazards, "source": source}


def canonical(smiles):
    p = subprocess.run([os.path.join(ROOT, "build", "chemlarp"), "canon", smiles], capture_output=True, text=True)
    lines = [l for l in p.stdout.split("\n") if l.strip()]
    return lines[-1].strip() if p.returncode == 0 and lines else ""


def cstr(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=os.path.join(ROOT, "build", "ghs_cache"))
    ap.add_argument("--out", default=os.path.join(ROOT, "src", "ghs_data.inc"))
    a = ap.parse_args()
    os.makedirs(a.cache, exist_ok=True)
    names = [n.strip() for n in NAMES.split(";") if n.strip()]
    rows, phrases, missing = [], {}, []
    for name in names:
        q = urllib.parse.quote(name)
        cids = json.loads(get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/cids/JSON" % q, a.cache, "cid_" + name) or b"{}")
        cid = (cids.get("IdentifierList", {}).get("CID") or [None])[0]
        if not cid:
            missing.append(name + " (no CID)")
            continue
        props = json.loads(get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%d/property/SMILES,IUPACName/JSON" % cid,
                               a.cache, "prop_%d" % cid) or b"{}")
        p = (props.get("PropertyTable", {}).get("Properties") or [{}])[0]
        smiles = p.get("SMILES") or p.get("IsomericSMILES") or p.get("CanonicalSMILES") or ""
        rec = json.loads(get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%d/JSON?heading=GHS%%20Classification" % cid,
                             a.cache, "ghs_%d" % cid) or b"{}")
        g = ghs(rec) if rec else None
        canon = canonical(smiles) if smiles else ""
        if not canon:
            missing.append("%s (CID %d: structure not read)" % (name, cid))
            continue
        if g is None:
            g = {"signal": "", "pictograms": [], "hazards": [], "source": ""}  # PubChem lists no GHS classification
        for code, _, phrase in g["hazards"]:
            if phrase and code not in phrases:
                phrases[code] = phrase
        rows.append((canon, name, cid, g))
        print("%-32s CID %-8d %-8s %s" % (name, cid, g["signal"] or "-", " ".join(h[0] for h in g["hazards"]) or "(no GHS data)"))
    seen = set()
    uniq = []
    for r in rows:
        if r[0] in seen:
            continue
        seen.add(r[0])
        uniq.append(r)
    with open(a.out, "w") as f:
        f.write("// Generated by scripts/gen_ghs_data.py from PubChem (GHS Classification section, first source,\n")
        f.write("// usually the ECHA C&L Inventory notification summary). %d compounds. Do not edit by hand.\n" % len(uniq))
        f.write("// {canonical SMILES, name, PubChem CID, signal word, pictograms, \"H225:100;H319:95.2;...\"}\n")
        for canon, name, cid, g in uniq:
            hz = ";".join("%s:%g" % (c, p) for c, p, _ in g["hazards"])
            f.write("{%s, %s, %d, %s, %s, %s},\n" % (cstr(canon), cstr(name), cid, cstr(g["signal"]),
                                                   cstr(" ".join(g["pictograms"])), cstr(hz)))
    with open(a.out.replace("ghs_data.inc", "ghs_phrases.inc"), "w") as f:
        f.write("// Generated by scripts/gen_ghs_data.py: UN GHS hazard statement texts as PubChem prints them.\n")
        for code in sorted(phrases):
            f.write("{%s, %s},\n" % (cstr(code), cstr(phrases[code])))
    print("\n%d compounds written, %d hazard statements; not found: %s" % (len(uniq), len(phrases), ", ".join(missing) or "none"))


if __name__ == "__main__":
    sys.exit(main())
