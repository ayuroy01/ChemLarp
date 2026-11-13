#!/usr/bin/env python3
"""Generates src/aqueous_data.inc from PHREEQC's phreeqc.dat (USGS, public domain; Parkhurst & Appelo,
U.S. Geological Survey Techniques and Methods 6-A43, 2013).

Usage: python3 scripts/gen_aqueous_data.py path/to/phreeqc.dat > src/aqueous_data.inc

What is taken over:
  - SOLUTION_MASTER_SPECIES: one component per element, or per valence state for redox elements
    ("Fe(2)", "Fe(3)", "C(4)", "C(-4)", "N(5)", "N(-3)", ...). H, O, e- and Alkalinity are not components.
  - SOLUTION_SPECIES: every species whose formation reaction, rewritten in terms of component master
    species, H+ and H2O, has no electrons (redox reactions between valence states are left out: each
    valence state is its own component). log K(T) is kept as the sum of the analytic expressions
    (-analytic A1..A6: log K = A1 + A2 T + A3/T + A4 log10 T + A5/T^2 + A6 T^2) and the van 't Hoff parts
    (-log_k at 25 C with -delta_h) of the reactions it was rewritten from, which PHREEQC uses in the same
    way (an analytic expression takes precedence over log_k and delta_h of the same reaction). The
    Debye-Hueckel parameters of -gamma (ion size a, b) are kept; the last -gamma line counts.
  - PHASES: minerals and gases whose dissolution reaction, rewritten the same way, has no electrons.
Molar volumes (-Vm), diffusion coefficients and Peng-Robinson constants (pressure effects) are not used.
"""
import re
import sys
from fractions import Fraction

R_LN10 = 8.314462618e-3 * 2.302585092994046  # kJ/(mol K)


def tokens_of_side(side):
    """'CO3-2 + 2 H+' -> [(Fraction(1), 'CO3-2'), (Fraction(2), 'H+')]."""
    out = []
    # Split on '+' that separates terms: a '+' preceded by whitespace and followed by whitespace.
    for term in re.split(r"\s\+\s", " " + side.strip() + " "):
        term = term.strip()
        if not term:
            continue
        m = re.match(r"^([0-9]*\.?[0-9]+)\s+(\S+)$", term) or re.match(r"^([0-9]*\.?[0-9]+)(\(?[A-Z].*)$", term)
        if m:
            out.append((Fraction(m.group(1)), m.group(2).strip()))
        else:
            out.append((Fraction(1), term))
    return out


def charge_of(species):
    if species == "e-":
        return -1
    m = re.search(r"([+-])(\d*)$", species)
    if not m:
        return 0
    n = int(m.group(2)) if m.group(2) else 1
    return n if m.group(1) == "+" else -n


def parse_logk_block(opts):
    """Returns (analytic[6] or None, log_k25, delta_h kJ/mol)."""
    analytic, logk, dh = None, 0.0, 0.0
    for key, vals in opts:
        k = key.lower()
        if k in ("-log_k", "log_k", "-logk"):
            logk = float(vals[0])
        elif k in ("-delta_h", "delta_h", "-deltah"):
            dh = float(vals[0])
            unit = vals[1].lower() if len(vals) > 1 and not vals[1].startswith("#") else "kj"
            if unit.startswith("kcal"):
                dh *= 4.184
            elif unit.startswith("cal"):
                dh *= 4.184e-3
            elif unit.startswith("j") :
                dh *= 1e-3
        elif k.startswith("-analytic") or k.startswith("-a_e") or k == "-analytical_expression":
            nums = []
            for v in vals:
                if v.startswith("#"):
                    break
                nums.append(float(v))
            analytic = (nums + [0.0] * 6)[:6]
    return analytic, logk, dh


def read_blocks(text):
    sections = {}
    current = None
    lines = []
    for raw in text.split("\n"):
        raw = raw.split("#", 1)[0]
        # ';' separates statements on one line; an indented line keeps its indentation for the pieces.
        lead = raw[: len(raw) - len(raw.lstrip())]
        parts = raw.split(";")
        lines.append(parts[0].rstrip())
        lines.extend(lead + "\t" + q.strip() for q in parts[1:] if q.strip())
    for line in lines:
        if not line.strip():
            continue
        if re.match(r"^[A-Z_]+\s*$", line):
            current = line.strip()
            sections.setdefault(current, [])
            continue
        if current:
            sections[current].append(line)
    return sections


def main():
    text = open(sys.argv[1], encoding="latin-1").read().replace("\r", "")
    sec = read_blocks(text)

    # ---- masters ----
    masters = []          # (component name, master species, alk, element, gfw)
    element_default = {}  # "Fe" -> "Fe(2)"
    gfw = {}
    raw_masters = []
    for line in sec["SOLUTION_MASTER_SPECIES"]:
        f = line.split()
        raw_masters.append(f)
    has_valences = {re.sub(r"\(.*", "", f[0]) for f in raw_masters if "(" in f[0]}
    skip = {"H", "O", "E", "Alkalinity"}
    for f in raw_masters:
        name, species, alk = f[0], f[1], float(f[2])
        el = re.sub(r"\(.*", "", name)
        if len(f) >= 5:
            gfw[el] = float(f[4])
        if el in skip:
            continue
        if "(" not in name and el in has_valences:
            continue  # the plain element line: its valence-state lines are the components
        comp = name
        m = re.match(r"^(\w+)\(\+?(-?\d+)\)$", name)
        if m:
            comp = f"{m.group(1)}({m.group(2)})"
        masters.append([comp, species, alk, el])
    for f in raw_masters:
        name, species = f[0], f[1]
        el = re.sub(r"\(.*", "", name)
        if "(" not in name and el in has_valences:
            for mm in masters:
                if mm[1] == species and mm[3] == el:
                    element_default[el] = mm[0]
    # De-duplicate (phreeqc.dat lists Sg twice).
    seen, uniq = set(), []
    for mm in masters:
        if mm[0] not in seen:
            seen.add(mm[0])
            uniq.append(mm)
    masters = uniq
    master_of_species = {mm[1]: mm[0] for mm in masters}

    # ---- species: parse equations and options ----
    def parse_equations(lines, phase_mode=False):
        items = []
        cur = None
        name = None
        for line in lines:
            s = line.strip()
            if phase_mode and not line.startswith((" ", "\t")) and "=" not in s:
                name = s.split()[0]
                continue
            if "=" in s and not s.startswith("-"):
                cur = {"eq": s, "opts": [], "name": name}
                items.append(cur)
            elif s.startswith("-") or s.split()[0].lower() in ("log_k", "delta_h"):
                if cur is not None:
                    p = s.split()
                    cur["opts"].append((p[0], p[1:]))
        return items

    species_raw = parse_equations(sec["SOLUTION_SPECIES"])
    # definitions: species -> (stoich dict over base symbols, analytic, logk, dh) with base symbols being
    # component masters, H+, H2O, e-.
    base = set(master_of_species) | {"H+", "H2O", "e-"}
    defs = {}
    gamma = {}
    order = []
    pending = []
    for it in species_raw:
        lhs, rhs = it["eq"].split("=", 1)
        L, Rr = tokens_of_side(lhs), tokens_of_side(rhs)
        if not Rr:
            continue
        coef0, sp = Rr[0]
        analytic, logk, dh = parse_logk_block(it["opts"])
        for key, vals in it["opts"]:
            if key.lower() == "-gamma":
                gamma[sp] = (float(vals[0]), float(vals[1]))
        if len(L) == 1 and len(Rr) == 1 and L[0][1] == sp:
            continue  # identity (a master species)
        stoich = {}
        for c, s in L:
            stoich[s] = stoich.get(s, 0) + c / coef0
        for c, s in Rr[1:]:
            stoich[s] = stoich.get(s, 0) - c / coef0
        scale = 1 / coef0
        pending.append((sp, stoich, analytic, logk * float(scale), dh * float(scale), float(scale)))
    # Resolve in order (a species may be defined from one defined earlier).
    def combine(terms, product_side=False):
        # A species s formed from masters with K_f (a_s = K_f prod a_m^nu) that appears c times among the
        # reactants of a formation reaction adds c log K_f; among the products of a dissolution reaction
        # (phases), it subtracts c log K_f. Its stoichiometry enters with +c either way.
        st, an, lk, dh = {}, None, 0.0, 0.0
        for i, (mult, (s2, a2, l2, d2)) in enumerate(terms):
            km = -mult if (product_side and i > 0) else mult
            for k, v in s2.items():
                st[k] = st.get(k, 0) + mult * v
            if a2 is not None:
                an = [x + float(km) * y for x, y in zip(an or [0.0] * 6, a2)]
            lk += float(km) * l2
            dh += float(km) * d2
        return {k: v for k, v in st.items() if v != 0}, an, lk, dh

    for sp, stoich, analytic, logk, dh, scale in pending:
        if sp in master_of_species and sp not in ("H+", "H2O"):
            continue  # a valence-state master (e.g. Fe+3 from Fe+2): its own component
        an = [x * scale for x in analytic] if analytic is not None else None
        terms = [(Fraction(1), ({}, an, 0.0 if an else logk, 0.0 if an else dh))]
        ok = True
        for s, c in stoich.items():
            if s in base:
                terms.append((Fraction(1), ({s: c}, None, 0.0, 0.0)))
            elif s in defs:
                terms.append((c, defs[s]))
            else:
                ok = False
        if not ok:
            continue
        st, a2, l2, d2 = combine(terms)
        if sp in defs:
            continue
        defs[sp] = (st, a2, l2, d2)
        order.append(sp)

    def usable(st):
        return "e-" not in st and all(k in base for k in st)

    # ---- phases ----
    phases = []
    for it in parse_equations(sec["PHASES"], phase_mode=True):
        lhs, rhs = it["eq"].split("=", 1)
        L, Rr = tokens_of_side(lhs), tokens_of_side(rhs)
        formula = L[0][1]
        stoich = {}
        for c, s in Rr:
            stoich[s] = stoich.get(s, 0) + c / L[0][0]
        for c, s in L[1:]:
            stoich[s] = stoich.get(s, 0) - c / L[0][0]
        analytic, logk, dh = parse_logk_block(it["opts"])
        terms = [(Fraction(1), ({}, analytic, 0.0 if analytic else logk, 0.0 if analytic else dh))]
        ok = True
        for s, c in stoich.items():
            if s in base:
                terms.append((Fraction(1), ({s: c}, None, 0.0, 0.0)))
            elif s in defs:
                terms.append((c, defs[s]))
            else:
                ok = False
        if not ok:
            continue
        st, a2, l2, d2 = combine(terms, product_side=True)
        if not usable(st):
            continue
        phases.append((it["name"], formula, st, a2, l2, d2))

    comp_index = {mm[0]: i for i, mm in enumerate(masters)}
    sym_index = lambda s: -1 if s == "H+" else -2 if s == "H2O" else comp_index[master_of_species[s]]

    def stoich_c(st):
        items = sorted(st.items(), key=lambda kv: sym_index(kv[0]))
        return "{" + ", ".join(f"{{{sym_index(k)}, {float(v)!r}}}" for k, v in items) + "}"

    def logk_c(an, lk, dh):
        a = an if an is not None else [0.0] * 6
        return "{" + ", ".join(repr(float(x)) for x in a) + "}, " + repr(lk) + ", " + repr(dh)

    out = []
    out.append("// Generated by scripts/gen_aqueous_data.py from PHREEQC's phreeqc.dat (USGS, public domain;")
    out.append("// Parkhurst & Appelo 2013). Do not edit by hand.")
    out.append("// Stoichiometry indices: >= 0 component, -1 H+, -2 H2O.")
    out.append("static const AqMaster kAqMasters[] = {")
    for comp, species, alk, el in masters:
        out.append(f'    {{"{comp}", "{species}", {charge_of(species)}, {alk!r}, "{el}", {gfw.get(el, 0.0)!r}}},')
    out.append("};")
    out.append("static const AqDefault kAqDefaults[] = {")
    for el, comp in sorted(element_default.items()):
        out.append(f'    {{"{el}", "{comp}"}},')
    out.append("};")
    out.append("static const AqSpeciesData kAqSpecies[] = {")
    n_sp = 0
    for sp in order:
        st, an, lk, dh = defs[sp]
        if not usable(st):
            continue
        a, b = gamma.get(sp, (0.0, 0.0))
        hg = "true" if sp in gamma else "false"
        out.append(f'    {{"{sp}", {charge_of(sp)}, {stoich_c(st)}, {logk_c(an, lk, dh)}, {hg}, {a!r}, {b!r}}},')
        n_sp += 1
    out.append("};")
    out.append("// Master species' own Debye-Hueckel parameters (index = component).")
    out.append("static const AqGamma kAqMasterGamma[] = {")
    for comp, species, alk, el in masters:
        a, b = gamma.get(species, (0.0, 0.0))
        out.append(f'    {{{"true" if species in gamma else "false"}, {a!r}, {b!r}}},  // {species}')
    out.append("};")
    a, b = gamma.get("H+", (9.0, 0.0))
    out.append(f"static const AqGamma kAqHydrogenGamma = {{true, {a!r}, {b!r}}};")
    out.append("static const AqPhaseData kAqPhases[] = {")
    for name, formula, st, an, lk, dh in phases:
        out.append(f'    {{"{name}", "{formula}", {stoich_c(st)}, {logk_c(an, lk, dh)}}},')
    out.append("};")
    nmax = max([len(d[0]) for d in defs.values()] + [len(p[2]) for p in phases])
    if nmax > 6:
        sys.exit(f"a reaction has {nmax} terms; raise AqTerms in src/aqueous.cpp")
    sys.stderr.write(f"{len(masters)} components, {n_sp} species, {len(phases)} phases\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
