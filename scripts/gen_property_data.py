"""Regenerates src/property_data.inc, the pure-component property table, from published compilations as
distributed with the `chemicals` Python package (MIT licence). Compounds are those of the ideal-gas heat
capacity databank of Poling, Prausnitz & O'Connell, The Properties of Gases and Liquids, 5th ed. (2001),
keyed by SMILES (resolved from CAS numbers with `chemicals`).

Field sources (first available wins; the source letter is stored with the critical constants):
  Tc, Pc, Vc   I: IUPAC critical-property series (Ambrose, Tsonopoulos, Young et al., J. Chem. Eng. Data
                  1995-2015), C: CRC Handbook of Chemistry and Physics, P: PSRK appendix (Horstmann et al.,
                  Fluid Phase Equilib. 227 (2005) 157)
  omega        by its definition, -log10(Psat(0.7 Tc)/Pc) - 1, with the Tc and Pc above and Psat from the Wagner
               equations of Poling et al. or McGarry (Ind. Eng. Chem. Process Des. Dev. 22 (1983) 313), or the
               Poling Antoine equation when 0.7 Tc is inside its range; otherwise a tabulated value (PSRK
               appendix, then Passut & Danner, Ind. Eng. Chem. Process Des. Dev. 12 (1973) 365) only when that
               table's own Tc and Pc agree with ours (1%, 5%), since omega depends on them
  Tm, Tb       CRC Handbook, Physical Constants of Organic (then Inorganic) Compounds; Tb then from the CRC
               enthalpy-of-vaporization table
  dHvap(Tb)    CRC Handbook, Enthalpy of Vaporization
  Cp ideal gas Poling et al. Appendix A: Cp/R = a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4 (Tmin..Tmax)
  Cp liquid    Poling et al. Appendix A, at 298.15 K
  Z_RA         Rackett parameter, COSTALD/Rackett tables (Hankinson & Thomson, AIChE J. 25 (1979) 653;
                  Spencer & Danner, J. Chem. Eng. Data 17 (1972) 236)
  dfH, dfG, S  ideal gas at 298.15 K and 1 bar, CRC Handbook, Standard Thermodynamic Properties

Two consistency checks guard against transcription errors in the compilations (reported on stderr):
critical constants come from IUPAC, then CRC, then PSRK, except that an IUPAC row disagreeing with CRC
(Tc 0.5%, Pc 5%) is replaced when another source (PSRK, Yaws, or the constants published with Wagner
equations) confirms CRC; unconfirmed rows are kept and listed. A Cp polynomial that misses the databank's own Cp(298.15 K) by more than 2% is kept only if it agrees
within 2% with a third source (NIST-JANAF tables, else the TRC fit) and dropped otherwise.
Missing values are written as NAN. Usage: python scripts/gen_property_data.py > src/property_data.inc
"""
import json
import math
import os
import sys

import chemicals
import chemicals.identifiers as I
import pandas as pd
from chemicals.heat_capacity import TRCCp
import chemicals.vapor_pressure as VP

D = os.path.dirname(chemicals.__file__)
# Antoine rows that miss the normal boiling point (see scripts/gen_antoine_data.py).
BAD_ANTOINE = {"2207-04-7", "7783-54-2", "10544-72-6"}


def table(rel, **kw):
    df = pd.read_csv(os.path.join(D, rel), sep=kw.pop("sep", "\t"), index_col=0, **kw)
    return df[~df.index.duplicated(keep="first")]


def num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return math.nan
    return v if math.isfinite(v) else math.nan


def first(*vals):
    for v in vals:
        if not math.isnan(v):
            return v
    return math.nan


def c(v):
    return "NAN" if math.isnan(v) else repr(v)


def main():
    poling = table("Heat Capacity/PolingDatabank.tsv")
    iupac = table("Critical Properties/IUPACOrganicCriticalProps.tsv")
    crc_crit = table("Critical Properties/CRCCriticalOrganics.tsv")
    psrk = table("Critical Properties/Appendix to PSRK Revision 4.tsv")
    pd73 = table("Critical Properties/PassutDanner1973.tsv")
    yaws = table("Critical Properties/Yaws Collection.tsv")
    phys = table("Misc/Physical Constants of Organic Compounds.csv")
    phys_inorg = table("Misc/Physical Constants of Inorganic Compounds.csv")
    hvap = table("Phase Change/CRC Handbook Heat of Vaporization.tsv")
    costald = table("Density/COSTALD Parameters.tsv")
    crc_thermo = table("Heat Capacity/CRC Standard Thermodynamic Properties of Chemical Substances.tsv")

    janaf = json.load(open(os.path.join(D, "Heat Capacity/JANAF_1998_gas_Cp.json")))
    trc = table("Heat Capacity/TRC Thermodynamics of Organic Compounds in the Gas State.tsv")

    def third_cp(cas):
        """Ideal-gas Cp at 298.15 K from the NIST-JANAF tables (1998), else the TRC fit."""
        if cas in janaf and 298.15 in janaf[cas][0]:
            return janaf[cas][1][janaf[cas][0].index(298.15)], "JANAF"
        if cas in trc.index and trc.at[cas, "Tmin"] <= 298.15 <= trc.at[cas, "Tmax"]:
            return TRCCp(298.15, *[float(trc.at[cas, f"a{i}"]) for i in range(8)]), "TRC"
        return math.nan, ""

    def psat_07(cas, T):
        """Vapour pressure from a published equation valid at T, or NaN."""
        w = VP.Psat_data_WagnerPoling
        if cas in w.index and w.at[cas, "Tmin"] <= T <= w.at[cas, "Tmax"]:
            r = w.loc[cas]
            return VP.Wagner(T, r.Tc, r.Pc, r.A, r.B, r.C, r.D)
        w = VP.Psat_data_WagnerMcGarry
        if cas in w.index and w.at[cas, "Tmin"] <= T < w.at[cas, "Tc"]:
            r = w.loc[cas]
            return VP.Wagner_original(T, r.Tc, r.Pc, r.A, r.B, r.C, r.D)
        a = VP.Psat_data_AntoinePoling
        if cas in a.index and a.at[cas, "Tmin"] <= T <= a.at[cas, "Tmax"] and cas not in BAD_ANTOINE:
            r = a.loc[cas]
            return 10 ** (r.A - r.B / (T + r.C))
        return math.nan

    def omega_definition(cas, Tc, Pc):
        if math.isnan(Tc) or math.isnan(Pc):
            return math.nan
        p = psat_07(cas, 0.7 * Tc)
        return -math.log10(p / Pc) - 1 if p > 0 else math.nan

    def get(df, cas, col):
        return num(df.at[cas, col]) if cas in df.index and col in df.columns else math.nan

    print("// Generated by scripts/gen_property_data.py from published compilations (chemicals "
          f"{chemicals.__version__}); see the script for the source of each field. Do not edit by hand.")
    print("// {SMILES, name, CAS, MW, Tc K, Pc Pa, Vc m3/mol, crit source, omega, Tm K, Tb K, dHvap(Tb) J/mol,")
    print("//  Cp/R a0..a4, Cp Tmin K, Cp Tmax K, Cp liquid 298 J/mol/K, Z_RA, dfH J/mol, dfG J/mol, S J/mol/K}")
    print("const PropertyRow kProperties[] = {")
    seen = set()
    for cas, row in poling.iterrows():
        try:
            meta = I.search_chemical(cas)
        except Exception:
            continue
        smiles = getattr(meta, "smiles", None)
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)
        cand = {}
        for letter, df in (("I", iupac), ("C", crc_crit), ("P", psrk)):
            Tc, Pc = get(df, cas, "Tc"), get(df, cas, "Pc")
            if not (math.isnan(Tc) or math.isnan(Pc)):
                cand[letter] = (Tc, Pc, get(df, cas, "Vc"))
        # Other published constants that vote but are not used directly: Yaws, and the Tc, Pc published
        # with the Wagner vapour-pressure equations of Poling et al. and McGarry.
        voters = dict(cand)
        for letter, df in (("Y", yaws), ("W", VP.Psat_data_WagnerPoling), ("M", VP.Psat_data_WagnerMcGarry)):
            Tc, Pc = get(df, cas, "Tc"), get(df, cas, "Pc")
            if not (math.isnan(Tc) or math.isnan(Pc)):
                voters[letter] = (Tc, Pc, math.nan)

        def agree(a, b):
            return abs(a[0] / b[0] - 1) <= 0.005 and abs(a[1] / b[1] - 1) <= 0.05

        # Critically evaluated compilations first: IUPAC, then CRC, then PSRK. An IUPAC row is replaced by
        # the CRC row only when the two disagree and another source confirms CRC: this catches
        # transcription errors (cyclohexane and methylcyclopentane exchanged, decane Tc = 611.7 K) without
        # letting older compilations, which often share one origin, outvote evaluated data.
        order = [k for k in "ICP" if k in cand]
        if "I" in cand and "C" in cand and not agree(cand["I"], cand["C"]) and \
                any(agree(cand["C"], v) for j, v in voters.items() if j not in ("I", "C")):
            print(f"# {cas} {row['Chemical'].strip()}: IUPAC critical constants {cand['I'][:2]} disagree with "
                  f"CRC {cand['C'][:2]}, which another source confirms; using CRC", file=sys.stderr)
            order.remove("I")
        elif order and not any(agree(cand[order[0]], v) for j, v in voters.items() if j != order[0]) and \
                len(voters) > 1:
            print(f"# {cas} {row['Chemical'].strip()}: critical constants of {order[0]} {cand[order[0]][:2]} "
                  f"not confirmed by {sorted(set(voters) - {order[0]})}; kept", file=sys.stderr)
        letter = order[0] if order else "-"
        Tc, Pc, Vc = cand[letter] if order else (math.nan, math.nan, math.nan)
        for k in order[1:]:
            if not math.isnan(Vc):
                break
            Vc = cand[k][2]
        omega = omega_definition(cas, Tc, Pc)
        for df in (psrk, pd73):
            if not math.isnan(omega) or math.isnan(Tc):
                break
            if abs(get(df, cas, "Tc") / Tc - 1) <= 0.01 and abs(get(df, cas, "Pc") / Pc - 1) <= 0.05:
                omega = get(df, cas, "omega")
        name = str(row["Chemical"]).strip().replace('"', "'")
        # The Cp polynomial must reproduce the databank's own Cp at 298.15 K (2%), which catches
        # mis-scaled coefficients (3-methylpentane).
        a = [num(row[f"a{i}"]) for i in range(5)]
        cp298 = 8.314462618 * sum(ai * 298.15 ** i for i, ai in enumerate(a))
        if not math.isnan(num(row["Cpg"])) and abs(cp298 / num(row["Cpg"]) - 1) > 0.02:
            # The tabulated value is sometimes the wrong one (acetylene, chlorine): a third source decides.
            ref, src = third_cp(cas)
            if not math.isnan(ref) and abs(cp298 / ref - 1) <= 0.02:
                print(f"# {cas} {name}: Cp polynomial {cp298:.1f} vs table {row['Cpg']} at 298 K; kept, {src} "
                      f"gives {ref:.1f}", file=sys.stderr)
            else:
                print(f"# {cas} {name}: Cp polynomial {cp298:.1f} vs table {row['Cpg']} at 298 K; dropped "
                      f"({src or 'no third source'} {ref:.1f})", file=sys.stderr)
                a = [math.nan] * 5
        f = [
            meta.MW, Tc, Pc, Vc, None, omega, first(get(phys, cas, "Tm"), get(phys_inorg, cas, "Tm")),
            first(get(phys, cas, "Tb"), get(phys_inorg, cas, "Tb"), get(hvap, cas, "Tb")), get(hvap, cas, "HvapTb"),
            *a, num(row["Tmin"]), num(row["Tmax"]), num(row["Cpl"]), get(costald, cas, "Z_RA"), get(crc_thermo, cas, "Hfg"),
            get(crc_thermo, cas, "Gfg"), get(crc_thermo, cas, "S0g"),
        ]
        vals = ", ".join(f"'{letter}'" if v is None else c(num(v)) for v in f)
        smiles = smiles.replace("\\", "\\\\")
        print(f'    {{"{smiles}", "{name}", "{cas}", {vals}}},')
    print("};")


if __name__ == "__main__":
    main()
