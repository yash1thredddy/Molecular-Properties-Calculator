"""Built-in formula presets + activity-column placeholder resolution."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

import pandas as pd

# Activity placeholder aliases (extends LEI DependencyChecker conventions).
ACTIVITY_ALIASES = {
    "pIC50": ["pIC50", "pic50", "pKi", "pki", "PKI", "p_ki"],
    "percent_inhibition": ["Percent_Inhibition", "percent_inhibition",
                           "%inhibition", "pct_inhibition", "inhibition",
                           "%_inhibition", "Inhibition_%"],
}

_PLACEHOLDER_RE = re.compile(r"\[activity:([a-zA-Z0-9_]+)\]")


@dataclass(frozen=True)
class Preset:
    name: str
    expression: str          # may contain [activity:...] placeholders
    description: str
    group: str               # "cited" | "extension" | "structure"
    citation: str
    needs_activity: tuple = ()   # placeholder keys required


def detect_activity_column(df: pd.DataFrame, key: str) -> Optional[str]:
    aliases = [a.lower() for a in ACTIVITY_ALIASES.get(key, [])]
    for col in df.columns:
        if col.lower() in aliases:
            return col
    return None


def resolve_placeholders(expression: str, mapping: dict[str, str],
                         raw_concentration: dict[str, bool] | None = None) -> str:
    """Replace [activity:key] with the mapped column. If raw_concentration[key]
    is True, the mapped column is a raw IC50/Ki in nM -> inject -log10(col/1e9)
    to convert to a p-scale value (resolves the unit-collision risk)."""
    raw_concentration = raw_concentration or {}
    def repl(m: re.Match) -> str:
        key = m.group(1)
        if key not in mapping:
            raise KeyError(f"No column mapped for activity:{key}")
        col = mapping[key]
        if raw_concentration.get(key):
            # nM -> M is /1e9; pX = -log10(M)
            return f"(0 - log10([{col}] / 1000000000))"
        return f"[{col}]"
    return _PLACEHOLDER_RE.sub(repl, expression)


PRESETS: list[Preset] = [
    Preset("PEI", "([activity:percent_inhibition]/100) / ([Molecular_Weight]/1000)",
           "Percentage Efficiency Index — % inhibition (fraction) per kDa.",
           "cited", "Abad-Zapatero & Metz, Drug Discov. Today 2005",
           ("percent_inhibition",)),
    Preset("PSEI", "([activity:percent_inhibition]/100) / ([TPSA]/100)",
           "PEI logic on the polar-surface (SEI) axis. Derived extension.",
           "extension", "Not a separately published metric",
           ("percent_inhibition",)),
    Preset("LEH", "1.37 * [activity:pIC50] / [Heavy_Atom_Count]",
           "Hopkins ligand efficiency (~kcal/mol per heavy atom). 1.37=2.303RT@298K.",
           "cited", "Hopkins et al. 2004", ("pIC50",)),
    Preset("LLE", "[activity:pIC50] - [LogP]",
           "Lipophilic ligand efficiency (LipE).",
           "cited", "Leeson & Springthorpe 2007", ("pIC50",)),
    Preset("LELP", "[LogP] / (1.37 * [activity:pIC50] / [Heavy_Atom_Count])",
           "LogP price paid per unit of ligand efficiency.",
           "cited", "Keseru & Makara 2009", ("pIC50",)),
    Preset("BEI", "[activity:pIC50] / ([Molecular_Weight]/1000)",
           "Binding Efficiency Index — potency per kDa.",
           "cited", "Abad-Zapatero & Metz 2005", ("pIC50",)),
    Preset("SEI", "[activity:pIC50] / ([TPSA]/100)",
           "Surface Efficiency Index — potency per 100 A^2 polar surface.",
           "cited", "Abad-Zapatero & Metz 2005", ("pIC50",)),
    Preset("SILE", "[activity:pIC50] / ([Heavy_Atom_Count] ^ 0.3)",
           "Size-independent ligand efficiency.",
           "cited", "Nissink, JCIM 2009", ("pIC50",)),
    Preset("Lipinski_Violation_Count",
           "if([Molecular_Weight]>500,1,0)+if([LogP]>5,1,0)"
           "+if([HB_Donors]>5,1,0)+if([HB_Acceptors]>10,1,0)",
           "Count of Rule-of-Five violations (0-4).",
           "structure", "Lipinski et al. 2001"),
    Preset("Veber_Pass", "if([Rotatable_Bonds]<=10 and [TPSA]<=140, 1, 0)",
           "1 if it passes Veber oral-bioavailability criteria.",
           "structure", "Veber et al. 2002"),
    Preset("TPSA_per_HeavyAtom", "[TPSA] / [Heavy_Atom_Count]",
           "Polar surface density per heavy atom.", "structure", "derived QC"),
    Preset("MW_per_HeavyAtom", "[Molecular_Weight] / [Heavy_Atom_Count]",
           "Mean atomic mass; flags heavy elements.", "structure", "derived QC"),
]
