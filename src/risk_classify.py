"""Phase 4 (part 2) — Transfer-risk classification.

Implements the four-tier framework from the research plan:

    High        on a conjugative plasmid  OR  IS elements flanking both sides
    Medium      on a non-conjugative plasmid  OR  a single flanking IS element
    Low         chromosomal, acquired, no mobility signal
    Negligible  intrinsic resistance, no mobility signal

"Intrinsic" resistance (e.g. the natural vancomycin resistance of most
Lactobacillus / Leuconostoc / Pediococcus, or vanC in Enterococcus gallinarum)
is configurable so the framework can be tuned per study.
"""
from __future__ import annotations

from typing import List, Optional

import pandas as pd

from .utils import get_logger

log = get_logger("risk")

RISK_ORDER = ["Negligible", "Low", "Medium", "High"]

# (species_substring or None, gene_prefix) pairs treated as intrinsic resistance.
DEFAULT_INTRINSIC: List[dict] = [
    {"species": "Lactobacillus", "gene_prefix": "van"},
    {"species": "Lactiplantibacillus", "gene_prefix": "van"},
    {"species": "Leuconostoc", "gene_prefix": "van"},
    {"species": "Pediococcus", "gene_prefix": "van"},
    {"species": None, "gene_prefix": "vanC"},   # E. gallinarum / casseliflavus
]


def _is_intrinsic(species: str, gene: str, rules: List[dict]) -> bool:
    g = gene.lower()
    sp = (species or "").lower()
    for r in rules:
        rsp = (r.get("species") or "").lower()
        gp = (r.get("gene_prefix") or "").lower()
        if rsp and rsp not in sp:
            continue
        if gp and not g.startswith(gp):
            continue
        if rsp or gp:
            return True
    return False


def classify(arg_df: pd.DataFrame,
             metadata: Optional[pd.DataFrame] = None,
             intrinsic_rules: Optional[List[dict]] = None) -> pd.DataFrame:
    """Add `risk_level` and `risk_reason` columns to the context-annotated ARG table."""
    rules = intrinsic_rules if intrinsic_rules is not None else DEFAULT_INTRINSIC
    out = arg_df.copy()

    species_lookup = {}
    if metadata is not None and "strain_id" in metadata and "species" in metadata:
        species_lookup = dict(zip(metadata.strain_id, metadata.species))

    levels, reasons = [], []
    for _, r in out.iterrows():
        species = species_lookup.get(r.strain_id, "")
        on_plasmid = bool(r.get("on_plasmid", False))
        conjugative = bool(r.get("conjugative", False))
        is_both = bool(r.get("is_flank_both", False))
        is_single = bool(r.get("is_upstream", False) or r.get("is_downstream", False))

        if (on_plasmid and conjugative) or is_both:
            level = "High"
            reason = ("conjugative plasmid" if (on_plasmid and conjugative)
                      else "IS elements flanking both sides")
        elif on_plasmid or is_single:
            level = "Medium"
            reason = ("non-conjugative plasmid" if on_plasmid
                      else "single flanking IS element")
        elif _is_intrinsic(species, r.gene, rules):
            level = "Negligible"
            reason = "intrinsic resistance, no mobility signal"
        else:
            level = "Low"
            reason = "chromosomal, acquired, no mobility signal"

        levels.append(level)
        reasons.append(reason)

    out["risk_level"] = levels
    out["risk_reason"] = reasons
    out["risk_level"] = pd.Categorical(out["risk_level"], categories=RISK_ORDER, ordered=True)
    return out


def risk_summary_by_species(classified: pd.DataFrame,
                            metadata: pd.DataFrame) -> pd.DataFrame:
    """Counts of ARGs per species per risk level (feeds Figure 2)."""
    df = classified.merge(metadata[["strain_id", "species", "group"]],
                          on="strain_id", how="left")
    tab = (df.groupby(["species", "risk_level"], observed=False)
             .size().unstack(fill_value=0)
             .reindex(columns=RISK_ORDER, fill_value=0))
    return tab
