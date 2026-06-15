"""Phase 2 - Curate ARG calls: separate acquired mobile cargo from core/
intrinsic / biocide / efflux determinants. Raw calls are NEVER deleted; every
original row is preserved with added curation columns.

Usage:
    py scripts\\curate_arg_calls.py
    py scripts\\curate_arg_calls.py --results results_scale50 --rules config\\intrinsic_core_gene_rules.yaml
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from _sota_common import df_to_md, read_csv_safe, sota_dir

INCLUDE_CATEGORIES = {"acquired_mobile_ARG", "ambiguous_requires_validation"}


def load_rules(path: Path) -> list[dict]:
    data = yaml.safe_load(Path(path).read_text())
    return data.get("rules", [])


def classify(gene: str, rules: list[dict]) -> dict:
    g = str(gene).lower()
    for r in rules:
        if r["pattern"].lower() in g:
            return {"curated_category": r["category"],
                    "curation_reason": r["reason"],
                    "confidence": r.get("confidence", "medium")}
    # default: treat as acquired mobile cargo, confidence depends on known families
    known = ("bla", "ctx", "kpc", "ndm", "oxa", "tem", "shv", "van", "erm",
             "tet", "aac", "aph", "ant", "aad", "sul", "dfr", "cat", "cml",
             "mph", "mef", "lnu", "lsa", "qnr", "mcr", "fos", "cfr", "rmt", "arm")
    conf = "medium" if g.startswith(known) else "low"
    return {"curated_category": "acquired_mobile_ARG",
            "curation_reason": "No core/intrinsic/biocide rule matched; treated as acquired cargo.",
            "confidence": conf}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_scale50")
    ap.add_argument("--rules", default="config/intrinsic_core_gene_rules.yaml")
    args = ap.parse_args(argv)
    out = sota_dir()

    arg = read_csv_safe(Path(args.results) / "table2_master_arg.csv")
    if arg is None:
        raise SystemExit(f"table2_master_arg.csv not found in {args.results}")
    rules = load_rules(Path(args.rules))

    arg = arg.copy()
    arg["raw_risk_level"] = arg["risk_level"]
    cls = arg["gene"].apply(lambda g: classify(g, rules))
    arg["curated_category"] = [c["curated_category"] for c in cls]
    arg["curation_reason"] = [c["curation_reason"] for c in cls]
    arg["confidence"] = [c["confidence"] for c in cls]
    arg["include_in_mobile_ARG_analysis"] = arg["curated_category"].isin(INCLUDE_CATEGORIES)
    # curated risk: keep tier if included, else mark Excluded (raw is preserved)
    arg["curated_risk_level"] = arg.apply(
        lambda r: r["raw_risk_level"] if r["include_in_mobile_ARG_analysis"]
        else "Excluded-curated", axis=1)

    arg.to_csv(out / "curated_ARG_calls_SCALE50.csv", index=False)

    # curated High/Medium = included rows still at High/Medium after curation
    inc = arg[arg.include_in_mobile_ARG_analysis]
    hm = inc[inc.curated_risk_level.isin(["High", "Medium"])]
    hm.to_csv(out / "curated_high_medium_SCALE50.csv", index=False)

    # ---- summary -------------------------------------------------------- #
    cat_counts = arg.curated_category.value_counts()
    excluded = arg[~arg.include_in_mobile_ARG_analysis]
    raw_hm = int(arg.raw_risk_level.isin(["High", "Medium"]).sum())
    cur_hm = len(hm)

    md = ["# ARG curation summary (SCALE50)\n",
          "Raw calls are preserved; curation only adds columns and flags which "
          "rows are eligible for the mobile-ARG analysis.\n",
          "## Category breakdown\n",
          df_to_md(cat_counts.rename_axis("curated_category").reset_index(name="count")),
          f"\nHigh/Medium loci before curation: **{raw_hm}**; after removing "
          f"core/intrinsic/biocide/efflux determinants: **{cur_hm}**.\n",
          f"\n## Excluded from mobile-ARG analysis ({len(excluded)} rows)\n"]
    if len(excluded):
        ecols = [c for c in ["strain_id", "gene", "raw_risk_level",
                             "curated_category", "confidence", "curation_reason"]
                 if c in excluded.columns]
        md.append(df_to_md(excluded[ecols].drop_duplicates(subset=["strain_id", "gene"]),
                           max_rows=80))
    else:
        md.append("_none_\n")
    md.append("\n## Curated High/Medium acquired cargo (kept)\n")
    kcols = [c for c in ["strain_id", "gene", "drug_class", "curated_risk_level",
                         "confidence"] if c in hm.columns]
    md.append(df_to_md(hm[kcols], max_rows=80))
    (out / "curation_summary.md").write_text("\n".join(md), encoding="utf-8")

    print(f"curated {len(arg)} rows -> {out/'curated_ARG_calls_SCALE50.csv'}")
    print(f"High/Medium: {raw_hm} raw -> {cur_hm} after curation "
          f"({len(excluded)} rows excluded as core/intrinsic/biocide/efflux)")


if __name__ == "__main__":
    main()
