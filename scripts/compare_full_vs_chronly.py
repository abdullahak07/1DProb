"""Phase 1 - Compare full-assembly SCALE50 against the MATCHED chromosome-only
SCALE50 baseline (NOT the old SCALE8 baseline).

Usage:
    py scripts\\compare_full_vs_chronly.py
    py scripts\\compare_full_vs_chronly.py --full results_scale50 --chronly results_scale50_chronly
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _sota_common import df_to_md, read_csv_safe, sota_dir

RISK = ["Low", "Medium", "High", "Negligible"]


def _summarise(results: Path) -> dict | None:
    arg = read_csv_safe(results / "table2_master_arg.csv")
    feat = read_csv_safe(results / "mge_features.csv")
    if arg is None:
        return None
    rc = arg["risk_level"].value_counts().to_dict() if "risk_level" in arg else {}
    d = {"arg_loci": len(arg)}
    for r in RISK:
        d[r] = int(rc.get(r, 0))
    d["insertion_sequence"] = int((feat.feature_type == "insertion_sequence").sum()) if feat is not None else 0
    d["plasmid_replicon"] = int((feat.feature_type == "plasmid_replicon").sum()) if feat is not None else 0
    d["conjugation_marker"] = int((feat.feature_type == "conjugation_marker").sum()) if feat is not None else 0
    d["_arg_df"] = arg
    return d


def _locus_key(df: pd.DataFrame) -> set:
    return set(zip(df.strain_id, df.gene))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", default="results_scale50")
    ap.add_argument("--chronly", default="results_scale50_chronly")
    args = ap.parse_args(argv)
    out = sota_dir()

    full = _summarise(Path(args.full))
    chr_ = _summarise(Path(args.chronly))

    if full is None:
        raise SystemExit(f"full-assembly results not found in {args.full}")

    if chr_ is None:
        msg = ("# Full-assembly vs chromosome-only comparison\n\n"
               "**Chromosome-only matched baseline NOT RUN YET.**\n\n"
               f"Run it first:\n\n"
               "    py -m src.pipeline --config config\\config_scale50_chronly.yaml\n\n"
               "Then re-run this script. (This deliberately does NOT fall back to "
               "the old SCALE8 baseline.)\n")
        (out / "full_vs_chronly_summary.md").write_text(msg, encoding="utf-8")
        # still write a one-column CSV with the full-assembly side
        rows = [{"metric": k, "full_assembly": v} for k, v in full.items() if k != "_arg_df"]
        pd.DataFrame(rows).to_csv(out / "full_vs_chronly_comparison.csv", index=False)
        print("chromosome-only baseline missing - wrote 'not run yet' summary.")
        return

    # ---- metric comparison table --------------------------------------- #
    metrics = ["arg_loci", "Low", "Medium", "High", "Negligible",
               "insertion_sequence", "plasmid_replicon", "conjugation_marker"]
    rows = []
    for m in metrics:
        f, c = full[m], chr_[m]
        rows.append({"metric": m, "chromosome_only": c, "full_assembly": f,
                     "delta": f - c})
    comp = pd.DataFrame(rows)
    comp.to_csv(out / "full_vs_chronly_comparison.csv", index=False)

    # ---- genes / loci gained & lost ------------------------------------ #
    f_arg, c_arg = full["_arg_df"], chr_["_arg_df"]
    f_keys, c_keys = _locus_key(f_arg), _locus_key(c_arg)
    gained = sorted(f_keys - c_keys)
    lost = sorted(c_keys - f_keys)

    # plasmid-associated ARGs gained in full assembly
    plas_gained = pd.DataFrame()
    if "on_plasmid" in f_arg.columns:
        fg = f_arg[f_arg.apply(lambda r: (r.strain_id, r.gene) in (f_keys - c_keys), axis=1)]
        plas_gained = fg[fg.on_plasmid == True]  # noqa: E712

    # ---- explicit reconciliation (A/B/C/D) ----------------------------- #
    net_row_delta = len(f_arg) - len(c_arg)
    accounting = pd.DataFrame([
        {"quantity": "A. net ARG delta (rows)", "value": net_row_delta,
         "unit": "ARG loci (rows)",
         "note": "total full rows minus total chromosome-only rows; equals the sum of risk-tier deltas"},
        {"quantity": "B. matched gained loci", "value": len(gained),
         "unit": "unique (strain,gene)",
         "note": "strain x gene combinations present in full assembly but not chromosome-only"},
        {"quantity": "C. matched lost loci", "value": len(lost),
         "unit": "unique (strain,gene)",
         "note": "combinations in chromosome-only but not full (expect ~0; >0 = coordinate/collapse differences)"},
        {"quantity": "D. plasmid-associated gained loci", "value": len(plas_gained),
         "unit": "ARG loci (rows)",
         "note": "rows among the gained set located on a plasmid contig"},
        {"quantity": "A - B (row vs key gap)", "value": net_row_delta - len(gained),
         "unit": "ARG loci (rows)",
         "note": "genes found at multiple loci in one strain (e.g. a gene x2 on two plasmids) add rows but not new keys"},
    ])
    accounting.to_csv(out / "full_vs_chronly_locus_accounting.csv", index=False)

    # ---- markdown summary ---------------------------------------------- #
    md = ["# Full-assembly vs chromosome-only comparison (matched SCALE50)\n",
          "_Same 50 strains; the only difference is `full_assembly` true vs false._\n",
          "## Metric comparison\n", df_to_md(comp),
          "\n## Locus accounting (A/B/C/D)\n",
          "These count different units, so they are reported separately:\n",
          df_to_md(accounting),
          f"\nThe net **+{net_row_delta} rows** equals the sum of risk-tier deltas; "
          f"the **{len(gained)} matched gained loci** count unique strain x gene "
          f"combinations, so the {net_row_delta - len(gained)}-row difference is "
          f"genes recovered at more than one plasmid locus in the same strain.\n",
          f"\n## ARG loci gained in full assembly ({len(gained)})\n",
          df_to_md(pd.DataFrame(gained, columns=["strain_id", "gene"]), max_rows=60),
          f"\n## ARG loci only seen in chromosome-only run ({len(lost)})\n",
          df_to_md(pd.DataFrame(lost, columns=["strain_id", "gene"]), max_rows=40),
          f"\n## Plasmid-associated ARGs gained in full assembly ({len(plas_gained)})\n"]
    if len(plas_gained):
        cols = [c for c in ["strain_id", "gene", "drug_class", "contig",
                            "risk_level"] if c in plas_gained.columns]
        md.append(df_to_md(plas_gained[cols], max_rows=60))
    else:
        md.append("_none_\n")
    (out / "full_vs_chronly_summary.md").write_text("\n".join(md), encoding="utf-8")
    print(f"wrote {out/'full_vs_chronly_comparison.csv'}, locus_accounting.csv and summary.md")
    print(f"  A net ARG delta (rows) = {net_row_delta:+d} | "
          f"B matched gained loci = {len(gained)} | "
          f"C matched lost loci = {len(lost)} | "
          f"D plasmid-associated gained = {len(plas_gained)}")


if __name__ == "__main__":
    main()
