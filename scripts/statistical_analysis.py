"""Phase 5 - Statistics & enrichment (pandas + stdlib; scipy optional).

If scipy is missing, Fisher's exact test is computed with math.comb and the
group burden comparison uses a pure-Python permutation test.

Usage:
    py scripts\\statistical_analysis.py
    py scripts\\statistical_analysis.py --curated results_sota\\curated_ARG_calls_SCALE50.csv
"""
from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import pandas as pd

from _sota_common import df_to_md, read_csv_safe, sota_dir

try:
    from scipy.stats import fisher_exact as _scipy_fisher  # type: ignore
    from scipy.stats import mannwhitneyu as _scipy_mwu      # type: ignore
    HAVE_SCIPY = True
except Exception:  # noqa: BLE001
    HAVE_SCIPY = False


# --------------------------------------------------------------------------- #
# pure-Python stats fallbacks                                                  #
# --------------------------------------------------------------------------- #
def fisher_exact(a, b, c, d):
    """Two-sided Fisher exact test on a 2x2 table. Returns (odds_ratio, p)."""
    if HAVE_SCIPY:
        try:
            orr, p = _scipy_fisher([[a, b], [c, d]])
            return float(orr), float(p)
        except Exception:  # noqa: BLE001
            pass
    n = a + b + c + d
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d

    def hyper(x):
        # P(X=x) given margins
        return (math.comb(row1, x) * math.comb(row2, col1 - x)) / math.comb(n, col1)

    lo = max(0, col1 - row2)
    hi = min(col1, row1)
    p_obs = hyper(a)
    p = sum(hyper(x) for x in range(lo, hi + 1) if hyper(x) <= p_obs + 1e-12)
    orr = float("inf") if (b * c) == 0 else (a * d) / (b * c)
    return orr, min(1.0, p)


def burden_test(x, y, n_perm=20000, seed=42):
    """Group difference in per-strain counts. scipy Mann-Whitney or permutation."""
    if HAVE_SCIPY and len(x) and len(y):
        try:
            u, p = _scipy_mwu(x, y, alternative="two-sided")
            return {"test": "Mann-Whitney U", "stat": float(u), "p": float(p)}
        except Exception:  # noqa: BLE001
            pass
    # permutation test on difference of means
    rng = random.Random(seed)
    obs = abs((sum(x) / len(x) if x else 0) - (sum(y) / len(y) if y else 0))
    pool = list(x) + list(y)
    nx = len(x)
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        dx, dy = pool[:nx], pool[nx:]
        d = abs((sum(dx) / len(dx) if dx else 0) - (sum(dy) / len(dy) if dy else 0))
        if d >= obs - 1e-12:
            ge += 1
    return {"test": "permutation (pure-Python)", "stat": obs, "p": (ge + 1) / (n_perm + 1)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", default="results_sota/curated_ARG_calls_SCALE50.csv")
    ap.add_argument("--results", default="results_scale50")
    args = ap.parse_args(argv)
    out = sota_dir()

    arg = read_csv_safe(args.curated)
    used_curated = arg is not None
    if arg is None:
        arg = read_csv_safe(Path(args.results) / "table2_master_arg.csv")
        if arg is None:
            raise SystemExit("Need curated_ARG_calls_SCALE50.csv or table2_master_arg.csv")
        arg["include_in_mobile_ARG_analysis"] = True
        arg["raw_risk_level"] = arg["risk_level"]
    meta = read_csv_safe(Path(args.results) / "table1_strain_metadata.csv")
    edges = read_csv_safe(Path(args.results) / "gene_sharing_edges.csv")

    risk_col = "curated_risk_level" if "curated_risk_level" in arg.columns else "risk_level"
    arg["_risk"] = arg[risk_col].where(arg.get("include_in_mobile_ARG_analysis", True),
                                       other="Excluded-curated")
    # group lookup
    grp = {}
    if meta is not None:
        grp = dict(zip(meta.strain_id, meta.group))
    arg["group"] = arg.strain_id.map(grp).fillna("unknown")

    # ---- 1. burden per strain ------------------------------------------ #
    all_strains = meta[["strain_id", "group"]].copy() if meta is not None else \
        arg[["strain_id", "group"]].drop_duplicates()
    inc = arg[arg.get("include_in_mobile_ARG_analysis", True)]
    burden = (inc.groupby("strain_id").size().rename("mobile_arg_count"))
    hi = inc[inc._risk == "High"].groupby("strain_id").size().rename("high")
    me = inc[inc._risk == "Medium"].groupby("strain_id").size().rename("medium")
    tot = arg.groupby("strain_id").size().rename("total_arg_loci")
    bs = all_strains.set_index("strain_id").join([tot, burden, hi, me]).fillna(0)
    for c in ["total_arg_loci", "mobile_arg_count", "high", "medium"]:
        bs[c] = bs[c].astype(int)
    bs = bs.reset_index()
    bs.to_csv(out / "group_burden_by_strain.csv", index=False)

    pro = bs[bs.group == "probiotic"].mobile_arg_count.tolist()
    pat = bs[bs.group == "pathogen"].mobile_arg_count.tolist()
    bt = burden_test(pat, pro)

    # ---- 2-4. enrichment via Fisher ------------------------------------ #
    def enrich(mask_name, mask):
        a = int(((arg.group == "pathogen") & mask).sum())
        b = int(((arg.group == "pathogen") & ~mask).sum())
        c = int(((arg.group == "probiotic") & mask).sum())
        d = int(((arg.group == "probiotic") & ~mask).sum())
        orr, p = fisher_exact(a, b, c, d)
        return {"comparison": mask_name, "pathogen_yes": a, "pathogen_no": b,
                "probiotic_yes": c, "probiotic_no": d,
                "odds_ratio": round(orr, 3) if orr != float("inf") else "inf",
                "p_value": p}

    enr_rows = [enrich("High/Medium risk", arg._risk.isin(["High", "Medium"]))]
    if "on_plasmid" in arg.columns:
        enr_rows.append(enrich("plasmid-associated", arg.on_plasmid == True))  # noqa: E712
    if "is_flank_both" in arg.columns:
        enr_rows.append(enrich("IS-flanked both sides", arg.is_flank_both == True))  # noqa: E712
    risk_enr = pd.DataFrame(enr_rows)
    risk_enr.to_csv(out / "risk_enrichment.csv", index=False)

    # ---- strain-level enrichment (presence/absence per strain) --------- #
    # Locus-level Fisher is underpowered (probiotics have few ARG loci); the
    # biologically meaningful test is "how many STRAINS carry >=1 such locus".
    strain_rows = []
    n_pat_str = int((bs.group == "pathogen").sum())
    n_pro_str = int((bs.group == "probiotic").sum())

    def strain_enrich(label, col):
        a = int(((bs.group == "pathogen") & (bs[col] > 0)).sum())   # pathogen strains with >=1
        c = int(((bs.group == "probiotic") & (bs[col] > 0)).sum())  # probiotic strains with >=1
        orr, p = fisher_exact(a, n_pat_str - a, c, n_pro_str - c)
        return {"comparison": label,
                "pathogen_strains_positive": f"{a}/{n_pat_str}",
                "probiotic_strains_positive": f"{c}/{n_pro_str}",
                "odds_ratio": round(orr, 3) if orr != float("inf") else "inf",
                "p_value": p}

    strain_rows.append(strain_enrich("strains with >=1 High ARG", "high"))
    bs["_hm"] = bs["high"] + bs["medium"]
    strain_rows.append(strain_enrich("strains with >=1 High/Medium ARG", "_hm"))
    strain_enr = pd.DataFrame(strain_rows)
    strain_enr.to_csv(out / "strain_level_enrichment.csv", index=False)

    # ---- 5. drug-class distribution by group --------------------------- #
    dc = (arg.groupby(["drug_class", "group"]).size().unstack(fill_value=0)
          if "drug_class" in arg.columns else pd.DataFrame())
    if not dc.empty:
        dc = dc.reset_index()
    dc.to_csv(out / "drug_class_enrichment.csv", index=False)

    # ---- 8. network summary -------------------------------------------- #
    if edges is not None and len(edges):
        within_pro = int(((edges.group_a == "probiotic") & (edges.group_b == "probiotic")).sum())
        within_pat = int(((edges.group_a == "pathogen") & (edges.group_b == "pathogen")).sum())
        cross = int(edges.get("cross_group", pd.Series([], dtype=bool)).sum())
        genes = ", ".join(sorted(set(edges.gene))[:25])
        net = pd.DataFrame([{"total_edges": len(edges), "within_probiotic": within_pro,
                             "within_pathogen": within_pat,
                             "probiotic_pathogen": cross,
                             "n_edge_genes": edges.gene.nunique(),
                             "edge_genes_sample": genes}])
    else:
        net = pd.DataFrame([{"total_edges": 0, "within_probiotic": 0,
                             "within_pathogen": 0, "probiotic_pathogen": 0,
                             "n_edge_genes": 0, "edge_genes_sample": ""}])
    net.to_csv(out / "network_summary.csv", index=False)

    # ---- stats_summary ------------------------------------------------- #
    stats_rows = [
        {"analysis": "mobile ARG burden pathogen vs probiotic", "test": bt["test"],
         "statistic": round(bt["stat"], 4), "p_value": bt["p"],
         "detail": f"pathogen mean={_mean(pat):.2f}, probiotic mean={_mean(pro):.2f}"},
    ]
    for r in enr_rows:
        stats_rows.append({"analysis": f"locus enrichment: {r['comparison']}",
                           "test": "Fisher exact (2-sided, locus-level)",
                           "statistic": r["odds_ratio"], "p_value": r["p_value"],
                           "detail": f"pathogen {r['pathogen_yes']}/{r['pathogen_yes']+r['pathogen_no']} "
                                     f"vs probiotic {r['probiotic_yes']}/{r['probiotic_yes']+r['probiotic_no']}"})
    for r in strain_rows:
        stats_rows.append({"analysis": f"strain enrichment: {r['comparison']}",
                           "test": "Fisher exact (2-sided, strain-level)",
                           "statistic": r["odds_ratio"], "p_value": r["p_value"],
                           "detail": f"pathogen {r['pathogen_strains_positive']} "
                                     f"vs probiotic {r['probiotic_strains_positive']} strains positive"})
    stats = pd.DataFrame(stats_rows)
    stats.to_csv(out / "stats_summary.csv", index=False)

    md = ["# Statistical analysis (SCALE50)\n",
          f"Engine: {'scipy' if HAVE_SCIPY else 'pure-Python fallback (scipy not installed)'}. "
          f"ARG set: {'curated (intrinsic/biocide excluded)' if used_curated else 'raw table2'}.\n",
          "## Key tests\n", df_to_md(stats),
          "\n## Enrichment detail\n", df_to_md(risk_enr),
          "\n## Strain-level enrichment (presence/absence per strain - better powered)\n",
          df_to_md(strain_enr),
          "\n## Network summary\n", df_to_md(net),
          "\n## Per-strain burden (top of table)\n",
          df_to_md(bs.sort_values('mobile_arg_count', ascending=False), max_rows=25),
          "\n_Note: with all High-risk loci in pathogens and 0 probiotic-pathogen "
          "sharing edges, enrichment p-values describe a pathogen-concentrated "
          "resistome, not probiotic risk._\n"]
    (out / "stats_summary.md").write_text("\n".join(md), encoding="utf-8")

    print(f"stats written to {out}/  (scipy={HAVE_SCIPY})")
    print("burden test:", bt)


def _mean(x):
    return sum(x) / len(x) if x else 0.0


if __name__ == "__main__":
    main()
