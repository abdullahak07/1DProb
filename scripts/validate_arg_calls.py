"""Independent ARG validation across multiple tools (AMRFinderPlus, ResFinder,
BLAST). Each locus gets a per-tool verdict and an aggregate status.

Status values (per locus):
  confirmed_multi_tool  - >=2 independent tools agree on the gene family
  confirmed_amrfinder   - only AMRFinderPlus confirms
  confirmed_resfinder   - only ResFinder confirms
  confirmed_blast_only  - only BLAST confirms
  discordant            - at least one tool reports a conflicting family, none confirming-by->=2
  not_checked           - no external tool covered this locus (builtin recheck only)

A builtin consistency recheck is recorded separately (builtin_recheck) and is
NEVER counted as external confirmation.

Inputs (drop in validation/, header-only template files are treated as absent):
  validation/amrfinderplus_results.tsv
  validation/resfinder_results.tsv
  validation/blast_results.tsv

Usage:
    py scripts\\validate_arg_calls.py
    py scripts\\validate_arg_calls.py --results results_scale50
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from _sota_common import df_to_md, locus_id, overlap, read_csv_safe, sota_dir

TOOLS = [("amrfinderplus_results.tsv", "AMRFinderPlus", "amrfinder"),
         ("resfinder_results.tsv", "ResFinder", "resfinder"),
         ("blast_results.tsv", "BLAST", "blast")]
WEAK_ID = 90.0


def _norm_token(name: str) -> str:
    """lowercase, drop leading 'bla', strip all non-alphanumerics (keep digits)."""
    g = str(name).strip().lower()
    g = re.sub(r"^bla", "", g)
    g = re.sub(r"[^a-z0-9]", "", g)
    return g


# Groups of names that denote the SAME gene/family across CARD / AMRFinderPlus /
# ResFinder nomenclature. Any member matches any other member.
SYNONYM_GROUPS = [
    ["CTX-M-14", "blaCTX-M-14"],
    ["KPC-2", "blaKPC-2"],
    ["ErmB", "erm(B)"],
    ["AAC6_Ie_APH2_Ia", "aac(6')-Ie/aph(2'')-Ia", "aac(6')-aph(2'')"],
    ["APH(3'')-Ib", "aph(3'')-Ib", "strA"],
    ["APH(6)-Id", "aph(6)-Id", "strB"],
    ["ANT(4')-Ia", "ant(4')-Ia", "aadD"],
    ["sul2"],
    ["catA1"],
    ["qacH"],
]
_ALIAS = {}
for _gi, _grp in enumerate(SYNONYM_GROUPS):
    for _n in _grp:
        _ALIAS[_norm_token(_n)] = _gi


def gene_stem(name: str) -> str:
    """Fallback stem (keeps trailing digits so sul1 != sul2)."""
    return _norm_token(name)


def families_match(a, b) -> bool:
    """True if a and b denote the same gene family.

    Empty / NaN / 'n/a' on either side is never a match (and never a conflict).
    Uses explicit synonym groups first, then a conservative token comparison.
    """
    ra, rb = str(a).strip().lower(), str(b).strip().lower()
    if ra in ("", "nan", "n/a", "na", "none") or rb in ("", "nan", "n/a", "na", "none"):
        return False
    ga, gb = _ALIAS.get(_norm_token(a)), _ALIAS.get(_norm_token(b))
    if ga is not None and gb is not None:
        return ga == gb
    sa, sb = _norm_token(a), _norm_token(b)
    if not sa or not sb:
        return False
    return sa == sb


def is_real_hit(gene) -> bool:
    return str(gene).strip().lower() not in ("", "nan", "n/a", "na", "none")


def _find_col(cols, *keys):
    low = {c.lower(): c for c in cols}
    for k in keys:
        for lc, orig in low.items():
            if k in lc:
                return orig
    return None


def has_data(path: Path) -> bool:
    """True only if the file exists AND has >=1 non-comment data row after a header."""
    if not path.exists():
        return False
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    return len(rows) >= 2


def load_external(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, comment="#")
    contig = _find_col(df.columns, "contig")
    start = _find_col(df.columns, "start", "position in contig")
    stop = _find_col(df.columns, "stop", "end")
    gene = _find_col(df.columns, "element symbol", "gene symbol", "resistance gene",
                     "best_hit_gene", "sseqid", "subject", "element", "gene")
    ident = _find_col(df.columns, "% identity", "pident", "identity")
    cov = _find_col(df.columns, "% coverage", "coverage", "qcov")
    lid = _find_col(df.columns, "locus_id", "qseqid", "query")
    # strain: avoid AMRFinder's 'Sequence name' description column
    strain = None
    for cand in ("strain_id", "strain", "isolate", "sample", "biosample"):
        strain = _find_col(df.columns, cand)
        if strain:
            break
    if strain is None:
        for c in df.columns:
            if c.strip().lower() == "name":  # AMRFinder --name column, exact match only
                strain = c
                break
    # coordinates: handle 'start..end' range strings (e.g. ResFinder 'Position in contig')
    if start and df[start].astype(str).str.contains(r"\.\.").any():
        parts = df[start].astype(str).str.extract(r"(\d+)\D+(\d+)")
        s_col = pd.to_numeric(parts[0], errors="coerce")
        e_col = pd.to_numeric(parts[1], errors="coerce")
    else:
        s_col = pd.to_numeric(df[start], errors="coerce") if start else 0
        e_col = pd.to_numeric(df[stop], errors="coerce") if stop else 0
    return pd.DataFrame({
        "locus_id": df[lid].astype(str) if lid else "",
        "strain_id": df[strain].astype(str) if strain else "",
        "contig": df[contig].astype(str) if contig else "",
        "start": s_col,
        "stop": e_col,
        "gene": df[gene] if gene else "",
        "identity": pd.to_numeric(df[ident], errors="coerce") if ident else float("nan"),
        "coverage": pd.to_numeric(df[cov], errors="coerce") if cov else float("nan"),
    })


def match(row, ext: pd.DataFrame):
    # ignore external rows without a real gene call
    ext = ext[ext.gene.apply(is_real_hit)]
    if not len(ext):
        return None
    rid = locus_id(row.strain_id, row.contig, row.start, row.end, row.gene)
    if "locus_id" in ext.columns:
        byid = ext[ext.locus_id == rid]
        if len(byid):
            return byid.iloc[0]
    sub = ext[ext.contig.astype(str) == str(row.contig)]
    # require same strain_id when the external file provides one
    if "strain_id" in ext.columns and (sub.strain_id.astype(str).str.len() > 0).any():
        same = sub[sub.strain_id.astype(str) == str(row.strain_id)]
        if len(same):
            sub = same
    best, best_ov = None, 0
    for _, e in sub.iterrows():
        ov = overlap(row.start, row.end, e.start, e.stop)
        if ov > best_ov:
            best, best_ov = e, ov
    return best if best_ov > 0 else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_scale50")
    ap.add_argument("--validation-dir", default="validation")
    args = ap.parse_args(argv)
    out = sota_dir()

    arg = read_csv_safe(Path(args.results) / "table2_master_arg.csv")
    if arg is None:
        raise SystemExit(f"table2_master_arg.csv not found in {args.results}")
    arg = arg.copy()

    vdir = Path(args.validation_dir)
    loaded = {}
    for fn, label, key in TOOLS:
        p = vdir / fn
        if has_data(p):
            try:
                loaded[key] = (label, load_external(p))
            except Exception as exc:  # noqa: BLE001
                print(f"  ! could not parse {p}: {exc}")
    present = list(loaded)
    print(f"external ARG tools imported: {[loaded[k][0] for k in present] or 'none'}")

    recs = []
    for _, row in arg.iterrows():
        rec = {"strain_id": row.strain_id, "builtin_gene": row.gene,
               "builtin_identity": row.get("pct_identity"),
               "builtin_coverage": row.get("pct_coverage"),
               "builtin_contig": row.contig, "builtin_start": row.start,
               "builtin_end": row.end, "drug_class": row.get("drug_class"),
               "risk_level": row.get("risk_level")}
        confirming, conflicting, notes = [], [], []
        for key in ("amrfinder", "resfinder", "blast"):
            hit_gene, hit_id, hit_cov, agree = "", float("nan"), float("nan"), False
            if key in loaded:
                e = match(row, loaded[key][1])
                if e is not None and is_real_hit(e.gene):
                    hit_gene = e.gene
                    hit_id = float(e.identity) if pd.notna(e.identity) else float("nan")
                    hit_cov = float(e.coverage) if pd.notna(e.coverage) else float("nan")
                    same_family = families_match(row.gene, e.gene)
                    # thresholds only apply if the columns exist (non-NaN)
                    id_ok = (pd.isna(hit_id) or hit_id >= 90.0)
                    cov_ok = (pd.isna(hit_cov) or hit_cov >= 80.0)
                    if same_family and id_ok and cov_ok:
                        agree = True
                        confirming.append(key)
                    elif same_family and not (id_ok and cov_ok):
                        # right gene but below identity/coverage: neither confirm nor conflict
                        notes.append(f"{key} same gene below threshold "
                                     f"(id={hit_id}, cov={hit_cov})")
                    else:
                        # a real, overlapping hit naming a DIFFERENT family = conflict
                        conflicting.append(key)
                # e is None or blank gene -> tool did not cover this locus: ignore
            rec[f"{key}_hit"] = hit_gene
            rec[f"{key}_identity"] = hit_id
            rec[f"{key}_agree"] = agree
        bi = float(row.get("pct_identity", 0) or 0)
        bc = float(row.get("pct_coverage", 0) or 0)
        rec["builtin_recheck"] = "confirmed" if (bi >= 95 and bc >= 90) else "weak"

        if len(confirming) >= 2:
            status = "confirmed_multi_tool"
        elif len(confirming) == 1:
            status = {"amrfinder": "confirmed_amrfinder",
                      "resfinder": "confirmed_resfinder",
                      "blast": "confirmed_blast_only"}[confirming[0]]
        elif conflicting:
            status = "discordant"
        else:
            status = "not_checked"
            notes.append("no external tool imported (builtin recheck only)"
                         if not present else "no external hit covered this locus")
        rec["n_tools_confirm"] = len(confirming)
        rec["n_tools_conflict"] = len(conflicting)
        rec["validation_status"] = status
        rec["validation_note"] = "; ".join(notes)
        best_tool = confirming[0] if confirming else (conflicting[0] if conflicting else "")
        rec["validation_tool"] = {"amrfinder": "AMRFinderPlus", "resfinder": "ResFinder",
                                  "blast": "BLAST"}.get(best_tool, "none")
        rec["validation_best_hit"] = rec.get(f"{best_tool}_hit", "") if best_tool else ""
        recs.append(rec)

    val = pd.DataFrame(recs)
    val.to_csv(out / "validated_ARG_calls_SCALE50.csv", index=False)
    val[val.validation_status == "discordant"].to_csv(
        out / "discordant_ARG_calls_SCALE50.csv", index=False)
    hr = val[val.risk_level == "High"].copy()
    hr.to_csv(out / "high_risk_validation_summary.csv", index=False)

    counts = {
        "amrfinder_present": "amrfinder" in present,
        "resfinder_present": "resfinder" in present,
        "blast_present": "blast" in present,
        "amrfinder_confirmed": int(val.amrfinder_agree.sum()),
        "resfinder_confirmed": int(val.resfinder_agree.sum()),
        "blast_confirmed": int(val.blast_agree.sum()),
        "multi_tool_confirmed": int((val.validation_status == "confirmed_multi_tool").sum()),
        "discordant": int((val.validation_status == "discordant").sum()),
        "not_checked": int((val.validation_status == "not_checked").sum()),
    }
    pd.DataFrame([counts]).to_csv(out / "validation_tool_counts.csv", index=False)

    sc = val.validation_status.value_counts().to_dict()
    md = ["# ARG validation summary (multi-tool)\n",
          f"External identity tools imported: "
          f"**{', '.join(loaded[k][0] for k in present) if present else 'none'}**.\n"]
    if not present:
        md.append("> No external tool imported. The builtin consistency recheck "
                  "confirms internal reproducibility only (NOT independent "
                  "validation); all loci are 'not_checked' externally.\n")
    md.append("## Status counts\n")
    md.append(df_to_md(pd.DataFrame([{"validation_status": k, "count": v}
                                     for k, v in sc.items()])))
    md.append("\n## Per-tool confirmed counts\n")
    md.append(df_to_md(pd.DataFrame([
        {"tool": "AMRFinderPlus", "imported": counts["amrfinder_present"], "confirmed_loci": counts["amrfinder_confirmed"]},
        {"tool": "ResFinder", "imported": counts["resfinder_present"], "confirmed_loci": counts["resfinder_confirmed"]},
        {"tool": "BLAST", "imported": counts["blast_present"], "confirmed_loci": counts["blast_confirmed"]},
        {"tool": "multi-tool (>=2 agree)", "imported": "-", "confirmed_loci": counts["multi_tool_confirmed"]},
    ])))
    md.append("\n## High-risk loci\n")
    hcols = ["strain_id", "builtin_gene", "amrfinder_hit", "resfinder_hit",
             "blast_hit", "n_tools_confirm", "validation_status"]
    md.append(df_to_md(hr[[c for c in hcols if c in hr.columns]]))
    disc = val[val.validation_status == "discordant"]
    md.append(f"\n## Discordant loci ({len(disc)})\n")
    md.append(df_to_md(disc[[c for c in hcols if c in disc.columns]], max_rows=50)
              if len(disc) else "_none_\n")
    md.append("\n_Allele names may be asserted only for loci confirmed by an "
              "independent tool; everything else stays 'best CARD match'._\n")
    (out / "validation_summary.md").write_text("\n".join(md), encoding="utf-8")

    print("status:", sc)
    print("per-tool confirmed:", {k: v for k, v in counts.items() if "confirmed" in k})


if __name__ == "__main__":
    main()
