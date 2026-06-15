"""Generate results/final_report.md from the pipeline outputs.

No admin, no extra dependencies (pandas + stdlib only; no `tabulate`). Run from
the project root AFTER a pipeline run:

    py src\\generate_report.py
    py src\\generate_report.py --results results --config config\\config.yaml

It reads the CSV tables in results/ and writes a paper-ready Markdown report
including the chromosome-only vs full-assembly comparison, the High/Medium
findings, a pathogen-sharing summary, a curation note on core/intrinsic genes,
and the limitations section.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
from pathlib import Path
from typing import Optional

import pandas as pd


DEFAULT_BASELINE = {
    "label": "chromosome-only",
    "args": 66,
    "insertion_sequence": 188,
    "plasmid_replicon": 0,
    "conjugation_marker": 0,
    "Low": 63,
    "Medium": 3,
    "High": 0,
}


INTRINSIC_WATCHLIST = {
    "eptb": (
        "phosphoethanolamine transferase; core LPS-modification locus "
        "(polymyxin context). Chromosomal, not typically mobile cargo."
    ),
    "ugd": (
        "UDP-glucose 6-dehydrogenase; core housekeeping gene in LPS/"
        "polymyxin modification. Chromosomal, not acquired cargo."
    ),
    "emre": (
        "small multidrug resistance (SMR) efflux pump; intrinsic and "
        "species-core in many Enterobacteriaceae. Chromosomal."
    ),
}


def _read(results: Path, name: str) -> Optional[pd.DataFrame]:
    fp = results / name
    if fp.exists():
        try:
            return pd.read_csv(fp)
        except Exception:
            return None
    return None


def _safe_cell(v) -> str:
    """Convert any table cell to safe Markdown text."""

    if v is None:
        return "NA"

    try:
        if pd.isna(v):
            return "NA"
    except Exception:
        pass

    if isinstance(v, bool):
        return "True" if v else "False"

    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return "NA"
        if v.is_integer():
            return str(int(v))
        return f"{v:.1f}"

    text = str(v)
    if text.lower() == "nan":
        return "NA"

    return text.replace("|", "\\|")


def df_to_md(df: Optional[pd.DataFrame], max_rows: int = 0) -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table."""

    if df is None or df.empty:
        return "_none_\n"

    work = df.copy()

    if max_rows and len(work) > max_rows:
        work = work.head(max_rows)

    cols = [str(c) for c in work.columns]

    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]

    for _, row in work.iterrows():
        cells = [_safe_cell(row[c]) for c in work.columns]
        lines.append("| " + " | ".join(cells) + " |")

    return "\n".join(lines) + "\n"


def _feature_count(feat: Optional[pd.DataFrame], ftype: str) -> int:
    if feat is None or feat.empty or "feature_type" not in feat:
        return 0
    return int((feat["feature_type"] == ftype).sum())


def _watchlist_hits(gene: str) -> Optional[str]:
    g = str(gene).lower()
    for key, note in INTRINSIC_WATCHLIST.items():
        if key in g:
            return note
    return None


def _genome_available_count(meta: Optional[pd.DataFrame]) -> int:
    if meta is None or meta.empty or "genome_available" not in meta.columns:
        return 0

    vals = meta["genome_available"].astype(str).str.lower()
    return int(vals.isin(["true", "1", "yes"]).sum())


def build_report(
    results: Path,
    config_path: Optional[Path],
    baseline: dict,
) -> str:
    arg = _read(results, "table2_master_arg.csv")
    meta = _read(results, "table1_strain_metadata.csv")
    feat = _read(results, "mge_features.csv")
    findings = _read(results, "table3_high_risk_findings.csv")
    edges = _read(results, "gene_sharing_edges.csv")

    if arg is None:
        raise SystemExit("table2_master_arg.csv not found - run the pipeline first.")

    n_args = len(arg)
    n_is = _feature_count(feat, "insertion_sequence")
    n_rep = _feature_count(feat, "plasmid_replicon")
    n_conj = _feature_count(feat, "conjugation_marker")

    risk_counts = (
        arg["risk_level"].value_counts().to_dict()
        if "risk_level" in arg.columns
        else {}
    )

    full = {
        "label": "full-assembly",
        "args": n_args,
        "insertion_sequence": n_is,
        "plasmid_replicon": n_rep,
        "conjugation_marker": n_conj,
        "Low": int(risk_counts.get("Low", 0)),
        "Medium": int(risk_counts.get("Medium", 0)),
        "High": int(risk_counts.get("High", 0)),
    }

    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    md: list[str] = []

    md.append("# AMR Gene Cargo Screen - Final Report\n")
    md.append(f"_Generated {now} from `{results}/`._\n")

    if config_path and config_path.exists():
        try:
            import yaml

            cfg = yaml.safe_load(config_path.read_text())
            scr = cfg.get("screen", {})
            mge = cfg.get("mge", {})
            fetch = cfg.get("fetch", {})

            md.append("## Methods snapshot\n")
            md.append(
                f"- ARG backend: `{scr.get('backend')}`; "
                f"ARG thresholds: identity >= {scr.get('min_identity')}%, "
                f"coverage >= {scr.get('min_coverage')}%; "
                f"overlap-collapse: {scr.get('collapse_overlaps')} "
                f"(reciprocal >= {scr.get('overlap_threshold')}).\n"
                f"- MGE thresholds: identity >= {mge.get('min_identity')}%, "
                f"coverage >= {mge.get('min_coverage')}%; flank window "
                f"{mge.get('flank_bp')} bp.\n"
                f"- Full assembly fetch: {fetch.get('full_assembly')}.\n"
            )
        except Exception:
            pass

    md.append("## Dataset\n")

    if meta is not None and not meta.empty:
        keep = [
            c
            for c in [
                "strain_id",
                "species",
                "group",
                "accession",
                "n_contigs",
                "genome_length_bp",
                "assembly_quality",
                "genome_available",
            ]
            if c in meta.columns
        ]

        md.append(df_to_md(meta[keep]))

        total_strains = len(meta)
        available = _genome_available_count(meta)

        md.append(
            f"\n{total_strains} strains; {available} with genomes. "
            "`n_contigs` reflects replicons per assembly "
            "(chromosome + plasmids). Missing genomes are retained in the "
            "metadata table as unavailable.\n"
        )

        if "genome_available" in meta.columns:
            missing = meta[
                ~meta["genome_available"].astype(str).str.lower().isin(
                    ["true", "1", "yes"]
                )
            ]
            if not missing.empty:
                miss_cols = [
                    c
                    for c in ["strain_id", "species", "accession", "genome_available"]
                    if c in missing.columns
                ]
                md.append("### Missing genomes\n")
                md.append(df_to_md(missing[miss_cols]))

    md.append("## Chromosome-only vs full-assembly\n")

    comp = pd.DataFrame([baseline, full]).set_index("label")
    comp = comp[
        [
            "args",
            "insertion_sequence",
            "plasmid_replicon",
            "conjugation_marker",
            "Low",
            "Medium",
            "High",
        ]
    ]

    comp = comp.rename(
        columns={
            "args": "ARG loci",
            "insertion_sequence": "IS features",
            "plasmid_replicon": "plasmid replicons",
            "conjugation_marker": "conjugation markers",
        }
    )

    comp_md = comp.reset_index().rename(columns={"label": "run"})
    md.append(df_to_md(comp_md))

    md.append(
        "\nMoving from chromosome-only accessions to full assemblies recovers "
        "the plasmid compartment, which is where clinically important "
        "transferable resistance may sit. In this run, plasmid replicons go "
        f"from {baseline['plasmid_replicon']} to {full['plasmid_replicon']}, "
        f"and High-risk calls from {baseline['High']} to {full['High']}. "
        "Chromosome-only screening can therefore under-call transfer risk.\n"
    )

    md.append("## Risk-level distribution (full assembly)\n")

    rc = pd.DataFrame(
        {
            "risk_level": ["High", "Medium", "Low", "Negligible"],
            "count": [
                int(risk_counts.get("High", 0)),
                int(risk_counts.get("Medium", 0)),
                int(risk_counts.get("Low", 0)),
                int(risk_counts.get("Negligible", 0)),
            ],
        }
    )

    md.append(df_to_md(rc))

    md.append("## High-risk ARG loci\n")

    high = arg[arg["risk_level"] == "High"] if "risk_level" in arg.columns else arg.iloc[0:0]

    hi_cols = [
        c
        for c in [
            "strain_id",
            "gene",
            "drug_class",
            "contig",
            "start",
            "end",
            "pct_identity",
            "pct_coverage",
            "on_plasmid",
            "conjugative",
            "is_flank_both",
            "risk_reason",
        ]
        if c in arg.columns
    ]

    md.append(df_to_md(high[hi_cols]))

    if len(high) and "gene" in high.columns:
        dup = high["gene"].value_counts()
        multi = dup[dup > 1]
        if len(multi):
            md.append(
                "\nGenes found at multiple independent High-risk loci "
                "(distinct contig/coordinates, not automatically duplicates): "
                + ", ".join(f"{g} (x{int(n)})" for g, n in multi.items())
                + ".\n"
            )

    md.append("## Medium-risk ARG loci\n")

    med = arg[arg["risk_level"] == "Medium"] if "risk_level" in arg.columns else arg.iloc[0:0]

    if "on_plasmid" in med.columns and "gene" in med.columns:
        plas = med[med["on_plasmid"].astype(str).str.lower().isin(["true", "1", "yes"])]
        chrom = med[~med["on_plasmid"].astype(str).str.lower().isin(["true", "1", "yes"])]

        md.append(
            "**Plasmid-associated (Medium):** "
            + (", ".join(sorted(plas["gene"].astype(str).unique())) if len(plas) else "none")
            + "\n"
        )

        md.append(
            "\n**IS-associated / chromosomal context (Medium):** "
            + (", ".join(sorted(chrom["gene"].astype(str).unique())) if len(chrom) else "none")
            + "\n"
        )

    md.append("\n")
    md.append(df_to_md(med[hi_cols]))

    md.append("## Probiotic-pathogen gene sharing\n")

    if findings is not None and not findings.empty:
        md.append(
            f"{len(findings)} cross-group probiotic-pathogen shared ARG edge(s):\n"
        )
        md.append(df_to_md(findings, max_rows=30))
    else:
        md.append(
            "No probiotic-pathogen shared-ARG edges above the identity "
            "threshold in this run.\n"
        )

    if edges is not None and not edges.empty:
        md.append(f"\nTotal gene-sharing edges, all pairs: {len(edges)}.\n")

    md.append("## Curation note: core/intrinsic genes, interpret with caution\n")

    md.append(
        "Some detected determinants are core chromosomal or intrinsic "
        "resistance-associated genes rather than horizontally acquired mobile "
        "cargo. These should not be over-interpreted as transferable ARG cargo, "
        "even when an IS happens to lie nearby. Detected here:\n"
    )

    flagged = []

    if "gene" in arg.columns:
        for g in sorted(set(arg["gene"].astype(str))):
            note = _watchlist_hits(g)
            if note:
                ctx = arg[arg["gene"].astype(str) == g]
                tiers = (
                    ", ".join(sorted(set(str(t) for t in ctx["risk_level"])))
                    if "risk_level" in ctx.columns
                    else "NA"
                )
                flagged.append(
                    {
                        "gene": g,
                        "current_risk": tiers,
                        "rationale": note,
                    }
                )

    if flagged:
        md.append(df_to_md(pd.DataFrame(flagged)))
        md.append(
            "\nRecommended handling: report these separately from acquired "
            "cargo, or add them to `risk_classify.DEFAULT_INTRINSIC` so they are "
            "scored Negligible. They are retained in this report to keep the "
            "screen transparent.\n"
        )
    else:
        md.append("_No watchlist genes detected in this run._\n")

    md.append("## Limitations\n")

    md.append(
        "- Transfer risk is inferred from genetic context, including plasmid "
        "replicon, IS flanking, and conjugation marker co-localisation. It is "
        "not demonstrated by conjugation assays.\n"
        "- Conjugation markers are detected by nucleotide similarity, a coarse "
        "proxy. Relaxase typing in the field uses protein HMMs such as "
        "MOB-suite, CONJscan, or oriTfinder.\n"
        "- The screening engine is ungapped. Hits with indels relative to the "
        "reference allele may score lower identity than with a gapped aligner "
        "such as BLAST.\n"
        "- Allele-level names reflect the best database match, not definitive "
        "clinical typing.\n"
        "- Reported gene sets depend on CARD version, thresholds, and the "
        "overlap-collapse setting.\n"
        "- Missing genomes in scaled runs should be re-fetched or replaced "
        "before final publication-grade interpretation.\n"
    )

    md.append("## Figures\n")

    figdir = results / "figures"
    figures = [
        ("fig1_heatmap.png", "ARG presence/absence heatmap"),
        ("fig2_risk_bars.png", "Risk profile by species"),
        ("fig3_gene_context.png", "Genetic context of elevated-risk ARGs"),
        ("fig4_network.png", "ARG gene-sharing network"),
    ]

    for fn, cap in figures:
        if (figdir / fn).exists():
            md.append(f"- `figures/{fn}` - {cap}\n")

    return "\n".join(md)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--baseline", help="JSON file overriding chromosome-only baseline")
    args = ap.parse_args(argv)

    results = Path(args.results)

    baseline = DEFAULT_BASELINE
    if args.baseline and Path(args.baseline).exists():
        baseline = {
            **DEFAULT_BASELINE,
            **json.loads(Path(args.baseline).read_text()),
        }

    report = build_report(
        results,
        Path(args.config) if args.config else None,
        baseline,
    )

    out = results / "final_report.md"
    out.write_text(report, encoding="utf-8")

    print(f"wrote {out}  ({len(report.splitlines())} lines)")


if __name__ == "__main__":
    main()