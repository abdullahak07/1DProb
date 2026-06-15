"""End-to-end orchestrator.

Run:
    python -m src.pipeline --config config/config.yaml

Reads a YAML config, resolves genome FASTAs (download from NCBI or use local
files), then runs Phases 2-6 and writes every table + figure to the output dir.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import pandas as pd
import yaml

from . import amr_screen, figures, genome_fetch, mge_analysis, pathogen_compare, risk_classify
from .utils import ensure_dir, get_logger

log = get_logger("pipeline")


def load_config(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def resolve_genomes(cfg: dict) -> Dict[str, Path]:
    """Return {strain_id: fasta_path}, downloading if fetch.enabled."""
    fetch = cfg.get("fetch", {})
    strains = cfg["strains"]
    full = fetch.get("full_assembly", False)
    if fetch.get("enabled"):
        if full:
            log.info("fetching %d FULL assemblies (chromosome + plasmids)", len(strains))
            return genome_fetch.fetch_all_assemblies(
                strains, fetch["genome_dir"], fetch["email"], fetch.get("api_key"))
        log.info("fetching %d genomes from NCBI", len(strains))
        return genome_fetch.fetch_all(
            strains, fetch["genome_dir"], fetch["email"], fetch.get("api_key"))

    genome_dir = Path(fetch.get("genome_dir", "data/genomes"))
    fasta_map: Dict[str, Path] = {}
    for s in strains:
        if s.get("fasta"):
            fp = Path(s["fasta"])
        else:
            sub = genome_dir / s["group"] / s["species"].replace(" ", "_")
            # prefer a full-assembly file if one was downloaded
            candidates = [sub / f"{s['accession']}_full.fna",
                          genome_dir / f"{s['accession']}_full.fna",
                          sub / f"{s['accession']}.fna",
                          genome_dir / f"{s['accession']}.fna"]
            fp = next((c for c in candidates if c.exists()), candidates[-1])
        if fp.exists():
            fasta_map[s["strain_id"]] = fp
        else:
            log.warning("missing genome for %s (%s)", s["strain_id"], fp)
    return fasta_map


def _load_cat(path):
    return amr_screen.load_catalogue(path) if path and Path(path).exists() else None


def run_pipeline(cfg: dict) -> dict:
    out_dir = ensure_dir(cfg.get("output_dir", "results"))
    fig_dir = ensure_dir(Path(out_dir) / "figures")
    strains = cfg["strains"]
    cats = cfg.get("catalogues", {})
    scr = cfg.get("screen", {})
    backend = scr.get("backend", "builtin")
    min_id = scr.get("min_identity", 80.0)
    min_cov = scr.get("min_coverage", 60.0)

    # ---- Phase 2: genomes + metadata ------------------------------------ #
    fasta_map = resolve_genomes(cfg)
    metadata = genome_fetch.build_metadata(strains, fasta_map)
    metadata.to_csv(Path(out_dir) / "table1_strain_metadata.csv", index=False)
    log.info("Table 1: %d strains (%d with genomes)",
             len(metadata), int(metadata.get("genome_available", pd.Series()).sum()))

    # ---- Phase 3: AMR screen -------------------------------------------- #
    amr_cat = _load_cat(cats.get("amr_fasta"))
    arg_df = amr_screen.screen_all(
        fasta_map, backend,
        catalogue=amr_cat, catalogue_fasta=cats.get("amr_fasta"),
        abricate_db=scr.get("abricate_db", "card"),
        min_cov=min_cov, min_id=min_id,
        collapse_overlaps=scr.get("collapse_overlaps", True),
        overlap_threshold=scr.get("overlap_threshold", 0.8))
    log.info("Phase 3: %d total ARG hits", len(arg_df))

    # ---- Phase 4: MGE context + risk ------------------------------------ #
    mcfg = cfg.get("mge", {})
    mge_min_id = mcfg.get("min_identity", 80.0)
    mge_min_cov = mcfg.get("min_coverage", 60.0)
    feat_df = mge_analysis.detect_all_features(
        fasta_map,
        _load_cat(cats.get("replicon_fasta")),
        _load_cat(cats.get("is_fasta")),
        _load_cat(cats.get("conjugation_fasta")),
        min_cov=mge_min_cov, min_id=mge_min_id)
    feat_df.to_csv(Path(out_dir) / "mge_features.csv", index=False)
    if not feat_df.empty:
        log.info("Phase 4: %d MGE features (%s)", len(feat_df),
                 ", ".join(f"{t}:{n}" for t, n in
                           feat_df.feature_type.value_counts().items()))

    arg_ctx = mge_analysis.annotate_arg_context(
        arg_df, feat_df, flank_bp=mcfg.get("flank_bp", 5000))
    classified = risk_classify.classify(arg_ctx, metadata)
    classified.to_csv(Path(out_dir) / "table2_master_arg.csv", index=False)

    risk_tab = risk_classify.risk_summary_by_species(classified, metadata)
    risk_tab.to_csv(Path(out_dir) / "risk_summary_by_species.csv")

    # ---- Phase 5: pathogen comparison ----------------------------------- #
    arg_seq = pathogen_compare.extract_arg_sequences(classified, fasta_map)
    edges = pathogen_compare.build_sharing_edges(
        arg_seq, metadata,
        identity_threshold=cfg.get("compare", {}).get("identity_threshold", 95.0))
    edges.to_csv(Path(out_dir) / "gene_sharing_edges.csv", index=False)
    network = pathogen_compare.build_network(edges, metadata)
    findings = pathogen_compare.key_findings(edges)
    findings.to_csv(Path(out_dir) / "table3_high_risk_findings.csv", index=False)
    log.info("Phase 5: %d sharing edges, %d cross-group (probiotic<->pathogen)",
             len(edges), len(findings))

    # ---- Phase 6: figures ----------------------------------------------- #
    f1 = figures.fig1_heatmap(classified, metadata, fig_dir / "fig1_heatmap.png")
    f2 = figures.fig2_risk_bars(risk_tab, fig_dir / "fig2_risk_bars.png")
    f3 = figures.fig3_gene_context(classified, feat_df, fig_dir / "fig3_gene_context.png")
    f4 = figures.fig4_network(network, fig_dir / "fig4_network.png")

    return {
        "metadata": metadata, "arg": classified, "features": feat_df,
        "risk_table": risk_tab, "edges": edges, "findings": findings,
        "network": network,
        "figures": [f1, f2, f3, f4],
        "output_dir": Path(out_dir),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="Probiotic AMR-cargo pipeline")
    ap.add_argument("--config", required=True)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    res = run_pipeline(cfg)
    print("\n=== DONE ===")
    print("Outputs in:", res["output_dir"])
    for f in res["figures"]:
        print("  figure:", f)
    return res


if __name__ == "__main__":
    main()
