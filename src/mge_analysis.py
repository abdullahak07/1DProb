"""Phase 4 (part 1) — Mobile genetic element context.

For every ARG hit, determine its mobility context:

* on_plasmid     : a plasmid replicon (PlasmidFinder-style) sits on the same contig
* conjugative    : a conjugation marker (relaxase / MOB / T4SS) sits on the same contig
* is_upstream    : an insertion sequence lies within `flank_bp` upstream
* is_downstream  : an insertion sequence lies within `flank_bp` downstream
* is_flank_both  : IS elements on both sides (composite-transposon signature)

Feature detection reuses the same sequence-search engine as the ARG screen, so a
plasmid replicon / IS element is simply "another catalogue". For real production
runs the same contexts come from `abricate --db plasmidfinder` and
MobileElementFinder; parsers for those formats are provided too.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

from .amr_screen import (_build_index, _reciprocal_overlap, _search_indexed,
                         parse_abricate)
from .utils import get_logger, read_fasta

log = get_logger("mge")

FEATURE_COLS = ["strain_id", "contig", "start", "end", "feature",
                "feature_type", "pct_identity"]


def load_feature_catalogue(fasta_path: str | Path, feature_type: str) -> List[dict]:
    """Header convention:  >NAME|...   (only NAME is used)."""
    from .amr_screen import load_catalogue
    cat = load_catalogue(fasta_path)
    for c in cat:
        c["feature_type"] = feature_type
    return cat


def _collapse_features(df: pd.DataFrame, threshold: float = 0.8) -> pd.DataFrame:
    """Collapse overlapping same-type features per contig to one best call.

    A locus where 40 near-identical ISfinder variants all match should count as
    one insertion sequence, not 40. Best-first by identity then length.
    """
    if df.empty:
        return df.reset_index(drop=True)
    work = df.copy()
    work["_len"] = (work["end"] - work["start"] + 1).abs()
    work = work.sort_values(["pct_identity", "_len"], ascending=[False, False])
    keep = []
    for _, grp in work.groupby(["strain_id", "contig", "feature_type"], sort=False):
        kept_spans = []
        for idx, r in grp.iterrows():
            s, e = int(r.start), int(r.end)
            if any(_reciprocal_overlap(s, e, ks, ke, threshold) for ks, ke in kept_spans):
                continue
            kept_spans.append((s, e))
            keep.append(idx)
    return (work.loc[keep].drop(columns="_len")
                .sort_values(["strain_id", "contig", "start"]).reset_index(drop=True))


def detect_features(strain_id: str, genome_path: str | Path,
                    catalogue: List[dict], feature_type: str,
                    min_cov: float = 60.0, min_id: float = 80.0,
                    k: int = 13, index_step: int = 3,
                    collapse: bool = True) -> pd.DataFrame:
    """Detect MGE features using the SAME k-mer-indexed engine as the AMR screen.

    The genome index is built once per contig, then every catalogue feature is
    queried against it — fast enough to screen the ~5,700-sequence ISfinder set
    on Windows with no external tools.
    """
    contigs = read_fasta(genome_path)
    rows: List[dict] = []
    for cid, cseq in contigs.items():
        if len(cseq) < k:
            continue
        index = _build_index(cseq, k, index_step)
        for ref in catalogue:
            for h in _search_indexed(ref["seq"], cseq, index, k, min_cov, min_id):
                rows.append({
                    "strain_id": strain_id, "contig": cid,
                    "start": h["start"], "end": h["end"],
                    "feature": ref["gene"], "feature_type": feature_type,
                    "pct_identity": h["pct_identity"],
                })
        del index
    df = pd.DataFrame(rows, columns=FEATURE_COLS)
    if collapse and not df.empty:
        df = _collapse_features(df)
    return df


def detect_all_features(fasta_map: Dict[str, Path],
                        replicon_cat: List[dict] | None,
                        is_cat: List[dict] | None,
                        conj_cat: List[dict] | None,
                        min_cov: float = 60.0, min_id: float = 80.0) -> pd.DataFrame:
    frames = []
    specs = [("plasmid_replicon", replicon_cat),
             ("insertion_sequence", is_cat),
             ("conjugation_marker", conj_cat)]
    for sid, fp in fasta_map.items():
        for ftype, cat in specs:
            if cat:
                frames.append(detect_features(sid, fp, cat, ftype, min_cov, min_id))
    if not frames:
        return pd.DataFrame(columns=FEATURE_COLS)
    return pd.concat(frames, ignore_index=True)


def parse_plasmidfinder_abricate(text: str, strain_id: str) -> pd.DataFrame:
    """Parse `abricate --db plasmidfinder` output into FEATURE_COLS."""
    df = parse_abricate(text, strain_id, "plasmidfinder")
    if df.empty:
        return pd.DataFrame(columns=FEATURE_COLS)
    out = df.rename(columns={"gene": "feature"})[["strain_id", "contig", "start", "end", "feature"]].copy()
    out["feature_type"] = "plasmid_replicon"
    return out


def annotate_arg_context(arg_df: pd.DataFrame, feat_df: pd.DataFrame,
                         flank_bp: int = 5000) -> pd.DataFrame:
    """Add mobility-context columns to the ARG table."""
    arg = arg_df.copy()
    for col in ["on_plasmid", "conjugative", "is_upstream", "is_downstream", "is_flank_both"]:
        arg[col] = False
    if arg.empty:
        return arg

    feat_df = feat_df if feat_df is not None else pd.DataFrame(columns=FEATURE_COLS)

    def feats(strain, contig, ftype):
        return feat_df[(feat_df.strain_id == strain) &
                       (feat_df.contig == contig) &
                       (feat_df.feature_type == ftype)]

    for i, row in arg.iterrows():
        reps = feats(row.strain_id, row.contig, "plasmid_replicon")
        conj = feats(row.strain_id, row.contig, "conjugation_marker")
        iss = feats(row.strain_id, row.contig, "insertion_sequence")

        arg.at[i, "on_plasmid"] = bool(len(reps))
        arg.at[i, "conjugative"] = bool(len(conj))

        up = ((iss.end >= row.start - flank_bp) & (iss.end <= row.start)).any() if len(iss) else False
        down = ((iss.start <= row.end + flank_bp) & (iss.start >= row.end)).any() if len(iss) else False
        arg.at[i, "is_upstream"] = bool(up)
        arg.at[i, "is_downstream"] = bool(down)
        arg.at[i, "is_flank_both"] = bool(up and down)
    return arg
