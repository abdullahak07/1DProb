"""Phase 5 — Pathogen comparison & gene-sharing network.

Extracts the actual nucleotide sequence of every detected ARG from its genome,
computes pairwise identity between strains carrying the same gene, and builds a
network whose edges are shared ARGs above an identity threshold. Probiotic<->
pathogen edges — especially for high-concern classes (glycopeptide, beta-lactam)
— are the publishable signal of cross-species transfer potential.
"""
from __future__ import annotations

from difflib import SequenceMatcher
from itertools import combinations
from pathlib import Path
from typing import Dict, List

import networkx as nx
import pandas as pd

from .utils import get_logger, read_fasta, revcomp

log = get_logger("compare")

HIGH_CONCERN_CLASSES = {"glycopeptide", "beta-lactam", "oxazolidinone",
                        "phenicol-oxazolidinone"}


def extract_arg_sequences(arg_df: pd.DataFrame,
                          fasta_map: Dict[str, Path]) -> pd.DataFrame:
    """Add the realised nucleotide sequence of each ARG hit to the table."""
    out = arg_df.copy()
    seqs: List[str] = []
    cache: Dict[str, Dict[str, str]] = {}
    for _, r in out.iterrows():
        fp = fasta_map.get(r.strain_id)
        if fp is None:
            seqs.append("")
            continue
        if r.strain_id not in cache:
            cache[r.strain_id] = read_fasta(fp)
        contig = cache[r.strain_id].get(r.contig, "")
        sub = contig[int(r.start) - 1:int(r.end)]
        if r.strand == "-":
            sub = revcomp(sub)
        seqs.append(sub)
    out["arg_sequence"] = seqs
    return out


def seq_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    if len(a) == len(b):
        m = sum(1 for x, y in zip(a, b) if x == y)
        return 100.0 * m / len(a)
    return 100.0 * SequenceMatcher(None, a, b).ratio()


def build_sharing_edges(arg_seq_df: pd.DataFrame, metadata: pd.DataFrame,
                        identity_threshold: float = 95.0) -> pd.DataFrame:
    """One row per pair of strains sharing the same gene above threshold."""
    grp = dict(zip(metadata.strain_id, metadata.group))
    rows: List[dict] = []
    for gene, sub in arg_seq_df.groupby("gene"):
        recs = sub.to_dict("records")
        for a, b in combinations(recs, 2):
            if a["strain_id"] == b["strain_id"]:
                continue
            ident = seq_identity(a["arg_sequence"], b["arg_sequence"])
            if ident < identity_threshold:
                continue
            ga, gb = grp.get(a["strain_id"], ""), grp.get(b["strain_id"], "")
            rows.append({
                "strain_a": a["strain_id"], "strain_b": b["strain_id"],
                "group_a": ga, "group_b": gb,
                "gene": gene, "drug_class": a["drug_class"],
                "pairwise_identity": round(ident, 2),
                "cross_group": ga != gb and "" not in (ga, gb),
                "high_concern": a["drug_class"] in HIGH_CONCERN_CLASSES,
            })
    return pd.DataFrame(rows, columns=[
        "strain_a", "strain_b", "group_a", "group_b", "gene", "drug_class",
        "pairwise_identity", "cross_group", "high_concern"])


def build_network(edges: pd.DataFrame, metadata: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    meta = metadata.set_index("strain_id")
    for sid, row in meta.iterrows():
        g.add_node(sid, species=row.get("species", ""), group=row.get("group", ""))
    for _, e in edges.iterrows():
        g.add_edge(e.strain_a, e.strain_b, gene=e.gene,
                   drug_class=e.drug_class, identity=e.pairwise_identity,
                   cross_group=bool(e.cross_group), high_concern=bool(e.high_concern))
    return g


def key_findings(edges: pd.DataFrame) -> pd.DataFrame:
    """Probiotic<->pathogen shared ARGs — sorted with the scariest first."""
    if edges.empty:
        return edges
    kf = edges[edges.cross_group].copy()
    kf = kf.sort_values(["high_concern", "pairwise_identity"],
                        ascending=[False, False])
    return kf.reset_index(drop=True)
