"""Phase 6 — Figures.

Figure 1  ARG presence/absence heatmap (strains x genes)
Figure 2  stacked risk-level bar chart per species
Figure 3  gene-neighbourhood maps for the top high-risk ARGs
Figure 4  probiotic<->pathogen gene-sharing network
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from .mge_analysis import FEATURE_COLS
from .risk_classify import RISK_ORDER
from .utils import ensure_dir, get_logger

log = get_logger("figures")
sns.set_theme(style="whitegrid")

_RISK_COLORS = {"High": "#c0392b", "Medium": "#e67e22",
                "Low": "#f1c40f", "Negligible": "#95a5a6"}


def fig1_heatmap(arg_df: pd.DataFrame, metadata: pd.DataFrame, out: str | Path) -> Path:
    out = Path(out)
    if arg_df.empty:
        return _empty(out, "No ARGs detected")
    pa = (arg_df.assign(present=1)
                .pivot_table(index="strain_id", columns="gene",
                             values="present", aggfunc="max", fill_value=0))
    order = (metadata.sort_values(["group", "species", "strain_id"])
                     .strain_id.tolist())
    pa = pa.reindex([s for s in order if s in pa.index])
    h = max(3.0, 0.32 * len(pa) + 1.5)
    w = max(5.0, 0.45 * pa.shape[1] + 3)
    fig, ax = plt.subplots(figsize=(w, h))
    sns.heatmap(pa, cmap=["#f7f9fb", "#2c3e50"], cbar=False, linewidths=0.5,
                linecolor="#dfe6e9", ax=ax)
    ax.set_title("Figure 1 — ARG presence/absence across strains", loc="left",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("Resistance gene"); ax.set_ylabel("Strain")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    log.info("wrote %s", out.name)
    return out


def fig2_risk_bars(risk_tab: pd.DataFrame, out: str | Path) -> Path:
    out = Path(out)
    if risk_tab.empty or risk_tab.values.sum() == 0:
        return _empty(out, "No ARGs to classify")
    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(risk_tab) + 3), 5))
    bottom = np.zeros(len(risk_tab))
    x = np.arange(len(risk_tab))
    for level in RISK_ORDER:
        vals = risk_tab.get(level, pd.Series(0, index=risk_tab.index)).values
        ax.bar(x, vals, bottom=bottom, label=level,
               color=_RISK_COLORS[level], edgecolor="white")
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(risk_tab.index, rotation=35, ha="right",
                       fontstyle="italic", fontsize=9)
    ax.set_ylabel("Number of ARGs")
    ax.set_title("Figure 2 — ARG transfer-risk profile by species", loc="left",
                 fontsize=12, fontweight="bold")
    ax.legend(title="Transfer risk", frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    log.info("wrote %s", out.name)
    return out


def fig3_gene_context(classified: pd.DataFrame, feat_df: pd.DataFrame,
                      out: str | Path, top_n: int = 4, flank_bp: int = 5000) -> Path:
    out = Path(out)
    high = classified[classified.risk_level == "High"]
    medium = classified[classified.risk_level == "Medium"]
    # show High first; if fewer than top_n, top up with Medium so the panel is useful
    sel = high.head(top_n)
    if len(sel) < top_n:
        sel = pd.concat([sel, medium.head(top_n - len(sel))], ignore_index=True)
    if sel.empty:
        return _empty(out, "No High/Medium-risk ARGs to map")

    # adaptive title that matches what is actually plotted
    tiers = set(str(t) for t in sel.risk_level.unique())
    if tiers == {"High"}:
        title = "Figure 3 — Genetic context of high-risk ARGs"
    elif "High" in tiers:
        title = "Figure 3 — Genetic context of high- and medium-risk ARGs"
    else:
        title = "Figure 3 — Genetic context of elevated-risk (Medium) ARGs"

    n = len(sel)
    fig, axes = plt.subplots(n, 1, figsize=(9, 1.9 * n + 0.5), squeeze=False)
    for ax, (_, arg) in zip(axes[:, 0], sel.iterrows()):
        win_start = arg.start - flank_bp
        win_end = arg.end + flank_bp
        ax.set_xlim(win_start, win_end); ax.set_ylim(0, 1)
        ax.axhline(0.5, color="#bdc3c7", lw=2, zorder=0)            # the contig
        _arrow(ax, arg.start, arg.end, arg.strand, "#c0392b", arg.gene, y=0.5)
        feats = feat_df[(feat_df.strain_id == arg.strain_id) &
                        (feat_df.contig == arg.contig)]
        for _, f in feats.iterrows():
            if f.end < win_start or f.start > win_end:
                continue
            color = {"insertion_sequence": "#2980b9",
                     "plasmid_replicon": "#27ae60",
                     "conjugation_marker": "#8e44ad"}.get(f.feature_type, "#7f8c8d")
            _arrow(ax, f.start, f.end, "+", color, f.feature, y=0.5)
        ax.set_yticks([])
        ax.set_title(f"{arg.strain_id}  |  {arg.gene}  |  {arg.contig}  "
                     f"({arg.risk_level}: {arg.risk_reason})", loc="left", fontsize=9)
        ax.set_xlabel("contig position (bp)", fontsize=8)
    fig.suptitle(title, x=0.01, ha="left",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(out, dpi=150); plt.close(fig)
    log.info("wrote %s", out.name)
    return out


def fig4_network(graph: nx.Graph, out: str | Path) -> Path:
    out = Path(out)
    if graph.number_of_edges() == 0:
        return _empty(out, "No shared-ARG edges")
    # show only connected nodes to keep it readable
    connected = [n for n in graph.nodes if graph.degree(n) > 0]
    g = graph.subgraph(connected)
    pos = nx.spring_layout(g, seed=42, k=0.9)
    fig, ax = plt.subplots(figsize=(9, 7))
    node_colors = ["#27ae60" if g.nodes[n].get("group") == "probiotic"
                   else "#c0392b" for n in g.nodes]
    nx.draw_networkx_nodes(g, pos, node_color=node_colors, node_size=650,
                           edgecolors="white", ax=ax)
    hi = [(u, v) for u, v, d in g.edges(data=True) if d.get("high_concern")]
    lo = [(u, v) for u, v, d in g.edges(data=True) if not d.get("high_concern")]
    nx.draw_networkx_edges(g, pos, edgelist=lo, edge_color="#95a5a6", width=1.5, ax=ax)
    nx.draw_networkx_edges(g, pos, edgelist=hi, edge_color="#c0392b", width=3.0, ax=ax)
    nx.draw_networkx_labels(g, pos, font_size=7, ax=ax)
    elabels = {(u, v): d.get("gene", "") for u, v, d in g.edges(data=True)}
    nx.draw_networkx_edge_labels(g, pos, edge_labels=elabels, font_size=6, ax=ax)
    ax.set_title("Figure 4 — Probiotic-pathogen ARG sharing network", loc="left",
                 fontsize=12, fontweight="bold")
    ax.text(0.01, -0.04, "green = probiotic   red = pathogen   "
            "thick red edge = high-concern shared gene",
            transform=ax.transAxes, fontsize=8, color="#555")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    log.info("wrote %s", out.name)
    return out


# --------------------------------------------------------------------------- #
def _arrow(ax, start, end, strand, color, label, y=0.5):
    width = end - start
    dx = width if strand == "+" else -width
    x0 = start if strand == "+" else end
    ax.annotate("", xy=(x0 + dx, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=6, shrinkA=0, shrinkB=0))
    ax.text((start + end) / 2, y + 0.18, label, ha="center", fontsize=7, color=color)


def _empty(out: Path, msg: str) -> Path:
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=12, color="#888")
    ax.axis("off"); fig.savefig(out, dpi=150); plt.close(fig)
    return out
