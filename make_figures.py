"""Generate the four figures that require per-locus data, directly from your
SCALE50 result CSVs. Run from the project root after the pipeline + SOTA steps:

    py make_figures.py --results results_scale50 --sota results_sota --out figs

Produces (vector PDF):
    figS4_strain_heatmap.pdf   strain x risk-tier ARG counts, grouped by group
    figS5_sharing_network.pdf  gene-sharing network coloured by group
    figS6_burden_boxplot.pdf   per-strain ARG burden by group (boxplot + points)
    figS7_drugclass.pdf        drug-class distribution by group

These are intentionally NOT pre-rendered from aggregate numbers: they depend on
the real per-locus tables, so they must be built from your data.
"""
import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "serif", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "figure.dpi": 200})
PRO, PAT = "#2a9d8f", "#e76f51"
TIERS = ["High", "Medium", "Low"]
TIER_C = {"High": "#b2182b", "Medium": "#ef8a62", "Low": "#f0f0f0"}


def load(results, sota):
    arg = pd.read_csv(Path(results) / "table2_master_arg.csv")
    meta = pd.read_csv(Path(results) / "table1_strain_metadata.csv")
    edges = None
    p = Path(results) / "gene_sharing_edges.csv"
    if p.exists():
        edges = pd.read_csv(p)
    grp = dict(zip(meta.strain_id, meta.group))
    arg["group"] = arg.strain_id.map(grp)
    return arg, meta, edges


def heatmap(arg, meta, out):
    order = (meta.sort_values(["group", "strain_id"]).strain_id.tolist())
    mat = (arg.pivot_table(index="strain_id", columns="risk_level",
                           values="gene", aggfunc="count", fill_value=0)
           .reindex(index=order, columns=TIERS, fill_value=0))
    fig, ax = plt.subplots(figsize=(6, max(6, len(order) * 0.16)))
    im = ax.imshow(mat.values, aspect="auto", cmap="magma_r")
    ax.set_xticks(range(len(TIERS))); ax.set_xticklabels(TIERS)
    ax.set_yticks(range(len(order)))
    gmap = dict(zip(meta.strain_id, meta.group))
    ax.set_yticklabels([f"{s}" for s in order], fontsize=5)
    for t, s in enumerate(order):
        ax.get_yticklabels()[t].set_color(PRO if gmap.get(s) == "probiotic" else PAT)
    fig.colorbar(im, ax=ax, shrink=0.5, label="ARG loci")
    ax.set_title("ARG loci by risk tier (rows: strains, coloured by group)")
    fig.savefig(Path(out) / "figS4_strain_heatmap.pdf"); plt.close(fig)


def network(edges, meta, out):
    if edges is None or not len(edges):
        print("  no edges file; skipping network"); return
    import itertools
    gmap = dict(zip(meta.strain_id, meta.group))
    nodes = sorted(set(edges.strain_a) | set(edges.strain_b))
    ang = {n: 2 * np.pi * i / len(nodes) for i, n in enumerate(nodes)}
    pos = {n: (np.cos(a), np.sin(a)) for n, a in ang.items()}
    fig, ax = plt.subplots(figsize=(6, 6)); ax.axis("off")
    for _, e in edges.iterrows():
        x1, y1 = pos[e.strain_a]; x2, y2 = pos[e.strain_b]
        cross = gmap.get(e.strain_a) != gmap.get(e.strain_b)
        ax.plot([x1, x2], [y1, y2], color=("#d62728" if cross else "#bbbbbb"),
                lw=(1.4 if cross else 0.5), alpha=0.8, zorder=1)
    for n in nodes:
        x, y = pos[n]
        ax.scatter(x, y, s=30, color=(PRO if gmap.get(n) == "probiotic" else PAT), zorder=2)
    ax.set_title("Gene-sharing network ($\\geq$95% identity); red = cross-group")
    fig.savefig(Path(out) / "figS5_sharing_network.pdf"); plt.close(fig)


def burden(arg, meta, out):
    counts = (arg.groupby("strain_id").size()
              .reindex(meta.strain_id, fill_value=0))
    g = meta.set_index("strain_id").group
    pro = counts[g == "probiotic"].values
    pat = counts[g == "pathogen"].values
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    ax.boxplot([pro, pat], labels=["Probiotic-\nassociated", "Pathogen/\ncomparator"],
               widths=0.5, showfliers=False)
    for i, d in enumerate([pro, pat], start=1):
        ax.scatter(np.random.normal(i, 0.05, len(d)), d, s=12,
                   color=(PRO if i == 1 else PAT), alpha=0.7, zorder=3)
    ax.set_ylabel("ARG loci per strain")
    ax.set_title("Per-strain ARG burden by group")
    fig.savefig(Path(out) / "figS6_burden_boxplot.pdf"); plt.close(fig)


def drugclass(arg, out):
    if "drug_class" not in arg.columns:
        print("  no drug_class column; skipping"); return
    tab = (arg.groupby(["drug_class", "group"]).size().unstack(fill_value=0))
    tab = tab.sort_values(by=list(tab.columns), ascending=False)
    fig, ax = plt.subplots(figsize=(6, max(3, len(tab) * 0.3)))
    y = np.arange(len(tab))
    if "probiotic" in tab.columns:
        ax.barh(y - 0.2, tab.get("probiotic", 0), height=0.4, color=PRO, label="probiotic")
    if "pathogen" in tab.columns:
        ax.barh(y + 0.2, tab.get("pathogen", 0), height=0.4, color=PAT, label="pathogen")
    ax.set_yticks(y); ax.set_yticklabels(tab.index, fontsize=7)
    ax.invert_yaxis(); ax.set_xlabel("Loci"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("Drug-class distribution by group")
    fig.savefig(Path(out) / "figS7_drugclass.pdf"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_scale50")
    ap.add_argument("--sota", default="results_sota")
    ap.add_argument("--out", default="figs")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    arg, meta, edges = load(a.results, a.sota)
    heatmap(arg, meta, a.out)
    network(edges, meta, a.out)
    burden(arg, meta, a.out)
    drugclass(arg, a.out)
    print("done; figures written to", a.out)


if __name__ == "__main__":
    main()
