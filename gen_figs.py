"""Generate the figures that are fully determined by the reported SCALE50
aggregate numbers (Fig 1, 2, 3, 4 and Supplementary S1-S3). Vector PDF output.
Figures requiring per-locus CSVs (strain heatmap, network, burden boxplot,
drug-class) are produced separately by make_figures.py against the result files.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.titlesize": 11,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.bbox": "tight",
    "figure.dpi": 200,
})
CHR, FULL = "#9ecae1", "#08519c"     # chromosome-only / full assembly
PRO, PAT = "#2a9d8f", "#e76f51"      # probiotic / pathogen
HIGH, MED, LOW = "#b2182b", "#ef8a62", "#bdbdbd"
OUT = "figs/"


# ---------------------------------------------------------------- Fig 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 4.6)); ax.axis("off")
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)

    def box(x, y, w, h, text, fc="#f0f0f0", ec="#333333"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                                    fc=fc, ec=ec, lw=1.1))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.6)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=11, lw=1.0, color="#555555"))

    box(3.2, 8.7, 3.6, 1.0, "50 reference genomes\n(25 probiotic-associated, 25 pathogen)", fc="#dbe9f6")
    box(0.5, 7.0, 4.0, 1.0, "Full assembly\n(chromosome + plasmids)", fc="#e8f0e3")
    box(5.5, 7.0, 4.0, 1.0, "Chromosome-only\n(matched panel)", fc="#e8f0e3")
    arrow(4.4, 8.7, 2.5, 8.0); arrow(5.6, 8.7, 7.5, 8.0)
    box(0.4, 5.3, 4.4, 1.1, "ARG detection (CARD)\n$\\geq$95% id, $\\geq$80% cov; overlap collapse", fc="#f7eedd")
    box(5.2, 5.3, 4.4, 1.1, "MGE context\nPlasmidFinder replicons; ISfinder IS", fc="#f7eedd")
    arrow(2.5, 7.0, 2.5, 6.4); arrow(7.5, 7.0, 7.4, 6.4); arrow(4.8, 5.85, 5.2, 5.85)
    box(2.4, 3.6, 5.2, 1.1, "Four-tier mobile-context risk + curation\n(acquired vs intrinsic/efflux/biocide)", fc="#efe3f0")
    arrow(2.6, 5.3, 4.0, 4.7); arrow(7.4, 5.3, 6.0, 4.7)
    box(0.4, 1.8, 4.4, 1.1, "External identity validation\nAMRFinderPlus, ResFinder, BLAST", fc="#f6dfe0")
    box(5.2, 1.8, 4.4, 1.1, "Mobility overlay\nMOB-suite (plasmid/MGE)", fc="#f6dfe0")
    arrow(4.0, 3.6, 2.6, 2.9); arrow(6.0, 3.6, 7.4, 2.9)
    box(2.7, 0.2, 4.6, 1.0, "Comparative outputs, statistics,\nfigures and supplementary tables", fc="#e0e0e0")
    arrow(2.6, 1.8, 4.4, 1.2); arrow(7.4, 1.8, 5.6, 1.2)
    fig.savefig(OUT + "fig1_workflow.pdf"); plt.close(fig)


# ---------------------------------------------------------------- Fig 2
def fig2():
    fig, axs = plt.subplots(1, 3, figsize=(8.4, 3.3))
    x = np.arange(2); w = 0.62

    cats = ["ARG loci", "High", "Medium", "Low"]
    chr_ = [338, 4, 21, 313]; full = [373, 15, 43, 315]
    xc = np.arange(len(cats))
    axs[0].bar(xc - w/4, chr_, width=w/2, color=CHR, label="Chromosome-only")
    axs[0].bar(xc + w/4, full, width=w/2, color=FULL, label="Full assembly")
    for i, (a, b) in enumerate(zip(chr_, full)):
        axs[0].text(i - w/4, a + 6, str(a), ha="center", fontsize=7)
        axs[0].text(i + w/4, b + 6, str(b), ha="center", fontsize=7)
    axs[0].set_xticks(xc); axs[0].set_xticklabels(cats, rotation=20, ha="right")
    axs[0].set_ylabel("Loci"); axs[0].set_title("(a) ARG loci by risk tier")
    axs[0].legend(fontsize=7, frameon=False)

    axs[1].bar(x, [0, 32], width=w, color=[CHR, FULL])
    for i, v in enumerate([0, 32]):
        axs[1].text(i, v + 0.6, str(v), ha="center", fontsize=8)
    axs[1].set_xticks(x); axs[1].set_xticklabels(["Chr-only", "Full"])
    axs[1].set_title("(b) Plasmid replicons"); axs[1].set_ylim(0, 36)

    axs[2].bar(x, [1003, 1116], width=w, color=[CHR, FULL])
    for i, v in enumerate([1003, 1116]):
        axs[2].text(i, v + 8, str(v), ha="center", fontsize=8)
    axs[2].set_xticks(x); axs[2].set_xticklabels(["Chr-only", "Full"])
    axs[2].set_title("(c) IS features"); axs[2].set_ylim(0, 1240)
    fig.tight_layout(); fig.savefig(OUT + "fig2_compare.pdf"); plt.close(fig)


# ---------------------------------------------------------------- Fig 3
def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
    axs[0].bar([0, 1], [0.52, 12.32], color=[PRO, PAT], width=0.6)
    for i, v in enumerate([0.52, 12.32]):
        axs[0].text(i, v + 0.3, f"{v:.2f}", ha="center", fontsize=9)
    axs[0].set_xticks([0, 1]); axs[0].set_xticklabels(["Probiotic-\nassociated", "Pathogen/\ncomparator"])
    axs[0].set_ylabel("Mean ARG loci per strain")
    axs[0].set_title("(a) Per-strain ARG burden")
    axs[0].set_ylim(0, 14)
    axs[0].text(0.5, 13.2, "Mann\u2013Whitney $U$=606.5, $P$<0.001",
                ha="center", fontsize=7.5, style="italic")

    axs[1].bar([0, 1], [0, 15], color=[PRO, PAT], width=0.6)
    for i, v in enumerate([0, 15]):
        axs[1].text(i, v + 0.3, str(v), ha="center", fontsize=9)
    axs[1].set_xticks([0, 1]); axs[1].set_xticklabels(["Probiotic-\nassociated", "Pathogen/\ncomparator"])
    axs[1].set_ylabel("High-risk ARG loci")
    axs[1].set_title("(b) High-risk localisation"); axs[1].set_ylim(0, 17)
    fig.tight_layout(); fig.savefig(OUT + "fig3_burden.pdf"); plt.close(fig)


# ---------------------------------------------------------------- Fig 4
def fig4():
    # 15 High-risk loci (all in pathogen genomes). status: 2=multi-tool,1=confirmed,0=flagged
    loci = [
        ("E. faecalis V583", "erm(B)", "MLS", 2),
        ("E. faecalis V583", "qacH \u2020", "biocide", 0),
        ("E. faecalis V583", "aac(6')-Ie/aph(2'')-Ia", "aminoglycoside", 1),
        ("K. pneumoniae HS11286", "bla$_{CTX-M-14}$ (p)", "beta-lactam", 2),
        ("K. pneumoniae HS11286", "bla$_{CTX-M-14}$ (p)", "beta-lactam", 2),
        ("K. pneumoniae MGH78578", "aph(3'')-Ib", "aminoglycoside", 1),
        ("K. pneumoniae MGH78578", "aph(6)-Id", "aminoglycoside", 1),
        ("K. pneumoniae MGH78578", "sul2", "sulfonamide", 2),
        ("K. pneumoniae MGH78578", "aph(3')-Ia", "aminoglycoside", 1),
        ("K. pneumoniae MGH78578", "catA1", "phenicol", 1),
        ("E. cloacae ATCC 13047", "sul2", "sulfonamide", 1),
        ("S. aureus N315", "ant(4')-Ia", "aminoglycoside", 1),
        ("S. aureus Mu50", "ant(4')-Ia", "aminoglycoside", 1),
        ("S. aureus Mu50", "aac(6')-Ie/aph(2'')-Ia", "aminoglycoside", 1),
        ("A. baumannii ATCC 17978", "sul2", "sulfonamide", 2),
    ]
    classes = {"beta-lactam": "#b2182b", "aminoglycoside": "#2166ac", "sulfonamide": "#1b7837",
               "MLS": "#762a83", "phenicol": "#e08214", "biocide": "#888888"}
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    y = np.arange(len(loci))[::-1]
    for yi, (strain, gene, cls, st) in zip(y, loci):
        ax.scatter(0, yi, s=70, color=classes[cls], zorder=3)
        marker = {2: "$\\bullet\\bullet$", 1: "$\\bullet$", 0: "$\\circ$"}[st]
        lab = {2: "multi-tool", 1: "confirmed", 0: "flagged"}[st]
        ax.text(0.32, yi, f"{gene}", va="center", fontsize=8.5)
        ax.text(3.1, yi, lab, va="center", fontsize=7.5,
                color=("#08519c" if st == 2 else ("#333333" if st == 1 else "#999999")))
    # strain labels grouped on the left
    ax.set_yticks(y); ax.set_yticklabels([l[0] for l in loci], fontsize=7.5)
    ax.set_xlim(-0.4, 4.0); ax.set_xticks([])
    ax.set_title("High-risk ARG loci (all in pathogen/comparator genomes)")
    # legend for drug classes
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
                          markersize=8, label=k) for k, c in classes.items()]
    ax.legend(handles=handles, fontsize=7, frameon=False, ncol=3,
              loc="lower center", bbox_to_anchor=(0.5, -0.13))
    fig.tight_layout(); fig.savefig(OUT + "fig4_highrisk.pdf"); plt.close(fig)


# ---------------------------------------------------------------- Supp bars
def supp_curation():
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    cats = ["acquired\nmobile ARG", "efflux/\nregulatory", "intrinsic/\ncore", "ambiguous", "biocide/\nstress"]
    vals = [316, 39, 9, 5, 4]
    ax.bar(range(5), vals, color=["#08519c", "#6baed6", "#9e9ac8", "#fdae6b", "#969696"])
    for i, v in enumerate(vals):
        ax.text(i, v + 4, str(v), ha="center", fontsize=8)
    ax.set_xticks(range(5)); ax.set_xticklabels(cats, fontsize=7.5)
    ax.set_ylabel("Loci"); ax.set_title("Curation categories")
    fig.tight_layout(); fig.savefig(OUT + "figS1_curation.pdf"); plt.close(fig)


def supp_validation():
    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    cats = ["AMRFinder\nPlus", "ResFinder", "BLAST", "multi-tool\n($\\geq$2)", "discordant", "not\nchecked"]
    vals = [94, 68, 46, 68, 33, 224]
    cols = ["#1b7837", "#1b7837", "#1b7837", "#08519c", "#b2182b", "#bdbdbd"]
    ax.bar(range(6), vals, color=cols)
    for i, v in enumerate(vals):
        ax.text(i, v + 3, str(v), ha="center", fontsize=8)
    ax.set_xticks(range(6)); ax.set_xticklabels(cats, fontsize=7.5)
    ax.set_ylabel("Loci"); ax.set_title("External identity-validation status")
    fig.tight_layout(); fig.savefig(OUT + "figS2_validation.pdf"); plt.close(fig)


def supp_mobility():
    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    cats = ["on plasmid", "IS both\nsides", "mobile-context\npositive", "MOB-suite\nconfirmed", "externally\nconjugative"]
    vals = [33, 15, 37, 33, 0]
    ax.bar(range(5), vals, color=["#4292c6", "#807dba", "#41ab5d", "#08519c", "#d9d9d9"])
    for i, v in enumerate(vals):
        ax.text(i, v + 0.6, str(v), ha="center", fontsize=8)
    ax.set_xticks(range(5)); ax.set_xticklabels(cats, fontsize=7.5)
    ax.set_ylabel("High/Medium loci"); ax.set_title("Mobility-context status (n=58)")
    ax.set_ylim(0, 42)
    fig.tight_layout(); fig.savefig(OUT + "figS3_mobility.pdf"); plt.close(fig)


for f in (fig1, fig2, fig3, fig4, supp_curation, supp_validation, supp_mobility):
    f()
print("generated:", *sorted(__import__("os").listdir(OUT)))
