"""Generate a small, fully deterministic synthetic dataset for the smoke test.

It plants ARGs / IS elements / plasmid replicons / conjugation markers into
random genomes so that EVERY code path and EVERY risk tier is exercised, with a
known expected outcome:

  strain                       gene     expected context        expected risk
  --------------------------------------------------------------------------
  Lacto_acidophilus_LA1  (pro) tetM     chromosomal, no MGE      Low
  Lacto_rhamnosus_GG     (pro) vanY     chromosomal, intrinsic   Negligible
  Bifido_longum_BL1      (pro) ermB     IS on both sides         High
  Entero_faecium_PRO     (pro) vanA     conjugative plasmid      High
  Saccharomyces_boul_SB  (pro) (none)   -                        -
  Entero_faecium_CLIN    (pat) vanA     conjugative plasmid      High   <- shares vanA with PRO
  Kleb_pneumoniae_KP     (pat) blaTEM   non-conjugative plasmid  Medium
  Ecoli_MDR              (pat) tetM     single IS upstream       Medium <- shares tetM with LA1
"""
from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, List

from src.utils import ensure_dir, write_fasta

RNG = random.Random(20260608)  # deterministic


def _rand(n: int) -> str:
    return "".join(RNG.choice("ACGT") for _ in range(n))


def _mutate(seq: str, k: int) -> str:
    s = list(seq)
    for _ in range(k):
        i = RNG.randrange(len(s))
        s[i] = RNG.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


# fixed reference sequences (the "truth") -------------------------------------
GENES = {
    "tetM":   _rand(900),
    "ermB":   _rand(750),
    "vanA":   _rand(1030),
    "vanY":   _rand(620),
    "blaTEM": _rand(860),
}
GENE_CLASS = {"tetM": "tetracycline", "ermB": "macrolide", "vanA": "glycopeptide",
              "vanY": "glycopeptide", "blaTEM": "beta-lactam"}

REPLICONS = {"rep_pVRE": _rand(450), "rep_pKP": _rand(430)}
IS_ELEMENTS = {"IS6": _rand(820), "IS256": _rand(790)}
CONJ_MARKERS = {"traA_relaxase": _rand(680)}


def _contig(parts: List[str], gaps: List[int]) -> str:
    """Stitch feature sequences together separated by random spacer DNA.

    len(gaps) == len(parts) + 1 (leading + between + trailing spacer).
    """
    out = [_rand(gaps[0])]
    for i, p in enumerate(parts):
        out.append(p)
        out.append(_rand(gaps[i + 1]))
    return "".join(out)


def write_catalogues(base: Path) -> Dict[str, Path]:
    cdir = ensure_dir(base / "catalogues")
    amr = {f"{g}|{GENE_CLASS[g]}|SYN:{g}|{g} resistance": s for g, s in GENES.items()}
    rep = {f"{n}|||{n}": s for n, s in REPLICONS.items()}
    iss = {f"{n}|||{n}": s for n, s in IS_ELEMENTS.items()}
    conj = {f"{n}|||{n}": s for n, s in CONJ_MARKERS.items()}
    paths = {
        "amr_fasta": cdir / "amr.fasta",
        "replicon_fasta": cdir / "replicons.fasta",
        "is_fasta": cdir / "insertion_sequences.fasta",
        "conjugation_fasta": cdir / "conjugation.fasta",
    }
    write_fasta(paths["amr_fasta"], amr)
    write_fasta(paths["replicon_fasta"], rep)
    write_fasta(paths["is_fasta"], iss)
    write_fasta(paths["conjugation_fasta"], conj)
    return paths


def write_genomes(base: Path) -> List[dict]:
    gdir = ensure_dir(base / "genomes")
    strains: List[dict] = []

    def emit(strain_id, species, group, accession, contigs, **meta):
        write_fasta(gdir / f"{accession}.fna", contigs)
        strains.append({"strain_id": strain_id, "species": species, "group": group,
                        "accession": accession,
                        "fasta": str(gdir / f"{accession}.fna"), **meta})

    # 1) Low: chromosomal tetM, no MGE
    emit("Lacto_acidophilus_LA1", "Lactobacillus acidophilus", "probiotic", "SYN_LA1",
         {"chromosome": _contig([GENES["tetM"]], [3000, 4000])},
         product_example="Culturelle", isolation_source="commercial product", qps_listed="yes")

    # 2) Negligible: intrinsic vanY in Lactobacillus, chromosomal, no MGE
    emit("Lacto_rhamnosus_GG", "Lactobacillus rhamnosus", "probiotic", "SYN_LGG",
         {"chromosome": _contig([GENES["vanY"]], [2500, 3500])},
         product_example="Culturelle", isolation_source="commercial product", qps_listed="yes")

    # 3) High: ermB flanked by IS on both sides (composite transposon)
    emit("Bifido_longum_BL1", "Bifidobacterium longum", "probiotic", "SYN_BL1",
         {"chromosome": _contig(
             [IS_ELEMENTS["IS256"], GENES["ermB"], IS_ELEMENTS["IS256"]],
             [2000, 1400, 1400, 2000])},
         product_example="Align", isolation_source="commercial product", qps_listed="yes")

    # 4) High: vanA on a conjugative plasmid (replicon + relaxase, same contig)
    emit("Entero_faecium_PRO", "Enterococcus faecium", "probiotic", "SYN_EFPRO",
         {"plasmid_pVRE": _contig(
             [REPLICONS["rep_pVRE"], CONJ_MARKERS["traA_relaxase"], GENES["vanA"]],
             [1500, 1200, 1200, 1500]),
          "chromosome": _rand(6000)},
         product_example="European probiotic", isolation_source="commercial product", qps_listed="no")

    # 5) yeast control: no bacterial ARGs at all
    emit("Saccharomyces_boul_SB", "Saccharomyces boulardii", "probiotic", "SYN_SB",
         {"chromosome": _rand(7000)},
         product_example="Florastor", isolation_source="commercial product", qps_listed="yes")

    # 6) pathogen High: clinical VRE, vanA on conjugative plasmid (SAME vanA seq as PRO)
    emit("Entero_faecium_CLIN", "Enterococcus faecium", "pathogen", "SYN_EFCLIN",
         {"plasmid_pVRE": _contig(
             [REPLICONS["rep_pVRE"], CONJ_MARKERS["traA_relaxase"], GENES["vanA"]],
             [1300, 1100, 1100, 1300]),
          "chromosome": _rand(6500)},
         isolation_source="clinical isolate")

    # 7) pathogen Medium: blaTEM on a non-conjugative plasmid (replicon only)
    emit("Kleb_pneumoniae_KP", "Klebsiella pneumoniae", "pathogen", "SYN_KP",
         {"plasmid_pKP": _contig([REPLICONS["rep_pKP"], GENES["blaTEM"]], [1500, 1300, 1500]),
          "chromosome": _rand(6000)},
         isolation_source="clinical isolate")

    # 8) pathogen Medium: tetM with a single upstream IS (SAME tetM seq as LA1)
    emit("Ecoli_MDR", "Escherichia coli", "pathogen", "SYN_ECOLI",
         {"chromosome": _contig([IS_ELEMENTS["IS6"], GENES["tetM"]], [2500, 1500, 8000])},
         isolation_source="clinical isolate")

    return strains


def generate(base: str | Path) -> dict:
    base = ensure_dir(base)
    cats = write_catalogues(Path(base))
    strains = write_genomes(Path(base))
    return {"catalogues": {k: str(v) for k, v in cats.items()},
            "strains": strains, "base": str(base)}


if __name__ == "__main__":
    info = generate("data/synthetic")
    print(f"generated {len(info['strains'])} synthetic genomes under {info['base']}")
