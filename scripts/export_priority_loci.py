"""External-validation prep: export the priority ARG loci as nucleotide FASTA,
write a metadata table, and drop blank import templates for AMRFinderPlus /
ResFinder / BLAST.

Priority sets:
  - high_risk_loci_SCALE50.fasta          : all raw High-risk loci
  - curated_high_medium_loci_SCALE50.fasta: curated acquired-cargo High/Medium loci

FASTA headers are the locus_id (strain|contig|start-end|gene) so results can be
matched back automatically by validate_arg_calls.py.

Usage:
    py scripts\\export_priority_loci.py
    py scripts\\export_priority_loci.py --results results_scale50 --genome-dir data\\genomes
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import pandas as pd

from _sota_common import locus_id, read_csv_safe

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                           "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def read_fasta(path: str) -> dict:
    """accession (first header token) -> sequence."""
    seqs, acc, buf = {}, None, []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def build_contig_index(genome_dir: str) -> dict:
    idx = {}
    pats = ["*.fna", "*.fasta", "*.fa"]
    files = []
    for p in pats:
        files += glob.glob(os.path.join(genome_dir, "**", p), recursive=True)
    for f in files:
        try:
            for acc, seq in read_fasta(f).items():
                idx.setdefault(acc, seq)
        except Exception:  # noqa: BLE001
            continue
    return idx


def extract(idx, contig, start, end, strand):
    seq = idx.get(str(contig))
    if seq is None:
        return None
    s, e = int(float(start)), int(float(end))
    if s > e:
        s, e = e, s
    sub = seq[max(0, s - 1):e]
    if str(strand) == "-":
        sub = revcomp(sub)
    return sub


TEMPLATES = {
    "amrfinderplus_results_template.tsv":
        "# Paste AMRFinderPlus output, or fill these columns. Keep the header row.\n"
        "# Match is by Contig id + coordinate overlap (whole-genome run) OR by locus_id.\n"
        "Contig id\tStart\tStop\tGene symbol\t% Identity to reference sequence\t% Coverage of reference sequence\tlocus_id\n",
    "resfinder_results_template.tsv":
        "# Paste ResFinder 'Resistance gene' table, or fill these columns. Keep the header row.\n"
        "Resistance gene\tIdentity\tCoverage\tContig\tPosition in contig\tlocus_id\n",
    "blast_results_template.tsv":
        "# For web/standalone BLAST of the exported loci FASTA: one row per query locus.\n"
        "# locus_id MUST equal the FASTA header you submitted (strain|contig|start-end|gene).\n"
        "locus_id\tbest_hit_gene\tpident\tcoverage\tevalue\tdatabase\tnotes\n",
}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_scale50")
    ap.add_argument("--sota", default="results_sota")
    ap.add_argument("--genome-dir", default="data/genomes")
    ap.add_argument("--out", default="validation")
    args = ap.parse_args(argv)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    arg = read_csv_safe(Path(args.results) / "table2_master_arg.csv")
    if arg is None:
        raise SystemExit(f"table2_master_arg.csv not found in {args.results}")
    cur_hm = read_csv_safe(Path(args.sota) / "curated_high_medium_SCALE50.csv")

    high = arg[arg.risk_level == "High"].copy()
    high["_set"] = "high_risk"
    if cur_hm is not None:
        chm = cur_hm.copy()
    else:
        chm = arg[arg.risk_level.isin(["High", "Medium"])].copy()
        print("  ! curated_high_medium_SCALE50.csv missing; using raw High/Medium")
    chm["_set"] = "curated_high_medium"

    idx = build_contig_index(args.genome_dir)
    print(f"indexed {len(idx)} contigs from {args.genome_dir}")

    def write_fasta(df, fname):
        n_ok = 0
        recs = []
        with open(outdir / fname, "w", encoding="utf-8") as fh:
            for _, r in df.iterrows():
                lid = locus_id(r.strain_id, r.contig, r.start, r.end, r.gene)
                seq = extract(idx, r.contig, r.start, r.end, r.get("strand", "+"))
                found = seq is not None and len(seq) > 0
                if found:
                    fh.write(f">{lid}\n")
                    for i in range(0, len(seq), 70):
                        fh.write(seq[i:i + 70] + "\n")
                    n_ok += 1
                recs.append({"locus_id": lid, "strain_id": r.strain_id, "gene": r.gene,
                             "drug_class": r.get("drug_class"), "contig": r.contig,
                             "start": r.start, "end": r.end, "strand": r.get("strand", "+"),
                             "raw_risk_level": r.get("raw_risk_level", r.get("risk_level")),
                             "curated_risk_level": r.get("curated_risk_level", ""),
                             "builtin_identity": r.get("pct_identity"),
                             "builtin_coverage": r.get("pct_coverage"),
                             "sequence_length": len(seq) if found else 0,
                             "sequence_found": found, "priority_set": df["_set"].iloc[0],
                             "fasta_file": fname})
        return recs, n_ok

    hr_recs, hr_ok = write_fasta(high, "high_risk_loci_SCALE50.fasta")
    chm_recs, chm_ok = write_fasta(chm, "curated_high_medium_loci_SCALE50.fasta")

    meta = pd.DataFrame(hr_recs + chm_recs)
    meta.to_csv(outdir / "priority_validation_metadata.csv", index=False)

    for fname, content in TEMPLATES.items():
        (outdir / fname).write_text(content, encoding="utf-8")

    miss = int((~meta.sequence_found).sum())
    print(f"high-risk loci: {len(high)} ({hr_ok} sequences written)")
    print(f"curated High/Medium loci: {len(chm)} ({chm_ok} sequences written)")
    print(f"metadata -> {outdir/'priority_validation_metadata.csv'} "
          f"({len(meta)} rows, {miss} without sequence)")
    print("templates written:", ", ".join(TEMPLATES))
    if miss:
        print(f"  ! {miss} loci had no genome sequence found - check --genome-dir "
              f"or that the contig accessions are present as *_full.fna")


if __name__ == "__main__":
    main()
