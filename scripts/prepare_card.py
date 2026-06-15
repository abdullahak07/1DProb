"""Extract CARD's nucleotide ARG FASTA from the downloaded broadstreet tarball.

No admin, no 7-Zip, no conda — pure Python stdlib (tarfile handles .tar.bz2).

Download in your browser from https://card.mcmaster.ca/download
  -> "Download CARD Data"  (a file like  broadstreet-vX.X.X.tar.bz2)
Save it anywhere, then run (from the project root):

    py scripts\\prepare_card.py "C:\\Users\\you\\Downloads\\broadstreet-v3.3.0.tar.bz2"

It writes:  data\\catalogues\\card_nucleotide.fasta

ResFinder alternative: download the resfinder_db zip from
https://bitbucket.org/genomicepidemiology/resfinder_db (or its GitHub mirror),
unzip it, and run this script with --resfinder pointing at the unzipped folder;
it concatenates every *.fsa into the same output file.
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path

TARGET_MEMBER = "nucleotide_fasta_protein_homolog_model.fasta"
OUT = Path("data") / "catalogues" / "card_nucleotide.fasta"


def from_card_tarball(tar_path: Path) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:bz2") as tar:
        member = next((m for m in tar.getmembers()
                       if m.name.endswith(TARGET_MEMBER)), None)
        if member is None:
            sys.exit(f"ERROR: {TARGET_MEMBER} not found inside {tar_path}")
        with tar.extractfile(member) as src, open(OUT, "wb") as dst:
            data = src.read()
            dst.write(data)
    n = data.count(b">")
    print(f"wrote {OUT}  ({n} reference genes)")
    return n


def from_resfinder_dir(folder: Path) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(OUT, "w", encoding="utf-8") as dst:
        for fsa in sorted(Path(folder).rglob("*.fsa")):
            text = fsa.read_text(encoding="utf-8", errors="ignore")
            n += text.count(">")
            dst.write(text)
            if not text.endswith("\n"):
                dst.write("\n")
    print(f"wrote {OUT}  ({n} reference genes from ResFinder)")
    return n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("tarball", nargs="?", help="path to CARD broadstreet .tar.bz2")
    ap.add_argument("--resfinder", help="path to an unzipped resfinder_db folder")
    args = ap.parse_args(argv)
    if args.resfinder:
        from_resfinder_dir(Path(args.resfinder))
    elif args.tarball:
        from_card_tarball(Path(args.tarball))
    else:
        ap.error("give a CARD .tar.bz2 path, or --resfinder <folder>")


if __name__ == "__main__":
    main()
