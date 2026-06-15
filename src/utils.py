"""Shared helpers: logging, FASTA I/O, sequence utilities."""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Dict, Iterator, Tuple

_LOG_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"


def get_logger(name: str = "amr") -> logging.Logger:
    """Return a configured logger (idempotent)."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(_LOG_FORMAT, datefmt="%H:%M:%S"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


def ensure_dir(path: str | Path) -> Path:
    """Create a directory (and parents) if missing; return it as Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# --------------------------------------------------------------------------- #
# FASTA                                                                        #
# --------------------------------------------------------------------------- #
def read_fasta(path: str | Path) -> Dict[str, str]:
    """Minimal, dependency-free multi-FASTA reader.

    Returns {record_id: sequence}. record_id is the first whitespace-delimited
    token of the header (matching how BLAST/abricate report SEQUENCE ids).
    """
    seqs: Dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        seqs[header] = "".join(chunks).upper()
    return seqs


def write_fasta(path: str | Path, records: Dict[str, str], width: int = 70) -> None:
    """Write {id: sequence} to a wrapped multi-FASTA file."""
    with open(path, "w") as fh:
        for rid, seq in records.items():
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def iter_fasta_files(folder: str | Path, suffixes=(".fna", ".fasta", ".fa")) -> Iterator[Path]:
    """Yield FASTA files under *folder*, sorted, by suffix."""
    folder = Path(folder)
    for p in sorted(folder.rglob("*")):
        if p.suffix.lower() in suffixes:
            yield p


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse complement (IUPAC-lite: handles N)."""
    return seq.translate(_COMP)[::-1]


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return 100.0 * gc / len(seq)


def n_contigs(seqs: Dict[str, str]) -> int:
    return len(seqs)


def genome_length(seqs: Dict[str, str]) -> int:
    return sum(len(s) for s in seqs.values())
