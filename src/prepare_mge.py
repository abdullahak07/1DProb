"""Build the three Phase-4 MGE catalogues with NO admin, NO conda, NO BLAST.

Pure Python standard library only (urllib / tarfile / zipfile / re). Run from
the project root:

    py scripts\\prepare_mge.py --email you@murdoch.edu.au

Creates:
    data/catalogues/replicons.fasta            (PlasmidFinder replicon types)
    data/catalogues/insertion_sequences.fasta  (ISfinder IS nucleotide seqs)
    data/catalogues/conjugation.fasta          (relaxase / transfer CDS, NCBI)

Sources (all fetched over plain HTTPS):
  * PlasmidFinder DB : bitbucket.org/genomicepidemiology/plasmidfinder_db
  * ISfinder mirror  : github.com/thanhleviet/ISfinder-sequences  (IS.fna)
  * Conjugation      : NCBI E-utilities, CDS whose product is a relaxase /
                       conjugal-transfer / MOB / T4SS protein

NOTE on conjugation: robust relaxase detection uses PROTEIN HMM profiles
(MOB-suite / CONJscan / oriTfinder), which need installs you can't do here.
Nucleotide matching of relaxases is therefore a COARSE proxy — treat the
`conjugative` flag as supporting evidence, not proof. The pipeline works fine
without this file (plasmid-borne ARGs then score Medium rather than High).

Flags:
  --email EMAIL     contact for NCBI (recommended; required for conjugation)
  --skip-conjugation   build only replicons + IS
  --conj-retmax N   how many NCBI records to mine for transfer CDS (default 40)
"""
from __future__ import annotations

import argparse
import io
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT_DIR = Path("data") / "catalogues"
UA = {"User-Agent": "Mozilla/5.0 (compatible; amr-pipeline/1.0)"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PLASMIDFINDER_BASE = "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/raw/master"
PLASMIDFINDER_FILES = ["enterobacteriales.fsa", "enterobacteriaceae.fsa", "gram_positive.fsa"]
ISFINDER_URL = "https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.fna"

CONJ_KEYWORDS = re.compile(
    r"relaxase|mobiliz|conjugal transfer|conjugative transfer|"
    r"\bmob[A-Z]?\b|\btra[A-Z]\b|\btrw[A-Z]\b|\btrb[A-Z]\b|"
    r"virb4|vird4|type iv (secretion|coupling)|t4cp|t4ss|"
    r"coupling protein|nickase",
    re.IGNORECASE)


def _get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _count(fasta_bytes: bytes) -> int:
    return fasta_bytes.count(b">")


# --------------------------------------------------------------------------- #
def build_replicons() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "replicons.fasta"
    total = 0
    with open(out, "wb") as fh:
        for name in PLASMIDFINDER_FILES:
            url = f"{PLASMIDFINDER_BASE}/{name}"
            try:
                data = _get(url)
            except Exception as exc:  # noqa: BLE001
                print(f"  - {name}: not available ({exc})")
                continue
            if not data.lstrip().startswith(b">"):
                print(f"  - {name}: not FASTA, skipped")
                continue
            fh.write(data)
            if not data.endswith(b"\n"):
                fh.write(b"\n")
            n = _count(data)
            total += n
            print(f"  + {name}: {n} replicon sequences")
    print(f"wrote {out}  ({total} total)")
    return total


def build_insertion_sequences() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "insertion_sequences.fasta"
    data = _get(ISFINDER_URL)
    out.write_bytes(data)
    n = _count(data)
    print(f"wrote {out}  ({n} IS sequences)")
    return n


# --------------------------------------------------------------------------- #
def _esearch(term: str, email: str, retmax: int) -> list[str]:
    params = {"db": "nuccore", "term": term, "retmax": str(retmax),
              "retmode": "json", "tool": "amr-pipeline"}
    if email:
        params["email"] = email
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    raw = _get(url).decode("utf-8", "ignore")
    return re.findall(r'"(\d+)"', raw.split('"idlist"', 1)[-1])


def _efetch_cds(ids: list[str], email: str) -> str:
    params = {"db": "nuccore", "id": ",".join(ids), "rettype": "fasta_cds_na",
              "retmode": "text", "tool": "amr-pipeline"}
    if email:
        params["email"] = email
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params)
    return _get(url).decode("utf-8", "ignore")


def _iter_fasta(text: str):
    header, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def build_conjugation(email: str, retmax: int = 40) -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "conjugation.fasta"
    # mine well-characterised conjugative plasmids for their transfer CDS
    term = ('(conjugative plasmid[Title] OR conjugal transfer[Title]) '
            'AND bacteria[Organism] AND plasmid[Title] '
            'AND 20000:250000[SLEN] AND refseq[Filter]')
    print("  querying NCBI for conjugative plasmid records ...")
    ids = _esearch(term, email, retmax)
    if not ids:
        print("  ! no records found; conjugation.fasta will be empty")
        out.write_text("", encoding="utf-8")
        return 0
    written = 0
    seen: set[str] = set()
    with open(out, "w", encoding="utf-8") as fh:
        # fetch in small batches to be gentle on E-utilities
        for i in range(0, len(ids), 5):
            batch = ids[i:i + 5]
            try:
                cds = _efetch_cds(batch, email)
            except Exception as exc:  # noqa: BLE001
                print(f"  ! efetch batch failed ({exc})")
                continue
            for header, seq in _iter_fasta(cds):
                if not seq or len(seq) < 300:
                    continue
                gene = ""
                gm = re.search(r"\[gene=([^\]]+)\]", header)
                pm = re.search(r"\[protein=([^\]]+)\]", header)
                label = (pm.group(1) if pm else "") + " " + (gm.group(1) if gm else "")
                if not CONJ_KEYWORDS.search(label):
                    continue
                gene = (gm.group(1) if gm else
                        re.sub(r"[^A-Za-z0-9]+", "_", (pm.group(1) if pm else "transfer")))[:30]
                key = seq[:60]
                if key in seen:
                    continue
                seen.add(key)
                drug = "conjugation_marker"
                fh.write(f">{gene}|{drug}|conj|{label.strip()}\n{seq}\n")
                written += 1
            time.sleep(0.4)
    print(f"wrote {out}  ({written} transfer/relaxase CDS)")
    if written == 0:
        print("  (note: empty is OK — plasmid ARGs will score Medium not High)")
    return written


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", default="", help="contact email for NCBI")
    ap.add_argument("--skip-conjugation", action="store_true")
    ap.add_argument("--conj-retmax", type=int, default=40)
    args = ap.parse_args(argv)

    print("[1/3] PlasmidFinder replicons")
    build_replicons()
    print("[2/3] ISfinder insertion sequences")
    build_insertion_sequences()
    if args.skip_conjugation:
        print("[3/3] conjugation: skipped")
    else:
        print("[3/3] conjugation markers (coarse nucleotide proxy)")
        if not args.email:
            print("  ! NCBI strongly recommends --email; continuing anyway")
        build_conjugation(args.email, args.conj_retmax)
    print("\nDone. Catalogues are in", OUT_DIR)


if __name__ == "__main__":
    sys.exit(main())
