"""Phase 3 — AMR gene screening.

Detects acquired antibiotic-resistance genes (ARGs) in each genome and returns
one normalised table regardless of which backend produced it.

Backends
--------
* ``builtin``   : pure-Python, stdlib-only, NO external binaries and NO admin.
                  Uses a k-mer index of each genome contig + ungapped extension,
                  so it tolerates point mutations and is fast enough to screen a
                  handful of real genomes against the full CARD catalogue.
                  This is the recommended backend on a locked-down Windows box.
* ``blast``     : runs ``blastn`` against a local nucleotide ARG catalogue.
                  Faster/gapped, but needs ncbi-blast+ installed.
* ``abricate``  : real production tool bundling CARD/ResFinder/AMRFinder/
                  PlasmidFinder. Needs conda/bioconda.

Every backend emits the SAME normalised schema, so Phases 4-6 never branch on
the backend:

    strain_id contig start end strand gene pct_coverage pct_identity
    database accession product drug_class

Catalogue header parsing is auto-detected and supports:
    CARD       >gb|ACC|+|s-e|ARO:xxxx|geneName [organism]   (or  Name:gene)
    ResFinder  >tet(M)_1_X92947
    pipeline   >GENE|DRUG_CLASS|ACCESSION|PRODUCT            (our own / synthetic)
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from .utils import get_logger, read_fasta, revcomp

log = get_logger("amr")

NORM_COLS = [
    "strain_id", "contig", "start", "end", "strand", "gene",
    "pct_coverage", "pct_identity", "database", "accession",
    "product", "drug_class",
]

# Gene-name prefix -> antibiotic class (used to fill drug_class when a backend
# does not provide it). Order matters: more specific prefixes first.
_DRUG_PREFIX: List[Tuple[str, str]] = [
    ("cfr", "phenicol-oxazolidinone"),
    ("optr", "oxazolidinone"),
    ("poxt", "oxazolidinone"),
    ("vanc", "glycopeptide"),
    ("van", "glycopeptide"),
    ("tet", "tetracycline"),
    ("erm", "macrolide"),
    ("mef", "macrolide"),
    ("mph", "macrolide"),
    ("msr", "macrolide"),
    ("ere", "macrolide"),
    ("aac", "aminoglycoside"),
    ("aph", "aminoglycoside"),
    ("ant", "aminoglycoside"),
    ("aad", "aminoglycoside"),
    ("str", "aminoglycoside"),
    ("arm", "aminoglycoside"),
    ("rmt", "aminoglycoside"),
    ("bla", "beta-lactam"),
    ("mec", "beta-lactam"),
    ("cmy", "beta-lactam"),
    ("ctx", "beta-lactam"),
    ("shv", "beta-lactam"),
    ("tem", "beta-lactam"),
    ("oxa", "beta-lactam"),
    ("kpc", "beta-lactam"),
    ("ndm", "beta-lactam"),
    ("imp", "beta-lactam"),
    ("vim", "beta-lactam"),
    ("cat", "phenicol"),
    ("cml", "phenicol"),
    ("flo", "phenicol"),
    ("fex", "phenicol"),
    ("lnu", "lincosamide"),
    ("lsa", "lincosamide"),
    ("vga", "lincosamide"),
    ("dfr", "trimethoprim"),
    ("sul", "sulfonamide"),
    ("qnr", "quinolone"),
    ("oqx", "quinolone"),
    ("mcr", "colistin"),
    ("fos", "fosfomycin"),
]


def drug_class_for(gene: str) -> str:
    g = gene.lower().lstrip("(")
    for prefix, klass in _DRUG_PREFIX:
        if g.startswith(prefix):
            return klass
    return "unclassified"


# --------------------------------------------------------------------------- #
# Reference catalogue (auto-detecting header parser)                          #
# --------------------------------------------------------------------------- #
def parse_header(header: str) -> Tuple[str, str, str, str]:
    """Return (gene, drug_class, accession, product) from a FASTA header.

    Auto-detects CARD, ResFinder and the pipeline's own pipe format.
    """
    h = header.lstrip(">").strip()

    # --- CARD ----------------------------------------------------------- #
    if "ARO:" in h or h.lower().startswith("name:") or "|name:" in h.lower():
        fields = [f.strip() for f in h.split("|")]
        gene = ""
        accession = ""
        product = h
        for f in fields:
            low = f.lower()
            if low.startswith("name:"):
                gene = f.split(":", 1)[1].strip()
            elif low.startswith("aro:"):
                pass  # marker only
        if not gene:
            # gb|ACC|+|range|ARO:xxx|GeneName [organism]
            for i, f in enumerate(fields):
                if f.startswith("ARO:") and i + 1 < len(fields):
                    product = fields[i + 1]
                    gene = product.split("[")[0].split()[0]
                    break
        if not gene:
            gene = fields[-1].split("[")[0].split()[0] if fields else h
        for f in fields:
            if f and f[0].isalnum() and any(c.isdigit() for c in f) and "." in f:
                accession = f
                break
        if not accession:
            accession = next((f for f in fields if f.startswith("ARO:")), gene)
        gene = gene.strip()
        return gene, drug_class_for(gene), accession, product.strip()

    # --- pipeline / synthetic:  GENE|CLASS|ACC|PRODUCT ------------------ #
    if "|" in h:
        parts = h.split("|")
        gene = parts[0].split()[0]
        drug = parts[1].strip() if len(parts) > 1 and parts[1].strip() else drug_class_for(gene)
        acc = parts[2].strip() if len(parts) > 2 else gene
        product = parts[3].strip() if len(parts) > 3 else f"{gene} resistance"
        return gene, drug, acc, product

    # --- ResFinder:  gene(_variant)?_<num>_<accession> ------------------ #
    token = h.split()[0]
    if "_" in token:
        gene = token.split("_")[0]
        accession = token.split("_")[-1]
        return gene, drug_class_for(gene), accession, gene

    return token, drug_class_for(token), token, token


def load_catalogue(fasta_path: str | Path) -> List[dict]:
    """Load an ARG reference FASTA into a list of catalogue entries."""
    seqs = read_fasta(fasta_path)
    headers: Dict[str, str] = {}
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                full = line[1:].strip()
                rid = full.split()[0]
                headers[rid] = full
    cat: List[dict] = []
    for rid, seq in seqs.items():
        gene, drug, acc, product = parse_header(headers.get(rid, rid))
        cat.append({"gene": gene, "drug_class": drug, "accession": acc,
                    "product": product, "seq": seq.upper(), "ref_id": rid})
    return cat


# --------------------------------------------------------------------------- #
# Backend: built-in k-mer-indexed aligner  (stdlib only, no admin)            #
# --------------------------------------------------------------------------- #
_BASE2BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


def _encode_positions(seq: str, k: int, step: int = 1):
    """Yield (pos, code) for every k-mer of *seq*.

    Uses a rolling 2-bit integer code. Any non-ACGT base resets the window so
    ambiguous bases never create false k-mers. With step>1 only positions whose
    index is divisible by step are yielded (used to shrink the genome index).
    """
    code = 0
    mask = (1 << (2 * k)) - 1
    valid = 0
    for i, b in enumerate(seq):
        bits = _BASE2BITS.get(b, -1)
        if bits < 0:
            valid = 0
            code = 0
            continue
        code = ((code << 2) | bits) & mask
        valid += 1
        if valid >= k:
            pos = i - k + 1
            if step == 1 or pos % step == 0:
                yield pos, code


def _build_index(contig: str, k: int, index_step: int) -> Dict[int, List[int]]:
    idx: Dict[int, List[int]] = defaultdict(list)
    for pos, code in _encode_positions(contig, k, index_step):
        idx[code].append(pos)
    return idx


def _hamming_identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    m = sum(1 for x, y in zip(a, b) if x == y)
    return 100.0 * m / n


def _extend_score(gene_seq: str, contig: str, base: int,
                  min_cov: float, min_id: float):
    """Ungapped extension of a full-length gene placement at contig pos *base*."""
    L = len(gene_seq)
    if base < 0:
        return None
    end = min(base + L, len(contig))
    window = contig[base:end]
    if not window:
        return None
    ident = _hamming_identity(gene_seq, window)
    cov = 100.0 * len(window) / L
    if ident >= min_id and cov >= min_cov:
        return {"start": base + 1, "end": base + len(window),
                "pct_identity": round(ident, 2), "pct_coverage": round(cov, 2)}
    return None


def _search_indexed(gene_seq: str, contig: str, index: Dict[int, List[int]],
                    k: int, min_cov: float, min_id: float,
                    min_seed_hits: int = 3) -> List[dict]:
    """Find near-matches of *gene_seq* in *contig* on both strands via the index."""
    hits: List[dict] = []
    L = len(gene_seq)
    if L < k:
        return hits
    for strand, gseq in (("+", gene_seq), ("-", revcomp(gene_seq))):
        votes: Counter = Counter()
        for offset, code in _encode_positions(gseq, k, 1):
            for gpos in index.get(code, ()):
                votes[gpos - offset] += 1          # anchor = genome start of gene
        if not votes:
            continue
        # evaluate the most-supported anchors first
        candidates = [a for a, v in votes.items() if v >= min_seed_hits]
        if not candidates:                          # short genes: relax threshold
            candidates = [a for a, v in votes.most_common(3)]
        scored: List[dict] = []
        for base in sorted(set(candidates)):
            res = _extend_score(gseq, contig, base, min_cov, min_id)
            if res:
                res["strand"] = strand
                scored.append(res)
        # keep best non-overlapping placements
        scored.sort(key=lambda h: (-h["pct_identity"], -h["pct_coverage"]))
        kept: List[dict] = []
        for h in scored:
            if not any(abs(h["start"] - x["start"]) < L * 0.5 for x in kept):
                kept.append(h)
        hits.extend(kept)
    return hits


# the small, index-free matcher is still used by mge_analysis for tiny POC
# catalogues and by tests; identical results, simpler code path.
def _search_gene(gene_seq: str, contig: str, seed_len: int = 15,
                 min_cov: float = 60.0, min_id: float = 80.0) -> List[dict]:
    hits: List[dict] = []
    L = len(gene_seq)
    if L < seed_len:
        return hits
    for strand, gseq in (("+", gene_seq), ("-", revcomp(gene_seq))):
        seed_positions = sorted({0, L // 4, L // 2, (3 * L) // 4})
        seen: set[int] = set()
        for sp in seed_positions:
            seed = gseq[sp:sp + seed_len]
            frm = 0
            while True:
                idx = contig.find(seed, frm)
                if idx == -1:
                    break
                frm = idx + 1
                base = idx - sp
                if base < 0 or base in seen:
                    continue
                seen.add(base)
                res = _extend_score(gseq, contig, base, min_cov, min_id)
                if res:
                    res["strand"] = strand
                    hits.append(res)
    hits.sort(key=lambda h: (-h["pct_identity"], -h["pct_coverage"]))
    kept: List[dict] = []
    for h in hits:
        if not any(abs(h["start"] - x["start"]) < L * 0.5 for x in kept):
            kept.append(h)
    return kept


def _reciprocal_overlap(s1: int, e1: int, s2: int, e2: int, t: float) -> bool:
    """True if [s1,e1] and [s2,e2] reciprocally overlap by >= t (both directions)."""
    ov = min(e1, e2) - max(s1, s2) + 1
    if ov <= 0:
        return False
    len1, len2 = e1 - s1 + 1, e2 - s2 + 1
    if len1 <= 0 or len2 <= 0:
        return False
    return (ov / len1) >= t and (ov / len2) >= t


def collapse_overlaps(df: pd.DataFrame, threshold: float = 0.8) -> pd.DataFrame:
    """Collapse overlapping hits to one best hit per strain+contig locus.

    Two hits on the SAME strain and SAME contig that reciprocally overlap by
    >= `threshold` are treated as the same locus (e.g. the many SHV/LEN/OKP
    beta-lactamase alleles that all map to one chromosomal locus). Within each
    overlapping cluster only the single best hit is kept, ranked by:
        1. highest % identity
        2. then highest % coverage
        3. then longest matched length (end - start + 1)

    All columns (contig, start, end, gene, pct_identity, pct_coverage, ...) are
    preserved; only rows are removed, so any column added later (e.g. risk_level
    during classification) is unaffected.
    """
    if df.empty:
        return df.reset_index(drop=True)

    work = df.copy()
    work["_matched_len"] = (work["end"] - work["start"] + 1).abs()
    # best-first ordering so the greedy keep picks the winner of each cluster
    work = work.sort_values(
        ["pct_identity", "pct_coverage", "_matched_len"],
        ascending=[False, False, False])

    keep_idx: List[int] = []
    for (_, _), grp in work.groupby(["strain_id", "contig"], sort=False):
        kept_rows: List[Tuple[int, int]] = []          # (start, end) of kept hits
        for idx, row in grp.iterrows():
            s, e = int(row.start), int(row.end)
            if any(_reciprocal_overlap(s, e, ks, ke, threshold) for ks, ke in kept_rows):
                continue                               # subsumed by a better hit
            kept_rows.append((s, e))
            keep_idx.append(idx)

    out = (work.loc[keep_idx]
               .drop(columns="_matched_len")
               .sort_values(["strain_id", "contig", "start"])
               .reset_index(drop=True))
    return out


def screen_builtin(strain_id: str, genome_path: str | Path,
                   catalogue: List[dict], min_cov: float = 60.0,
                   min_id: float = 80.0, k: int = 13,
                   index_step: int = 3, collapse: bool = True,
                   overlap_threshold: float = 0.8) -> pd.DataFrame:
    """Index each contig once, then screen every catalogue gene against it."""
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
                    "start": h["start"], "end": h["end"], "strand": h["strand"],
                    "gene": ref["gene"], "pct_coverage": h["pct_coverage"],
                    "pct_identity": h["pct_identity"], "database": "builtin",
                    "accession": ref["accession"], "product": ref["product"],
                    "drug_class": ref["drug_class"],
                })
        del index  # free memory before next contig
    df = pd.DataFrame(rows, columns=NORM_COLS)
    if df.empty:
        return df
    # drop exact same-gene/same-start duplicates first (keep best identity)
    df = (df.sort_values("pct_identity", ascending=False)
            .drop_duplicates(subset=["strain_id", "contig", "gene", "start"])
            .reset_index(drop=True))
    # then collapse overlapping multi-allele hits down to one best per locus
    if collapse:
        before = len(df)
        df = collapse_overlaps(df, overlap_threshold)
        log.info("%-22s collapsed %d -> %d hits (overlap>=%.2f)",
                 strain_id, before, len(df), overlap_threshold)
    return df


# --------------------------------------------------------------------------- #
# Backend: BLAST                                                               #
# --------------------------------------------------------------------------- #
def screen_blast(strain_id: str, genome_path: str | Path,
                 catalogue_fasta: str | Path, min_cov: float = 60.0,
                 min_id: float = 80.0) -> pd.DataFrame:
    if shutil.which("blastn") is None:
        raise RuntimeError("blastn not on PATH (install ncbi-blast+).")
    cat = {c["ref_id"]: c for c in load_catalogue(catalogue_fasta)}
    with tempfile.TemporaryDirectory() as tmp:
        db = str(Path(tmp) / "genome")
        subprocess.run(["makeblastdb", "-in", str(genome_path), "-dbtype", "nucl",
                        "-out", db], check=True, capture_output=True)
        fmt = "6 qseqid sseqid sstart send sstrand pident length qlen"
        proc = subprocess.run(
            ["blastn", "-query", str(catalogue_fasta), "-db", db,
             "-outfmt", fmt, "-perc_identity", str(min_id)],
            check=True, capture_output=True, text=True)
    rows: List[dict] = []
    for line in proc.stdout.splitlines():
        q, s, sstart, send, sstrand, pident, length, qlen = line.split("\t")
        cov = 100.0 * int(length) / int(qlen)
        if cov < min_cov:
            continue
        ref = cat.get(q, {})
        start, end = sorted((int(sstart), int(send)))
        rows.append({
            "strain_id": strain_id, "contig": s, "start": start, "end": end,
            "strand": "+" if sstrand == "plus" else "-",
            "gene": ref.get("gene", q.split("|")[0]),
            "pct_coverage": round(cov, 2), "pct_identity": float(pident),
            "database": "blast", "accession": ref.get("accession", q),
            "product": ref.get("product", ""),
            "drug_class": ref.get("drug_class", drug_class_for(q)),
        })
    return pd.DataFrame(rows, columns=NORM_COLS)


# --------------------------------------------------------------------------- #
# Backend: abricate                                                            #
# --------------------------------------------------------------------------- #
def screen_abricate(strain_id: str, genome_path: str | Path,
                    db: str = "card", min_cov: float = 60.0,
                    min_id: float = 80.0) -> pd.DataFrame:
    if shutil.which("abricate") is None:
        raise RuntimeError("abricate not on PATH.")
    proc = subprocess.run(
        ["abricate", "--db", db, "--minid", str(min_id),
         "--mincov", str(min_cov), str(genome_path)],
        check=True, capture_output=True, text=True)
    return parse_abricate(proc.stdout, strain_id, db)


def parse_abricate(text: str, strain_id: str, db: str) -> pd.DataFrame:
    """Parse abricate TSV into the normalised schema."""
    rows: List[dict] = []
    header: List[str] | None = None
    for line in text.splitlines():
        cols = line.split("\t")
        if line.startswith("#FILE"):
            header = [c.lstrip("#").strip() for c in cols]
            continue
        if header is None or len(cols) < len(header):
            continue
        rec = dict(zip(header, cols))
        gene = rec.get("GENE", "")
        rows.append({
            "strain_id": strain_id, "contig": rec.get("SEQUENCE", ""),
            "start": int(float(rec.get("START", 0))),
            "end": int(float(rec.get("END", 0))),
            "strand": rec.get("STRAND", "+"),
            "gene": gene,
            "pct_coverage": float(rec.get("%COVERAGE", 0) or 0),
            "pct_identity": float(rec.get("%IDENTITY", 0) or 0),
            "database": rec.get("DATABASE", db),
            "accession": rec.get("ACCESSION", ""),
            "product": rec.get("PRODUCT", ""),
            "drug_class": rec.get("RESISTANCE", "") or drug_class_for(gene),
        })
    return pd.DataFrame(rows, columns=NORM_COLS)


# --------------------------------------------------------------------------- #
# Dispatcher                                                                   #
# --------------------------------------------------------------------------- #
def screen_genome(strain_id: str, genome_path: str | Path, backend: str,
                  catalogue: List[dict] | None = None,
                  catalogue_fasta: str | Path | None = None,
                  abricate_db: str = "card",
                  min_cov: float = 60.0, min_id: float = 80.0,
                  collapse_overlaps: bool = True,
                  overlap_threshold: float = 0.8) -> pd.DataFrame:
    if backend == "abricate":
        df = screen_abricate(strain_id, genome_path, abricate_db, min_cov, min_id)
    elif backend == "blast":
        if catalogue_fasta is None:
            raise ValueError("blast backend needs catalogue_fasta")
        df = screen_blast(strain_id, genome_path, catalogue_fasta, min_cov, min_id)
    elif backend == "builtin":
        if catalogue is None:
            if catalogue_fasta is None:
                raise ValueError("builtin backend needs catalogue or catalogue_fasta")
            catalogue = load_catalogue(catalogue_fasta)
        # builtin collapses internally so the locus stats log per-strain
        return screen_builtin(strain_id, genome_path, catalogue, min_cov, min_id,
                              collapse=collapse_overlaps,
                              overlap_threshold=overlap_threshold)
    else:
        raise ValueError(f"unknown backend: {backend}")
    # blast / abricate: collapse here so behaviour is identical across backends
    if collapse_overlaps and not df.empty:
        df = collapse_overlaps_df(df, overlap_threshold)
    return df


# alias kept so the dispatcher's local arg name doesn't shadow the function
collapse_overlaps_df = collapse_overlaps


def screen_all(fasta_map: Dict[str, Path], backend: str, **kwargs) -> pd.DataFrame:
    # pre-load the catalogue once for the builtin backend (big speed win)
    if backend == "builtin" and kwargs.get("catalogue") is None and kwargs.get("catalogue_fasta"):
        kwargs["catalogue"] = load_catalogue(kwargs["catalogue_fasta"])
        log.info("loaded %d reference genes", len(kwargs["catalogue"]))
    frames = []
    for sid, fp in fasta_map.items():
        df = screen_genome(sid, fp, backend, **kwargs)
        log.info("%-22s %d ARG hit(s)", sid, len(df))
        frames.append(df)
    frames = [f for f in frames if not f.empty]   # avoid pandas concat FutureWarning
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=NORM_COLS)
