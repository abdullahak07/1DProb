"""Phase 2 — Strain selection & genome collection.

Downloads genome assemblies for the configured probiotic and pathogen
accessions from NCBI and builds the Table-1 strain-metadata sheet.

Network access (NCBI) is only used by `fetch_accession`. Everything else
(`build_metadata`) works offline on already-downloaded FASTA files, which is
why the smoke test can exercise the metadata logic without touching NCBI.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .utils import (ensure_dir, gc_content, genome_length, get_logger,
                    n_contigs, read_fasta, write_fasta)

log = get_logger("fetch")


def fetch_accession(
    accession: str,
    out_dir: str | Path,
    email: str,
    api_key: Optional[str] = None,
    db: str = "nuccore",
    retries: int = 3,
) -> Path:
    """Download one genome (by nuccore/assembly accession) to *out_dir* as FASTA.

    Requires internet + Biopython. Used on the researcher's machine, not in the
    offline smoke test.
    """
    from Bio import Entrez, SeqIO  # imported lazily so offline runs don't need it

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    out_dir = ensure_dir(out_dir)
    out_path = Path(out_dir) / f"{accession}.fna"
    if out_path.exists() and out_path.stat().st_size > 0:
        log.info("cached %s", accession)
        return out_path

    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as h:
                records = list(SeqIO.parse(h, "fasta"))
            if not records:
                raise ValueError(f"no sequence returned for {accession}")
            seqs = {r.id: str(r.seq) for r in records}
            write_fasta(out_path, seqs)
            log.info("downloaded %s (%d contigs)", accession, len(seqs))
            return out_path
        except Exception as exc:  # noqa: BLE001 - we retry on any transient failure
            last_err = exc
            log.warning("fetch %s failed (attempt %d/%d): %s",
                        accession, attempt, retries, exc)
            time.sleep(2 * attempt)
    raise RuntimeError(f"could not download {accession}: {last_err}")


def fetch_all(strains: List[dict], out_dir: str | Path, email: str,
              api_key: Optional[str] = None) -> Dict[str, Path]:
    """Download every strain in the config list. Returns {strain_id: fasta_path}."""
    paths: Dict[str, Path] = {}
    for s in strains:
        sub = ensure_dir(Path(out_dir) / s["group"] / s["species"].replace(" ", "_"))
        try:
            paths[s["strain_id"]] = fetch_accession(s["accession"], sub, email, api_key)
        except Exception as exc:  # noqa: BLE001
            log.error("skipping %s: %s", s["strain_id"], exc)
    return paths


# --------------------------------------------------------------------------- #
# Full-assembly fetching (chromosome + ALL plasmid replicons)                 #
# --------------------------------------------------------------------------- #
def _extract_link_ids(elink_result, linkname: Optional[str] = None) -> List[str]:
    """Pull linked UIDs out of a Biopython Entrez.read(elink) result.

    Pure function (no network) so the link-walking logic is unit-testable.
    """
    ids: List[str] = []
    if not elink_result:
        return ids
    for ls in elink_result[0].get("LinkSetDb", []):
        if linkname and ls.get("LinkName") != linkname:
            continue
        for lnk in ls.get("Link", []):
            ids.append(lnk["Id"])
    return ids


def fetch_assembly(accession: str, out_dir: str | Path, email: str,
                   api_key: Optional[str] = None, retries: int = 3) -> Path:
    """Download the FULL assembly that *accession* belongs to.

    Walks  nuccore(accession) -> assembly -> all RefSeq replicons  and writes a
    single combined FASTA (chromosome + every plasmid). Falls back to fetching
    the single accession if no assembly link exists.
    """
    from Bio import Entrez, SeqIO  # lazy import; offline runs don't need it

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    out_dir = ensure_dir(out_dir)
    out_path = Path(out_dir) / f"{accession}_full.fna"
    if out_path.exists() and out_path.stat().st_size > 0:
        log.info("cached assembly %s", accession)
        return out_path

    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with Entrez.esearch(db="nuccore", term=accession) as h:
                rec = Entrez.read(h)
            if not rec["IdList"]:
                raise ValueError(f"no nuccore record for {accession}")
            nuc_uid = rec["IdList"][0]

            with Entrez.elink(dbfrom="nuccore", db="assembly", id=nuc_uid) as h:
                asm_link = Entrez.read(h)
            asm_ids = _extract_link_ids(asm_link)
            if not asm_ids:
                log.warning("no assembly link for %s; fetching single accession", accession)
                return fetch_accession(accession, out_dir, email, api_key)
            asm_uid = asm_ids[0]

            nuc_ids: List[str] = []
            for ln in ("assembly_nuccore_refseq", "assembly_nuccore_insdc",
                       "assembly_nuccore"):
                with Entrez.elink(dbfrom="assembly", db="nuccore",
                                  id=asm_uid, linkname=ln) as h:
                    rep_link = Entrez.read(h)
                nuc_ids = _extract_link_ids(rep_link, ln)
                if nuc_ids:
                    break
            if not nuc_ids:
                nuc_ids = [nuc_uid]

            with Entrez.efetch(db="nuccore", id=",".join(nuc_ids),
                               rettype="fasta", retmode="text") as h:
                records = list(SeqIO.parse(h, "fasta"))
            seqs = {r.id: str(r.seq) for r in records if str(r.seq)}
            if not seqs:
                raise ValueError("no sequences fetched")
            write_fasta(out_path, seqs)
            log.info("assembly %s -> %d replicon(s): %s", accession, len(seqs),
                     ", ".join(list(seqs)[:6]))
            return out_path
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            log.warning("assembly fetch %s failed (attempt %d/%d): %s",
                        accession, attempt, retries, exc)
            time.sleep(2 * attempt)
    raise RuntimeError(f"could not fetch assembly for {accession}: {last_err}")


def fetch_all_assemblies(strains: List[dict], out_dir: str | Path, email: str,
                         api_key: Optional[str] = None) -> Dict[str, Path]:
    """Full-assembly version of fetch_all (chromosome + plasmids per strain)."""
    paths: Dict[str, Path] = {}
    for s in strains:
        sub = ensure_dir(Path(out_dir) / s["group"] / s["species"].replace(" ", "_"))
        try:
            paths[s["strain_id"]] = fetch_assembly(s["accession"], sub, email, api_key)
        except Exception as exc:  # noqa: BLE001
            log.error("skipping %s: %s", s["strain_id"], exc)
    return paths


def build_metadata(strains: List[dict], fasta_map: Dict[str, Path]) -> pd.DataFrame:
    """Phase 2 Table 1: one row per strain with genome-quality metrics.

    Works purely on local FASTA files, so it is fully testable offline.
    """
    rows = []
    for s in strains:
        sid = s["strain_id"]
        fp = fasta_map.get(sid)
        if fp is None or not Path(fp).exists():
            log.warning("no FASTA for %s; recording as missing", sid)
            rows.append({**_meta_skeleton(s), "genome_available": False})
            continue
        seqs = read_fasta(fp)
        rows.append({
            **_meta_skeleton(s),
            "genome_available": True,
            "genome_length_bp": genome_length(seqs),
            "n_contigs": n_contigs(seqs),
            "gc_percent": round(
                sum(gc_content(x) * len(x) for x in seqs.values())
                / max(genome_length(seqs), 1), 2),
            "fasta_path": str(fp),
        })
    df = pd.DataFrame(rows)
    # quality flag: complete genomes (few contigs) are higher confidence
    if "n_contigs" in df:
        df["assembly_quality"] = df["n_contigs"].apply(
            lambda n: "complete" if pd.notna(n) and n <= 3
            else ("draft" if pd.notna(n) else "missing"))
    return df


def _meta_skeleton(s: dict) -> dict:
    return {
        "strain_id": s["strain_id"],
        "species": s["species"],
        "group": s["group"],                 # 'probiotic' | 'pathogen'
        "accession": s["accession"],
        "product_example": s.get("product_example", ""),
        "isolation_source": s.get("isolation_source", ""),
        "qps_listed": s.get("qps_listed", ""),
    }
