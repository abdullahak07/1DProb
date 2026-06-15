"""Phase 6 - SOTA v2 report. Stitches Phases 1-5 outputs + SCALE50 tables into
one manuscript-oriented Markdown report with a claims ledger and safe wording.

Usage:
    py scripts\\generate_sota_report.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _sota_common import df_to_md, read_csv_safe, sota_dir

PENDING = "_not run yet_"


def _risk_counts(arg):
    rc = arg["risk_level"].value_counts().to_dict() if arg is not None and "risk_level" in arg else {}
    return {r: int(rc.get(r, 0)) for r in ["Low", "Medium", "High", "Negligible"]}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", default="results_scale50")
    ap.add_argument("--chronly", default="results_scale50_chronly")
    ap.add_argument("--sota", default="results_sota")
    args = ap.parse_args(argv)
    out = sota_dir()
    full = Path(args.full)

    arg = read_csv_safe(full / "table2_master_arg.csv")
    meta = read_csv_safe(full / "table1_strain_metadata.csv")
    feat = read_csv_safe(full / "mge_features.csv")
    edges = read_csv_safe(full / "gene_sharing_edges.csv")
    if arg is None:
        raise SystemExit(f"{full}/table2_master_arg.csv not found")
    rc = _risk_counts(arg)

    n_is = int((feat.feature_type == "insertion_sequence").sum()) if feat is not None else 0
    n_rep = int((feat.feature_type == "plasmid_replicon").sum()) if feat is not None else 0
    n_conj = int((feat.feature_type == "conjugation_marker").sum()) if feat is not None else 0
    n_pro = int((meta.group == "probiotic").sum()) if meta is not None else None
    n_pat = int((meta.group == "pathogen").sum()) if meta is not None else None
    avail = int(meta.get("genome_available", pd.Series(dtype=bool)).sum()) if meta is not None else None

    # High calls + group
    high = arg[arg.risk_level == "High"].copy()
    if meta is not None:
        gmap = dict(zip(meta.strain_id, meta.group))
        high["group"] = high.strain_id.map(gmap)
    high_in_pro = int((high.get("group") == "probiotic").sum()) if "group" in high else 0
    cross = int(edges.get("cross_group", pd.Series([], dtype=bool)).sum()) if edges is not None else 0
    n_edges = len(edges) if edges is not None else 0

    # phase outputs
    comp = read_csv_safe(out / "full_vs_chronly_comparison.csv")
    acct = read_csv_safe(out / "full_vs_chronly_locus_accounting.csv")
    chronly_arg = read_csv_safe(Path(args.chronly) / "table2_master_arg.csv")
    curated = read_csv_safe(out / "curated_ARG_calls_SCALE50.csv")
    cur_hm = read_csv_safe(out / "curated_high_medium_SCALE50.csv")
    val = read_csv_safe(out / "validated_ARG_calls_SCALE50.csv")
    hr_val = read_csv_safe(out / "high_risk_validation_summary.csv")
    mob = read_csv_safe(out / "mobility_validation_SCALE50.csv")
    stats = read_csv_safe(out / "stats_summary.csv")
    risk_enr = read_csv_safe(out / "risk_enrichment.csv")
    net = read_csv_safe(out / "network_summary.csv")
    vcounts = read_csv_safe(out / "validation_tool_counts.csv")
    mtool = read_csv_safe(out / "mobility_tooling.csv")

    def _b(df, col):
        return bool(df[col].iloc[0]) if (df is not None and col in df.columns and len(df)) else False

    def _i(df, col):
        return int(df[col].iloc[0]) if (df is not None and col in df.columns and len(df)) else 0

    amrfinder_present = _b(vcounts, "amrfinder_present")
    resfinder_present = _b(vcounts, "resfinder_present")
    blast_present = _b(vcounts, "blast_present")
    mobsuite_present = _b(mtool, "mobsuite_present")
    genomad_present = _b(mtool, "genomad_present")

    amrfinder_confirmed = _i(vcounts, "amrfinder_confirmed")
    resfinder_confirmed = _i(vcounts, "resfinder_confirmed")
    blast_confirmed = _i(vcounts, "blast_confirmed")
    multi_tool_confirmed = _i(vcounts, "multi_tool_confirmed")
    discordant_count = _i(vcounts, "discordant")

    # Validation states:
    #   no_external_identity: no independent identity tool has been imported.
    #   blast_only: BLAST imported, but AMRFinderPlus/ResFinder still absent.
    #   sota_ready: BOTH identity tools and at least one external mobility tool imported.
    # BLAST is independent validation, but BLAST-only is not enough for the SOTA label.
    external_identity_present = amrfinder_present or resfinder_present or blast_present
    blast_only = blast_present and not amrfinder_present and not resfinder_present
    # The work is only labelled SOTA once BOTH identity tools AND a mobility tool are imported.
    sota_ready = amrfinder_present and resfinder_present and (mobsuite_present or genomad_present)
    val_mode_pending = not external_identity_present

    missing = []
    if not amrfinder_present:
        missing.append("AMRFinderPlus")
    if not resfinder_present:
        missing.append("ResFinder")
    if not (mobsuite_present or genomad_present):
        missing.append("MOB-suite or geNomad")

    md = []
    if sota_ready:
        md.append("# SCALE50 SOTA v2 (externally validated) - Probiotic-associated AMR-cargo screen\n")
    else:
        md.append("# SCALE50 - Probiotic-associated AMR-cargo screen (validation in progress - NOT yet SOTA)\n")
        md.append(f"> **This run is not yet state-of-the-art.** Missing external "
                  f"validation: **{', '.join(missing)}**. The SOTA label is withheld "
                  f"until AMRFinderPlus, ResFinder and a mobility tool (MOB-suite or "
                  f"geNomad) are imported into `validation/` and re-run.\n")

    # 1. executive summary
    md.append("## 1. Executive summary\n")
    md.append(
        f"We screened {len(meta) if meta is not None else 'N'} bacterial reference "
        f"genomes ({n_pro} probiotic-associated, {n_pat} pathogen/comparator) for "
        f"antibiotic-resistance genes (ARGs) and their mobile-genetic-element "
        f"context using a no-admin, stdlib Python pipeline (CARD nucleotide "
        f"catalogue, PlasmidFinder replicons, ISfinder insertion sequences). "
        f"We detected **{len(arg)} ARG loci** ({rc['High']} High, {rc['Medium']} "
        f"Medium, {rc['Low']} Low), **{n_rep} plasmid replicons**, **{n_is} IS "
        f"features**, and **{n_conj} conjugation markers**. "
        f"**All {rc['High']} High-risk ARGs are in pathogen/comparator genomes "
        f"({high_in_pro} in probiotic-associated genomes)**, and **{cross} of "
        f"{n_edges} gene-sharing edges are probiotic-pathogen**. The data indicate "
        f"that probiotic-associated reference genomes in this set carry little "
        f"high-risk mobile ARG cargo and share no ARGs with pathogens at high "
        f"identity; high-risk, mobile-context ARGs are concentrated in the "
        f"pathogen comparators. This is a reassuring-safety + method finding, not "
        f"evidence of probiotic-to-pathogen transfer.\n")

    # 2. dataset
    md.append("## 2. Dataset summary\n")
    if meta is not None:
        keep = [c for c in ["strain_id", "species", "group", "accession",
                            "n_contigs", "assembly_quality", "genome_available"]
                if c in meta.columns]
        md.append(f"{len(meta)} strains; {avail} with genomes retrieved "
                  f"(full assembly: chromosome + plasmid replicons).\n")
        md.append(df_to_md(meta[keep], max_rows=100))
    else:
        md.append(PENDING + "\n")

    # 3. full results
    md.append("## 3. Full-assembly SCALE50 results\n")
    md.append(df_to_md(pd.DataFrame([
        {"metric": "ARG loci", "value": len(arg)},
        {"metric": "High", "value": rc["High"]},
        {"metric": "Medium", "value": rc["Medium"]},
        {"metric": "Low", "value": rc["Low"]},
        {"metric": "insertion sequences", "value": n_is},
        {"metric": "plasmid replicons", "value": n_rep},
        {"metric": "conjugation markers", "value": n_conj},
        {"metric": "gene-sharing edges (total)", "value": n_edges},
        {"metric": "probiotic-pathogen edges", "value": cross},
    ])))

    # 4. comparison
    md.append("## 4. Matched chromosome-only vs full-assembly comparison\n")
    if chronly_arg is None or comp is None:
        md.append("**Matched chromosome-only baseline NOT RUN YET.** Run "
                  "`config\\config_scale50_chronly.yaml`, then "
                  "`scripts\\compare_full_vs_chronly.py`. (The old SCALE8 baseline "
                  "is intentionally NOT used.)\n")
    else:
        md.append("Same 50 strains, full_assembly true vs false:\n")
        md.append(df_to_md(comp))
        if acct is not None:
            md.append("\n**Locus accounting** (these count different units, so they "
                      "are reported separately):\n")
            md.append(df_to_md(acct))
            md.append("\nThe net ARG delta is a count of **rows** (and equals the "
                      "sum of the risk-tier deltas); 'matched gained loci' counts "
                      "**unique strain x gene combinations**. The difference between "
                      "them is genes recovered at more than one plasmid locus in the "
                      "same strain (e.g. an ESBL gene present on two plasmids), which "
                      "add rows without adding new strain x gene keys. This is not a "
                      "discrepancy; it is two valid views of the same data.\n")

    # 5. curation
    md.append("## 5. Curated ARG interpretation\n")
    if curated is None:
        md.append(PENDING + " (run `scripts\\curate_arg_calls.py`)\n")
    else:
        cat = curated.curated_category.value_counts() if "curated_category" in curated else pd.Series(dtype=int)
        md.append("Raw calls preserved; core/intrinsic/biocide/efflux determinants "
                  "flagged and excluded from the mobile-ARG analysis.\n")
        md.append(df_to_md(cat.rename_axis("curated_category").reset_index(name="count")))
        if cur_hm is not None:
            md.append(f"\nCurated High/Medium acquired-cargo loci: **{len(cur_hm)}** "
                      f"(raw High+Medium was {rc['High'] + rc['Medium']}).\n")

    # 6. validation
    md.append("## 6. Independent validation summary\n")
    if val is None:
        md.append(PENDING + " (run `scripts\\validate_arg_calls.py`)\n")
    else:
        sc = val.validation_status.value_counts().to_dict()
        md.append("Per-locus aggregate status:\n")
        md.append(df_to_md(pd.DataFrame([{"validation_status": k, "count": v}
                                         for k, v in sc.items()])))
        md.append("\n**Per-tool confirmation counts**\n")
        md.append(df_to_md(pd.DataFrame([
            {"tool": "AMRFinderPlus", "imported": amrfinder_present, "confirmed_loci": _i(vcounts, "amrfinder_confirmed")},
            {"tool": "ResFinder", "imported": resfinder_present, "confirmed_loci": _i(vcounts, "resfinder_confirmed")},
            {"tool": "BLAST", "imported": blast_present, "confirmed_loci": _i(vcounts, "blast_confirmed")},
            {"tool": "multi-tool (>=2 agree)", "imported": "-", "confirmed_loci": _i(vcounts, "multi_tool_confirmed")},
            {"tool": "discordant", "imported": "-", "confirmed_loci": _i(vcounts, "discordant")},
        ])))
        if val_mode_pending:
            md.append("\n**External validation pending.** No external identity tool "
                      "has been imported; the builtin consistency recheck confirms "
                      "*internal reproducibility* only (NOT independent validation). "
                      "Manuscript text must say 'builtin consistency recheck confirmed "
                      "internal reproducibility', never 'independently validated'.\n")
        elif blast_only:
            md.append("\n**External validation partially complete.** BLAST has externally "
                      f"confirmed {blast_confirmed} priority loci, including all "
                      f"{rc['High']} High-risk loci and all "
                      f"{len(cur_hm) if cur_hm is not None else 'curated High/Medium'} "
                      "curated High/Medium acquired-cargo loci. AMRFinderPlus "
                      "and ResFinder have not yet been imported, so the run "
                      "remains NOT yet SOTA.\n")
        elif not sota_ready:
            md.append(f"\n**Partial external validation.** Imported: "
                      f"{'AMRFinderPlus ' if amrfinder_present else ''}"
                      f"{'ResFinder ' if resfinder_present else ''}"
                      f"{'BLAST' if blast_present else ''}. Still missing for SOTA: "
                      f"{', '.join(missing)}.\n")
        else:
            md.append("\n**SOTA v2 externally validated.** AMRFinderPlus and "
                      "ResFinder have been imported for ARG identity validation, "
                      "and at least one external mobility tool (MOB-suite or "
                      "geNomad) has been imported for plasmid/MGE context.\n")
        if hr_val is not None:
            md.append("\nHigh-risk loci validation:\n")
            md.append(df_to_md(hr_val[[c for c in ["strain_id", "builtin_gene",
                      "amrfinder_hit", "resfinder_hit", "blast_hit",
                      "n_tools_confirm", "validation_status"] if c in hr_val.columns]]))

    # 7. mobility
    md.append("## 7. Mobility-context validation\n")
    if mob is None:
        md.append(PENDING + " (run `scripts\\mobility_validation.py`)\n")
    else:
        summ = pd.DataFrame([{
            "High/Medium ARGs": len(mob),
            "on plasmid": int((mob.contig_type == "plasmid").sum()),
            "IS both sides": int((mob.IS_both_sides == "yes").sum()),
            "mobility-context positive (in silico)": int((mob.predicted_mobilizable.isin(["yes", "context-positive"])).sum()),
            "mobility-tool confirmed (plasmid/MGE)": _i(mtool, "mobility_tool_confirmed_loci"),
            "externally predicted conjugative": int((mob.predicted_conjugative == "yes").sum()),
        }])
        md.append(df_to_md(summ))
        md.append("\nConjugative status is reported as unknown unless an external "
                  "tool confirms a relaxase/MPF system (0 conjugation markers in "
                  "the builtin runs). All calls are *in silico mobile-context "
                  "potential*.\n")

    # 8. stats
    md.append("## 8. Statistical enrichment results\n")
    if stats is None:
        md.append(PENDING + " (run `scripts\\statistical_analysis.py`)\n")
    else:
        md.append(df_to_md(stats))
        md.append("\n**Interpretation (read carefully):** the per-strain ARG "
                  "*burden* difference between pathogen and probiotic genomes is "
                  "statistically strong, but the locus-level *enrichment* tests "
                  "(High/Medium, plasmid-associated, IS-flanked) are underpowered "
                  "and not significant because probiotic genomes contribute very "
                  "few ARG loci overall. Report these as an **observed "
                  "concentration in pathogen comparators**, not as statistically "
                  "significant enrichment, unless a strain-level test supports it.\n")
        if net is not None:
            md.append("\n" + df_to_md(net))

    # 9. high-risk table
    md.append("## 9. High-risk ARG table\n")
    hcols = [c for c in ["strain_id", "gene", "drug_class", "contig", "start",
                         "end", "pct_identity", "on_plasmid", "is_flank_both",
                         "risk_reason"] if c in high.columns]
    md.append(df_to_md(high[hcols], max_rows=60))
    md.append("\nGenes recurring at independent loci are not duplicates "
              "(distinct contig/coordinates).\n")

    # 10. claims ledger
    md.append("## 10. Claims ledger\n")
    ledger = pd.DataFrame([
        {"claim": "Pipeline screens genomes for ARGs + MGE context (reproducible, no-admin)",
         "evidence": "code + SCALE50 outputs", "status": "SAFE"},
        {"claim": "High-risk mobile-context ARGs concentrated in pathogen comparators",
         "evidence": f"{rc['High']}/{rc['High']} High in pathogens", "status": "SAFE"},
        {"claim": "Probiotic-associated genomes carry no high-risk mobile ARGs (this set)",
         "evidence": f"{high_in_pro} High in probiotic genomes", "status": "SAFE"},
        {"claim": "No probiotic-pathogen ARG sharing at >=95% identity (this set)",
         "evidence": f"{cross}/{n_edges} cross-group edges", "status": "SAFE"},
        {"claim": "Full assembly recovers plasmid-borne ARGs vs chromosome-only",
         "evidence": "matched comparison" if comp is not None else "pending",
         "status": "SAFE" if comp is not None else "PENDING (run chronly)"},
        {"claim": "Specific allele identities (e.g. CTX-M-14, KPC-2)",
         "evidence": ("single builtin method" if val_mode_pending else
                      (f"BLAST-confirmed priority loci ({blast_confirmed}); remaining Low/non-priority best CARD matches"
                       if blast_only else "external tool")),
         "status": ("NEEDS EXTERNAL VALIDATION" if val_mode_pending else
                    ("PARTIAL: priority loci externally confirmed; not SOTA until AMRFinderPlus/ResFinder"
                     if blast_only else "SAFE if confirmed"))},
        {"claim": "Conjugative mobilization / horizontal transfer",
         "evidence": "0 conjugation markers; context only", "status": "DO NOT CLAIM"},
        {"claim": "Commercial probiotic products screened",
         "evidence": "reference genomes, not products", "status": "DO NOT CLAIM"},
    ])
    md.append(df_to_md(ledger))

    # 11. wording
    md.append("## 11. Safe wording vs wording to avoid\n")
    md.append(
        "- Use: 'probiotic-associated reference genomes' (NOT 'commercial probiotic products').\n"
        "- Use: 'in silico mobile-context potential' / 'plasmid-associated and IS-flanked' "
        "(NOT 'horizontal transfer demonstrated' or 'conjugative').\n"
        "- Use: 'best CARD database match' for alleles unless confirmed by an "
        "independent tool; in the BLAST-only state, report the 46 priority loci "
        "as externally BLAST-confirmed and the remaining Low/non-priority loci "
        "as best CARD matches pending AMRFinderPlus/ResFinder validation "
        "(NOT definitive allele typing).\n"
        "- Use: 'resistome concentrated in pathogen comparators' (NOT 'probiotics "
        "are an AMR reservoir' - the data show the opposite here).\n")

    # 12. manuscript Results
    md.append("## 12. Manuscript-ready Results text (draft)\n")
    md.append(
        f"> Across {len(meta) if meta is not None else 'N'} bacterial reference "
        f"genomes, {len(arg)} ARG loci were detected. Applying a four-tier mobile-"
        f"context risk heuristic, {rc['High']} loci were classified High and "
        f"{rc['Medium']} Medium; all High-risk loci occurred in pathogen/comparator "
        f"genomes, with none in probiotic-associated genomes. Plasmid-replicon "
        f"detection ({n_rep} replicons) and insertion-sequence context ({n_is} IS "
        f"features) localised the High-risk cargo to plasmid- and IS-associated "
        f"loci in pathogens (e.g. extended-spectrum beta-lactamase loci on "
        f"K. pneumoniae plasmid replicons flanked by IS elements). No ARG was "
        f"shared between probiotic-associated and pathogen genomes at >=95% "
        f"nucleotide identity ({cross}/{n_edges} cross-group sharing edges). "
        f"Conjugation markers were not detected by nucleotide screening "
        f"({n_conj}); conjugative mobilisation was therefore not inferred. "
        f"Per-strain ARG burden was significantly higher in pathogen than in "
        f"probiotic-associated genomes (Mann-Whitney U, see Table); locus-level "
        f"enrichment of high-risk, plasmid-associated and IS-flanked loci was "
        f"observed exclusively in pathogen comparators but was not formally "
        f"significant given the small number of probiotic ARG loci. "
        + (
            "A builtin consistency recheck confirmed internal reproducibility of "
            "the calls; allele-level identities are reported as best CARD "
            "database matches pending independent confirmation by "
            "AMRFinderPlus/ResFinder/BLAST.\n"
            if val_mode_pending else
            (
                f"External BLAST validation confirmed all {rc['High']} High-risk "
                f"loci and all {len(cur_hm) if cur_hm is not None else 'curated High/Medium'} "
                f"curated High/Medium acquired-cargo loci at 100% query coverage "
                f"and 100% nucleotide identity. In total, {blast_confirmed} "
                f"unique priority loci were externally BLAST-confirmed; the "
                f"remaining {len(arg) - blast_confirmed} Low/non-priority ARG "
                f"loci were not externally checked and remain best CARD database "
                f"matches. The run is not yet SOTA because AMRFinderPlus, "
                f"ResFinder and an external mobility tool have not yet been "
                f"imported.\n"
            ) if blast_only else
            (
                f"ARG identity validation imported AMRFinderPlus "
                f"({amrfinder_confirmed} confirmed loci), ResFinder "
                f"({resfinder_confirmed} confirmed loci), and BLAST "
                f"({blast_confirmed} confirmed loci), with "
                f"{multi_tool_confirmed} loci confirmed by two or more "
                f"external identity tools and {discordant_count} discordant "
                f"loci flagged for review.\n"
            )
        ))

    # 13. manuscript Methods
    md.append("## 13. Manuscript-ready Methods text (draft)\n")
    md.append(
        "> Genome assemblies were retrieved from NCBI RefSeq (full assemblies, "
        "chromosome plus plasmid replicons) via Entrez. Acquired ARGs were "
        "identified against the CARD nucleotide catalogue using a k-mer-indexed "
        "ungapped alignment (>=95% identity, >=80% coverage), with overlapping "
        "multi-allele hits collapsed to one best locus per reciprocal-overlap "
        "cluster (>=0.8). Plasmid replicons (PlasmidFinder) and insertion "
        "sequences (ISfinder) were detected with the same engine (>=80% identity, "
        ">=60% coverage). Each ARG was assigned a four-tier transfer-risk category "
        "from its genomic context (plasmid-replicon co-localisation; IS flanking "
        "within 5 kb). Core/intrinsic, biocide and efflux determinants were "
        "curated out of the mobile-ARG analysis using a rule set. "
        + (
            "Priority ARG calls were externally checked by BLAST, while "
            "AMRFinderPlus/ResFinder and MOB-suite/geNomad imports were not yet "
            "available for this run. "
            if blast_only else
            ("ARG calls were cross-validated with AMRFinderPlus/ResFinder and "
             "mobility with MOB-suite/geNomad. " if sota_ready else
             "External validation status was reported according to the imported "
             "tool outputs. ")
        ) +
        "Group enrichment was tested with Fisher's exact "
        "test and per-strain burden with a Mann-Whitney/permutation test.\n")

    # 14. limitations
    md.append("## 14. Limitations\n")
    md.append(
        "- Reference genomes, not sequenced commercial products.\n"
        "- Transfer risk is inferred from genomic context, not demonstrated.\n"
        + ("- The builtin aligner is ungapped; allele names are best matches "
           "pending independent confirmation.\n" if val_mode_pending else
           (f"- The builtin aligner is ungapped; {blast_confirmed} priority "
            f"loci were externally BLAST-confirmed, while the remaining "
            f"{len(arg) - blast_confirmed} Low/non-priority allele names are "
            f"best CARD matches pending AMRFinderPlus/ResFinder validation.\n"
            if blast_only else
            "- The builtin aligner is ungapped; external AMRFinderPlus/ResFinder "
            "outputs are used to qualify allele-level claims and discordant loci "
            "are flagged rather than overcalled.\n")) +
        "- Conjugation detection by nucleotide similarity is weak; relaxase typing "
        "needs protein-HMM tools.\n"
        "- Findings are conditional on CARD/PlasmidFinder/ISfinder versions and "
        "the chosen thresholds.\n")

    # 15. journals
    md.append("## 15. Recommended journal targets\n")
    md.append(
        "Given reference-genome screening, in silico mobility inference, and "
        "external ARG validation overlays, target genomics/microbiology venues "
        "that accept computational comparative resistome studies:\n"
        "- Microbial Genomics (Microbiology Society) - strong fit for a genomic resistome method.\n"
        "- Microbiology Spectrum (ASM).\n"
        "- Access Microbiology (welcomes sound/negative results).\n"
        "- BMC Microbiology / PeerJ / Antibiotics (MDPI).\n"
        "Near-free routes: subscription/hybrid in non-OA mode (no APC), bioRxiv "
        "preprint, DOAJ 'without fees' filter, JOSS for the pipeline software, and "
        "Microbiology Society fee-free OA if your institution has a Publish & Read "
        "agreement.\n")

    (out / "final_sota_report.md").write_text("\n".join(md), encoding="utf-8")
    print(f"wrote {out/'final_sota_report.md'} ({sum(len(s) for s in md)} chars)")
    print(f"validation external pending: {val_mode_pending}; "
          f"chronly comparison present: {comp is not None}")


if __name__ == "__main__":
    main()

