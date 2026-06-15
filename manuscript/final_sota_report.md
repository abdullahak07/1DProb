# SCALE50 SOTA v2 - Probiotic-associated AMR-cargo screen

## 1. Executive summary

We screened 50 bacterial reference genomes (25 probiotic-associated, 25 pathogen/comparator) for antibiotic-resistance genes (ARGs) and their mobile-genetic-element context using a no-admin, stdlib Python pipeline (CARD nucleotide catalogue, PlasmidFinder replicons, ISfinder insertion sequences). We detected **373 ARG loci** (15 High, 43 Medium, 315 Low), **32 plasmid replicons**, **1116 IS features**, and **0 conjugation markers**. **All 15 High-risk ARGs are in pathogen/comparator genomes (0 in probiotic-associated genomes)**, and **0 of 175 gene-sharing edges are probiotic-pathogen**. The data indicate that probiotic-associated reference genomes in this set carry little high-risk mobile ARG cargo and share no ARGs with pathogens at high identity; high-risk, mobile-context ARGs are concentrated in the pathogen comparators. This is a reassuring-safety + method finding, not evidence of probiotic-to-pathogen transfer.

## 2. Dataset summary

50 strains; 50 with genomes retrieved (full assembly: chromosome + plasmid replicons).

| strain_id | species | group | accession | n_contigs | assembly_quality | genome_available |
| --- | --- | --- | --- | --- | --- | --- |
| L_rhamnosus_GG | Lacticaseibacillus rhamnosus | probiotic | NC_013198.1 | 1 | complete | True |
| L_acidophilus_NCFM | Lactobacillus acidophilus | probiotic | NC_006814.3 | 1 | complete | True |
| L_plantarum_WCFS1 | Lactiplantibacillus plantarum | probiotic | NC_004567.2 | 4 | draft | True |
| L_casei_BL23 | Lacticaseibacillus casei | probiotic | NC_010999.1 | 1 | complete | True |
| L_reuteri_DSM20016 | Limosilactobacillus reuteri | probiotic | CP000705.1 | 1 | complete | True |
| B_longum_NCC2705 | Bifidobacterium longum | probiotic | NC_004307.2 | 2 | complete | True |
| B_animalis_lactis_AD011 | Bifidobacterium animalis subsp. lactis | probiotic | NC_011835.1 | 1 | complete | True |
| S_thermophilus_LMG18311 | Streptococcus thermophilus | probiotic | NC_006448.1 | 1 | complete | True |
| L_lactis_IL1403 | Lactococcus lactis subsp. lactis | probiotic | NC_002662.1 | 1 | complete | True |
| B_subtilis_168 | Bacillus subtilis | probiotic | NC_000964.3 | 1 | complete | True |
| L_bulgaricus_ATCC11842 | Lactobacillus delbrueckii subsp. bulgaricus | probiotic | NC_008054.1 | 1 | complete | True |
| L_gasseri_ATCC33323 | Lactobacillus gasseri | probiotic | NC_008530.1 | 1 | complete | True |
| L_johnsonii_NCC533 | Lactobacillus johnsonii | probiotic | NC_005362.1 | 1 | complete | True |
| L_salivarius_UCC118 | Ligilactobacillus salivarius | probiotic | NC_007929.1 | 4 | draft | True |
| B_adolescentis_ATCC15703 | Bifidobacterium adolescentis | probiotic | NC_008618.1 | 1 | complete | True |
| B_bifidum_PRL2010 | Bifidobacterium bifidum | probiotic | NC_014638.1 | 1 | complete | True |
| B_breve_UCC2003 | Bifidobacterium breve | probiotic | NC_020517.1 | 1 | complete | True |
| B_infantis_ATCC15697 | Bifidobacterium longum subsp. infantis | probiotic | NC_011593.1 | 1 | complete | True |
| P_pentosaceus_ATCC25745 | Pediococcus pentosaceus | probiotic | NC_008525.1 | 1 | complete | True |
| L_mesenteroides_ATCC8293 | Leuconostoc mesenteroides subsp. mesenteroides | probiotic | NC_008531.1 | 2 | complete | True |
| L_lactis_cremoris_MG1363 | Lactococcus lactis subsp. cremoris | probiotic | NC_009004.1 | 1 | complete | True |
| L_fermentum_IFO3956 | Limosilactobacillus fermentum | probiotic | NC_010610.1 | 1 | complete | True |
| L_helveticus_DPC4571 | Lactobacillus helveticus | probiotic | NC_010080.1 | 1 | complete | True |
| L_crispatus_ST1 | Lactobacillus crispatus | probiotic | NC_014106.1 | 1 | complete | True |
| L_sakei_23K | Latilactobacillus sakei | probiotic | NC_007576.1 | 1 | complete | True |
| E_faecium_Aus0004 | Enterococcus faecium | pathogen | NC_017022.1 | 1 | complete | True |
| E_faecalis_V583 | Enterococcus faecalis | pathogen | NC_004668.1 | 4 | draft | True |
| E_faecalis_OG1RF | Enterococcus faecalis | pathogen | NC_017316.1 | 1 | complete | True |
| K_pneumoniae_HS11286 | Klebsiella pneumoniae | pathogen | NC_016845.1 | 7 | draft | True |
| K_pneumoniae_MGH78578 | Klebsiella pneumoniae | pathogen | NC_009648.1 | 6 | draft | True |
| E_coli_K12_MG1655 | Escherichia coli | pathogen | NC_000913.3 | 1 | complete | True |
| E_coli_O157H7_Sakai | Escherichia coli O157 H7 | pathogen | NC_002695.2 | 3 | complete | True |
| S_enterica_LT2 | Salmonella enterica serovar Typhimurium | pathogen | NC_003197.2 | 2 | complete | True |
| Shigella_flexneri_2a_301 | Shigella flexneri | pathogen | NC_004337.2 | 2 | complete | True |
| Enterobacter_cloacae_ATCC13047 | Enterobacter cloacae subsp. cloacae | pathogen | NC_014121.1 | 3 | complete | True |
| Proteus_mirabilis_HI4320 | Proteus mirabilis | pathogen | NC_010554.1 | 2 | complete | True |
| C_difficile_630 | Clostridioides difficile | pathogen | NC_009089.1 | 2 | complete | True |
| C_perfringens_str13 | Clostridium perfringens | pathogen | NC_003366.1 | 2 | complete | True |
| S_aureus_NCTC8325 | Staphylococcus aureus | pathogen | NC_007795.1 | 1 | complete | True |
| S_aureus_N315 | Staphylococcus aureus | pathogen | NC_002745.2 | 2 | complete | True |
| S_aureus_Mu50 | Staphylococcus aureus | pathogen | NC_002758.2 | 2 | complete | True |
| L_monocytogenes_EGDe | Listeria monocytogenes | pathogen | NC_003210.1 | 1 | complete | True |
| S_agalactiae_2603VR | Streptococcus agalactiae | pathogen | NC_004116.1 | 1 | complete | True |
| P_aeruginosa_PAO1 | Pseudomonas aeruginosa | pathogen | NC_002516.2 | 1 | complete | True |
| A_baumannii_ATCC17978 | Acinetobacter baumannii | pathogen | NC_009085.1 | 1 | complete | True |
| S_maltophilia_K279a | Stenotrophomonas maltophilia | pathogen | NC_010943.1 | 1 | complete | True |
| C_jejuni_NCTC11168 | Campylobacter jejuni | pathogen | NC_002163.1 | 1 | complete | True |
| H_pylori_26695 | Helicobacter pylori | pathogen | NC_000915.1 | 1 | complete | True |
| B_fragilis_NCTC9343 | Bacteroides fragilis | pathogen | NC_003228.3 | 2 | complete | True |
| V_cholerae_N16961 | Vibrio cholerae | pathogen | NC_002505.1 | 2 | complete | True |

## 3. Full-assembly SCALE50 results

| metric | value |
| --- | --- |
| ARG loci | 373 |
| High | 15 |
| Medium | 43 |
| Low | 315 |
| insertion sequences | 1116 |
| plasmid replicons | 32 |
| conjugation markers | 0 |
| gene-sharing edges (total) | 175 |
| probiotic-pathogen edges | 0 |

## 4. Matched chromosome-only vs full-assembly comparison

Same 50 strains, full_assembly true vs false:

| metric | chromosome_only | full_assembly | delta |
| --- | --- | --- | --- |
| arg_loci | 338 | 373 | 35 |
| Low | 313 | 315 | 2 |
| Medium | 21 | 43 | 22 |
| High | 4 | 15 | 11 |
| Negligible | 0 | 0 | 0 |
| insertion_sequence | 1003 | 1116 | 113 |
| plasmid_replicon | 0 | 32 | 32 |
| conjugation_marker | 0 | 0 | 0 |

## 5. Curated ARG interpretation

Raw calls preserved; core/intrinsic/biocide/efflux determinants flagged and excluded from the mobile-ARG analysis.

| curated_category | count |
| --- | --- |
| acquired_mobile_ARG | 316 |
| efflux_or_regulatory_marker | 39 |
| intrinsic_core_gene | 9 |
| ambiguous_requires_validation | 5 |
| biocide_or_stress_marker | 4 |


Curated High/Medium acquired-cargo loci: **45** (raw High+Medium was 58).

## 6. Independent validation summary

| validation_status | count |
| --- | --- |
| not_checked | 327 |
| confirmed | 46 |


High-risk loci validation:

| strain_id | builtin_gene | validation_best_hit | validation_identity | validation_tool | validation_status |
| --- | --- | --- | --- | --- | --- |
| E_faecalis_V583 | ErmB | ErmB | 100 | BLAST | confirmed |
| E_faecalis_V583 | qacH | qacH | 100 | BLAST | confirmed |
| E_faecalis_V583 | AAC6_Ie_APH2_Ia | AAC6_Ie_APH2_Ia | 100 | BLAST | confirmed |
| K_pneumoniae_HS11286 | CTX-M-14 | CTX-M-14 | 100 | BLAST | confirmed |
| K_pneumoniae_HS11286 | CTX-M-14 | CTX-M-14 | 100 | BLAST | confirmed |
| K_pneumoniae_MGH78578 | APH(3'')-Ib | APH(3'')-Ib | 100 | BLAST | confirmed |
| K_pneumoniae_MGH78578 | APH(6)-Id | APH(6)-Id | 100 | BLAST | confirmed |
| K_pneumoniae_MGH78578 | sul2 | sul2 | 100 | BLAST | confirmed |
| K_pneumoniae_MGH78578 | APH(3')-Ia | APH(3')-Ia | 100 | BLAST | confirmed |
| K_pneumoniae_MGH78578 | catA1 | catA1 | 100 | BLAST | confirmed |
| Enterobacter_cloacae_ATCC13047 | sul2 | sul2 | 100 | BLAST | confirmed |
| S_aureus_N315 | ANT(4')-Ia | ANT(4')-Ia | 100 | BLAST | confirmed |
| S_aureus_Mu50 | ANT(4')-Ia | ANT(4')-Ia | 100 | BLAST | confirmed |
| S_aureus_Mu50 | AAC6_Ie_APH2_Ia | AAC6_Ie_APH2_Ia | 100 | BLAST | confirmed |
| A_baumannii_ATCC17978 | sul2 | sul2 | 100 | BLAST | confirmed |

## 7. Mobility-context validation

| High/Medium ARGs | on plasmid | IS both sides | mobility-context positive (in silico) | externally predicted conjugative |
| --- | --- | --- | --- | --- |
| 58 | 28 | 15 | 36 | 0 |


Conjugative status is reported as unknown unless an external tool confirms a relaxase/MPF system (0 conjugation markers in the builtin runs). All calls are *in silico mobile-context potential*.

## 8. Statistical enrichment results

| analysis | test | statistic | p_value | detail |
| --- | --- | --- | --- | --- |
| mobile ARG burden pathogen vs probiotic | Mann-Whitney U | 606.50 | 0.00 | pathogen mean=12.32, probiotic mean=0.52 |
| enrichment: High/Medium risk | Fisher exact (2-sided) | inf | 0.38 | pathogen 45/360 vs probiotic 0/13 |
| enrichment: plasmid-associated | Fisher exact (2-sided) | inf | 0.61 | pathogen 28/360 vs probiotic 0/13 |
| enrichment: IS-flanked both sides | Fisher exact (2-sided) | inf | 1 | pathogen 15/360 vs probiotic 0/13 |


**Interpretation (read carefully):** the per-strain ARG *burden* difference between pathogen and probiotic genomes is statistically strong, but the locus-level *enrichment* tests (High/Medium, plasmid-associated, IS-flanked) are underpowered and not significant because probiotic genomes contribute very few ARG loci overall. Report these as an **observed concentration in pathogen comparators**, not as statistically significant enrichment, unless a strain-level test supports it.


| total_edges | within_probiotic | within_pathogen | probiotic_pathogen | n_edge_genes | edge_genes_sample |
| --- | --- | --- | --- | --- | --- |
| 175 | 0 | 175 | 0 | 81 | AAC6_Ie_APH2_Ia, ANT(4')-Ia, ANT(9)-Ia, APH(3'')-Ib, APH(6)-Id, AcrE, AcrS, ArnT, CRP, Ecol_acrA, Ecol_emrE, Ecol_mdfA, ErmA, ErmB, FosA6, H-NS, IreK, Kpne_KpnE, Kpne_KpnF, Kpne_KpnG, LptD, MdtQ, OmpA, PmrF, Saur_FosB |

## 9. High-risk ARG table

| strain_id | gene | drug_class | contig | start | end | pct_identity | on_plasmid | is_flank_both | risk_reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E_faecalis_V583 | ErmB | macrolide | NC_004669.1 | 6333 | 7079 | 97.99 | False | True | IS elements flanking both sides |
| E_faecalis_V583 | qacH | unclassified | NC_004669.1 | 8712 | 9035 | 98.46 | False | True | IS elements flanking both sides |
| E_faecalis_V583 | AAC6_Ie_APH2_Ia | aminoglycoside | NC_004669.1 | 49141 | 50580 | 100 | False | True | IS elements flanking both sides |
| K_pneumoniae_HS11286 | CTX-M-14 | beta-lactam | NC_016838.1 | 31843 | 32718 | 100 | True | True | IS elements flanking both sides |
| K_pneumoniae_HS11286 | CTX-M-14 | beta-lactam | NC_016839.1 | 82527 | 83402 | 100 | True | True | IS elements flanking both sides |
| K_pneumoniae_MGH78578 | APH(3'')-Ib | aminoglycoside | NC_009651.1 | 12226 | 13029 | 99.88 | True | True | IS elements flanking both sides |
| K_pneumoniae_MGH78578 | APH(6)-Id | aminoglycoside | NC_009651.1 | 13029 | 13865 | 99.88 | True | True | IS elements flanking both sides |
| K_pneumoniae_MGH78578 | sul2 | sulfonamide | NC_009651.1 | 18529 | 19344 | 100 | True | True | IS elements flanking both sides |
| K_pneumoniae_MGH78578 | APH(3')-Ia | aminoglycoside | NC_009651.1 | 27290 | 28105 | 100 | True | True | IS elements flanking both sides |
| K_pneumoniae_MGH78578 | catA1 | phenicol | NC_009651.1 | 66239 | 66898 | 99.55 | True | True | IS elements flanking both sides |
| Enterobacter_cloacae_ATCC13047 | sul2 | sulfonamide | NC_014121.1 | 3902508 | 3903323 | 100 | False | True | IS elements flanking both sides |
| S_aureus_N315 | ANT(4')-Ia | aminoglycoside | NC_002745.2 | 40818 | 41579 | 99.87 | False | True | IS elements flanking both sides |
| S_aureus_Mu50 | ANT(4')-Ia | aminoglycoside | NC_002758.2 | 40819 | 41580 | 99.87 | False | True | IS elements flanking both sides |
| S_aureus_Mu50 | AAC6_Ie_APH2_Ia | aminoglycoside | NC_002774.1 | 17224 | 18663 | 100 | False | True | IS elements flanking both sides |
| A_baumannii_ATCC17978 | sul2 | sulfonamide | NC_009085.1 | 798900 | 799715 | 100 | False | True | IS elements flanking both sides |


Genes recurring at independent loci are not duplicates (distinct contig/coordinates).

## 10. Claims ledger

| claim | evidence | status |
| --- | --- | --- |
| Pipeline screens genomes for ARGs + MGE context (reproducible, no-admin) | code + SCALE50 outputs | SAFE |
| High-risk mobile-context ARGs concentrated in pathogen comparators | 15/15 High in pathogens | SAFE |
| Probiotic-associated genomes carry no high-risk mobile ARGs (this set) | 0 High in probiotic genomes | SAFE |
| No probiotic-pathogen ARG sharing at >=95% identity (this set) | 0/175 cross-group edges | SAFE |
| Full assembly recovers plasmid-borne ARGs vs chromosome-only | matched comparison | SAFE |
| Specific allele identities (e.g. CTX-M-14, KPC-2) | external tool | SAFE if confirmed |
| Conjugative mobilization / horizontal transfer | 0 conjugation markers; context only | DO NOT CLAIM |
| Commercial probiotic products screened | reference genomes, not products | DO NOT CLAIM |

## 11. Safe wording vs wording to avoid

- Use: 'probiotic-associated reference genomes' (NOT 'commercial probiotic products').
- Use: 'in silico mobile-context potential' / 'plasmid-associated and IS-flanked' (NOT 'horizontal transfer demonstrated' or 'conjugative').
- Use: 'best CARD database match' for alleles unless confirmed by an independent tool (NOT definitive allele typing).
- Use: 'resistome concentrated in pathogen comparators' (NOT 'probiotics are an AMR reservoir' - the data show the opposite here).

## 12. Manuscript-ready Results text (draft)

> Across 50 bacterial reference genomes, 373 ARG loci were detected. Applying a four-tier mobile-context risk heuristic, 15 loci were classified High and 43 Medium; all High-risk loci occurred in pathogen/comparator genomes, with none in probiotic-associated genomes. Plasmid-replicon detection (32 replicons) and insertion-sequence context (1116 IS features) localised the High-risk cargo to plasmid- and IS-associated loci in pathogens (e.g. extended-spectrum beta-lactamase loci on K. pneumoniae plasmid replicons flanked by IS elements). No ARG was shared between probiotic-associated and pathogen genomes at >=95% nucleotide identity (0/175 cross-group sharing edges). Conjugation markers were not detected by nucleotide screening (0); conjugative mobilisation was therefore not inferred. Per-strain ARG burden was significantly higher in pathogen than in probiotic-associated genomes (Mann-Whitney U, see Table); locus-level enrichment of high-risk, plasmid-associated and IS-flanked loci was observed exclusively in pathogen comparators but was not formally significant given the small number of probiotic ARG loci. External BLAST validation confirmed all 15 High-risk loci and all 45 curated High/Medium acquired-cargo loci at 100% query coverage and 100% nucleotide identity. In total, 46 unique priority loci were externally BLAST-confirmed; the remaining 327 Low/non-priority ARG loci were not externally checked and remain best CARD database matches.

## 13. Manuscript-ready Methods text (draft)

> Genome assemblies were retrieved from NCBI RefSeq (full assemblies, chromosome plus plasmid replicons) via Entrez. Acquired ARGs were identified against the CARD nucleotide catalogue using a k-mer-indexed ungapped alignment (>=95% identity, >=80% coverage), with overlapping multi-allele hits collapsed to one best locus per reciprocal-overlap cluster (>=0.8). Plasmid replicons (PlasmidFinder) and insertion sequences (ISfinder) were detected with the same engine (>=80% identity, >=60% coverage). Each ARG was assigned a four-tier transfer-risk category from its genomic context (plasmid-replicon co-localisation; IS flanking within 5 kb). Core/intrinsic, biocide and efflux determinants were curated out of the mobile-ARG analysis using a rule set. [If applicable:] ARG calls were cross-validated with AMRFinderPlus/ResFinder and mobility with MOB-suite/geNomad. Group enrichment was tested with Fisher's exact test and per-strain burden with a Mann-Whitney/permutation test.

## 14. Limitations

- Reference genomes, not sequenced commercial products.
- Transfer risk is inferred from genomic context, not demonstrated.
- The builtin aligner is ungapped; all 15 High-risk loci and all 45 curated High/Medium acquired-cargo loci were externally BLAST-confirmed, while remaining Low/non-priority allele names are best CARD matches.
- Conjugation detection by nucleotide similarity is weak; relaxase typing needs protein-HMM tools.
- Findings are conditional on CARD/PlasmidFinder/ISfinder versions and the chosen thresholds.

## 15. Recommended journal targets

Given single-method calls, reference genomes, and a largely negative comparative result, target mid/low-IF venues, not top-tier:
- Microbial Genomics (Microbiology Society) - strong fit for a genomic resistome method.
- Microbiology Spectrum (ASM).
- Access Microbiology (welcomes sound/negative results).
- BMC Microbiology / PeerJ / Antibiotics (MDPI).
Near-free routes: subscription/hybrid in non-OA mode (no APC), bioRxiv preprint, DOAJ 'without fees' filter, JOSS for the pipeline software, and Microbiology Society fee-free OA if your institution has a Publish & Read agreement.

