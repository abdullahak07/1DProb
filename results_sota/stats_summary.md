# Statistical analysis (SCALE50)

Engine: scipy. ARG set: curated (intrinsic/biocide excluded).

## Key tests

| analysis | test | statistic | p_value | detail |
| --- | --- | --- | --- | --- |
| mobile ARG burden pathogen vs probiotic | Mann-Whitney U | 606.50 | 0.00 | pathogen mean=12.32, probiotic mean=0.52 |
| enrichment: High/Medium risk | Fisher exact (2-sided) | inf | 0.38 | pathogen 45/360 vs probiotic 0/13 |
| enrichment: plasmid-associated | Fisher exact (2-sided) | inf | 0.61 | pathogen 28/360 vs probiotic 0/13 |
| enrichment: IS-flanked both sides | Fisher exact (2-sided) | inf | 1 | pathogen 15/360 vs probiotic 0/13 |


## Enrichment detail

| comparison | pathogen_yes | pathogen_no | probiotic_yes | probiotic_no | odds_ratio | p_value |
| --- | --- | --- | --- | --- | --- | --- |
| High/Medium risk | 45 | 315 | 0 | 13 | inf | 0.38 |
| plasmid-associated | 28 | 332 | 0 | 13 | inf | 0.61 |
| IS-flanked both sides | 15 | 345 | 0 | 13 | inf | 1 |


## Network summary

| total_edges | within_probiotic | within_pathogen | probiotic_pathogen | n_edge_genes | edge_genes_sample |
| --- | --- | --- | --- | --- | --- |
| 175 | 0 | 175 | 0 | 81 | AAC6_Ie_APH2_Ia, ANT(4')-Ia, ANT(9)-Ia, APH(3'')-Ib, APH(6)-Id, AcrE, AcrS, ArnT, CRP, Ecol_acrA, Ecol_emrE, Ecol_mdfA, ErmA, ErmB, FosA6, H-NS, IreK, Kpne_KpnE, Kpne_KpnF, Kpne_KpnG, LptD, MdtQ, OmpA, PmrF, Saur_FosB |


## Per-strain burden (top of table)

| strain_id | group | total_arg_loci | mobile_arg_count | high | medium |
| --- | --- | --- | --- | --- | --- |
| P_aeruginosa_PAO1 | pathogen | 49 | 48 | 0 | 0 |
| S_aureus_N315 | pathogen | 26 | 26 | 1 | 3 |
| E_coli_K12_MG1655 | pathogen | 41 | 25 | 0 | 0 |
| E_coli_O157H7_Sakai | pathogen | 40 | 24 | 0 | 0 |
| K_pneumoniae_MGH78578 | pathogen | 27 | 23 | 5 | 9 |
| S_aureus_Mu50 | pathogen | 22 | 21 | 2 | 3 |
| K_pneumoniae_HS11286 | pathogen | 23 | 19 | 2 | 10 |
| Shigella_flexneri_2a_301 | pathogen | 27 | 19 | 0 | 6 |
| A_baumannii_ATCC17978 | pathogen | 19 | 19 | 1 | 0 |
| E_faecalis_V583 | pathogen | 15 | 14 | 2 | 0 |
| S_aureus_NCTC8325 | pathogen | 12 | 12 | 0 | 0 |
| B_subtilis_168 | probiotic | 10 | 10 | 0 | 0 |
| S_maltophilia_K279a | pathogen | 9 | 9 | 0 | 0 |
| E_faecium_Aus0004 | pathogen | 9 | 9 | 0 | 0 |
| C_difficile_630 | pathogen | 6 | 6 | 0 | 0 |
| S_enterica_LT2 | pathogen | 7 | 6 | 0 | 0 |
| C_jejuni_NCTC11168 | pathogen | 5 | 5 | 0 | 0 |
| V_cholerae_N16961 | pathogen | 5 | 5 | 0 | 0 |
| E_faecalis_OG1RF | pathogen | 4 | 4 | 0 | 0 |
| Enterobacter_cloacae_ATCC13047 | pathogen | 3 | 3 | 1 | 0 |
| L_monocytogenes_EGDe | pathogen | 3 | 3 | 0 | 0 |
| S_agalactiae_2603VR | pathogen | 3 | 3 | 0 | 0 |
| Proteus_mirabilis_HI4320 | pathogen | 2 | 2 | 0 | 0 |
| L_lactis_IL1403 | probiotic | 1 | 1 | 0 | 0 |
| H_pylori_26695 | pathogen | 1 | 1 | 0 | 0 |


_Note: with all High-risk loci in pathogens and 0 probiotic-pathogen sharing edges, enrichment p-values describe a pathogen-concentrated resistome, not probiotic risk._
