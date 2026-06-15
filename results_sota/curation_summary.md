# ARG curation summary (SCALE50)

Raw calls are preserved; curation only adds columns and flags which rows are eligible for the mobile-ARG analysis.

## Category breakdown

| curated_category | count |
| --- | --- |
| acquired_mobile_ARG | 316 |
| efflux_or_regulatory_marker | 39 |
| intrinsic_core_gene | 9 |
| ambiguous_requires_validation | 5 |
| biocide_or_stress_marker | 4 |


High/Medium loci before curation: **58**; after removing core/intrinsic/biocide/efflux determinants: **45**.


## Excluded from mobile-ARG analysis (52 rows)

| strain_id | gene | raw_risk_level | curated_category | confidence | curation_reason |
| --- | --- | --- | --- | --- | --- |
| E_faecalis_V583 | qacH | High | biocide_or_stress_marker | high | Quaternary-ammonium (biocide) efflux; not a classic antibiotic ARG. |
| K_pneumoniae_HS11286 | qacEdelta1 | Medium | biocide_or_stress_marker | high | qacE-delta1; class 1 integron 3'-CS marker. Integron-associated biocide marker, not an antibiotic ARG itself. |
| K_pneumoniae_HS11286 | MdtQ | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| K_pneumoniae_HS11286 | ArnT | Low | intrinsic_core_gene | medium | 4-amino-4-deoxy-L-arabinose transferase; core LPS modification (polymyxin). Chromosomal. |
| K_pneumoniae_HS11286 | eptB | Medium | intrinsic_core_gene | high | Phosphoethanolamine transferase; core LPS-modification locus (polymyxin context). Chromosomal, not mobile cargo. |
| K_pneumoniae_MGH78578 | MdtQ | Medium | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| K_pneumoniae_MGH78578 | ArnT | Low | intrinsic_core_gene | medium | 4-amino-4-deoxy-L-arabinose transferase; core LPS modification (polymyxin). Chromosomal. |
| K_pneumoniae_MGH78578 | eptB | Medium | intrinsic_core_gene | high | Phosphoethanolamine transferase; core LPS-modification locus (polymyxin context). Chromosomal, not mobile cargo. |
| K_pneumoniae_MGH78578 | qacEdelta1 | Medium | biocide_or_stress_marker | high | qacE-delta1; class 1 integron 3'-CS marker. Integron-associated biocide marker, not an antibiotic ARG itself. |
| E_coli_K12_MG1655 | acrB | Low | efflux_or_regulatory_marker | medium | AcrAB-TolC RND efflux; intrinsic chromosomal multidrug efflux. |
| E_coli_K12_MG1655 | Ecol_acrA | Low | efflux_or_regulatory_marker | medium | AcrAB-TolC RND efflux component; intrinsic chromosomal. |
| E_coli_K12_MG1655 | Ecol_emrE | Medium | efflux_or_regulatory_marker | high | Small multidrug resistance (SMR) efflux pump; intrinsic/species-core in many taxa. Chromosomal. |
| E_coli_K12_MG1655 | mdtG | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtH | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | ugd | Medium | intrinsic_core_gene | high | UDP-glucose 6-dehydrogenase; core housekeeping gene in LPS/polymyxin modification. Chromosomal. |
| E_coli_K12_MG1655 | mdtA | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtB | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtC | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | PmrF | Low | intrinsic_core_gene | medium | arn/pmr operon LPS modification; core chromosomal polymyxin-associated locus. |
| E_coli_K12_MG1655 | mdtE | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtF | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtP | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtO | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtN | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_K12_MG1655 | mdtM | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | acrB | Low | efflux_or_regulatory_marker | medium | AcrAB-TolC RND efflux; intrinsic chromosomal multidrug efflux. |
| E_coli_O157H7_Sakai | Ecol_acrA | Low | efflux_or_regulatory_marker | medium | AcrAB-TolC RND efflux component; intrinsic chromosomal. |
| E_coli_O157H7_Sakai | mdtG | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtH | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | Ecol_emrE | Low | efflux_or_regulatory_marker | high | Small multidrug resistance (SMR) efflux pump; intrinsic/species-core in many taxa. Chromosomal. |
| E_coli_O157H7_Sakai | ugd | Low | intrinsic_core_gene | high | UDP-glucose 6-dehydrogenase; core housekeeping gene in LPS/polymyxin modification. Chromosomal. |
| E_coli_O157H7_Sakai | mdtA | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtB | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtC | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | PmrF | Low | intrinsic_core_gene | medium | arn/pmr operon LPS modification; core chromosomal polymyxin-associated locus. |
| E_coli_O157H7_Sakai | mdtE | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtF | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtP | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtO | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtN | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| E_coli_O157H7_Sakai | mdtM | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| S_enterica_LT2 | MdtK | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| Shigella_flexneri_2a_301 | acrB | Low | efflux_or_regulatory_marker | medium | AcrAB-TolC RND efflux; intrinsic chromosomal multidrug efflux. |
| Shigella_flexneri_2a_301 | Sfle_acrA | Low | efflux_or_regulatory_marker | medium | AcrAB-TolC RND efflux component; intrinsic chromosomal. |
| Shigella_flexneri_2a_301 | mdtH | Low | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| Shigella_flexneri_2a_301 | mdtB | Medium | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| Shigella_flexneri_2a_301 | PmrF | Low | intrinsic_core_gene | medium | arn/pmr operon LPS modification; core chromosomal polymyxin-associated locus. |
| Shigella_flexneri_2a_301 | mdtF | Medium | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| Shigella_flexneri_2a_301 | mdtP | Medium | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| Shigella_flexneri_2a_301 | mdtN | Medium | efflux_or_regulatory_marker | medium | Multidrug efflux transporter (mdt family); commonly intrinsic/chromosomal. |
| S_aureus_Mu50 | qacA | Medium | biocide_or_stress_marker | medium | Quaternary-ammonium-compound resistance (biocide/disinfectant marker). |
| P_aeruginosa_PAO1 | Paer_emrE | Low | efflux_or_regulatory_marker | high | Small multidrug resistance (SMR) efflux pump; intrinsic/species-core in many taxa. Chromosomal. |


## Curated High/Medium acquired cargo (kept)

| strain_id | gene | drug_class | curated_risk_level | confidence |
| --- | --- | --- | --- | --- |
| E_faecalis_V583 | ErmB | macrolide | High | medium |
| E_faecalis_V583 | AAC6_Ie_APH2_Ia | aminoglycoside | High | medium |
| K_pneumoniae_HS11286 | CTX-M-14 | beta-lactam | High | medium |
| K_pneumoniae_HS11286 | TEM-1 | beta-lactam | Medium | medium |
| K_pneumoniae_HS11286 | rmtB | aminoglycoside | Medium | medium |
| K_pneumoniae_HS11286 | cmlA9 | phenicol | Medium | medium |
| K_pneumoniae_HS11286 | aadA2 | aminoglycoside | Medium | medium |
| K_pneumoniae_HS11286 | TEM-1 | beta-lactam | Medium | medium |
| K_pneumoniae_HS11286 | APH(6)-Id | aminoglycoside | Medium | medium |
| K_pneumoniae_HS11286 | APH(3'')-Ib | aminoglycoside | Medium | medium |
| K_pneumoniae_HS11286 | sul2 | sulfonamide | Medium | medium |
| K_pneumoniae_HS11286 | CTX-M-14 | beta-lactam | High | medium |
| K_pneumoniae_HS11286 | KPC-2 | beta-lactam | Medium | medium |
| K_pneumoniae_HS11286 | TEM-1 | beta-lactam | Medium | medium |
| K_pneumoniae_MGH78578 | SHV-12 | beta-lactam | Medium | medium |
| K_pneumoniae_MGH78578 | AAC(6')-Ib' | aminoglycoside | Medium | medium |
| K_pneumoniae_MGH78578 | OXA-9 | beta-lactam | Medium | medium |
| K_pneumoniae_MGH78578 | TEM-90 | beta-lactam | Medium | medium |
| K_pneumoniae_MGH78578 | APH(3'')-Ib | aminoglycoside | High | medium |
| K_pneumoniae_MGH78578 | APH(6)-Id | aminoglycoside | High | medium |
| K_pneumoniae_MGH78578 | sul2 | sulfonamide | High | medium |
| K_pneumoniae_MGH78578 | APH(3')-Ia | aminoglycoside | High | medium |
| K_pneumoniae_MGH78578 | sul1 | sulfonamide | Medium | medium |
| K_pneumoniae_MGH78578 | cmlA5 | phenicol | Medium | medium |
| K_pneumoniae_MGH78578 | ANT(2'')-Ia | aminoglycoside | Medium | medium |
| K_pneumoniae_MGH78578 | catA1 | phenicol | High | medium |
| K_pneumoniae_MGH78578 | tet(D) | tetracycline | Medium | medium |
| K_pneumoniae_MGH78578 | TEM-166 | beta-lactam | Medium | medium |
| Shigella_flexneri_2a_301 | leuO | unclassified | Medium | low |
| Shigella_flexneri_2a_301 | H-NS | unclassified | Medium | low |
| Shigella_flexneri_2a_301 | marA | unclassified | Medium | low |
| Shigella_flexneri_2a_301 | gadX | unclassified | Medium | low |
| Shigella_flexneri_2a_301 | cpxA | unclassified | Medium | low |
| Shigella_flexneri_2a_301 | EC-8 | unclassified | Medium | low |
| Enterobacter_cloacae_ATCC13047 | sul2 | sulfonamide | High | medium |
| S_aureus_N315 | ANT(4')-Ia | aminoglycoside | High | medium |
| S_aureus_N315 | mecA | beta-lactam | Medium | low |
| S_aureus_N315 | mecR1 | beta-lactam | Medium | low |
| S_aureus_N315 | kdpD | unclassified | Medium | low |
| S_aureus_Mu50 | ANT(4')-Ia | aminoglycoside | High | medium |
| S_aureus_Mu50 | mecA | beta-lactam | Medium | low |
| S_aureus_Mu50 | mecR1 | beta-lactam | Medium | low |
| S_aureus_Mu50 | kdpD | unclassified | Medium | low |
| S_aureus_Mu50 | AAC6_Ie_APH2_Ia | aminoglycoside | High | medium |
| A_baumannii_ATCC17978 | sul2 | sulfonamide | High | medium |
