# ARG validation summary (multi-tool)

External identity tools imported: **AMRFinderPlus, ResFinder, BLAST**.

## Status counts

| validation_status | count |
| --- | --- |
| not_checked | 224 |
| confirmed_multi_tool | 68 |
| discordant | 33 |
| confirmed_amrfinder | 31 |
| confirmed_blast_only | 11 |
| confirmed_resfinder | 6 |


## Per-tool confirmed counts

| tool | imported | confirmed_loci |
| --- | --- | --- |
| AMRFinderPlus | True | 94 |
| ResFinder | True | 68 |
| BLAST | True | 46 |
| multi-tool (>=2 agree) | - | 68 |


## High-risk loci

| strain_id | builtin_gene | amrfinder_hit | resfinder_hit | blast_hit | n_tools_confirm | validation_status |
| --- | --- | --- | --- | --- | --- | --- |
| E_faecalis_V583 | ErmB | erm(B) | erm(B) | ErmB | 3 | confirmed_multi_tool |
| E_faecalis_V583 | qacH | qacZ |  | qacH | 1 | confirmed_blast_only |
| E_faecalis_V583 | AAC6_Ie_APH2_Ia | aac(6')-Ie/aph(2'')-Ia | aac(6')-aph(2'') | AAC6_Ie_APH2_Ia | 3 | confirmed_multi_tool |
| K_pneumoniae_HS11286 | CTX-M-14 | blaCTX-M-14 | blaCTX-M-14 | CTX-M-14 | 3 | confirmed_multi_tool |
| K_pneumoniae_HS11286 | CTX-M-14 | blaCTX-M-14 | blaCTX-M-14 | CTX-M-14 | 3 | confirmed_multi_tool |
| K_pneumoniae_MGH78578 | APH(3'')-Ib | aph(3'')-Ib | aph(3'')-Ib | APH(3'')-Ib | 3 | confirmed_multi_tool |
| K_pneumoniae_MGH78578 | APH(6)-Id | aph(6)-Id | aph(6)-Id | APH(6)-Id | 3 | confirmed_multi_tool |
| K_pneumoniae_MGH78578 | sul2 | sul2 | sul2 | sul2 | 3 | confirmed_multi_tool |
| K_pneumoniae_MGH78578 | APH(3')-Ia | aph(3')-Ia | aph(3')-Ia | APH(3')-Ia | 3 | confirmed_multi_tool |
| K_pneumoniae_MGH78578 | catA1 | catA1 | catA1 | catA1 | 3 | confirmed_multi_tool |
| Enterobacter_cloacae_ATCC13047 | sul2 | sul2 | sul2 | sul2 | 3 | confirmed_multi_tool |
| S_aureus_N315 | ANT(4')-Ia | aadD1 | aadD | ANT(4')-Ia | 2 | confirmed_multi_tool |
| S_aureus_Mu50 | ANT(4')-Ia | aadD1 | aadD | ANT(4')-Ia | 2 | confirmed_multi_tool |
| S_aureus_Mu50 | AAC6_Ie_APH2_Ia | aac(6')-Ie/aph(2'')-Ia | aac(6')-aph(2'') | AAC6_Ie_APH2_Ia | 3 | confirmed_multi_tool |
| A_baumannii_ATCC17978 | sul2 | sul2 | sul2 | sul2 | 3 | confirmed_multi_tool |


## Discordant loci (33)

| strain_id | builtin_gene | amrfinder_hit | resfinder_hit | blast_hit | n_tools_confirm | validation_status |
| --- | --- | --- | --- | --- | --- | --- |
| E_faecium_Aus0004 | mel | mef(A) | mef(A) |  | 0 | discordant |
| E_faecium_Aus0004 | vanX_in_vanB_cl | vanX-B | VanHBX |  | 0 | discordant |
| E_faecium_Aus0004 | vanW_in_vanB_cl | vanW-B | VanHBX |  | 0 | discordant |
| E_faecium_Aus0004 | vanS_in_vanB_cl | vanS-B |  |  | 0 | discordant |
| E_faecalis_V583 | vanX_in_vanB_cl | vanX-B | VanHBX |  | 0 | discordant |
| E_faecalis_V583 | vanH_in_vanB_cl | vanH-B | VanHBX |  | 0 | discordant |
| E_faecalis_V583 | vanW_in_vanB_cl | vanW-B | VanHBX |  | 0 | discordant |
| E_faecalis_V583 | vanY_in_vanB_cl | vanY-B |  |  | 0 | discordant |
| E_faecalis_V583 | vanS_in_vanB_cl | vanS-B |  |  | 0 | discordant |
| E_faecalis_V583 | vanR_in_vanB_cl | vanR-B |  |  | 0 | discordant |
| K_pneumoniae_MGH78578 | SHV-40 | blaSHV-11 | blaSHV-89 |  | 0 | discordant |
| E_coli_K12_MG1655 | Ecol_emrE | emrE |  |  | 0 | discordant |
| E_coli_K12_MG1655 | Ecol_ampC_BLA | blaEC |  |  | 0 | discordant |
| E_coli_O157H7_Sakai | Ecol_emrE | emrE |  |  | 0 | discordant |
| E_coli_O157H7_Sakai | EC-15 | blaEC |  |  | 0 | discordant |
| S_enterica_LT2 | mdsC | mdsB |  |  | 0 | discordant |
| C_difficile_630 | CdnimB | nimB-Cd |  |  | 0 | discordant |
| S_aureus_NCTC8325 | Saur_LmrS | lmrS |  |  | 0 | discordant |
| S_aureus_NCTC8325 | Saur_FosB | fosB |  |  | 0 | discordant |
| S_aureus_N315 | Saur_LmrS | lmrS |  |  | 0 | discordant |
| S_aureus_N315 | Saur_FosB | fosB |  |  | 0 | discordant |
| S_aureus_N315 | PC1_blaZ | blaZ | blaZ |  | 0 | discordant |
| S_aureus_Mu50 | Saur_LmrS | lmrS |  |  | 0 | discordant |
| S_aureus_Mu50 | Saur_FosB | fosB |  |  | 0 | discordant |
| L_monocytogenes_EGDe | lin | vga(G) |  |  | 0 | discordant |
| P_aeruginosa_PAO1 | Paer_catB7 | catB7 | catB7 |  | 0 | discordant |
| P_aeruginosa_PAO1 | MexD |  | tmexD2 |  | 0 | discordant |
| P_aeruginosa_PAO1 | OXA-904 | blaOXA-50 | blaOXA-50 |  | 0 | discordant |
| A_baumannii_ATCC17978 | Abau_AbaF | abaF |  |  | 0 | discordant |
| A_baumannii_ATCC17978 | Abau_AmvA | amvA |  |  | 0 | discordant |
| A_baumannii_ATCC17978 | ADC-66 | blaADC-26 | blaADC-25 |  | 0 | discordant |
| B_fragilis_NCTC9343 | CepA-29 | cepA | cepA |  | 0 | discordant |
| V_cholerae_N16961 | Vcho_varG | varG |  |  | 0 | discordant |


_Allele names may be asserted only for loci confirmed by an independent tool; everything else stays 'best CARD match'._
