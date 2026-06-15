# Full-assembly vs chromosome-only comparison (matched SCALE50)

_Same 50 strains; the only difference is `full_assembly` true vs false._

## Metric comparison

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


Full assembly recovers the plasmid compartment: plasmid replicons 0 -> 32, High-risk ARGs 4 -> 15, total ARG loci 338 -> 373.


## ARG loci gained in full assembly (33)

| strain_id | gene |
| --- | --- |
| E_faecalis_V583 | AAC6_Ie_APH2_Ia |
| E_faecalis_V583 | ErmB |
| E_faecalis_V583 | qacH |
| K_pneumoniae_HS11286 | APH(3'')-Ib |
| K_pneumoniae_HS11286 | APH(6)-Id |
| K_pneumoniae_HS11286 | CTX-M-14 |
| K_pneumoniae_HS11286 | KPC-2 |
| K_pneumoniae_HS11286 | TEM-1 |
| K_pneumoniae_HS11286 | aadA2 |
| K_pneumoniae_HS11286 | cmlA9 |
| K_pneumoniae_HS11286 | qacEdelta1 |
| K_pneumoniae_HS11286 | rmtB |
| K_pneumoniae_HS11286 | sul2 |
| K_pneumoniae_MGH78578 | AAC(6')-Ib' |
| K_pneumoniae_MGH78578 | ANT(2'')-Ia |
| K_pneumoniae_MGH78578 | APH(3'')-Ib |
| K_pneumoniae_MGH78578 | APH(3')-Ia |
| K_pneumoniae_MGH78578 | APH(6)-Id |
| K_pneumoniae_MGH78578 | OXA-9 |
| K_pneumoniae_MGH78578 | SHV-12 |
| K_pneumoniae_MGH78578 | SHV-40 |
| K_pneumoniae_MGH78578 | TEM-166 |
| K_pneumoniae_MGH78578 | TEM-90 |
| K_pneumoniae_MGH78578 | catA1 |
| K_pneumoniae_MGH78578 | cmlA5 |
| K_pneumoniae_MGH78578 | qacEdelta1 |
| K_pneumoniae_MGH78578 | sul1 |
| K_pneumoniae_MGH78578 | sul2 |
| K_pneumoniae_MGH78578 | tet(D) |
| S_aureus_Mu50 | AAC6_Ie_APH2_Ia |
| S_aureus_Mu50 | qacA |
| S_aureus_N315 | PC1_blaZ |
| V_cholerae_N16961 | catB9 |


## ARG loci only seen in chromosome-only run (1)

| strain_id | gene |
| --- | --- |
| K_pneumoniae_MGH78578 | SHV-79 |


## Plasmid-associated ARGs gained in full assembly (28)

| strain_id | gene | drug_class | contig | risk_level |
| --- | --- | --- | --- | --- |
| K_pneumoniae_HS11286 | CTX-M-14 | beta-lactam | NC_016838.1 | High |
| K_pneumoniae_HS11286 | TEM-1 | beta-lactam | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | rmtB | aminoglycoside | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | cmlA9 | phenicol | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | qacEdelta1 | unclassified | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | aadA2 | aminoglycoside | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | TEM-1 | beta-lactam | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | APH(6)-Id | aminoglycoside | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | APH(3'')-Ib | aminoglycoside | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | sul2 | sulfonamide | NC_016839.1 | Medium |
| K_pneumoniae_HS11286 | CTX-M-14 | beta-lactam | NC_016839.1 | High |
| K_pneumoniae_HS11286 | KPC-2 | beta-lactam | NC_016846.1 | Medium |
| K_pneumoniae_HS11286 | TEM-1 | beta-lactam | NC_016846.1 | Medium |
| K_pneumoniae_MGH78578 | SHV-12 | beta-lactam | NC_009650.1 | Medium |
| K_pneumoniae_MGH78578 | AAC(6')-Ib' | aminoglycoside | NC_009650.1 | Medium |
| K_pneumoniae_MGH78578 | OXA-9 | beta-lactam | NC_009650.1 | Medium |
| K_pneumoniae_MGH78578 | TEM-90 | beta-lactam | NC_009650.1 | Medium |
| K_pneumoniae_MGH78578 | APH(3'')-Ib | aminoglycoside | NC_009651.1 | High |
| K_pneumoniae_MGH78578 | APH(6)-Id | aminoglycoside | NC_009651.1 | High |
| K_pneumoniae_MGH78578 | sul2 | sulfonamide | NC_009651.1 | High |
| K_pneumoniae_MGH78578 | APH(3')-Ia | aminoglycoside | NC_009651.1 | High |
| K_pneumoniae_MGH78578 | sul1 | sulfonamide | NC_009651.1 | Medium |
| K_pneumoniae_MGH78578 | qacEdelta1 | unclassified | NC_009651.1 | Medium |
| K_pneumoniae_MGH78578 | cmlA5 | phenicol | NC_009651.1 | Medium |
| K_pneumoniae_MGH78578 | ANT(2'')-Ia | aminoglycoside | NC_009651.1 | Medium |
| K_pneumoniae_MGH78578 | catA1 | phenicol | NC_009651.1 | High |
| K_pneumoniae_MGH78578 | tet(D) | tetracycline | NC_009651.1 | Medium |
| K_pneumoniae_MGH78578 | TEM-166 | beta-lactam | NC_009651.1 | Medium |
