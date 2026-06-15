import pandas as pd

builtin = pd.read_csv("results_scale50/table2_master_arg.csv")
amr = pd.read_csv("validation/amrfinderplus_results.tsv", sep="\t")

b = builtin[
    (builtin["strain_id"] == "A_baumannii_ATCC17978") &
    (builtin["gene"].astype(str).str.contains("sul2", case=False, na=False))
]

a = amr[
    (amr["strain_id"] == "A_baumannii_ATCC17978") &
    (amr["Element symbol"].astype(str).str.contains("sul2", case=False, na=False))
]

print("BUILTIN")
print(b[["strain_id","gene","contig","start","end","pct_identity","pct_coverage"]].to_string(index=False))

print("\nAMRFINDER")
print(a[["strain_id","Element symbol","Contig id","Start","Stop","% Identity to reference","% Coverage of reference"]].to_string(index=False))
