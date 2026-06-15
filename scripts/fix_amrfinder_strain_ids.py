import pandas as pd
from pathlib import Path

meta_path = Path("results_scale50/table1_strain_metadata.csv")
amr_path = Path("validation/amrfinderplus_results.tsv")
out_path = Path("validation/amrfinderplus_results.fixed.tsv")

meta = pd.read_csv(meta_path)
amr = pd.read_csv(amr_path, sep="\t")

# Build accession -> strain_id map
acc_map = {}
for _, r in meta.iterrows():
    strain = str(r.get("strain_id", "")).strip()
    acc = str(r.get("accession", "")).strip()
    if acc:
        acc_map[acc] = strain

def map_strain(row):
    # Prefer Contig id because AMRFinder output has accession there
    contig = str(row.get("Contig id", "")).strip()
    raw = str(row.get("strain_id", "")).strip()
    
    if contig in acc_map:
        return acc_map[contig]
    
    # Fallback: accession appears inside path-like strain_id
    for acc, strain in acc_map.items():
        if acc and acc in raw:
            return strain
    
    return raw

amr["original_strain_id"] = amr["strain_id"]
amr["strain_id"] = amr.apply(map_strain, axis=1)

amr.to_csv(out_path, sep="\t", index=False)

print("Wrote:", out_path)
print("Rows:", len(amr))
print("Unique original strain_id:", amr["original_strain_id"].nunique())
print("Unique fixed strain_id:", amr["strain_id"].nunique())
print(amr[["original_strain_id", "Contig id", "strain_id", "Element symbol"]].head(20).to_string(index=False))
