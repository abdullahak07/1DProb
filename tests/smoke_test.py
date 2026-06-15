"""End-to-end smoke test.

Generates synthetic genomes, runs the full pipeline with the dependency-free
`builtin` backend, and asserts that every phase produced the expected result.

Run from the project root:
    py tests/smoke_test.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# ---------------------------------------------------------------------
# Make imports work on Windows, PowerShell, PyCharm, and direct script runs
# ---------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
TESTS_DIR = ROOT / "tests"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(TESTS_DIR))

from src.pipeline import run_pipeline  # noqa: E402
from make_synthetic_data import generate  # noqa: E402


PASS, FAIL = "PASS", "FAIL"
results: list[tuple[str, str, str]] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    results.append((PASS if cond else FAIL, name, detail))


def arg_rows(arg_df, strain, gene=None):
    m = arg_df.strain_id == strain
    if gene is not None:
        m &= arg_df.gene == gene
    return arg_df[m]


def risk_of(arg_df, strain, gene):
    r = arg_rows(arg_df, strain, gene)
    return None if r.empty else str(r.iloc[0].risk_level)


def has_edge(edges, s1, s2, gene):
    if edges.empty:
        return False

    pair = {s1, s2}

    m = edges.apply(
        lambda e: {e.strain_a, e.strain_b} == pair and e.gene == gene,
        axis=1,
    )

    return bool(m.any())


def main() -> int:
    base = ROOT / "data" / "synthetic"
    info = generate(base)
    out_dir = ROOT / "results"

    cfg = {
        "project_name": "smoke",
        "output_dir": str(out_dir),
        "fetch": {
            "enabled": False,
            "genome_dir": str(base / "genomes"),
        },
        "catalogues": info["catalogues"],
        "screen": {
            "backend": "builtin",
            "min_identity": 80,
            "min_coverage": 60,
        },
        "mge": {
            "flank_bp": 5000,
        },
        "compare": {
            "identity_threshold": 95,
        },
        "strains": info["strains"],
    }

    res = run_pipeline(cfg)

    meta = res["metadata"]
    arg = res["arg"]
    edges = res["edges"]
    findings = res["findings"]

    # ---- Phase 2 -------------------------------------------------------- #
    check("metadata has 8 strains", len(meta) == 8, f"got {len(meta)}")
    check("all genomes available", bool(meta["genome_available"].all()))
    check("metadata has quality column", "assembly_quality" in meta.columns)

    # ---- Phase 3: ARG detections --------------------------------------- #
    check("LA1 detects tetM", not arg_rows(arg, "Lacto_acidophilus_LA1", "tetM").empty)
    check("LGG detects vanY", not arg_rows(arg, "Lacto_rhamnosus_GG", "vanY").empty)
    check("BL1 detects ermB", not arg_rows(arg, "Bifido_longum_BL1", "ermB").empty)
    check("EF_PRO detects vanA", not arg_rows(arg, "Entero_faecium_PRO", "vanA").empty)
    check("EF_CLIN detects vanA", not arg_rows(arg, "Entero_faecium_CLIN", "vanA").empty)
    check("KP detects blaTEM", not arg_rows(arg, "Kleb_pneumoniae_KP", "blaTEM").empty)
    check("Ecoli detects tetM", not arg_rows(arg, "Ecoli_MDR", "tetM").empty)
    check("yeast control has NO ARGs", arg_rows(arg, "Saccharomyces_boul_SB").empty)

    # ---- Phase 4: risk tiers ------------------------------------------ #
    check(
        "LA1 tetM = Low",
        risk_of(arg, "Lacto_acidophilus_LA1", "tetM") == "Low",
        f"got {risk_of(arg, 'Lacto_acidophilus_LA1', 'tetM')}",
    )

    check(
        "LGG vanY = Negligible",
        risk_of(arg, "Lacto_rhamnosus_GG", "vanY") == "Negligible",
        f"got {risk_of(arg, 'Lacto_rhamnosus_GG', 'vanY')}",
    )

    check(
        "BL1 ermB = High",
        risk_of(arg, "Bifido_longum_BL1", "ermB") == "High",
        f"got {risk_of(arg, 'Bifido_longum_BL1', 'ermB')}",
    )

    check(
        "EF_PRO vanA = High",
        risk_of(arg, "Entero_faecium_PRO", "vanA") == "High",
        f"got {risk_of(arg, 'Entero_faecium_PRO', 'vanA')}",
    )

    check(
        "KP blaTEM = Medium",
        risk_of(arg, "Kleb_pneumoniae_KP", "blaTEM") == "Medium",
        f"got {risk_of(arg, 'Kleb_pneumoniae_KP', 'blaTEM')}",
    )

    check(
        "Ecoli tetM = Medium",
        risk_of(arg, "Ecoli_MDR", "tetM") == "Medium",
        f"got {risk_of(arg, 'Ecoli_MDR', 'tetM')}",
    )

    # ---- Mobility flags ------------------------------------------------ #
    bl1 = arg_rows(arg, "Bifido_longum_BL1", "ermB").iloc[0]
    check("BL1 ermB flagged IS-both", bool(bl1.is_flank_both))

    efp = arg_rows(arg, "Entero_faecium_PRO", "vanA").iloc[0]
    check(
        "EF_PRO vanA flagged on conjugative plasmid",
        bool(efp.on_plasmid) and bool(efp.conjugative),
    )

    # ---- Phase 5: gene sharing ---------------------------------------- #
    check(
        "vanA shared EF_PRO<->EF_CLIN",
        has_edge(edges, "Entero_faecium_PRO", "Entero_faecium_CLIN", "vanA"),
    )

    check(
        "tetM shared LA1<->Ecoli",
        has_edge(edges, "Lacto_acidophilus_LA1", "Ecoli_MDR", "tetM"),
    )

    check(
        "headline finding is high-concern cross-group vanA",
        (not findings.empty)
        and bool(
            (
                (findings.gene == "vanA")
                & findings.high_concern
                & findings.cross_group
            ).any()
        ),
    )

    # ---- Phase 6: figures --------------------------------------------- #
    for fp in res["figures"]:
        check(
            f"figure exists: {Path(fp).name}",
            Path(fp).exists() and Path(fp).stat().st_size > 0,
        )

    # ---- Output tables ------------------------------------------------- #
    expected_tables = [
        "table1_strain_metadata.csv",
        "table2_master_arg.csv",
        "table3_high_risk_findings.csv",
        "gene_sharing_edges.csv",
        "risk_summary_by_species.csv",
        "mge_features.csv",
    ]

    for tbl in expected_tables:
        check(f"table written: {tbl}", (out_dir / tbl).exists())

    # ---- Report -------------------------------------------------------- #
    print("\n" + "=" * 66)
    print("SMOKE TEST RESULTS")
    print("=" * 66)

    n_fail = 0

    for status, name, detail in results:
        line = f"[{status}] {name}"

        if status == FAIL:
            n_fail += 1
            if detail:
                line += f"   <-- {detail}"

        print(line)

    print("=" * 66)

    total = len(results)
    print(f"{total - n_fail}/{total} checks passed")

    print("=" * 66)

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())