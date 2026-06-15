"""Phase 4 - Mobility-context validation for High/Medium ARGs.

Uses the pipeline's own context (plasmid replicon, IS flanking) and OPTIONALLY
overlays external predictions (MOB-suite / geNomad / MobileElementFinder) if you
drop their result files in validation/. Because the builtin runs have 0
conjugation markers, nothing is called conjugative unless an external tool
confirms it. Wording is 'mobile-context potential', never 'HGT demonstrated'.

Usage:
    py scripts\\mobility_validation.py
    py scripts\\mobility_validation.py --results results_scale50
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _sota_common import HIGH_MED, df_to_md, read_csv_safe, sota_dir


def _find_col(cols, *keys):
    low = {c.lower(): c for c in cols}
    for k in keys:
        for lc, orig in low.items():
            if k in lc:
                return orig
    return None


def load_mobsuite(path: Path) -> dict:
    """contig -> {'molecule_type','mobility'} from a mob_recon contig report."""
    df = pd.read_csv(path, sep="\t")
    c = _find_col(df.columns, "contig")
    mol = _find_col(df.columns, "molecule_type", "molecule")
    mob = _find_col(df.columns, "mobility")
    res = {}
    if c:
        for _, r in df.iterrows():
            res[str(r[c])] = {"molecule_type": str(r[mol]).lower() if mol else "",
                              "mobility": str(r[mob]).lower() if mob else ""}
    return res


def load_genomad(path: Path) -> dict:
    """Return {'plasmid': set(contigs), 'virus': set(contigs)} from geNomad summary."""
    df = pd.read_csv(path, sep="\t", comment="#")
    c = _find_col(df.columns, "seq_name", "contig", "scaffold")
    cls = _find_col(df.columns, "classification", "class", "type", "topology")
    res = {"plasmid": set(), "virus": set()}
    if not c:
        return res
    for _, r in df.iterrows():
        name = str(r[c])
        label = str(r[cls]).lower() if cls else ""
        if "virus" in label or "provirus" in label or "phage" in label:
            res["virus"].add(name)
        elif "plasmid" in label:
            res["plasmid"].add(name)
        elif not cls:
            # no class column -> assume the file lists plasmid contigs
            res["plasmid"].add(name)
    return res


def load_mef(path: Path) -> set:
    """Set of contigs carrying a mobile element per MobileElementFinder."""
    df = pd.read_csv(path, sep="\t", comment="#") if path.suffix in (".tsv", ".txt") \
        else pd.read_csv(path, comment="#")
    c = _find_col(df.columns, "contig", "seq", "scaffold")
    return set(df[c].astype(str)) if c else set()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_scale50")
    ap.add_argument("--validation-dir", default="validation")
    args = ap.parse_args(argv)
    out = sota_dir()

    arg = read_csv_safe(Path(args.results) / "table2_master_arg.csv")
    feat = read_csv_safe(Path(args.results) / "mge_features.csv")
    if arg is None:
        raise SystemExit(f"table2_master_arg.csv not found in {args.results}")

    vdir = Path(args.validation_dir)

    def has_data(p):
        if not p.exists():
            return False
        with open(p, encoding="utf-8", errors="replace") as fh:
            rows = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
        return len(rows) >= 2

    mobsuite = load_mobsuite(vdir / "mobsuite_results.tsv") if has_data(vdir / "mobsuite_results.tsv") else {}
    genomad = load_genomad(vdir / "genomad_results.tsv") if has_data(vdir / "genomad_results.tsv") else {"plasmid": set(), "virus": set()}
    mef = load_mef(vdir / "mobileelementfinder_results.tsv") if has_data(vdir / "mobileelementfinder_results.tsv") else set()
    ext_note = []
    if mobsuite:
        ext_note.append("MOB-suite")
    if genomad["plasmid"] or genomad["virus"]:
        ext_note.append("geNomad")
    if mef:
        ext_note.append("MobileElementFinder")

    # replicon contigs from our own MGE features
    replicon_contigs = set()
    if feat is not None:
        replicon_contigs = set(feat[feat.feature_type == "plasmid_replicon"].contig.astype(str))

    hm = arg[arg.risk_level.isin(HIGH_MED)].copy()
    rows = []
    for _, r in hm.iterrows():
        contig = str(r.contig)
        on_plasmid = bool(r.get("on_plasmid", False))
        repl_near = on_plasmid or (contig in replicon_contigs)
        is_l = bool(r.get("is_upstream", False))
        is_r = bool(r.get("is_downstream", False))
        is_both = bool(r.get("is_flank_both", False))

        # contig type (external tools take precedence)
        if contig in mobsuite and mobsuite[contig]["molecule_type"]:
            ctype = mobsuite[contig]["molecule_type"]  # 'plasmid'/'chromosome'
        elif contig in genomad["plasmid"]:
            ctype = "plasmid"
        elif repl_near:
            ctype = "plasmid"
        else:
            ctype = "chromosome"

        # conjugative ONLY from external evidence
        conj = "unknown"
        if contig in mobsuite and mobsuite[contig]["mobility"]:
            mob = mobsuite[contig]["mobility"]
            conj = "yes" if "conjugative" in mob else ("no" if "non" in mob else "unknown")

        phage = "yes" if contig in genomad["virus"] else "unknown"

        # external mobility/MGE confirmation (a tool, not just our context)
        mob_tool_confirmed = (
            (contig in mobsuite and mobsuite[contig]["molecule_type"] == "plasmid")
            or contig in genomad["plasmid"] or contig in mef)

        # mobilizable potential
        ext_mob = (contig in mobsuite and mobsuite[contig].get("mobility"))
        if ext_mob and "conjugative" in mobsuite[contig]["mobility"]:
            mobiliz = "yes"  # external tool asserts conjugative/mobilizable
        elif mob_tool_confirmed:
            mobiliz = "yes"  # external tool confirms plasmid/MGE context
        elif repl_near or is_both:
            mobiliz = "context-positive"  # in silico context only, no external proof
        elif ctype == "chromosome" and not (is_l or is_r):
            mobiliz = "no"
        else:
            mobiliz = "unknown"

        bits = []
        if ctype == "plasmid":
            bits.append("plasmid-located")
        if is_both:
            bits.append("IS-flanked both sides")
        elif is_l or is_r:
            bits.append("single flanking IS")
        if conj == "yes":
            bits.append("externally predicted conjugative")
        elif conj == "unknown":
            bits.append("conjugative status unknown (no relaxase evidence)")
        if bits:
            interp = ", ".join(bits) + " - in silico mobile-context potential"
        else:
            interp = "no mobile-context signal"

        rows.append({
            "strain_id": r.strain_id, "gene": r.gene, "drug_class": r.get("drug_class"),
            "contig": contig, "risk_level": r.risk_level,
            "contig_type": ctype, "plasmid_replicon_nearby": "yes" if repl_near else "no",
            "IS_left": "yes" if is_l else "no", "IS_right": "yes" if is_r else "no",
            "IS_both_sides": "yes" if is_both else "no",
            "predicted_mobilizable": mobiliz, "predicted_conjugative": conj,
            "phage_or_viral_context": phage,
            "mobility_tool_confirmed": "yes" if mob_tool_confirmed else "no",
            "final_mobility_interpretation": interp,
        })

    mob = pd.DataFrame(rows)
    mob.to_csv(out / "mobility_validation_SCALE50.csv", index=False)

    md = ["# Mobility-context validation (SCALE50 High/Medium ARGs)\n",
          f"External overlays applied: **{', '.join(ext_note) if ext_note else 'none - external validation pending'}**.\n",
          "All calls are *in silico mobile-context potential*, not demonstrated HGT. "
          "Conjugative status is 'unknown' unless an external tool (e.g. MOB-suite) "
          "confirms a relaxase/MPF system - the builtin runs detected 0 conjugation "
          "markers.\n",
          "## Summary\n"]
    summ = pd.DataFrame([{
        "High/Medium ARGs": len(mob),
        "on plasmid": int((mob.contig_type == "plasmid").sum()),
        "IS both sides": int((mob.IS_both_sides == "yes").sum()),
        "mobility-context positive (in silico)": int(mob.predicted_mobilizable.isin(["yes", "context-positive"]).sum()),
        "mobility-tool confirmed (plasmid/MGE)": int((mob.mobility_tool_confirmed == "yes").sum()),
        "externally predicted conjugative": int((mob.predicted_conjugative == "yes").sum()),
    }])
    md.append(df_to_md(summ))

    # machine-readable tooling tally for the report + SOTA gate
    pd.DataFrame([{
        "mobsuite_present": bool(mobsuite),
        "genomad_present": bool(genomad["plasmid"] or genomad["virus"]),
        "mef_present": bool(mef),
        "mobility_tool_confirmed_loci": int((mob.mobility_tool_confirmed == "yes").sum()),
        "externally_conjugative_loci": int((mob.predicted_conjugative == "yes").sum()),
    }]).to_csv(out / "mobility_tooling.csv", index=False)
    md.append("\n## Per-locus interpretation\n")
    cols = ["strain_id", "gene", "contig_type", "plasmid_replicon_nearby",
            "IS_both_sides", "predicted_mobilizable", "predicted_conjugative",
            "final_mobility_interpretation"]
    md.append(df_to_md(mob[cols], max_rows=80))
    (out / "mobility_summary.md").write_text("\n".join(md), encoding="utf-8")

    print(f"mobility-validated {len(mob)} High/Medium loci -> "
          f"{out/'mobility_validation_SCALE50.csv'}")
    print("external overlays:", ext_note or "none")


if __name__ == "__main__":
    main()
