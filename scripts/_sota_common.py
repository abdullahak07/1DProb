"""Shared helpers for the SOTA v2 scripts (stdlib + pandas only)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SOTA_DIR = ROOT / "results_sota"


def sota_dir() -> Path:
    SOTA_DIR.mkdir(parents=True, exist_ok=True)
    return SOTA_DIR


def read_csv_safe(path) -> pd.DataFrame | None:
    p = Path(path)
    if not p.exists():
        return None
    try:
        return pd.read_csv(p)
    except Exception:  # noqa: BLE001
        return None


def df_to_md(df: pd.DataFrame | None, max_rows: int = 0, floatfmt: int = 2) -> str:
    """GitHub-flavoured Markdown table without the `tabulate` dependency."""
    if df is None or len(df) == 0:
        return "_none_\n"
    if max_rows and len(df) > max_rows:
        df = df.head(max_rows)
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |",
             "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        cells = []
        for c in df.columns:
            v = row[c]
            if isinstance(v, float):
                if math.isinf(v) or math.isnan(v):
                    v = "inf" if (isinstance(v, float) and math.isinf(v) and v > 0) else \
                        ("-inf" if isinstance(v, float) and math.isinf(v) else "n/a")
                else:
                    v = (f"{v:.{floatfmt}f}" if abs(v - int(v)) > 1e-9 else str(int(v)))
            cells.append(str(v).replace("|", "\\|"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def overlap(a_start, a_end, b_start, b_end) -> int:
    try:
        vals = [float(a_start), float(a_end), float(b_start), float(b_end)]
    except (TypeError, ValueError):
        return 0
    if any(math.isnan(v) for v in vals):
        return 0
    return max(0, min(int(vals[1]), int(vals[3])) - max(int(vals[0]), int(vals[2])) + 1)


HIGH_MED = ["High", "Medium"]


def locus_id(strain_id, contig, start, end, gene) -> str:
    """Stable identifier for one ARG locus; used by export + validation matching."""
    try:
        s, e = int(float(start)), int(float(end))
    except (TypeError, ValueError):
        s, e = start, end
    return f"{strain_id}|{contig}|{s}-{e}|{gene}"
