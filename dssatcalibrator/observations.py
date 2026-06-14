"""Read observed data into one tidy long-format ``Observations`` table.

Three sources, normalised to the same schema (one row per measurement):

    exp_id | treatment | variable | kind | date | value | sigma | weight

* ``kind`` is ``timeseries`` (a value on a date), ``scalar`` (end-of-season,
  ``date`` = NaT) or ``phenology`` (an event whose ``value`` is the event date,
  also surfaced in ``date``).
* DSSAT **FileA** (``*.??A``) -> end-of-season scalars + phenology dates.
* DSSAT **FileT** (``*.??T``) -> in-season time-series (replicates kept as rows).
* A user **CSV** in the §5 long format -> passed through with light normalisation.

``-99`` is treated as missing throughout.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .dssat_io import MISSING, yyddd_to_date

SCHEMA = ["exp_id", "treatment", "variable", "kind", "date", "value", "sigma", "weight"]

# Columns whose values are DSSAT date codes (YYDDD) rather than magnitudes.
_DATE_COLS = {
    "EDAT", "ADAT", "MDAT", "IDAT", "DRAT", "GDAT", "PD1T", "PDFT",
    "R1", "R3", "R5", "R7", "R8", "TSAT", "HDAT",
}


def _read_abt_blocks(path: str | Path) -> list[pd.DataFrame]:
    """Read a DSSAT A/T file into a list of wide DataFrames (one per @ header)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    blocks: list[pd.DataFrame] = []
    header: list[str] | None = None
    rows: list[list[str]] = []

    def flush():
        nonlocal header, rows
        if header and rows:
            ncol = len(header)
            norm = [r[:ncol] + [str(MISSING)] * (ncol - len(r)) for r in rows]
            blocks.append(pd.DataFrame(norm, columns=header))
        header, rows = None, []

    for ln in lines:
        s = ln.rstrip()
        if not s or s.startswith("!"):
            continue
        if s.startswith("*"):
            flush()
            continue
        if s.lstrip().startswith("@"):
            flush()
            header = s.lstrip().lstrip("@").split()
            # the first header token is the treatment id ("TRNO" / "TRT")
            if header and header[0] in ("TRNO", "TRT", "TR"):
                header[0] = "TRNO"
            rows = []
        elif header is not None and re.match(r"\s*[0-9]", s):
            rows.append(s.split())
    flush()
    return blocks


def read_filea(path: str | Path, exp_id: str | None = None) -> pd.DataFrame:
    """Read a DSSAT FileA (end-of-season averages) into the long schema."""
    exp_id = exp_id or Path(path).stem
    out = []
    for wide in _read_abt_blocks(path):
        if "TRNO" not in wide.columns:
            continue
        wide = wide.apply(pd.to_numeric, errors="coerce")
        for var in [c for c in wide.columns if c != "TRNO"]:
            for _, r in wide.iterrows():
                val = r[var]
                if pd.isna(val) or np.isclose(val, MISSING):
                    continue
                if var in _DATE_COLS:
                    d = yyddd_to_date(val)
                    out.append((exp_id, int(r["TRNO"]), var, "phenology", d, float(val), np.nan, 1.0))
                else:
                    out.append((exp_id, int(r["TRNO"]), var, "scalar", pd.NaT, float(val), np.nan, 1.0))
    return pd.DataFrame(out, columns=SCHEMA)


def read_filet(path: str | Path, exp_id: str | None = None) -> pd.DataFrame:
    """Read a DSSAT FileT (in-season time-series; replicate rows preserved)."""
    exp_id = exp_id or Path(path).stem
    out = []
    for wide in _read_abt_blocks(path):
        if "TRNO" not in wide.columns or "DATE" not in wide.columns:
            continue
        nums = wide.apply(pd.to_numeric, errors="coerce")
        value_cols = [c for c in wide.columns if c not in ("TRNO", "DATE")]
        for _, r in nums.iterrows():
            d = yyddd_to_date(r["DATE"])
            if pd.isna(d):
                continue
            for var in value_cols:
                val = r[var]
                if pd.isna(val) or np.isclose(val, MISSING):
                    continue
                out.append((exp_id, int(r["TRNO"]), var, "timeseries", d, float(val), np.nan, 1.0))
    return pd.DataFrame(out, columns=SCHEMA)


def read_csv(path: str | Path) -> pd.DataFrame:
    """Read a user long-format observations CSV and normalise to the schema."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    ren = {"experiment": "exp_id", "exp": "exp_id", "trt": "treatment",
           "var": "variable", "obs": "value", "val": "value"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    if "kind" not in df.columns:
        df["kind"] = np.where(df.get("date").notna() if "date" in df.columns else False,
                              "timeseries", "scalar")
    for col, default in (("sigma", np.nan), ("weight", 1.0), ("date", pd.NaT)):
        if col not in df.columns:
            df[col] = default
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.replace(MISSING, np.nan)
    keep = [c for c in SCHEMA if c in df.columns]
    return df[keep].copy()


@dataclass
class Observations:
    """A bundle of observations across one or more experiments."""

    table: pd.DataFrame

    @classmethod
    def from_dssat(cls, hemp_dir: str | Path, experiments: list[str], crop_ext: str = "HM") -> "Observations":
        """Load FileA/FileT observations for a list of DSSAT experiment codes.

        Missing files are skipped silently; an experiment with neither A nor T
        contributes no rows (and is reported by :meth:`coverage`).
        """
        hemp_dir = Path(hemp_dir)
        frames = []
        for exp in experiments:
            fa = hemp_dir / f"{exp}.{crop_ext}A"
            ft = hemp_dir / f"{exp}.{crop_ext}T"
            if fa.exists():
                frames.append(read_filea(fa, exp))
            if ft.exists():
                frames.append(read_filet(ft, exp))
        table = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=SCHEMA)
        return cls(table)

    @classmethod
    def from_csv(cls, path: str | Path) -> "Observations":
        return cls(read_csv(path))

    def coverage(self) -> pd.DataFrame:
        """Per-experiment × variable count of observations (for sanity checks)."""
        if self.table.empty:
            return pd.DataFrame()
        return (self.table.groupby(["exp_id", "kind", "variable"])
                .size().rename("n").reset_index())

    def experiments(self) -> list[str]:
        return sorted(self.table["exp_id"].unique().tolist()) if not self.table.empty else []

    def variables(self) -> list[str]:
        return sorted(self.table["variable"].unique().tolist()) if not self.table.empty else []
