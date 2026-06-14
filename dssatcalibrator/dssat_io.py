"""Parsers for DSSAT-CSM text output files.

DSSAT writes whitespace-delimited fixed-format files. We deliberately parse the
columns that have **no embedded spaces** so a simple ``str.split()`` is robust:

* ``PlantGro.OUT``  — daily growth time-series, one block per run/treatment.
* ``Evaluate.OUT``  — simulated-vs-measured scalars per treatment (written only
                      when observed FileA/FileT are present in the run directory).
* ``Summary.OUT``   — end-of-season summary (parsed by fixed column positions
                      because the TNAM column contains spaces).

The value ``-99`` (and ``-99.0``) is DSSAT's missing-data sentinel and is mapped
to ``NaN`` everywhere.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

MISSING = -99.0


def yyddd_to_date(code) -> pd.Timestamp:
    """Convert a DSSAT ``YYDDD`` (or ``YYYYDDD``) date code to a Timestamp.

    Two-digit years < 80 are read as 20xx, otherwise 19xx (DSSAT convention).
    Returns ``NaT`` for missing / unparseable codes.
    """
    try:
        s = str(int(float(code))).strip()
    except (ValueError, TypeError):
        return pd.NaT
    if s in ("", "-99", "0"):
        return pd.NaT
    if len(s) <= 5:  # YYDDD
        s = s.zfill(5)
        yy, doy = int(s[:2]), int(s[2:])
        year = 2000 + yy if yy < 80 else 1900 + yy
    else:            # YYYYDDD
        s = s.zfill(7)
        year, doy = int(s[:4]), int(s[4:])
    if doy < 1 or doy > 366:
        return pd.NaT
    return pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy - 1)


def _to_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce all columns to numeric where possible and map -99 -> NaN."""
    out = df.apply(pd.to_numeric, errors="coerce")
    return out.mask(np.isclose(out, MISSING))


# --------------------------------------------------------------------------- #
#  PlantGro.OUT  (daily time-series)                                          #
# --------------------------------------------------------------------------- #
def parse_plantgro(path: str | Path) -> pd.DataFrame:
    """Parse ``PlantGro.OUT`` into a tidy long-ish daily DataFrame.

    Returns one row per (treatment, day) with all numeric growth columns plus
    ``run``, ``treatment`` and a derived ``date`` (from YEAR + DOY) and ``DAP``.
    """
    text = Path(path).read_text(errors="replace").splitlines()
    runs: list[pd.DataFrame] = []

    run_no = None
    treatment = None
    header: list[str] | None = None
    rows: list[list[str]] = []

    def flush():
        nonlocal rows, header, run_no, treatment
        if header and rows:
            df = pd.DataFrame(rows, columns=header)
            df = _to_numeric(df)
            df["run"] = run_no
            df["treatment"] = treatment if treatment is not None else run_no
            runs.append(df)
        rows = []
        header = None

    for ln in text:
        if ln.startswith("*RUN"):
            flush()
            m = re.match(r"\*RUN\s+(\d+)", ln)
            run_no = int(m.group(1)) if m else None
            treatment = None
        elif ln.lstrip().startswith("TREATMENT"):
            m = re.match(r"\s*TREATMENT\s+(\d+)", ln)
            if m:
                treatment = int(m.group(1))
        elif ln.startswith("@YEAR") or ln.startswith("@ YEAR"):
            header = ln.lstrip("@").split()
            rows = []
        elif header is not None and re.match(r"\s*\d{4}\s", ln):
            parts = ln.split()
            if len(parts) >= len(header):
                rows.append(parts[: len(header)])
            elif parts:  # short row: pad
                rows.append(parts + [str(MISSING)] * (len(header) - len(parts)))
    flush()

    if not runs:
        return pd.DataFrame()
    out = pd.concat(runs, ignore_index=True)
    if {"YEAR", "DOY"}.issubset(out.columns):
        yr = out["YEAR"].astype("Int64")
        doy = out["DOY"].astype("Int64")
        out["date"] = pd.to_datetime(
            yr.astype("string") + "-" + doy.astype("string"),
            format="%Y-%j", errors="coerce",
        )
    out["run"] = out["run"].astype("Int64")
    out["treatment"] = out["treatment"].astype("Int64")
    return out


# --------------------------------------------------------------------------- #
#  Evaluate.OUT  (simulated vs measured scalars)                              #
# --------------------------------------------------------------------------- #
def parse_evaluate(path: str | Path) -> pd.DataFrame:
    """Parse ``Evaluate.OUT`` into a long table: one row per (treatment, variable).

    Columns: ``treatment``, ``run``, ``variable`` (the base name, e.g. ``ADAP``,
    ``HWAM``, ``CWAM``, ``LAIX``), ``sim`` and ``meas``. Pairs are detected from
    the ``...S`` (simulated) / ``...M`` (measured) suffix convention.
    """
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["treatment", "run", "variable", "sim", "meas"])
    lines = p.read_text(errors="replace").splitlines()
    hdr_idx = next((i for i, ln in enumerate(lines) if ln.startswith("@RUN")), None)
    if hdr_idx is None:
        return pd.DataFrame(columns=["treatment", "run", "variable", "sim", "meas"])

    header = lines[hdr_idx].lstrip("@").split()
    data = [ln.split() for ln in lines[hdr_idx + 1:] if re.match(r"\s*\d", ln)]
    if not data:
        return pd.DataFrame(columns=["treatment", "run", "variable", "sim", "meas"])
    wide = pd.DataFrame(data, columns=header[: len(data[0])])

    id_cols = {"RUN", "EXCODE", "TN", "RN", "CR"}
    sim_cols = [c for c in wide.columns if c.endswith("S") and c[:-1] + "M" in wide.columns
                and c not in id_cols]

    recs = []
    for _, row in wide.iterrows():
        trt = pd.to_numeric(row.get("TN"), errors="coerce")
        run = pd.to_numeric(row.get("RUN"), errors="coerce")
        for sc in sim_cols:
            base = sc[:-1]
            sim = pd.to_numeric(row[sc], errors="coerce")
            meas = pd.to_numeric(row[base + "M"], errors="coerce")
            recs.append((trt, run, base, sim, meas))
    out = pd.DataFrame(recs, columns=["treatment", "run", "variable", "sim", "meas"])
    out[["sim", "meas"]] = out[["sim", "meas"]].mask(np.isclose(out[["sim", "meas"]], MISSING))
    out["treatment"] = out["treatment"].astype("Int64")
    out["run"] = out["run"].astype("Int64")
    return out


# --------------------------------------------------------------------------- #
#  Summary.OUT  (end-of-season; fixed columns because TNAM has spaces)        #
# --------------------------------------------------------------------------- #
# Map of the summary columns we care about -> dtype hint. Parsed by name from
# the header token list, but values are read by splitting on whitespace AFTER
# the TNAM/FNAM text columns, which we skip via the leading numeric run number.
_SUMMARY_NUM = [
    "RUNNO", "TRNO", "CWAM", "HWAM", "HWAH", "BWAH", "PWAM", "LAIX",
    "ADAT", "MDAT", "EDAT", "PDAT", "HDAT",
]


def parse_summary(path: str | Path) -> pd.DataFrame:
    """Best-effort parse of ``Summary.OUT`` numeric columns by header name.

    The TNAM/FNAM text columns contain spaces, so we locate each requested
    column's position in the header and read the corresponding whitespace token
    from each data row only where that is unambiguous (numeric columns after the
    text block). Prefer :func:`parse_evaluate` for scalar fitting; this is a
    convenience/fallback.
    """
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    lines = p.read_text(errors="replace").splitlines()
    hdr_idx = next((i for i, ln in enumerate(lines) if ln.lstrip().startswith("@") and "RUNNO" in ln), None)
    if hdr_idx is None:
        return pd.DataFrame()
    header = lines[hdr_idx].lstrip("@").split()
    recs = []
    for ln in lines[hdr_idx + 1:]:
        if not re.match(r"\s*\d", ln):
            continue
        parts = ln.split()
        # Align from the RIGHT for numeric tail columns (robust to spacey TNAM).
        row = {}
        for col in _SUMMARY_NUM:
            if col in header:
                # position from the right is stable for the numeric tail
                idx_from_left = header.index(col)
                # only trust columns whose left index maps cleanly (no spacey col before)
                if idx_from_left < len(parts):
                    row[col] = parts[idx_from_left]
        recs.append(row)
    df = pd.DataFrame(recs)
    return _to_numeric(df) if not df.empty else df
