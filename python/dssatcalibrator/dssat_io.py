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

ID_COLUMNS = {
    "RUN", "RUNNO", "TRNO", "TRT", "TRTNO", "TN", "RN", "SN", "ON", "REP",
    "YEAR", "DOY", "DATE", "DAP", "DAS", "EXCODE", "CR", "MODEL", "WSTA",
    "FILEX", "FNAM", "TNAM", "CNAME", "ENAME", "CODE", "ID_FIELD", "ID_SOIL",
}


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


def _coerce_numeric_values(df: pd.DataFrame, skip: set[str] | None = None) -> pd.DataFrame:
    """Coerce numeric-looking object columns while preserving textual id columns."""
    skip = skip or set()
    out = df.copy()
    for col in out.columns:
        if col in skip:
            continue
        nums = pd.to_numeric(out[col], errors="coerce")
        if nums.notna().any():
            out[col] = nums.mask(np.isclose(nums, MISSING))
    return out


def _looks_like_data_row(line: str) -> bool:
    s = line.strip()
    if not s or s[0] in "*!@#$":
        return False
    return True


def _derive_date_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if {"YEAR", "DOY"}.issubset(out.columns):
        yr = pd.to_numeric(out["YEAR"], errors="coerce").astype("Int64")
        doy = pd.to_numeric(out["DOY"], errors="coerce").astype("Int64")
        out["date"] = pd.to_datetime(
            yr.astype("string") + "-" + doy.astype("string"),
            format="%Y-%j", errors="coerce",
        )
    elif "DATE" in out.columns:
        out["date"] = out["DATE"].map(yyddd_to_date)
    return out


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


# --------------------------------------------------------------------------- #
#  Generic DSSAT OUT table collector                                          #
# --------------------------------------------------------------------------- #
def parse_dssat_output_tables(path: str | Path) -> pd.DataFrame:
    """Parse table-like ``*.OUT`` sections that use DSSAT ``@`` headers.

    This is intentionally broad and best-effort. It preserves one row per output
    table row, adds ``source_file`` and ``section`` metadata, and coerces numeric
    columns where possible. Dedicated parsers above remain the source of truth for
    scoring-specific files such as ``PlantGro.OUT`` and ``Evaluate.OUT``.
    """
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    lines = p.read_text(errors="replace").splitlines()
    frames: list[pd.DataFrame] = []
    header: list[str] | None = None
    rows: list[list[str]] = []
    section = -1

    def flush():
        nonlocal rows, header
        if header and rows:
            ncol = len(header)
            norm = [r[:ncol] + [str(MISSING)] * max(0, ncol - len(r)) for r in rows]
            df = pd.DataFrame(norm, columns=header)
            df.insert(0, "row_index", range(len(df)))
            df.insert(0, "section", section)
            df.insert(0, "source_file", p.name)
            frames.append(_derive_date_columns(_coerce_numeric_values(df, skip={"source_file"})))
        rows = []
        header = None

    for ln in lines:
        stripped = ln.strip()
        if not stripped:
            continue
        if stripped.startswith("@"):
            flush()
            section += 1
            header = stripped.lstrip("@").split()
            rows = []
            continue
        if header is not None and _looks_like_data_row(ln):
            parts = ln.split()
            if parts:
                rows.append(parts)
    flush()
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def dssat_output_long(wide: pd.DataFrame) -> pd.DataFrame:
    """Convert a generic wide DSSAT output table into ``variable``/``value`` rows."""
    if wide.empty:
        return pd.DataFrame(columns=[
            "source_file", "section", "row_index", "run", "treatment", "date",
            "dap", "das", "variable", "value",
        ])
    id_map = {
        "RUN": "run", "RUNNO": "run",
        "TRNO": "treatment", "TRT": "treatment", "TRTNO": "treatment", "TN": "treatment",
        "DAP": "dap", "DAS": "das",
    }
    meta = pd.DataFrame(index=wide.index)
    for src, dst in id_map.items():
        if src in wide.columns and dst not in meta:
            meta[dst] = wide[src]
    for col in wide.columns:
        if col in ID_COLUMNS and col not in id_map:
            dst = col.lower()
            if dst not in meta:
                meta[dst] = wide[col]
    for col in ("source_file", "section", "row_index", "date"):
        if col in wide.columns:
            meta[col] = wide[col]

    value_cols = []
    for col in wide.columns:
        if col in ID_COLUMNS or col in {"source_file", "section", "row_index", "date"}:
            continue
        if pd.api.types.is_numeric_dtype(wide[col]):
            value_cols.append(col)
    if not value_cols:
        return pd.DataFrame(columns=[
            "source_file", "section", "row_index", "run", "treatment", "date",
            "dap", "das", "variable", "value",
        ])
    long = wide[value_cols].copy()
    long = long.where(np.isfinite(long), np.nan)
    long = long.join(meta)
    preferred = ["source_file", "section", "row_index", "run", "treatment", "date", "dap", "das"]
    id_vars = [c for c in preferred if c in long.columns]
    id_vars.extend([c for c in meta.columns if c not in id_vars and c in long.columns])
    out = long.melt(id_vars=id_vars, value_vars=value_cols, var_name="variable", value_name="value")
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    out = out.dropna(subset=["value"]).reset_index(drop=True)
    return out


DEFAULT_OUTPUT_FILES = [
    "PlantGro.OUT", "Evaluate.OUT", "Summary.OUT", "OVERVIEW.OUT", "Weather.OUT",
    "ET.OUT", "SoilWat.OUT", "SoilWater.OUT", "SoilNi.OUT", "SoilTemp.OUT",
    "PlantN.OUT", "PlantNBal.OUT", "PlantC.OUT", "SoilNBalSum.OUT", "SoilWatBal.OUT",
]


def collect_run_outputs(run_dir: str | Path, output_files: list[str] | None = None) -> dict[str, pd.DataFrame]:
    """Collect broad DSSAT outputs from one run directory.

    Returns ``{"wide": ..., "long": ..., "plantgro": ..., "evaluate": ...,
    "summary": ..., "manifest": ...}``. Missing files are represented in the
    manifest and otherwise skipped.
    """
    run_dir = Path(run_dir)
    files = output_files or DEFAULT_OUTPUT_FILES
    manifest = []
    wide_frames = []
    long_frames = []

    for name in files:
        p = run_dir / name
        manifest.append({
            "source_file": name,
            "exists": p.exists(),
            "size_bytes": int(p.stat().st_size) if p.exists() else 0,
        })
        if not p.exists() or p.stat().st_size == 0:
            continue
        generic = parse_dssat_output_tables(p)
        if not generic.empty:
            wide_frames.append(generic)
            long = dssat_output_long(generic)
            if not long.empty:
                long_frames.append(long)

    plantgro = parse_plantgro(run_dir / "PlantGro.OUT") if (run_dir / "PlantGro.OUT").exists() else pd.DataFrame()
    evaluate = parse_evaluate(run_dir / "Evaluate.OUT")
    summary = parse_summary(run_dir / "Summary.OUT") if (run_dir / "Summary.OUT").exists() else pd.DataFrame()

    return {
        "wide": pd.concat(wide_frames, ignore_index=True, sort=False) if wide_frames else pd.DataFrame(),
        "long": pd.concat(long_frames, ignore_index=True, sort=False) if long_frames else pd.DataFrame(),
        "plantgro": plantgro,
        "evaluate": evaluate,
        "summary": summary,
        "manifest": pd.DataFrame(manifest),
    }
