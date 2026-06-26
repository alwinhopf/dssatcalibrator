"""Pluggable weather drivers for the simulation inputs (optional, default OFF).

By default the calibrator uses the ``.WTH`` files DSSAT already resolves for each
experiment (``provider: file`` — a no-op acquisition layer). When a site has no
shipped weather, or you want an operational in-season nowcast, switch on a
provider that *acquires* daily weather and writes a DSSAT ``.WTH``:

    weather:
      provider: file | nasa_power | dssatutils      # default: file
      cache_dir: "weather_cache"
      gap_fill: none | persistence | climatology   # fill days past the last record
      horizon: 0                       # extra days to extend past the last obs (forecast)

Design notes
------------
* Everything here is *optional*. With the default ``provider: file`` nothing in
  this module runs, so existing runs are unaffected.
* ``NasaPowerProvider`` hits the keyless NASA POWER daily-point API at call time
  (needs network); the parsing, ``.WTH`` writing and gap-fill are pure functions
  and are unit-tested offline.
* NASA POWER is reanalysis with ~1–2 week latency, so an in-season forecast needs
  the ``gap_fill`` step to extend the record up to (and a ``horizon`` past) today.
  ``climatology`` uses the day-of-year mean of the available record; ``persistence``
  repeats a trailing window. Neither is a true NWP forecast — they are documented
  stand-ins until a forecast provider is wired in.

Layering note (workspace AGENTS.md principle 1 — documented divergence)
-----------------------------------------------------------------------
The foundation library ``dssatutils`` already owns weather acquisition
(``process_weather_nasapower`` and AgERA5 / Daymet / ERA5-Land / … providers). This
module is a **deliberately lightweight, single-point provider for the in-season
nowcast**: ``dssatutils.process_weather_nasapower`` is a *batch, whole-year,
geopandas/shapefile* downloader, whereas the nowcast needs a *single point, partial
season, gap-filled to "today" plus a horizon*. We intentionally do not fork the rich
multi-source download layer — for bulk/site acquisition use ``dssatutils`` directly;
this provider exists only for the operational forward-simulation path. Consuming
``provider: dssatutils`` is the shared-stack path when the optional acquisition
extra is installed; ``provider: nasa_power`` remains the zero-dependency fallback.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# canonical daily columns we manage in a .WTH
WTH_COLS = ["SRAD", "TMAX", "TMIN", "RAIN"]

# NASA POWER parameter -> our column
_POWER_MAP = {
    "ALLSKY_SFC_SW_DWN": "SRAD",
    "T2M_MAX": "TMAX",
    "T2M_MIN": "TMIN",
    "PRECTOTCORR": "RAIN",
}


def _yyddd(d: pd.Timestamp) -> str:
    return f"{d.year % 100:02d}{d.dayofyear:03d}"


def write_wth(path: str | Path, station: str, lat: float, lon: float,
              df: pd.DataFrame, elev: float = -99.0) -> Path:
    """Write a minimal but valid DSSAT ``.WTH`` from a daily DataFrame.

    ``df`` must have a ``date`` column plus ``SRAD/TMAX/TMIN/RAIN`` (missing
    columns are written as -99). ``TAV``/``AMP`` are derived from the data.
    """
    path = Path(path)
    df = df.sort_values("date").reset_index(drop=True)
    tmean = ((df.get("TMAX", pd.Series(dtype=float)) + df.get("TMIN", pd.Series(dtype=float))) / 2)
    tav = float(np.nanmean(tmean)) if len(tmean) else -99.0
    monthly = tmean.groupby(pd.to_datetime(df["date"]).dt.month).mean() if len(tmean) else pd.Series(dtype=float)
    amp = float((monthly.max() - monthly.min())) if len(monthly) else -99.0

    insi = (station or "CALB")[:4].upper().ljust(4)
    head = [
        f"*WEATHER DATA : {station}",
        "",
        "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
        f"  {insi} {lat:8.3f} {lon:8.3f} {elev:5.0f} {tav:5.1f} {amp:5.1f}   2.0  10.0",
        "@DATE  SRAD  TMAX  TMIN  RAIN",
    ]
    rows = []
    for _, r in df.iterrows():
        vals = []
        for c in WTH_COLS:
            v = r.get(c, np.nan)
            vals.append(f"{float(v):6.1f}" if pd.notna(v) else f"{-99.0:6.1f}")
        rows.append(f"{_yyddd(pd.Timestamp(r['date']))}" + "".join(vals))
    path.write_text("\n".join(head + rows) + "\n", encoding="utf-8")
    return path


def _date_from_dssat_code(code: str) -> pd.Timestamp:
    code = str(code).strip()
    if len(code) <= 5:
        yy = int(code[:-3])
        year = 2000 + yy if yy < 80 else 1900 + yy
    else:
        year = int(code[:-3])
    doy = int(code[-3:])
    return pd.to_datetime(f"{year}{doy:03d}", format="%Y%j")


def read_wth(path: str | Path) -> pd.DataFrame:
    """Read a DSSAT `.WTH` file into ``date`` + ``SRAD/TMAX/TMIN/RAIN``."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if not parts or not parts[0].isdigit():
            continue
        values = parts[1:]
        row = {"date": _date_from_dssat_code(parts[0])}
        for idx, col in enumerate(WTH_COLS):
            try:
                val = float(values[idx])
            except (IndexError, ValueError):
                val = np.nan
            row[col] = np.nan if val in (-99, -99.0, -999, -999.0) else val
        rows.append(row)
    return pd.DataFrame(rows, columns=["date", *WTH_COLS])


def fill_gap(df: pd.DataFrame, end_date, *, method: str = "none",
             window: int = 7) -> pd.DataFrame:
    """Extend a daily weather record forward to ``end_date``.

    ``persistence`` repeats the mean of the trailing ``window`` days; ``climatology``
    uses the day-of-year mean over the whole record. Returns the record with the
    appended (synthetic) days flagged in a ``filled`` boolean column.
    """
    if method in (None, "none") or df.empty:
        out = df.copy()
        out["filled"] = False
        return out

    df = df.sort_values("date").reset_index(drop=True)
    df["date"] = pd.to_datetime(df["date"])
    last = df["date"].max()
    end_date = pd.Timestamp(end_date)
    if end_date <= last:
        out = df.copy()
        out["filled"] = False
        return out

    future = pd.date_range(last + pd.Timedelta(days=1), end_date, freq="D")
    rows = []
    if method == "climatology":
        doy_mean = df.groupby(df["date"].dt.dayofyear)[WTH_COLS].mean()
        for d in future:
            base = doy_mean.reindex([d.dayofyear]).iloc[0] if d.dayofyear in doy_mean.index else df[WTH_COLS].mean()
            rows.append({"date": d, **{c: base.get(c, np.nan) for c in WTH_COLS}})
    elif method == "persistence":
        tail = df.tail(window)[WTH_COLS].mean()
        for d in future:
            rows.append({"date": d, **{c: tail.get(c, np.nan) for c in WTH_COLS}})
    else:
        raise ValueError(f"Unknown gap_fill method '{method}'")

    filled = pd.DataFrame(rows)
    df["filled"] = False
    filled["filled"] = True
    return pd.concat([df, filled], ignore_index=True)


@dataclass
class WeatherProvider:
    """Base class. ``fetch`` returns a daily DataFrame (date + WTH_COLS)."""
    cfg: dict

    def fetch(self, station: str, lat: float, lon: float, start, end) -> pd.DataFrame:  # pragma: no cover
        raise NotImplementedError


class FileWeatherProvider(WeatherProvider):
    """Default no-op: rely on the ``.WTH`` DSSAT already resolves. ``fetch`` reads
    an existing station file under ``dssat_dir/Weather`` if asked, else returns empty."""

    def fetch(self, station, lat, lon, start, end) -> pd.DataFrame:
        return pd.DataFrame(columns=["date", *WTH_COLS])


class NasaPowerProvider(WeatherProvider):
    """Acquire daily weather from the keyless NASA POWER point API (needs network)."""

    URL = ("https://power.larc.nasa.gov/api/temporal/daily/point"
           "?parameters={params}&community=AG&longitude={lon}&latitude={lat}"
           "&start={start}&end={end}&format=JSON")

    def fetch(self, station, lat, lon, start, end) -> pd.DataFrame:
        import json
        import urllib.request

        start = pd.Timestamp(start)
        end = pd.Timestamp(end)
        cache_dir = Path(self.cfg.get("weather", {}).get("cache_dir", "weather_cache"))
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache = cache_dir / f"power_{lat:.3f}_{lon:.3f}_{start:%Y%m%d}_{end:%Y%m%d}.json"

        if cache.exists():
            payload = json.loads(cache.read_text(encoding="utf-8"))
        else:
            url = self.URL.format(params=",".join(_POWER_MAP), lon=lon, lat=lat,
                                  start=start.strftime("%Y%m%d"), end=end.strftime("%Y%m%d"))
            logger.info("Fetching NASA POWER: %s", url)
            with urllib.request.urlopen(url, timeout=60) as resp:
                payload = json.loads(resp.read().decode())
            cache.write_text(json.dumps(payload), encoding="utf-8")
        return self._parse_power(payload)

    @staticmethod
    def _parse_power(payload: dict) -> pd.DataFrame:
        """Pure parse of a NASA POWER JSON payload into the daily schema."""
        params = payload["properties"]["parameter"]
        dates = sorted(next(iter(params.values())).keys())
        rows = []
        for d in dates:
            row = {"date": pd.to_datetime(d, format="%Y%m%d")}
            for p, col in _POWER_MAP.items():
                v = params.get(p, {}).get(d, -999)
                row[col] = np.nan if v in (-999, -99) else float(v)
            rows.append(row)
        return pd.DataFrame(rows)


class DssatutilsWeatherProvider(WeatherProvider):
    """Acquire whole-year weather through the shared ``dssatutils`` API."""

    def fetch(self, station, lat, lon, start, end) -> pd.DataFrame:
        from .acquisition import _single_point_gdf, _valid_coordinate

        lat_f, lon_f = _valid_coordinate(lat, lon)
        start = pd.Timestamp(start)
        end = pd.Timestamp(end)
        wcfg = self.cfg.get("weather", {}) or {}
        cache_dir = Path(wcfg.get("cache_dir", "weather_cache"))
        cache_dir.mkdir(parents=True, exist_ok=True)

        provider = str(wcfg.get("dssatutils_provider", "nasapower")).lower()
        provider = "nasapower" if provider in ("nasa_power", "nasa-power") else provider
        fn_name = f"process_weather_{provider}"

        try:
            import dssatutils
        except ImportError as exc:
            raise ImportError(
                "weather.provider: dssatutils requires the optional 'acquire' extra."
            ) from exc
        if not hasattr(dssatutils, fn_name):
            raise ValueError(f"dssatutils has no weather provider {fn_name!r}")

        gdf = _single_point_gdf(station, lat_f, lon_f)
        getattr(dssatutils, fn_name)(
            gdf,
            start_year=start.year,
            end_year=end.year,
            output_dir=str(cache_dir),
            id_col="ID",
            lat_col="LAT",
            lon_col="LONG",
            n_cores=int(wcfg.get("n_cores", 1) or 1),
            log_file=str(cache_dir / "dssatutils_weather.log"),
        )

        wth = cache_dir / f"{station}.WTH"
        if not wth.exists():
            raise FileNotFoundError(f"dssatutils did not produce {wth}")
        df = read_wth(wth)
        if df.empty:
            return df
        return df[(df["date"] >= start) & (df["date"] <= end)].reset_index(drop=True)


PROVIDERS = {
    "file": FileWeatherProvider,
    "nasa_power": NasaPowerProvider,
    "dssatutils": DssatutilsWeatherProvider,
}


def build_provider(cfg: dict) -> WeatherProvider:
    name = str(cfg.get("weather", {}).get("provider", "file")).lower()
    if name not in PROVIDERS:
        raise ValueError(f"Unknown weather provider '{name}'. Available: {list(PROVIDERS)}")
    return PROVIDERS[name](cfg)


def acquire_wth(cfg: dict, *, station: str, lat: float, lon: float, start, end,
                out_path: str | Path) -> Path | None:
    """Top-level helper: acquire weather, gap-fill to ``end`` (+horizon), write a ``.WTH``.

    Returns the written path, or ``None`` for the default ``file`` provider (which
    leaves DSSAT's own resolution untouched). Optional — only called when a caller
    explicitly opts into weather acquisition.
    """
    wcfg = cfg.get("weather", {})
    provider = build_provider(cfg)
    if isinstance(provider, FileWeatherProvider):
        return None
    horizon = int(wcfg.get("horizon", 0))
    end = pd.Timestamp(end) + pd.Timedelta(days=horizon)
    df = provider.fetch(station, lat, lon, start, end)
    if df.empty:
        logger.warning("Weather provider returned no data for %s", station)
        return None
    df = fill_gap(df, end, method=wcfg.get("gap_fill", "none"))
    return write_wth(out_path, station, lat, lon, df)
