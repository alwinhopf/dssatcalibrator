"""Optional shared-stack acquisition helpers for synthesized/new-site runs.

Real experiment calibration keeps using the FileX, weather and soil shipped
with the DSSAT experiment. These helpers are only for cases where a site has
coordinates but no local `.WTH` / `.SOL`, and they delegate acquisition to the
workspace foundation package `dssatutils`.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def _valid_coordinate(lat: float | int | str | None,
                      lon: float | int | str | None) -> tuple[float, float]:
    try:
        lat_f = float(lat)
        lon_f = float(lon)
    except (TypeError, ValueError) as exc:
        raise ValueError("Weather/soil acquisition needs numeric latitude and longitude.") from exc
    if abs(lat_f) > 90 or abs(lon_f) > 180 or lat_f <= -98 or lon_f <= -998:
        raise ValueError(f"Invalid or missing site coordinates: lat={lat!r}, lon={lon!r}")
    return lat_f, lon_f


def _single_point_gdf(site_id: str, lat: float, lon: float):
    try:
        import geopandas as gpd
        from shapely.geometry import Point
    except ImportError as exc:
        raise ImportError(
            "dssatutils acquisition requires the optional geospatial stack. "
            "Install with 'pip install -e .[acquire]'."
        ) from exc

    return gpd.GeoDataFrame(
        [{"ID": str(site_id), "LAT": float(lat), "LONG": float(lon)}],
        geometry=[Point(float(lon), float(lat))],
        crs="EPSG:4326",
    )


def _copy_mapped_soil(sol_dir: Path, map_csv: Path, site_id: str, out_path: Path) -> Path:
    direct = sol_dir / f"{site_id}.SOL"
    src = direct if direct.exists() else None

    if src is None and map_csv.exists():
        mapping = pd.read_csv(map_csv)
        if "ID" in mapping.columns:
            row = mapping[mapping["ID"].astype(str) == str(site_id)]
            if not row.empty:
                for col in ("SOIL_ID", "soil_id", "SOURCE_SOIL_ID"):
                    if col in row.columns and pd.notna(row.iloc[0][col]):
                        candidate = sol_dir / f"{row.iloc[0][col]}.SOL"
                        if candidate.exists():
                            src = candidate
                            break

    if src is None:
        raise FileNotFoundError(
            f"dssatutils did not produce a .SOL for site {site_id!r} in {sol_dir}"
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, out_path)
    return out_path


def acquire_soil_profile(cfg: dict, *, site_id: str, lat: float, lon: float,
                         out_path: str | Path) -> Path | None:
    """Acquire one DSSAT soil profile with `dssatutils` and copy it to `out_path`.

    `soil.provider: file` is a no-op. `soil.provider: dssatutils` supports
    `soil.source: ssurgo`, `soilgrids`, and `soilgrids_online`.
    """
    scfg = cfg.get("soil", {}) or {}
    provider = str(scfg.get("provider", "file")).lower()
    if provider in ("", "file", "none"):
        return None
    if provider != "dssatutils":
        raise ValueError("soil.provider must be 'file' or 'dssatutils'.")

    lat_f, lon_f = _valid_coordinate(lat, lon)
    gdf = _single_point_gdf(site_id, lat_f, lon_f)
    cache_dir = Path(scfg.get("cache_dir", "soil_cache"))
    sol_dir = cache_dir / f"{str(scfg.get('source', 'ssurgo')).lower()}_individual_SOL"
    sol_dir.mkdir(parents=True, exist_ok=True)
    map_csv = cache_dir / f"{str(scfg.get('source', 'ssurgo')).lower()}_soil_map.csv"
    source = str(scfg.get("source", "ssurgo")).lower()
    n_cores = int(scfg.get("n_cores", 1) or 1)

    try:
        import dssatutils
    except ImportError as exc:
        raise ImportError(
            "soil.provider: dssatutils requires the optional 'acquire' extra."
        ) from exc

    if source == "ssurgo":
        dssatutils.process_soils_ssurgo(
            gdf,
            output_dir_csv=str(map_csv),
            output_dir_individual=str(sol_dir),
            n_cores=n_cores,
            id_col="ID",
            lat_col="LAT",
            long_col="LONG",
        )
    elif source == "soilgrids":
        source_sol = scfg.get("source_sol_file") or scfg.get("external_soil_file")
        if not source_sol:
            raise ValueError(
                "soil.source: soilgrids requires soil.source_sol_file "
                "or soil.external_soil_file."
            )
        dssatutils.process_soils_soilgrids(
            gdf,
            source_sol_file=str(source_sol),
            output_csv_path=str(map_csv),
            output_sol_dir=str(sol_dir),
            id_col="ID",
        )
    elif source == "soilgrids_online":
        soilgrids = __import__("dssatutils.soil_soilgrids_online",
                               fromlist=["process_soils_soilgrids_online"])
        mode = str(scfg.get("soilgrids_mode", "REST")).upper()
        if hasattr(soilgrids, "USE_REST_API"):
            soilgrids.USE_REST_API = mode != "VRT"
        soilgrids.process_soils_soilgrids_online(
            gdf,
            soilfile_csv_path=str(map_csv),
            output_sol_dir=str(sol_dir),
            id_col="ID",
        )
    else:
        raise ValueError(
            "soil.source must be one of: ssurgo, soilgrids, soilgrids_online."
        )

    acquired = _copy_mapped_soil(sol_dir, map_csv, str(site_id), Path(out_path))
    logger.info("Acquired soil profile for %s via dssatutils:%s -> %s",
                site_id, source, acquired)
    return acquired
