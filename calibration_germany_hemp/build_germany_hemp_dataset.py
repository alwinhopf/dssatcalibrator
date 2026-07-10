"""Build an approximate Germany hemp calibration dataset from the two papers.

The papers do not publish raw plot-level data. This script creates an auditable
derived dataset from:
  * manual visual digitization of Agronomy 2020 Figure 3,
  * scalar yield/stand/soil/weather aggregates reported in the two PDFs,
  * DWD station weather via the local dssatutils implementation,
  * SoilGrids REST profile via the local dssatutils implementation.

Outputs are written below calibration_germany_hemp/derived/.
"""

from __future__ import annotations

import math
import os
import json
import shutil
import subprocess
import sys
from datetime import date, timedelta
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parents[1]
DSSATUTILS_PYTHON = WORKSPACE / "dssatutils" / "python"
DERIVED = ROOT / "derived"
WEATHER_DIR = DERIVED / "weather_dwd"
DWD_CACHE_DIR = DERIVED / "dwd_cache"
SOIL_DIR = DERIVED / "soil"
FIGURE_DIR = DERIVED / "source_figures"

SITE_ID = "GEMQ2018"
LAT = 52 + 28 / 60 + 2 / 3600
LON = 12 + 57 / 60 + 39 / 3600
SITE_NAME = "Marquardt/Potsdam field site"

SOWING_DATES = {
    "Santhica 27": date(2018, 5, 4),
    "Ivory": date(2018, 5, 22),
}

PDF_AGRONOMY = ROOT / "agronomy_10_01361_v2.pdf"
PDF_WATER = ROOT / "water_12_02982_v2.pdf"


def ensure_import_paths() -> None:
    path = str(DSSATUTILS_PYTHON)
    if path not in sys.path:
        sys.path.insert(0, path)
    install_requests_fallback()


def install_requests_fallback() -> None:
    """Provide the tiny requests.get surface dssatutils needs if requests is absent."""
    try:
        import requests  # noqa: F401

        return
    except ImportError:
        pass

    import types
    import ssl
    import urllib.error
    import urllib.parse
    import urllib.request

    class SimpleResponse:
        def __init__(self, status_code: int, content: bytes, url: str):
            self.status_code = status_code
            self.content = content
            self.url = url
            self.ok = 200 <= status_code < 300

        @property
        def text(self) -> str:
            try:
                return self.content.decode("utf-8")
            except UnicodeDecodeError:
                return self.content.decode("latin-1")

        def json(self):
            return json.loads(self.text)

        def raise_for_status(self) -> None:
            if not self.ok:
                raise RuntimeError(f"HTTP {self.status_code} for {self.url}")

    def get(url: str, params=None, timeout=None):
        if params:
            query = urllib.parse.urlencode(params, doseq=True)
            sep = "&" if "?" in url else "?"
            url = f"{url}{sep}{query}"
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return SimpleResponse(response.status, response.read(), url)
        except urllib.error.HTTPError as exc:
            return SimpleResponse(exc.code, exc.read(), url)
        except urllib.error.URLError as exc:
            if "CERTIFICATE_VERIFY_FAILED" not in str(exc):
                raise
            context = ssl._create_unverified_context()
            with urllib.request.urlopen(url, timeout=timeout, context=context) as response:
                return SimpleResponse(response.status, response.read(), url)

    module = types.ModuleType("requests")
    module.get = get
    sys.modules["requests"] = module


def render_figure3_source() -> None:
    """Render and crop Agronomy Figure 3 for audit of digitized points."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out_png = FIGURE_DIR / "agronomy_figure3_digitization_source.png"
    if out_png.exists():
        return

    bundled = Path.home() / ".cache" / "codex-runtimes" / "codex-primary-runtime" / "dependencies" / "native" / "poppler" / "Library" / "bin" / "pdftoppm.exe"
    pdftoppm = str(bundled) if bundled.exists() else shutil.which("pdftoppm")
    if pdftoppm is None:
        print("Could not find pdftoppm; skipping Figure 3 source crop render.")
        return

    try:
        from PIL import Image
    except ImportError:
        print("Pillow is unavailable; skipping Figure 3 source crop render.")
        return

    tmp_prefix = FIGURE_DIR / "agronomy_figure3_page"
    subprocess.run(
        [
            pdftoppm,
            "-png",
            "-r",
            "300",
            "-f",
            "7",
            "-l",
            "7",
            str(PDF_AGRONOMY),
            str(tmp_prefix),
        ],
        check=True,
    )
    rendered = FIGURE_DIR / "agronomy_figure3_page-07.png"
    if not rendered.exists():
        matches = sorted(FIGURE_DIR.glob("agronomy_figure3_page*.png"))
        if not matches:
            return
        rendered = matches[0]

    with Image.open(rendered) as img:
        # Crop includes all four panels, axis labels, legend, and caption lead.
        crop = img.crop((430, 880, 2100, 2050))
        crop.save(out_png)


def obs_date(cultivar: str, das: float) -> str:
    return (SOWING_DATES[cultivar] + timedelta(days=int(round(das)))).isoformat()


def yyddd(value: date) -> str:
    return f"{value.year % 100:02d}{value.timetuple().tm_yday:03d}"


def digitized_growth_observations() -> pd.DataFrame:
    """Approximate points digitized visually from Agronomy 2020 Figure 3."""
    rows: list[dict] = []

    def add(cultivar: str, variable: str, unit: str, points: Iterable[tuple[float, float]], panel: str) -> None:
        for das, value in points:
            rows.append(
                {
                    "site_id": SITE_ID,
                    "site_name": SITE_NAME,
                    "cultivar": cultivar,
                    "sowing_date": SOWING_DATES[cultivar].isoformat(),
                    "date": obs_date(cultivar, das),
                    "das": das,
                    "variable": variable,
                    "value": value,
                    "unit": unit,
                    "source_pdf": PDF_AGRONOMY.name,
                    "source_detail": f"Figure 3{panel}",
                    "digitization_method": "manual visual digitization from 300 dpi rendered PDF; axis-calibrated; approximate",
                    "notes": "Plot mean point digitized; printed error bars are not encoded as observation uncertainty.",
                }
            )

    # Panel a: stem/plant height. Axis: days after sowing vs plant height (cm).
    add(
        "Santhica 27",
        "plant_height",
        "cm",
        [
            (0, 0), (8, 4), (15, 13), (22, 28), (29, 46), (36, 74),
            (43, 104), (50, 136), (57, 151), (65, 163), (75, 190),
            (83, 205), (90, 207), (98, 206), (105, 217), (113, 218),
            (121, 219), (128, 222), (134, 223),
        ],
        "a",
    )
    add(
        "Ivory",
        "plant_height",
        "cm",
        [
            (0, 0), (8, 2), (16, 5), (23, 14), (30, 34), (37, 55),
            (44, 83), (52, 104), (58, 132), (65, 160), (72, 164),
            (80, 176), (88, 176), (96, 179), (104, 180), (112, 181),
        ],
        "a",
    )

    # Panel b: individual plant leaf area.
    add(
        "Santhica 27",
        "leaf_area_per_plant",
        "m2 plant-1",
        [
            (43, 0.040), (50, 0.050), (64, 0.062), (70, 0.068),
            (77, 0.080), (86, 0.075), (92, 0.075), (100, 0.064),
            (113, 0.070), (121, 0.025), (128, 0.027), (134, 0.022),
        ],
        "b",
    )
    add(
        "Ivory",
        "leaf_area_per_plant",
        "m2 plant-1",
        [
            (22, 0.010), (29, 0.040), (43, 0.040), (50, 0.096),
            (59, 0.140), (64, 0.118), (70, 0.128), (78, 0.080),
            (92, 0.073), (99, 0.064), (106, 0.035), (113, 0.032),
        ],
        "b",
    )

    # Panel c: plant density.
    add(
        "Santhica 27",
        "plant_density",
        "plants m-2",
        [
            (43, 149), (50, 131), (56, 128), (64, 119), (72, 112),
            (80, 106), (86, 102), (93, 98), (100, 103), (106, 91),
            (113, 89), (121, 82), (128, 78), (134, 83),
        ],
        "c",
    )
    add(
        "Ivory",
        "plant_density",
        "plants m-2",
        [
            (22, 58), (29, 54), (36, 57), (43, 59), (51, 52),
            (58, 53), (65, 50), (72, 49), (79, 47), (86, 43),
            (93, 44), (100, 44), (106, 44), (113, 42),
        ],
        "c",
    )

    # Panel d: LAI. Values are cross-checked to Water 2020 text maxima on 2018-07-20.
    add(
        "Santhica 27",
        "leaf_area_index",
        "m2 m-2",
        [
            (43, 5.6), (50, 6.1), (64, 7.4), (70, 7.8), (77, 8.8),
            (86, 7.8), (92, 7.4), (100, 6.4), (113, 6.1),
            (121, 2.0), (128, 2.0), (134, 1.8),
        ],
        "d",
    )
    add(
        "Ivory",
        "leaf_area_index",
        "m2 m-2",
        [
            (22, 0.6), (29, 2.1), (43, 2.1), (50, 5.0), (59, 7.5),
            (64, 6.0), (70, 6.4), (78, 3.4), (92, 2.7),
            (99, 2.5), (106, 1.5), (113, 1.2),
        ],
        "d",
    )

    return pd.DataFrame(rows)


def scalar_observations() -> pd.DataFrame:
    rows = [
        ("Ivory", "fresh_mass_whole_plant", 14.4, 6.3, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Ivory", "dry_mass_whole_plant", 10.0, 3.9, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Ivory", "dry_mass_straw", 8.8, 3.4, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Ivory", "fresh_mass_bast", 4.2, 1.9, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Ivory", "dry_mass_bast", 2.1, 0.9, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Ivory", "field_emergence", 30.0, np.nan, "%", "Water Section 2.1"),
        ("Ivory", "plant_density_start_growing_period", 60.0, np.nan, "plants m-2", "Water Table 1"),
        ("Ivory", "plant_density_harvest", 43.0, np.nan, "plants m-2", "Water Table 1"),
        ("Santhica 27", "fresh_mass_whole_plant", 24.5, 10.3, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Santhica 27", "dry_mass_whole_plant", 17.9, 6.3, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Santhica 27", "dry_mass_straw", 16.0, 5.1, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Santhica 27", "fresh_mass_bast", 6.3, 2.6, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Santhica 27", "dry_mass_bast", 2.8, 0.8, "t ha-1", "Agronomy Table 1; Water Table 1"),
        ("Santhica 27", "field_emergence", 75.0, np.nan, "%", "Water Section 2.1"),
        ("Santhica 27", "plant_density_start_growing_period", 149.0, np.nan, "plants m-2", "Water Table 1"),
        ("Santhica 27", "plant_density_harvest", 79.0, np.nan, "plants m-2", "Water Table 1"),
    ]
    out = []
    for cultivar, variable, value, sd, unit, source in rows:
        out.append(
            {
                "site_id": SITE_ID,
                "cultivar": cultivar,
                "sowing_date": SOWING_DATES[cultivar].isoformat(),
                "observation_date": date(2018, 9, 19).isoformat()
                if "mass" in variable or variable.endswith("harvest")
                else "",
                "variable": variable,
                "value": value,
                "sd": sd,
                "unit": unit,
                "source": source,
                "notes": "Scalar paper value; not plot-level raw data.",
            }
        )
    return pd.DataFrame(out)


def management_summary() -> pd.DataFrame:
    rows = [
        ("site_id", SITE_ID, ""),
        ("site_name", SITE_NAME, ""),
        ("latitude", LAT, "decimal degrees"),
        ("longitude", LON, "decimal degrees"),
        ("soil_description", "Eutric Lamellic Brunic Arenosol (Aric), low loamy sand SL2", "reported in Water Figure 1 text/caption"),
        ("groundwater_depth", 10, "m below surface"),
        ("santhica_sowing_date", "2018-05-04", ""),
        ("ivory_sowing_date", "2018-05-22", ""),
        ("planting_density", 200, "seeds m-2"),
        ("row_spacing", 12.5, "cm"),
        ("n_fertilizer", 70, "kg N ha-1 calcium ammonium nitrate"),
        ("sprinkler_irrigation", 10, "mm on 2018-05-30"),
        ("harvest_or_oven_dry_date", "2018-09-19", "Water paper; Agronomy has minor 2018-09-18/19 ambiguity"),
    ]
    return pd.DataFrame(rows, columns=["field", "value", "unit_or_note"])


def filex_management_inputs() -> pd.DataFrame:
    """Paper-derived management fields useful for a DSSAT FileX scaffold."""
    cultivar_rows = [
        {
            "exp_id": "GEMQ18S1",
            "cultivar": "Santhica 27",
            "ingeno_placeholder": "DE0001",
            "pdate": SOWING_DATES["Santhica 27"],
            "area_ha": 0.12,
            "field_emergence_pct": 75,
            "ppoe_plants_m2": 149,
            "plant_density_harvest_plants_m2": 79,
        },
        {
            "exp_id": "GEMQ18I1",
            "cultivar": "Ivory",
            "ingeno_placeholder": "DE0002",
            "pdate": SOWING_DATES["Ivory"],
            "area_ha": 0.08,
            "field_emergence_pct": 30,
            "ppoe_plants_m2": 60,
            "plant_density_harvest_plants_m2": 43,
        },
    ]
    rows = []
    for item in cultivar_rows:
        rows.append(
            {
                "exp_id": item["exp_id"],
                "cultivar": item["cultivar"],
                "treatment": 1,
                "cname": item["cultivar"],
                "ingeno": item["ingeno_placeholder"],
                "wsta": SITE_ID,
                "id_soil": SITE_ID,
                "field_area_ha": item["area_ha"],
                "pdate": item["pdate"].isoformat(),
                "pdate_yydoy": yyddd(item["pdate"]),
                "ppop_planned_plants_m2": 200,
                "ppop_source_value": "200 seeds m-2",
                "ppoe_emerged_plants_m2": item["ppoe_plants_m2"],
                "field_emergence_pct": item["field_emergence_pct"],
                "plrs_row_spacing_cm": 12.5,
                "pldp_planting_depth_cm": 2.0,
                "plme": "S",
                "plds": "R",
                "initial_relative_soil_moisture_pct": 75,
                "initial_soil_water_method": "75 percent relative soil moisture; convert per layer as SLLL + 0.75 * (SDUL - SLLL)",
                "initial_residual_mineral_n_kg_ha": 100,
                "initial_n_distribution": "user-specified residual N total; layer distribution still needs a DSSAT convention",
                "residue_management": "none",
                "tillage_management": "none",
                "chemical_applications": "none assumed",
                "irrigation_date": date(2018, 5, 30).isoformat(),
                "irrigation_yydoy": yyddd(date(2018, 5, 30)),
                "irrigation_amount_mm": 10,
                "irrigation_method_hint": "sprinkler",
                "n_fertilizer_date_assumed": (item["pdate"] - timedelta(days=1)).isoformat(),
                "n_fertilizer_yydoy_assumed": yyddd(item["pdate"] - timedelta(days=1)),
                "n_fertilizer_amount_kg_ha": 70,
                "n_fertilizer_material": "calcium ammonium nitrate",
                "n_fertilizer_date_confidence": "assumed as the day before each cultivar planting date; PDF gives rate but not exact date",
                "harvest_date": date(2018, 9, 19).isoformat(),
                "harvest_yydoy": yyddd(date(2018, 9, 19)),
                "plot_replicates": 3,
                "plot_area_m2": 0.25,
                "directly_reported_fields": "pdate, ppop seed rate, ppoe/start density, row spacing, irrigation, N rate, harvest date",
                "user_specified_assumptions": "PLDP=2 cm; initial soil water=75 percent relative; residual N=100 kg ha-1; no tillage; no residue",
                "fields_needing_assumption": "INGENO/ecotype, exact N date, initial N layer distribution convention, WSTA/header harmonization",
                "source": "Agronomy 2020 Section 2.3; Water 2020 Section 2.1 and Table 1",
            }
        )
    return pd.DataFrame(rows)


def write_filex_management_audit(df: pd.DataFrame) -> None:
    path = DERIVED / "filex_management_extractability.md"
    lines = [
        "# FileX Management Extractability",
        "",
        "The papers provide enough management information to create a useful FileX scaffold,",
        "but not enough to create a fully assumption-free DSSAT experiment file.",
        "",
        "## Directly Reported",
        "",
        "- Site coordinates: 52.467222 N, 12.960833 E.",
        "- Cultivars: Santhica 27 and Ivory.",
        "- Sowing dates: Santhica 27 on 2018-05-04; Ivory on 2018-05-22.",
        "- Sowing density: 200 seeds m-2.",
        "- Inter-row distance: 12.5 cm.",
        "- N fertilization rate: 70 kg N ha-1 as calcium ammonium nitrate.",
        "- Irrigation: 10 mm sprinkler irrigation on 2018-05-30.",
        "- Harvest/oven-dry date: 2018-09-19 in the Water paper; Agronomy has an 18/19 September ambiguity for yield handling.",
        "- Emergence/stand: Santhica 27 75 percent emergence, 149 plants m-2 start, 79 plants m-2 harvest; Ivory 30 percent emergence, 60 plants m-2 start, 43 plants m-2 harvest.",
        "- Plot geometry: three 0.25 m2 plots per cultivar; total cultivation area 0.12 ha Santhica 27 and 0.08 ha Ivory.",
        "",
        "## User-Specified Scaffold Assumptions",
        "",
        "- Planting depth: 2 cm.",
        "- Initial soil water: 75 percent relative soil moisture, to be converted per layer as `SLLL + 0.75 * (SDUL - SLLL)`.",
        "- Initial residual mineral N: 100 kg ha-1 total; layer distribution still needs a DSSAT convention.",
        "- Tillage: none.",
        "- Residue/organic amendments: none.",
        "- Kieserite/MgSO4 application: omitted from the runnable scaffold.",
        "",
        "## FileX Assumptions Still Needed",
        "",
        "- INGENO/genotype rows for Santhica 27 and Ivory are not present in the local hemp genotype file.",
        "- Exact N fertilizer application date is not reported; the scaffold assumes the day before each cultivar planting date.",
        "- Initial residual N is now specified as a total, but DSSAT still needs NH4/NO3 values by layer. Use an explicit distribution rule before writing FileX.",
        "- Weather station naming should be harmonized if a runnable FileX is generated: the derived WTH is `GEMQ2018.WTH` from DWD Potsdam station 03987.",
        "",
        "## Scaffold Rows",
        "",
        "```text",
        df.to_string(index=False),
        "```",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def calibrator_observation_table(growth: pd.DataFrame, scalars: pd.DataFrame) -> pd.DataFrame:
    """Create the native dssatcalibrator long observation schema."""
    exp_ids = {"Santhica 27": "GEMQ18S1", "Ivory": "GEMQ18I1"}
    rows = []

    growth_map = {
        "plant_height": ("CHTD", 0.01, "cm converted to m"),
        "leaf_area_index": ("LAID", 1.0, "unitless LAI"),
        # The digitized stand-density series is retained in
        # digitized_growth_observations.csv, but not mapped into the native
        # calibration table: CRGRO PlantGro P#AD is not a safe plant-density
        # analogue for hemp stand counts.
    }
    for _, row in growth.iterrows():
        mapped = growth_map.get(row["variable"])
        if mapped is None:
            continue
        variable, factor, note = mapped
        rows.append(
            {
                "exp_id": exp_ids[row["cultivar"]],
                "treatment": 1,
                "variable": variable,
                "kind": "timeseries",
                "date": row["date"],
                "value": float(row["value"]) * factor,
                "sigma": np.nan,
                "weight": 1.0,
                "cultivar": row["cultivar"],
                "source_variable": row["variable"],
                "source_detail": row["source_detail"],
                "notes": note,
            }
        )

    scalar_map = {
        "dry_mass_whole_plant": ("CWAD", 1000.0, "t ha-1 converted to kg ha-1"),
        "dry_mass_straw": ("SWAD", 1000.0, "t ha-1 converted to kg ha-1"),
    }
    for _, row in scalars.iterrows():
        mapped = scalar_map.get(row["variable"])
        if mapped is None:
            continue
        variable, factor, note = mapped
        rows.append(
            {
                "exp_id": exp_ids[row["cultivar"]],
                "treatment": 1,
                "variable": variable,
                "kind": "timeseries",
                "date": row["observation_date"],
                "value": float(row["value"]) * factor,
                "sigma": float(row["sd"]) * factor if pd.notna(row["sd"]) else np.nan,
                "weight": 1.0,
                "cultivar": row["cultivar"],
                "source_variable": row["variable"],
                "source_detail": row["source"],
                "notes": note,
            }
        )

    return pd.DataFrame(rows)


def phenology_targets() -> pd.DataFrame:
    rows = [
        {
            "site_id": SITE_ID,
            "cultivar": "Ivory",
            "target": "end_flowering_begin_senescence",
            "date": obs_date("Ivory", 52),
            "das": 52,
            "confidence": "low",
            "source_pdf": PDF_AGRONOMY.name,
            "source_location": "Section 4.3, page 11",
            "interpretation": "The paper interprets maximum photosynthesis around 52 DAS as possibly marking the end of flowering and beginning of senescence for early-mature Ivory, by analogy to Chameleon.",
            "usable_for_calibration": "weak phenology prior, not a direct observed flowering date",
        },
        {
            "site_id": SITE_ID,
            "cultivar": "Santhica 27",
            "target": "begin_senescence",
            "date": obs_date("Santhica 27", 70),
            "das": 70,
            "confidence": "low",
            "source_pdf": PDF_AGRONOMY.name,
            "source_location": "Section 4.3, page 11",
            "interpretation": "Seasonal gas exchange is described as roughly constant until at least 70 DAS; the authors state this timing fits earlier reported beginning of senescence for Santhica 27.",
            "usable_for_calibration": "weak senescence prior, not a direct observed stage date",
        },
        {
            "site_id": SITE_ID,
            "cultivar": "Ivory",
            "target": "lai_decline_after_peak",
            "date": obs_date("Ivory", 70),
            "das": 70,
            "confidence": "low",
            "source_pdf": PDF_AGRONOMY.name,
            "source_location": "Figure 3 and results/discussion text",
            "interpretation": "Digitized LAI and leaf area decline after the mid-season peak; this is an indirect development/senescence signal.",
            "usable_for_calibration": "growth-shape constraint only",
        },
        {
            "site_id": SITE_ID,
            "cultivar": "Santhica 27",
            "target": "late_lai_collapse",
            "date": obs_date("Santhica 27", 121),
            "das": 121,
            "confidence": "low",
            "source_pdf": PDF_AGRONOMY.name,
            "source_location": "Figure 3",
            "interpretation": "Digitized LAI drops strongly after about 113-121 DAS, indicating late canopy senescence.",
            "usable_for_calibration": "growth-shape constraint only",
        },
        {
            "site_id": SITE_ID,
            "cultivar": "both",
            "target": "flowering_date_reported_directly",
            "date": "",
            "das": np.nan,
            "confidence": "none",
            "source_pdf": f"{PDF_AGRONOMY.name}; {PDF_WATER.name}",
            "source_location": "full-text keyword audit",
            "interpretation": "No direct observed flowering dates or BBCH-style flowering stage dates were found in either PDF.",
            "usable_for_calibration": "no direct flowering target available",
        },
    ]
    return pd.DataFrame(rows)


def write_phenology_markdown(df: pd.DataFrame) -> None:
    path = DERIVED / "phenology_targets_from_text.md"
    lines = [
        "# Germany Hemp Phenology Targets From Paper Text",
        "",
        "The PDFs do not report a direct measured flowering date for Santhica 27 or Ivory.",
        "The entries below are weak inferred targets from development and senescence statements.",
        "",
    ]
    for _, row in df.iterrows():
        lines.extend(
            [
                f"## {row['cultivar']} - {row['target']}",
                f"- Date/DAS: {row['date'] or 'not reported'} / {'' if pd.isna(row['das']) else row['das']}",
                f"- Confidence: {row['confidence']}",
                f"- Source: {row['source_pdf']}, {row['source_location']}",
                f"- Interpretation: {row['interpretation']}",
                f"- Calibration use: {row['usable_for_calibration']}",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_weather_from_dwd() -> Path | None:
    ensure_import_paths()
    from dssatutils.weather_dwd import (
        _build_point_frame,
        _dwd_stations,
        _fetch_station,
        _historical_index,
        _write_wth,
        process_weather_dwd,
    )

    WEATHER_DIR.mkdir(parents=True, exist_ok=True)
    DWD_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    log_file = DERIVED / "dwd_weather_errors.log"
    paper_station_id = "03987"

    try:
        stations = _dwd_stations(str(DWD_CACHE_DIR))
        hist_idx = _historical_index(str(DWD_CACHE_DIR))
        station = stations[stations["station_id"] == paper_station_id].iloc[0]
        daily = _fetch_station(paper_station_id, hist_idx.get(paper_station_id), str(DWD_CACHE_DIR))
        frame = _build_point_frame(daily, float(station["lat"]), 2018, 2018)
        if not frame.empty:
            optional = ["SRAD", "RAIN", "TDEW", "RH2M", "WIND"]
            frame[optional] = frame[optional].fillna(-99)
            _write_wth(
                frame,
                SITE_ID,
                float(station["lat"]),
                float(station["lon"]),
                float(station["elev"]),
                str(WEATHER_DIR),
            )
            pd.DataFrame(
                [
                    {
                        "site_id": SITE_ID,
                        "selected_station_id": paper_station_id,
                        "selected_station_name": "Potsdam",
                        "station_latitude": float(station["lat"]),
                        "station_longitude": float(station["lon"]),
                        "station_elevation_m": float(station["elev"]),
                        "selection_reason": "Paper aggregate cites DWD Potsdam station 03987; forced to match published weather aggregate.",
                    }
                ]
            ).to_csv(DERIVED / "weather_station_selection.csv", index=False)
            wth = WEATHER_DIR / f"{SITE_ID}.WTH"
            return wth if wth.exists() else None
    except Exception as exc:  # noqa: BLE001
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write(f"Forced DWD station {paper_station_id} failed: {exc}\n")

    point = pd.DataFrame(
        [
            {
                "ID": SITE_ID,
                "lat": LAT,
                "lon": LON,
            }
        ]
    )
    process_weather_dwd(
        point,
        start_year=2018,
        end_year=2018,
        output_dir=str(WEATHER_DIR),
        id_col="ID",
        lat_col="lat",
        lon_col="lon",
        n_cores=1,
        log_file=str(log_file),
        dwd_cache_dir=str(DWD_CACHE_DIR),
        max_station_km=70.0,
    )
    wth = WEATHER_DIR / f"{SITE_ID}.WTH"
    return wth if wth.exists() else None


def read_wth(path: Path) -> pd.DataFrame:
    rows = []
    header = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("@  DATE"):
            header = line.split()
            continue
        if header and line.strip() and not line.startswith("*") and not line.startswith("$"):
            parts = line.split()
            if len(parts) >= 8:
                yyyydoy = parts[0]
                yr = int(yyyydoy[:4])
                doy = int(yyyydoy[4:])
                dt = date(yr, 1, 1) + timedelta(days=doy - 1)
                rows.append(
                    {
                        "date": dt,
                        "SRAD": float(parts[1]),
                        "TMAX": float(parts[2]),
                        "TMIN": float(parts[3]),
                        "RAIN": float(parts[4]),
                        "TDEW": float(parts[5]),
                        "RH2M": float(parts[6]),
                        "WIND": float(parts[7]),
                    }
                )
    return pd.DataFrame(rows)


def weather_comparison(wth: Path | None) -> pd.DataFrame:
    paper_rows = [
        ("annual_2018", date(2018, 1, 1), date(2018, 12, 31), "RAIN_sum_mm", 346.0, "Water Table A2, DWD Potsdam station annual precipitation 2018"),
        ("may_sep_2018", date(2018, 5, 1), date(2018, 9, 30), "RAIN_sum_mm", 128.0, "Water Table A2, DWD Potsdam station May-Sep precipitation 2018"),
        ("annual_2018", date(2018, 1, 1), date(2018, 12, 31), "TAVG_mean_C", 11.3, "Water Table A2, DWD Potsdam station annual mean air temperature 2018"),
        ("annual_2018", date(2018, 1, 1), date(2018, 12, 31), "WIND_mean_m_s", 4.1, "Water Table A2, DWD Potsdam station annual mean wind 2018"),
        ("santhica_sowing_to_harvest", date(2018, 5, 4), date(2018, 9, 19), "RAIN_sum_mm", 56.0, "Agronomy text, on-site precipitation during entire growing season"),
        ("water_flux_period", date(2018, 6, 15), date(2018, 9, 9), "RAIN_sum_mm", 44.0, "Water Table 6/abstract, measured plot precipitation during flux period"),
    ]
    rows = []
    if wth is None:
        for period, start, end, metric, paper_value, source in paper_rows:
            rows.append(
                {
                    "period": period,
                    "start_date": start.isoformat(),
                    "end_date": end.isoformat(),
                    "metric": metric,
                    "dwd_value": np.nan,
                    "paper_value": paper_value,
                    "difference_dwd_minus_paper": np.nan,
                    "source": source,
                    "notes": "DWD WTH was not created.",
                }
            )
        return pd.DataFrame(rows)

    df = read_wth(wth)
    df["TAVG"] = (df["TMAX"] + df["TMIN"]) / 2.0
    for period, start, end, metric, paper_value, source in paper_rows:
        sub = df[(df["date"] >= start) & (df["date"] <= end)].copy()
        if metric == "RAIN_sum_mm":
            value = sub.loc[sub["RAIN"] > -90, "RAIN"].sum()
        elif metric == "TAVG_mean_C":
            value = sub.loc[sub["TAVG"] > -90, "TAVG"].mean()
        elif metric == "WIND_mean_m_s":
            value = sub.loc[sub["WIND"] > -90, "WIND"].mean()
        else:
            value = np.nan
        rows.append(
            {
                "period": period,
                "start_date": start.isoformat(),
                "end_date": end.isoformat(),
                "metric": metric,
                "dwd_value": round(float(value), 2) if pd.notna(value) else np.nan,
                "paper_value": paper_value,
                "difference_dwd_minus_paper": round(float(value - paper_value), 2)
                if pd.notna(value)
                else np.nan,
                "source": source,
                "notes": "DWD source may use station coordinates and rain gauge; on-site plot aggregates can differ.",
            }
        )
    return pd.DataFrame(rows)


def fetch_soilgrids_profile() -> pd.DataFrame | None:
    ensure_import_paths()
    from dssatutils.soil_soilgrids_online import _calculate_soil_physics, _fetch_soilgrids_rest

    raw = _fetch_soilgrids_rest(LAT, LON)
    if raw is None or raw.empty:
        return None

    wide = raw.pivot_table(
        index=["depth_label", "depth_bottom", "depth_center"],
        columns="prop",
        values="value",
        aggfunc="first",
    ).reset_index()

    for col in ["clay", "sand", "silt", "soc", "bdod", "cfvo"]:
        if col not in wide:
            wide[col] = np.nan

    wide["ID"] = SITE_ID
    wide["latitude"] = LAT
    wide["longitude"] = LON
    wide["clay"] = wide["clay"] / 10.0
    wide["sand"] = wide["sand"] / 10.0
    wide["silt"] = wide["silt"] / 10.0
    wide["soc_pct"] = wide["soc"] / 100.0
    wide["om_pct"] = wide["soc_pct"] * 1.724
    wide["bdod"] = wide["bdod"] / 100.0
    wide["cfvo"] = wide["cfvo"] / 10.0
    physics = wide.apply(
        lambda r: pd.Series(
            _calculate_soil_physics(
                float(r["sand"]) if pd.notna(r["sand"]) else 79.5,
                float(r["clay"]) if pd.notna(r["clay"]) else 9.0,
                float(r["om_pct"]) if pd.notna(r["om_pct"]) else 1.8,
            )
        ),
        axis=1,
    )
    return pd.concat([wide, physics], axis=1).sort_values("depth_bottom")


def texture_ksat_cm_h(sand: float, clay: float) -> float:
    return min(999.0, 60.96 * (10 ** (0.0126 * sand - 0.0064 * clay - 0.6)))


def reported_approx_soil_profile() -> pd.DataFrame:
    depths = [
        ("0-5cm", 5, 2.5),
        ("5-15cm", 15, 10.0),
        ("15-30cm", 30, 22.5),
        ("30-60cm", 60, 45.0),
        ("60-100cm", 100, 80.0),
        ("100-200cm", 200, 150.0),
    ]
    # Reported ranges: sand 71-88%, silt 8-15%, clay 4-14%, Corg 0.4-1.7%,
    # pH 6.6-8.4, bulk density 1.44 +/- 0.09 g cm-3, PWP 8.2%, FC 26.5%.
    rows = []
    for depth_label, bottom, center in depths:
        sand = 79.5
        silt = 11.5
        clay = 9.0
        sbdm = 1.44 if bottom <= 30 else 1.50
        soc_pct = 1.05 if bottom <= 30 else 0.55
        ssat = max(0.35, min(0.55, 1.0 - sbdm / 2.65))
        rows.append(
            {
                "ID": SITE_ID,
                "latitude": LAT,
                "longitude": LON,
                "depth_label": depth_label,
                "depth_bottom": bottom,
                "depth_center": center,
                "SLLL": 0.082,
                "SDUL": 0.265,
                "SSAT": ssat,
                "SRGF": max(0.0, math.exp(-0.02 * center)),
                "SSKS": texture_ksat_cm_h(sand, clay),
                "bdod": sbdm,
                "soc_pct": soc_pct,
                "om_pct": soc_pct * 1.724,
                "clay": clay,
                "silt": silt,
                "sand": sand,
                "cfvo": 0.0,
                "ph": 7.5,
                "source": "Approximation from reported soil properties in the PDFs; homogeneous weakly loamy sand profile.",
            }
        )
    return pd.DataFrame(rows)


def write_sol(df: pd.DataFrame, path: Path, source_name: str) -> None:
    lines = [
        f"*SOILS: {source_name}",
        f"! Site: {SITE_NAME}; coordinates {LAT:.6f}, {LON:.6f}",
        "",
        f"*{SITE_ID[:10]:<10s}  GERMANYHEMP  {LAT:9.3f} {LON:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {SITE_ID[:11]:<11s} Germany       {LAT:9.3f} {LON:9.3f} ",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    ]
    for _, layer in df.sort_values("depth_bottom").iterrows():
        center = float(layer["depth_center"])
        srgf = float(layer["SRGF"]) if "SRGF" in layer and pd.notna(layer["SRGF"]) else max(0.0, math.exp(-0.02 * center))
        ssks = float(layer["SSKS"]) if "SSKS" in layer and pd.notna(layer["SSKS"]) else texture_ksat_cm_h(float(layer["sand"]), float(layer["clay"]))
        lines.append(
            f"{int(layer['depth_bottom']):6d}   -99"
            f" {float(layer['SLLL']):5.3f} {float(layer['SDUL']):5.3f} {float(layer['SSAT']):5.3f}"
            f" {srgf:5.2f} {ssks:5.1f} {float(layer['bdod']):5.2f} {float(layer['soc_pct']):5.2f}"
            f" {float(layer['clay']):5.1f} {float(layer['silt']):5.1f} {float(layer['cfvo']):5.1f}"
            "   -99   -99   -99   -99   -99"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def soil_outputs() -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame]:
    SOIL_DIR.mkdir(parents=True, exist_ok=True)
    reported = reported_approx_soil_profile()
    reported.to_csv(SOIL_DIR / "reported_approx_soil_profile.csv", index=False)
    write_sol(reported, SOIL_DIR / f"{SITE_ID}_reported_approx.SOL", "Reported approximate German hemp soil")

    soilgrids = fetch_soilgrids_profile()
    if soilgrids is not None:
        soilgrids.to_csv(SOIL_DIR / "soilgrids_profile.csv", index=False)
        write_sol(soilgrids, SOIL_DIR / f"{SITE_ID}_soilgrids_rest.SOL", "ISRIC SoilGrids 2.0 REST")

    if soilgrids is None:
        comparison = reported.copy()
        comparison["soilgrids_available"] = False
    else:
        keep = ["depth_bottom", "SLLL", "SDUL", "SSAT", "bdod", "soc_pct", "clay", "silt", "sand", "cfvo"]
        comparison = reported[keep].merge(
            soilgrids[keep],
            on="depth_bottom",
            how="outer",
            suffixes=("_reported_approx", "_soilgrids"),
        )
        for col in ["SLLL", "SDUL", "SSAT", "bdod", "soc_pct", "clay", "silt", "sand", "cfvo"]:
            comparison[f"{col}_difference_reported_minus_soilgrids"] = (
                comparison[f"{col}_reported_approx"] - comparison[f"{col}_soilgrids"]
            )
    comparison.to_csv(SOIL_DIR / "soil_profile_comparison.csv", index=False)
    return reported, soilgrids, comparison


def write_dataset_readme(wth: Path | None, soilgrids: pd.DataFrame | None) -> None:
    readme = DERIVED / "README_germany_hemp_derived_dataset.md"
    lines = [
        "# Germany Hemp Derived Dataset",
        "",
        "This directory contains an approximate calibration dataset derived from the two local PDFs.",
        "It is suitable for scaffold/testing and weak cultivar calibration, not as a full raw experimental dataset.",
        "",
        "## Key limits",
        "",
        "- Growth curves are manual visual digitizations from printed Figure 3 means.",
        "- Error bars from the figure are not encoded as per-observation uncertainty.",
        "- Weather is reconstructed from the nearest suitable DWD daily climate station, not the field logger.",
        "- Soil is an approximate DSSAT profile from reported bulk soil properties and a SoilGrids comparison.",
        "- No direct observed flowering dates were found in either PDF.",
        "",
        "## Main files",
        "",
        "- `digitized_growth_observations.csv`: height, leaf area per plant, plant density, and LAI.",
        "- `dssatcalibrator_observations_long.csv`: DSSAT-mappable observations in the native long schema.",
        "- `scalar_observations.csv`: yield, emergence percent, and stand-density scalar values.",
        "- `management_summary.csv`: sowing, fertilization, irrigation, harvest, coordinates.",
        "- `filex_management_inputs.csv`: FileX-oriented management scaffold with reported values and assumptions.",
        "- `filex_management_extractability.md`: audit of what is directly reported vs still assumed for FileX.",
        "- `weather_summary_comparison.csv`: DWD reconstructed aggregates vs paper aggregates.",
        "- `soil/reported_approx_soil_profile.csv` and `.SOL`: approximate profile from reported properties.",
        "- `soil/soilgrids_profile.csv` and `.SOL`: SoilGrids REST profile, when network extraction succeeds.",
        "- `phenology_targets_from_text.csv` and `.md`: weak inferred phenology/senescence targets.",
        "",
        "## Extraction status",
        "",
        f"- DWD weather file created: {'yes' if wth is not None else 'no'}",
        f"- SoilGrids REST profile created: {'yes' if soilgrids is not None else 'no'}",
        "",
    ]
    readme.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    DERIVED.mkdir(parents=True, exist_ok=True)
    render_figure3_source()

    growth = digitized_growth_observations()
    growth.to_csv(DERIVED / "digitized_growth_observations.csv", index=False)

    scalars = scalar_observations()
    scalars.to_csv(DERIVED / "scalar_observations.csv", index=False)
    calibrator_observation_table(growth, scalars).to_csv(
        DERIVED / "dssatcalibrator_observations_long.csv",
        index=False,
    )

    management_summary().to_csv(DERIVED / "management_summary.csv", index=False)
    filex_inputs = filex_management_inputs()
    filex_inputs.to_csv(DERIVED / "filex_management_inputs.csv", index=False)
    write_filex_management_audit(filex_inputs)

    phen = phenology_targets()
    phen.to_csv(DERIVED / "phenology_targets_from_text.csv", index=False)
    write_phenology_markdown(phen)

    wth = build_weather_from_dwd()
    weather_comparison(wth).to_csv(DERIVED / "weather_summary_comparison.csv", index=False)

    _, soilgrids, _ = soil_outputs()
    write_dataset_readme(wth, soilgrids)

    print(f"Wrote derived Germany hemp dataset to: {DERIVED}")


if __name__ == "__main__":
    main()
