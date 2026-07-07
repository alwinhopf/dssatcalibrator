from __future__ import annotations

import json
import os
import re
import shutil
import sys
from copy import deepcopy
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from dssatcalibrator import objective as obj, orchestrator, viz  # noqa: E402
from dssatcalibrator.config import load_config, resolve_dssat_paths  # noqa: E402
from dssatcalibrator.orchestrator import CalibrationResult  # noqa: E402
from dssatcalibrator.spaces import ParameterSpace  # noqa: E402
from dssatcalibrator.writers import edit_species, read_cultivar_values, read_ecotype_values  # noqa: E402


DATE_TAG = os.environ.get("DSSATCAL_RUN_TAG", "20260704_lai_initial")
RESULTS_DIR = Path("results/china_hemp_calibration")
FIGURES_DIR = Path("figures/china_hemp_calibration")
CULTIVARS = ["IB0002", "IB0008"]
SOURCE_ECOTYPES = {"IB0002": "HM0002", "IB0008": "HM0003"}
NEW_IDS = {
    "IB0002": ("IBCL02", "HMLA02", "NWG2730_CN_LAI"),
    "IB0008": ("IBCL08", "HMLA08", "YUNMA8_CN_LAI"),
}
FALLBACK_GENOTYPE_VALUES = {
    "XFRT": 0.6,
    "WTPSD": 0.035,
    "SFDUR": 20.0,
    "SDPDV": 1.0,
    "PODUR": 30.0,
    "THRSH": 80.0,
    "SDPRO": 0.256,
    "SDLIP": 0.319,
    "PM06": 0.0,
    "PM09": 0.8,
    "LNGSH": 5.0,
    "R7-R8": 20.0,
}

DEFAULT_STAGE1_THETA = (
    RESULTS_DIR
    / "china_hemp_long_stage1_phenology_node_20260704_split3_allcores"
    / "stage1_best_theta.json"
)
DEFAULT_STAGE2_THETA = (
    RESULTS_DIR
    / "china_hemp_long_stage2_biomass_protected_20260704_split3_allcores"
    / "stage2_best_theta.json"
)

STAGE1_PHENO_PARAMS = {
    "genetic_cultivar": ["CSDL", "PPSEN", "EM-FL", "FL-SH", "FL-SD", "SD-PM"],
    "genetic_ecotype": ["THVAR", "PL-EM", "EM-V1", "V1-JU", "JU-R0", "R1PPO", "OPTBI", "SLOBI"],
}
LAI_PARAMS = {
    "genetic_cultivar": ["FL-LF", "LFMAX", "SLAVR", "SIZLF"],
    "genetic_ecotype": ["FL-VS", "TRIFL", "RWDTH", "RHGHT"],
}
CAPACITY_PARAMS = {
    "genetic_cultivar": ["XFRT", "WTPSD", "SFDUR", "SDPDV", "PODUR", "THRSH", "SDPRO"],
    "genetic_ecotype": ["PM09", "LNGSH", "R7-R8"],
}
SPECIES_PARAMS = {
    "genetic_species": [
        "PGEFF", "KDIF", "SLWREF", "LNREF", "PGREF",
        "FINREF", "SLAREF", "SIZREF", "VSSINK",
        "SLAMAX", "SLAMIN", "TURSLA", "NSLA",
        "YVREF_V3", "YVREF_V4", "YVREF_V5",
        "YLEAF_MID", "YLEAF_LATE", "YSTEM_MID", "YSTEM_LATE",
        "FRSTMF", "FRLFF", "FRLFMX",
        "YVSHT_7", "YVSHT_10", "YVSHT_13",
        "YVSWH_7", "YVSWH_10", "YVSWH_13",
        "NHGT",
    ],
}


def _merged_param_groups(*groups: dict[str, list[str]]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for group in groups:
        for name, values in group.items():
            out.setdefault(name, [])
            for value in values:
                if value not in out[name]:
                    out[name].append(value)
    return out


def _read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def _theta_value(theta: dict, name: str, cultivar: str, fallback: float) -> float:
    return float(theta.get(f"{name}__{cultivar}", theta.get(name, fallback)))


def _common(cfg: dict, name: str) -> dict:
    cfg = deepcopy(cfg)
    cfg["calibrator"]["name"] = name
    cfg["calibrator"]["results_dir"] = str(RESULTS_DIR)
    cfg["calibrator"]["figures_dir"] = str(FIGURES_DIR)
    cfg["calibrator"]["workdir"] = "results/_workdir_lai_initial"
    cfg["calibrator"]["cache_spawns"] = False
    cfg["calibrator"]["keep_run_dirs"] = False
    cfg["calibrator"]["num_cores"] = int(os.environ.get("DSSATCAL_LONG_CORES", os.cpu_count() or 1))
    cfg["calibrator"]["batch_size"] = int(os.environ.get("DSSATCAL_BATCH_SIZE", 24))
    cfg["calibrator"]["spawn_timeout"] = max(90, int(cfg["calibrator"].get("spawn_timeout", 60)))
    cfg["method"].setdefault("sensitivity", {})["active"] = False
    cfg["method"].setdefault("select", {})["active"] = False
    cfg["method"].setdefault("surrogate", {})["active"] = False
    cfg["method"].setdefault("multiobjective", {})["engine"] = "none"
    cfg["method"]["preset"] = "B"
    cfg["method"].setdefault("bayesian", {})["engine"] = "none"
    cfg["crops"][0]["calibration_cultivars"] = list(CULTIVARS)
    return cfg


def _deactivate_all_parameters(cfg: dict) -> None:
    for params in (cfg.get("parameters") or {}).values():
        if not isinstance(params, dict):
            continue
        for spec in params.values():
            if isinstance(spec, dict):
                spec["active"] = False


def _set_cultivar_fixed(cfg: dict, theta: dict, names_by_group: dict[str, list[str]]) -> None:
    for group, names in names_by_group.items():
        params = cfg["parameters"].get(group, {})
        for name in names:
            if name not in params:
                continue
            spec = params[name]
            spec["active"] = False
            spec["fixed"] = True
            spec["scope"] = "cultivar"
            spec.pop("cultivars", None)
            fallback = float(spec.get("start", 0.5 * (float(spec["min"]) + float(spec["max"]))))
            starts = {cultivar: _theta_value(theta, name, cultivar, fallback) for cultivar in CULTIVARS}
            spec["start_by_cultivar"] = starts
            spec["start"] = starts.get("IB0008", next(iter(starts.values())))


def _activate_lai_parameters(cfg: dict, lai_theta: dict | None = None) -> None:
    bounds = {
        "FL-LF": (5.0, 110.0),
        "LFMAX": (0.3, 3.0),
        "SLAVR": (80.0, 800.0),
        "SIZLF": (50.0, 900.0),
        "FL-VS": (10.0, 110.0),
        "TRIFL": (0.08, 0.55),
        "RWDTH": (0.4, 5.0),
        "RHGHT": (0.4, 8.0),
    }
    for group, names in LAI_PARAMS.items():
        for name in names:
            spec = cfg["parameters"][group][name]
            lo, hi = bounds[name]
            spec.update({"active": True, "fixed": False, "scope": "cultivar", "min": lo, "max": hi})
            spec.pop("cultivars", None)
            if lai_theta:
                fallback = float(spec.get("start", 0.5 * (lo + hi)))
                spec["start_by_cultivar"] = {
                    cultivar: _theta_value(lai_theta, name, cultivar, fallback)
                    for cultivar in CULTIVARS
                }
                spec["start"] = spec["start_by_cultivar"].get("IB0008", fallback)


def _lai_objective(cfg: dict, *, include_biomass: bool = False) -> None:
    cfg["engine"].setdefault("timeseries_outputs", {}).update({
        "biomass": "CWAD",
        "leaf": "LWAD",
        "stem": "SWAD",
        "height": "CHTD",
        "width": "CWID",
        "node_stage": "L#SD",
        "LAI": "LAID",
    })
    cfg["engine"].setdefault("scalar_outputs", {}).update({
        "emergence": "EDAP",
        "anthesis": "ADAP",
    })
    cfg["objective"]["obs_autocorr"] = True
    cfg["objective"]["ignore_zero_observations"] = ["width"]
    cfg["objective"]["weights"] = {
        "LAI": 3.0,
        "leaf": 1.3,
        "node_stage": 1.2,
        "height": 1.0,
        "width": 1.0,
        "emergence": 0.6,
        "anthesis": 0.8,
        "biomass": 0.8 if include_biomass else 0.35,
        "stem": 0.7 if include_biomass else 0.25,
    }
    cfg["objective"]["error_model"] = {
        "LAI": {"type": "absolute", "value": 0.8},
        "leaf": {"type": "relative", "value": 0.30},
        "node_stage": {"type": "absolute", "value": 1.0},
        "height": {"type": "absolute", "value": 0.35},
        "width": {"type": "absolute", "value": 0.20},
        "emergence": {"type": "absolute", "value": 3.0},
        "anthesis": {"type": "absolute", "value": 3.0},
        "biomass": {"type": "relative", "value": 0.30},
        "stem": {"type": "relative", "value": 0.35},
    }


def _all_observed_output_mapping(cfg: dict) -> None:
    cfg["engine"].setdefault("timeseries_outputs", {}).update({
        "CHTD": "CHTD",
        "CNAD": "CNAD",
        "CWAD": "CWAD",
        "CWID": "CWID",
        "G#AD": "G#AD",
        "GN%D": "GN%D",
        "GNAD": "GNAD",
        "GWAD": "GWAD",
        "GWGD": "GWGD",
        "HIAD": "HIAD",
        "HIPD": "HIPD",
        "L#SD": "L#SD",
        "LAID": "LAID",
        "LI%N": "LI%N",
        "LN%D": "LN%D",
        "LNAD": "LNAD",
        "LWAD": "LWAD",
        "NFXD": "NFXD",
        "NWAD": "NWAD",
        "P#AD": "P#AD",
        "PWAD": "PWAD",
        "SH%D": "SH%D",
        "SHAD": "SHAD",
        "SHND": "SHND",
        "SLAD": "SLAD",
        "SN%D": "SN%D",
        "SNAD": "SNAD",
        "SNHD": "SNHD",
        "SWAD": "SWAD",
        "VNAD": "VNAD",
        "grain": "GWAD",
    })


def _capacity_objective(cfg: dict) -> None:
    _lai_objective(cfg, include_biomass=True)
    _all_observed_output_mapping(cfg)
    weights = cfg["objective"]["weights"]
    weights.update({
        "biomass": 1.4,
        "CWAD": 1.4,
        "stem": 1.3,
        "SWAD": 1.3,
        "leaf": 1.0,
        "LWAD": 1.0,
        "LAI": 1.4,
        "LAID": 1.4,
        "height": 0.9,
        "CHTD": 0.9,
        "width": 0.6,
        "CWID": 0.6,
        "grain": 0.6,
        "GWAD": 0.6,
        "GWGD": 0.4,
        "GNAD": 0.35,
        "GN%D": 0.35,
        "LN%D": 0.35,
        "LI%N": 0.35,
        "SN%D": 0.35,
        "SLAD": 0.35,
        "SHAD": 0.35,
        "PWAD": 0.35,
        "HIAD": 0.25,
        "HIPD": 0.25,
        "node_stage": 0.8,
        "L#SD": 0.8,
    })
    errors = cfg["objective"]["error_model"]
    for var in (
        "CWAD", "SWAD", "LWAD", "GWAD", "GWGD", "GNAD", "CNAD", "LNAD", "SNAD",
        "PWAD", "SHAD", "NWAD", "VNAD", "SLAD",
    ):
        errors[var] = {"type": "relative", "value": 0.35}
    for var in ("GN%D", "LN%D", "LI%N", "SN%D", "SH%D", "HIAD", "HIPD"):
        errors[var] = {"type": "absolute", "value": 1.0}
    errors["G#AD"] = {"type": "absolute", "value": 80.0}
    errors["P#AD"] = {"type": "absolute", "value": 80.0}
    errors["SHND"] = {"type": "absolute", "value": 80.0}


def _activate_capacity_parameters(cfg: dict, genotype_theta: dict | None = None) -> None:
    bounds = {
        "XFRT": (0.15, 0.85),
        "WTPSD": (0.01, 0.8),
        "SFDUR": (5.0, 60.0),
        "SDPDV": (1.0, 80.0),
        "PODUR": (5.0, 60.0),
        "THRSH": (55.0, 90.0),
        "SDPRO": (0.12, 0.35),
        "PM09": (0.2, 1.0),
        "LNGSH": (2.0, 15.0),
        "R7-R8": (1.0, 30.0),
    }
    for group, names in CAPACITY_PARAMS.items():
        for name in names:
            spec = cfg["parameters"][group][name]
            lo, hi = bounds[name]
            spec.update({"active": True, "fixed": False, "scope": "cultivar", "min": lo, "max": hi})
            spec.pop("cultivars", None)
            if genotype_theta:
                fallback = float(spec.get("start", 0.5 * (lo + hi)))
                spec["start_by_cultivar"] = {
                    cultivar: _theta_value(genotype_theta, name, cultivar, fallback)
                    for cultivar in CULTIVARS
                }
                spec["start"] = spec["start_by_cultivar"].get("IB0008", fallback)


def _add_species_parameter(
    cfg: dict,
    name: str,
    *,
    lo: float,
    hi: float,
    start: float,
    spe_key: str,
    spe_index: int,
) -> None:
    params = cfg["parameters"].setdefault("genetic_species", {})
    params[name] = {
        "active": True,
        "min": lo,
        "max": hi,
        "start": start,
        "scope": "global",
        "spe_key": spe_key,
        "spe_index": int(spe_index),
        "prior": {"dist": "uniform"},
    }


def _activate_species_canopy_source_parameters(cfg: dict) -> None:
    """Curated .SPE expansion for canopy/source diagnostics.

    These are global species parameters by design. Cultivar/ecotype parameters
    remain cultivar-specific; species edits test whether the hemp module's shared
    leaf expansion, source capacity, partitioning, and geometry curves are
    structurally limiting all China trials.
    """
    _add_species_parameter(cfg, "PGEFF", lo=0.045, hi=0.070, start=0.0541,
                           spe_key="PGEFF SCV KDIF, LFANGB", spe_index=0)
    _add_species_parameter(cfg, "KDIF", lo=0.65, hi=1.00, start=0.80,
                           spe_key="PGEFF SCV KDIF, LFANGB", spe_index=2)
    _add_species_parameter(cfg, "SLWREF", lo=0.0032, hi=0.0058, start=0.0046,
                           spe_key="SLWREF,SLWSLO,NSLOPE,LNREF,PGREF", spe_index=0)
    _add_species_parameter(cfg, "LNREF", lo=4.0, hi=5.8, start=4.90,
                           spe_key="SLWREF,SLWSLO,NSLOPE,LNREF,PGREF", spe_index=3)
    _add_species_parameter(cfg, "PGREF", lo=0.85, hi=1.35, start=1.03,
                           spe_key="SLWREF,SLWSLO,NSLOPE,LNREF,PGREF", spe_index=4)

    _add_species_parameter(cfg, "FINREF", lo=180.0, hi=320.0, start=220.0,
                           spe_key="FINREF,SLAREF,SIZREF,VSSINK,EVMODC", spe_index=0)
    _add_species_parameter(cfg, "SLAREF", lo=180.0, hi=420.0, start=240.0,
                           spe_key="FINREF,SLAREF,SIZREF,VSSINK,EVMODC", spe_index=1)
    _add_species_parameter(cfg, "SIZREF", lo=170.0, hi=480.0, start=220.0,
                           spe_key="FINREF,SLAREF,SIZREF,VSSINK,EVMODC", spe_index=2)
    _add_species_parameter(cfg, "VSSINK", lo=2.4, hi=5.5, start=3.0,
                           spe_key="FINREF,SLAREF,SIZREF,VSSINK,EVMODC", spe_index=3)
    _add_species_parameter(cfg, "SLAMAX", lo=300.0, hi=650.0, start=350.0,
                           spe_key="SLAMAX,SLAMIN,SLAPAR,TURSLA,NSLA", spe_index=0)
    _add_species_parameter(cfg, "SLAMIN", lo=150.0, hi=300.0, start=200.0,
                           spe_key="SLAMAX,SLAMIN,SLAPAR,TURSLA,NSLA", spe_index=1)
    _add_species_parameter(cfg, "TURSLA", lo=0.8, hi=2.4, start=1.5,
                           spe_key="SLAMAX,SLAMIN,SLAPAR,TURSLA,NSLA", spe_index=3)
    _add_species_parameter(cfg, "NSLA", lo=0.2, hi=2.2, start=1.0,
                           spe_key="SLAMAX,SLAMIN,SLAPAR,TURSLA,NSLA", spe_index=4)

    _add_species_parameter(cfg, "YVREF_V3", lo=90.0, hi=220.0, start=110.0,
                           spe_key="YVREF(1-6), LEAF AREA VALUES", spe_index=3)
    _add_species_parameter(cfg, "YVREF_V4", lo=170.0, hi=500.0, start=200.0,
                           spe_key="YVREF(1-6), LEAF AREA VALUES", spe_index=4)
    _add_species_parameter(cfg, "YVREF_V5", lo=260.0, hi=850.0, start=320.0,
                           spe_key="YVREF(1-6), LEAF AREA VALUES", spe_index=5)

    _add_species_parameter(cfg, "YLEAF_MID", lo=0.30, hi=0.52, start=0.40,
                           spe_key="YLEAF VALUES", spe_index=3)
    _add_species_parameter(cfg, "YLEAF_LATE", lo=0.25, hi=0.50, start=0.40,
                           spe_key="YLEAF VALUES", spe_index=6)
    _add_species_parameter(cfg, "YSTEM_MID", lo=0.25, hi=0.50, start=0.35,
                           spe_key="YSTEM VALUES", spe_index=3)
    _add_species_parameter(cfg, "YSTEM_LATE", lo=0.25, hi=0.55, start=0.35,
                           spe_key="YSTEM VALUES", spe_index=6)
    _add_species_parameter(cfg, "FRSTMF", lo=0.35, hi=0.70, start=0.45,
                           spe_key="WTFSD,PORPT,FRSTMF,FRLFF,ATOP,FRCNOD", spe_index=2)
    _add_species_parameter(cfg, "FRLFF", lo=0.18, hi=0.50, start=0.35,
                           spe_key="WTFSD,PORPT,FRSTMF,FRLFF,ATOP,FRCNOD", spe_index=3)
    _add_species_parameter(cfg, "FRLFMX", lo=0.60, hi=0.90, start=0.80,
                           spe_key="FRLFMX", spe_index=0)

    _add_species_parameter(cfg, "YVSHT_7", lo=0.08, hi=0.20, start=0.13,
                           spe_key="YVSHT(1-10)", spe_index=3)
    _add_species_parameter(cfg, "YVSHT_10", lo=0.09, hi=0.23, start=0.14,
                           spe_key="YVSHT(1-10)", spe_index=4)
    _add_species_parameter(cfg, "YVSHT_13", lo=0.10, hi=0.26, start=0.16,
                           spe_key="YVSHT(1-10)", spe_index=5)
    _add_species_parameter(cfg, "YVSWH_7", lo=0.05, hi=0.12, start=0.06,
                           spe_key="YVSWH(1-10)", spe_index=3)
    _add_species_parameter(cfg, "YVSWH_10", lo=0.05, hi=0.14, start=0.07,
                           spe_key="YVSWH(1-10)", spe_index=4)
    _add_species_parameter(cfg, "YVSWH_13", lo=0.05, hi=0.14, start=0.06,
                           spe_key="YVSWH(1-10)", spe_index=5)
    _add_species_parameter(cfg, "NHGT", lo=2.0, hi=7.0, start=5.0,
                           spe_key="NHGT", spe_index=0)


def configure_lai_stage(stage1_theta: dict, stage2_theta: dict | None = None) -> dict:
    cfg = load_config("calibration_china_hemp/stage2_biomass_after_phenology_cultivar.yaml")
    cfg = _common(cfg, f"china_hemp_lai_leafexp_protected_{DATE_TAG}")
    _deactivate_all_parameters(cfg)
    _set_cultivar_fixed(cfg, stage1_theta, STAGE1_PHENO_PARAMS)
    _activate_lai_parameters(cfg, stage2_theta)
    _lai_objective(cfg, include_biomass=False)
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_LAI_MAXITER", 8)),
        "popsize": int(os.environ.get("DSSATCAL_LAI_POPSIZE", 12)),
        "restarts": 1,
    })
    return cfg


def configure_resource_diagnostic(stage1_theta: dict, stage2_theta: dict, *, rich: bool) -> tuple[dict, dict]:
    cfg = configure_lai_stage(stage1_theta, stage2_theta)
    cfg["calibrator"]["name"] = f"china_hemp_resource_{'rich' if rich else 'baseline'}_{DATE_TAG}"
    fixed = {
        **stage2_theta,
    }
    if rich:
        init = cfg["parameters"].setdefault("initial_conditions", {})
        init["initial_soil_water_mult"] = {
            "active": False,
            "fixed": True,
            "min": 1.0,
            "max": 1.6,
            "start": 1.45,
            "dssat": "SH2O",
            "op": "mult",
            "clip_01": True,
        }
        init["initial_nh4_mult"] = {
            "active": False,
            "fixed": True,
            "min": 1.0,
            "max": 5.0,
            "start": 3.0,
            "dssat": "SNH4",
            "op": "mult",
        }
        init["initial_no3_mult"] = {
            "active": False,
            "fixed": True,
            "min": 1.0,
            "max": 5.0,
            "start": 3.0,
            "dssat": "SNO3",
            "op": "mult",
        }
    return cfg, fixed


def configure_initial_stage(stage1_theta: dict, lai_theta: dict) -> dict:
    cfg = configure_lai_stage(stage1_theta, lai_theta)
    cfg["calibrator"]["name"] = f"china_hemp_initial_water_no3_after_lai_{DATE_TAG}"
    genotype_theta = {**stage1_theta, **lai_theta}
    _set_cultivar_fixed(cfg, genotype_theta, _merged_param_groups(STAGE1_PHENO_PARAMS, LAI_PARAMS))
    for group in ("genetic_cultivar", "genetic_ecotype"):
        for spec in cfg["parameters"].get(group, {}).values():
            if isinstance(spec, dict):
                spec["active"] = False
    init = cfg["parameters"].setdefault("initial_conditions", {})
    init["initial_soil_water_mult"] = {
        "active": True,
        "min": 1.0,
        "max": 1.8,
        "start": 1.1,
        "scope": "experiment",
        "dssat": "SH2O",
        "op": "mult",
        "clip_01": True,
    }
    init["initial_no3_mult"] = {
        "active": True,
        "min": 1.0,
        "max": 6.0,
        "start": 1.5,
        "scope": "experiment",
        "dssat": "SNO3",
        "op": "mult",
    }
    _lai_objective(cfg, include_biomass=True)
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_INITIAL_MAXITER", 6)),
        "popsize": int(os.environ.get("DSSATCAL_INITIAL_POPSIZE", 12)),
        "restarts": 1,
    })
    return cfg


def configure_capacity_stage(stage1_theta: dict, lai_theta: dict) -> dict:
    cfg = configure_lai_stage(stage1_theta, lai_theta)
    cfg["calibrator"]["name"] = f"china_hemp_capacity_partition_expanded_{DATE_TAG}"
    genotype_theta = {**stage1_theta, **lai_theta}
    _set_cultivar_fixed(cfg, stage1_theta, STAGE1_PHENO_PARAMS)
    _activate_lai_parameters(cfg, lai_theta)
    _activate_capacity_parameters(cfg, genotype_theta)
    _capacity_objective(cfg)
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_CAPACITY_MAXITER", 6)),
        "popsize": int(os.environ.get("DSSATCAL_CAPACITY_POPSIZE", 12)),
        "restarts": 1,
    })
    return cfg


def configure_capacity_locked_stage(stage1_theta: dict, lai_theta: dict) -> dict:
    cfg = configure_lai_stage(stage1_theta, lai_theta)
    cfg["calibrator"]["name"] = f"china_hemp_capacity_partition_locked_lai_{DATE_TAG}"
    genotype_theta = {**stage1_theta, **lai_theta}
    _set_cultivar_fixed(
        cfg,
        genotype_theta,
        _merged_param_groups(STAGE1_PHENO_PARAMS, LAI_PARAMS),
    )
    for group in ("genetic_cultivar", "genetic_ecotype"):
        for spec in cfg["parameters"].get(group, {}).values():
            if isinstance(spec, dict):
                spec["active"] = False
    _activate_capacity_parameters(cfg, genotype_theta)
    _lai_objective(cfg, include_biomass=True)
    cfg["engine"].setdefault("timeseries_outputs", {})["grain"] = "GWAD"
    cfg["objective"]["weights"].update({
        "biomass": 1.4,
        "stem": 1.4,
        "leaf": 0.9,
        "grain": 0.7,
        "LAI": 1.2,
        "node_stage": 0.8,
        "height": 0.7,
        "width": 0.5,
    })
    cfg["objective"]["error_model"]["grain"] = {"type": "relative", "value": 0.35}
    cfg["objective"].setdefault("likelihood", {})["type"] = "huber"
    cfg["objective"]["likelihood"]["delta"] = 2.0
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_CAPACITY_LOCKED_MAXITER", 6)),
        "popsize": int(os.environ.get("DSSATCAL_CAPACITY_LOCKED_POPSIZE", 12)),
        "restarts": 1,
    })
    return cfg


def configure_species_large_stage(stage1_theta: dict, lai_theta: dict) -> dict:
    cfg = configure_lai_stage(stage1_theta, lai_theta)
    cfg["calibrator"]["name"] = f"china_hemp_species_canopy_source_large_{DATE_TAG}"
    _set_cultivar_fixed(cfg, stage1_theta, STAGE1_PHENO_PARAMS)
    _activate_lai_parameters(cfg, lai_theta)
    protect_lai = os.environ.get("DSSATCAL_SPECIES_PROTECT_LAI", "1").lower() in {"1", "true", "yes"}
    if protect_lai:
        _set_cultivar_fixed(cfg, lai_theta, LAI_PARAMS)
    _activate_capacity_parameters(cfg, {**stage1_theta, **lai_theta})
    _activate_species_canopy_source_parameters(cfg)
    cfg["gating"]["species"] = "free"
    cfg["engine"].setdefault("timeseries_outputs", {})["grain"] = "GWAD"
    cfg["objective"]["weights"].update({
        "LAI": 2.4,
        "biomass": 1.4,
        "stem": 1.4,
        "leaf": 1.1,
        "height": 1.1,
        "width": 1.0,
        "node_stage": 0.8,
        "grain": 0.5,
        "anthesis": 0.7,
        "emergence": 0.4,
    })
    cfg["objective"]["error_model"].update({
        "LAI": {"type": "absolute", "value": 0.9},
        "biomass": {"type": "relative", "value": 0.35},
        "stem": {"type": "relative", "value": 0.35},
        "leaf": {"type": "relative", "value": 0.35},
        "height": {"type": "absolute", "value": 0.45},
        "width": {"type": "absolute", "value": 0.25},
        "grain": {"type": "relative", "value": 0.45},
    })
    strict_biomass = os.environ.get("DSSATCAL_SPECIES_STRICT", "0").lower() in {"1", "true", "yes"}
    if strict_biomass:
        cfg["objective"]["weights"].update({
            "LAI": 2.0,
            "biomass": 2.2,
            "stem": 2.0,
            "leaf": 1.7,
            "height": 2.0,
            "width": 1.2,
            "node_stage": 1.2,
            "grain": 0.8,
        })
        cfg["objective"]["error_model"].update({
            "LAI": {"type": "absolute", "value": 0.8},
            "biomass": {"type": "relative", "value": 0.25},
            "stem": {"type": "relative", "value": 0.25},
            "leaf": {"type": "relative", "value": 0.25},
            "height": {"type": "absolute", "value": 0.25},
            "width": {"type": "absolute", "value": 0.20},
            "grain": {"type": "relative", "value": 0.35},
        })
        cfg["objective"].setdefault("likelihood", {})["type"] = "gaussian"
    else:
        cfg["objective"].setdefault("likelihood", {})["type"] = "huber"
        cfg["objective"]["likelihood"]["delta"] = 2.0
    cfg["method"].setdefault("optimizer", {}).update({
        "engine": "cmaes",
        "maxiter": int(os.environ.get("DSSATCAL_SPECIES_MAXITER", 20)),
        "popsize": int(os.environ.get("DSSATCAL_SPECIES_POPSIZE", 16)),
        "restarts": 1,
    })
    return cfg


def configure_species_lhs_stage(stage1_theta: dict, lai_theta: dict) -> dict:
    cfg = configure_species_large_stage(stage1_theta, lai_theta)
    cfg["calibrator"]["name"] = f"china_hemp_species_lhs_protected_lai_{DATE_TAG}"
    cfg["method"].setdefault("optimizer", {})["engine"] = "none"
    cfg["method"].setdefault("bayesian", {})["engine"] = "glue"
    cfg["method"].setdefault("sample", {}).update({
        "engine": "lhs",
        "n": int(os.environ.get("DSSATCAL_LHS_N", 160)),
    })
    return cfg


def _write_run_tables(result, outdir: Path, *, label: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{label}_best_theta.json").write_text(json.dumps(result.best_theta, indent=2, sort_keys=True))
    if result.best.residuals is None or result.best.residuals.empty:
        return
    resid = result.best.residuals.copy()
    resid.to_csv(outdir / f"{label}_residuals.csv", index=False)
    summary = resid.groupby(["exp_id", "user_var"], as_index=False).agg(
        n=("resid", "size"),
        rmse=("resid", lambda x: float((x.pow(2).mean()) ** 0.5)),
        mbe=("resid", "mean"),
        mean_obs=("obs", "mean"),
        mean_sim=("sim", "mean"),
        last_date=("date", "max"),
    )
    summary.to_csv(outdir / f"{label}_residual_summary_by_exp_var.csv", index=False)
    print(summary.round(3).to_string(index=False), flush=True)


def _run_and_report(cfg: dict, *, label: str, extra_json: dict | None = None):
    name = cfg["calibrator"]["name"]
    print(f"\n=== Running {name} ===", flush=True)
    result = orchestrator.calibrate(cfg, progress=True)
    spawns = orchestrator.spawn_results_for(cfg, result.best_theta, result.experiments)
    outdir = RESULTS_DIR / name
    figdir = FIGURES_DIR / name
    paths = viz.make_report(result, outdir, best_spawns=spawns, figdir=figdir)
    _write_run_tables(result, outdir, label=label)
    if extra_json:
        for filename, payload in extra_json.items():
            (outdir / filename).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(viz.summary_fit_table(result).round(3).to_string(index=False), flush=True)
    print(f"Tables:  {outdir.resolve()}", flush=True)
    print(f"Figures: {figdir.resolve()}", flush=True)
    print(f"Experiment panels: {len(paths.get('experiment_diagnostics', []))}", flush=True)
    return result


def _fixed_report(cfg: dict, theta: dict, *, label: str):
    print(f"\n=== Running fixed diagnostic {cfg['calibrator']['name']} ===", flush=True)
    setup = orchestrator._setup(cfg)
    space, _crop, _exe, _specs, _run_root, obs, experiments, _treatments = setup
    spawns = orchestrator.spawn_results_for(cfg, theta, experiments)
    best = obj.score(spawns, obs.table, cfg)
    row = {"sample_id": 0, "score": best.score, "loglik": best.loglik, "n_obs": len(best.residuals)}
    for name in space.names:
        row[name] = float(theta.get(name, space.start[space.names.index(name)]))
    design = pd.DataFrame([row])
    result = CalibrationResult(
        cfg=cfg,
        space=space,
        obs=obs,
        experiments=experiments,
        design=design,
        obj_results={0: best},
        best_theta={name: row[name] for name in space.names},
        best=best,
        glue=None,
        nsga2=None,
        extras={},
    )
    outdir = RESULTS_DIR / cfg["calibrator"]["name"]
    figdir = FIGURES_DIR / cfg["calibrator"]["name"]
    paths = viz.make_report(result, outdir, best_spawns=spawns, figdir=figdir)
    _write_run_tables(result, outdir, label=label)
    print(viz.summary_fit_table(result).round(3).to_string(index=False), flush=True)
    print(f"Tables:  {outdir.resolve()}", flush=True)
    print(f"Figures: {figdir.resolve()}", flush=True)
    print(f"Experiment panels: {len(paths.get('experiment_diagnostics', []))}", flush=True)
    return result


def _section(lines: list[str], title: str) -> tuple[int, int] | None:
    start = next((i for i, ln in enumerate(lines) if ln.upper().startswith(f"*{title.upper()}")), None)
    if start is None:
        return None
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("*")), len(lines))
    return start, end


def _split_header_row(header: str, row: str) -> dict[str, str]:
    names = header.lstrip("@").split()
    values = row.split(maxsplit=max(0, len(names) - 1))
    return dict(zip(names, values))


def _input_qa_tables(cfg: dict) -> None:
    outdir = RESULTS_DIR / f"china_hemp_input_qa_{DATE_TAG}"
    outdir.mkdir(parents=True, exist_ok=True)
    hemp_dir = Path(cfg["source"]["hemp_dir"])
    rows = []
    init_rows = []
    for exp_id in cfg["experiments"]:
        path = hemp_dir / f"{exp_id}.HMX"
        lines = path.read_text(errors="replace").splitlines()
        sec = _section(lines, "PLANTING DETAILS")
        if sec:
            start, end = sec
            header = next((lines[i] for i in range(start, end) if lines[i].lstrip().startswith("@P")), None)
            if header:
                for ln in lines[lines.index(header) + 1:end]:
                    if not re.match(r"\s*\d", ln) or ln.lstrip().startswith("!"):
                        continue
                    rec = _split_header_row(header, ln)
                    rows.append({
                        "exp_id": exp_id,
                        "treatment": rec.get("P"),
                        "PDATE": rec.get("PDATE"),
                        "EDATE": rec.get("EDATE"),
                        "PPOP": rec.get("PPOP"),
                        "PPOE": rec.get("PPOE"),
                        "PLRS": rec.get("PLRS"),
                        "PLDP": rec.get("PLDP"),
                    })
        sec = _section(lines, "INITIAL CONDITIONS")
        if sec:
            start, end = sec
            header_idx = next((i for i in range(start, end) if lines[i].lstrip().startswith("@C  ICBL")), None)
            if header_idx:
                header = lines[header_idx]
                for ln in lines[header_idx + 1:end]:
                    if not re.match(r"\s*\d", ln) or ln.lstrip().startswith("!"):
                        continue
                    rec = _split_header_row(header, ln)
                    init_rows.append({
                        "exp_id": exp_id,
                        "ICBL": rec.get("ICBL"),
                        "SH2O": rec.get("SH2O"),
                        "SNH4": rec.get("SNH4"),
                        "SNO3": rec.get("SNO3"),
                    })
    pd.DataFrame(rows).to_csv(outdir / "translated_planting_inputs.csv", index=False)
    pd.DataFrame(init_rows).to_csv(outdir / "translated_initial_conditions.csv", index=False)
    coverage = pd.read_csv("calibration_china_hemp/observation_coverage.csv")
    coverage.to_csv(outdir / "observation_coverage.csv", index=False)
    summary = [
        "# China Hemp Input QA",
        "",
        f"Experiments: {', '.join(cfg['experiments'])}",
        "YUKU2202 and YUKU2203 are not included.",
        "Planting date, plant population, emergence population, row spacing, planting depth, fertilizer, and irrigation remain sourced from the original HMX files.",
        "CNKU2101 keeps the configured spawn-time weather/soil override.",
        "",
        "Tables written:",
        "- translated_planting_inputs.csv",
        "- translated_initial_conditions.csv",
        "- observation_coverage.csv",
    ]
    (outdir / "input_qa_summary.md").write_text("\n".join(summary) + "\n")
    print(f"Input QA tables: {outdir.resolve()}", flush=True)


def _combined_genotype_theta(stage1_theta: dict, lai_theta: dict) -> dict:
    return {**stage1_theta, **lai_theta}


def _value_for(values: dict, theta: dict, name: str, cultivar: str) -> float:
    value = theta.get(f"{name}__{cultivar}", values.get(name, FALLBACK_GENOTYPE_VALUES.get(name, 0.0)))
    if value is None or (isinstance(value, float) and pd.isna(value)):
        value = FALLBACK_GENOTYPE_VALUES.get(name, 0.0)
    return float(value)


def _cultivar_row(code: str, name: str, eco: str, values: dict, theta: dict, cultivar: str) -> str:
    order = [
        "CSDL", "PPSEN", "EM-FL", "FL-SH", "FL-SD", "SD-PM", "FL-LF", "LFMAX",
        "SLAVR", "SIZLF", "XFRT", "WTPSD", "SFDUR", "SDPDV", "PODUR", "THRSH",
        "SDPRO", "SDLIP",
    ]
    vals = {k: _value_for(values, theta, k, cultivar) for k in order}
    return (
        f"{code:<6} {name:<18} . {eco:<6}"
        f" {vals['CSDL']:5.2f} {vals['PPSEN']:5.3f} {vals['EM-FL']:5.2f}"
        f" {vals['FL-SH']:5.2f} {vals['FL-SD']:5.2f} {vals['SD-PM']:5.2f}"
        f" {vals['FL-LF']:5.2f} {vals['LFMAX']:5.2f} {vals['SLAVR']:5.1f}"
        f" {vals['SIZLF']:5.1f} {vals['XFRT']:5.3f} {vals['WTPSD']:5.3f}"
        f" {vals['SFDUR']:5.2f} {vals['SDPDV']:5.2f} {vals['PODUR']:5.2f}"
        f" {vals['THRSH']:5.2f} {vals['SDPRO']:5.3f} {vals['SDLIP']:5.3f}"
    )


def _ecotype_row(code: str, name: str, values: dict, theta: dict, cultivar: str) -> str:
    order = [
        "THVAR", "PL-EM", "EM-V1", "V1-JU", "JU-R0", "PM06", "PM09", "LNGSH",
        "R7-R8", "FL-VS", "TRIFL", "RWDTH", "RHGHT", "R1PPO", "OPTBI", "SLOBI",
    ]
    vals = {k: _value_for(values, theta, k, cultivar) for k in order}
    return (
        f"{code:<6} {name:<16} 07 01"
        f" {vals['THVAR']:5.3f} {vals['PL-EM']:5.2f} {vals['EM-V1']:5.2f}"
        f" {vals['V1-JU']:5.2f} {vals['JU-R0']:5.2f} {vals['PM06']:5.3f}"
        f" {vals['PM09']:5.3f} {vals['LNGSH']:5.2f} {vals['R7-R8']:5.2f}"
        f" {vals['FL-VS']:5.2f} {vals['TRIFL']:5.3f} {vals['RWDTH']:5.3f}"
        f" {vals['RHGHT']:5.3f} {vals['R1PPO']:5.3f} {vals['OPTBI']:5.2f}"
        f" {vals['SLOBI']:5.3f}"
    )


def _export_genotypes(cfg: dict, theta: dict, outdir: Path, *, value_source: str) -> None:
    paths = resolve_dssat_paths(cfg)
    stem = cfg["crops"][0]["genotype_stem"]
    src_cul = paths["genotype"] / f"{stem}.CUL"
    src_eco = paths["genotype"] / f"{stem}.ECO"
    src_spe = paths["genotype"] / f"{stem}.SPE"
    dst_cul = outdir / f"{stem}.CUL"
    dst_eco = outdir / f"{stem}.ECO"
    dst_spe = outdir / f"{stem}.SPE"
    shutil.copy(src_cul, dst_cul)
    shutil.copy(src_eco, dst_eco)
    if src_spe.exists():
        shutil.copy(src_spe, dst_spe)
    cul_lines = dst_cul.read_text(errors="replace").splitlines()
    eco_lines = dst_eco.read_text(errors="replace").splitlines()
    cul_lines.extend(["", "! China LAI/resource calibration generated by dssatcalibrator"])
    eco_lines.extend(["", "! China LAI/resource calibration generated by dssatcalibrator"])
    records = []
    for cultivar in CULTIVARS:
        new_cul, new_eco, label = NEW_IDS[cultivar]
        source_eco = SOURCE_ECOTYPES[cultivar]
        cul_values = read_cultivar_values(src_cul, cultivar)
        eco_values = read_ecotype_values(src_eco, source_eco)
        cul_lines.append(_cultivar_row(new_cul, label, new_eco, cul_values, theta, cultivar))
        eco_lines.append(_ecotype_row(new_eco, f"{label[:12]} cal", eco_values, theta, cultivar))
        for group, names in _merged_param_groups(STAGE1_PHENO_PARAMS, LAI_PARAMS, CAPACITY_PARAMS).items():
            for name in names:
                value = _theta_value(theta, name, cultivar, float("nan"))
                if pd.isna(value):
                    continue
                records.append({
                    "source_cultivar": cultivar,
                    "new_cultivar": new_cul,
                    "source_ecotype": source_eco,
                    "new_ecotype": new_eco,
                    "group": group,
                    "parameter": name,
                    "value": value,
                    "value_source": value_source,
                })
    dst_cul.write_text("\n".join(cul_lines) + "\n")
    dst_eco.write_text("\n".join(eco_lines) + "\n")
    pd.DataFrame(records).to_csv(outdir / "calibrated_genotype_parameters.csv", index=False)
    species_records = []
    if dst_spe.exists():
        for name, spec in (cfg.get("parameters", {}).get("genetic_species", {}) or {}).items():
            if not isinstance(spec, dict) or name not in theta:
                continue
            update = {
                "value": float(theta[name]),
                "index": int(spec.get("spe_index", spec.get("token_index", 0))),
            }
            edit_species(dst_spe, {spec.get("spe_key", name): update})
            species_records.append({
                "group": "genetic_species",
                "parameter": name,
                "value": float(theta[name]),
                "spe_key": spec.get("spe_key", name),
                "spe_index": int(spec.get("spe_index", spec.get("token_index", 0))),
                "value_source": value_source,
            })
        if species_records:
            pd.DataFrame(species_records).to_csv(outdir / "calibrated_species_parameters.csv", index=False)


def main() -> None:
    part = os.environ.get("DSSATCAL_RUN_PART", "all").lower()
    stage1_theta = _read_json(Path(os.environ.get("DSSATCAL_STAGE1_THETA", DEFAULT_STAGE1_THETA)))
    stage2_theta = _read_json(Path(os.environ.get("DSSATCAL_STAGE2_THETA", DEFAULT_STAGE2_THETA)))
    base_cfg = configure_lai_stage(stage1_theta, stage2_theta)

    if part in {"all", "audit"}:
        _input_qa_tables(base_cfg)
        if part == "audit":
            return

    if part in {"all", "diagnostic"}:
        baseline_cfg, baseline_theta = configure_resource_diagnostic(stage1_theta, stage2_theta, rich=False)
        rich_cfg, rich_theta = configure_resource_diagnostic(stage1_theta, stage2_theta, rich=True)
        _fixed_report(baseline_cfg, baseline_theta, label="resource_baseline")
        _fixed_report(rich_cfg, rich_theta, label="resource_rich")
        if part == "diagnostic":
            return

    lai_result = None
    if part in {"all", "lai"}:
        lai_result = _run_and_report(
            configure_lai_stage(stage1_theta, stage2_theta),
            label="lai",
            extra_json={
                "protected_stage1_phenology.json": stage1_theta,
                "stage2_theta_used_as_lai_start.json": stage2_theta,
            },
        )
        _export_genotypes(
            lai_result.cfg,
            _combined_genotype_theta(stage1_theta, lai_result.best_theta),
            RESULTS_DIR / lai_result.cfg["calibrator"]["name"],
            value_source="stage1+lai",
        )
        if part == "lai":
            return

    if part in {"all", "initial"}:
        if lai_result is not None:
            lai_theta = lai_result.best_theta
        else:
            lai_theta = _read_json(Path(os.environ["DSSATCAL_LAI_THETA"]))
        initial_result = _run_and_report(
            configure_initial_stage(stage1_theta, lai_theta),
            label="initial",
            extra_json={
                "protected_genotype_theta.json": _combined_genotype_theta(stage1_theta, lai_theta),
            },
        )
        _export_genotypes(
            initial_result.cfg,
            _combined_genotype_theta(stage1_theta, lai_theta),
            RESULTS_DIR / initial_result.cfg["calibrator"]["name"],
            value_source="stage1+lai",
        )

    if part in {"capacity"}:
        lai_theta = _read_json(Path(os.environ["DSSATCAL_LAI_THETA"]))
        capacity_result = _run_and_report(
            configure_capacity_stage(stage1_theta, lai_theta),
            label="capacity",
            extra_json={
                "protected_stage1_phenology.json": stage1_theta,
                "lai_theta_used_as_start.json": lai_theta,
            },
        )
        _export_genotypes(
            capacity_result.cfg,
            _combined_genotype_theta(stage1_theta, capacity_result.best_theta),
            RESULTS_DIR / capacity_result.cfg["calibrator"]["name"],
            value_source="stage1+capacity",
        )

    if part in {"capacity_locked"}:
        lai_theta = _read_json(Path(os.environ["DSSATCAL_LAI_THETA"]))
        capacity_result = _run_and_report(
            configure_capacity_locked_stage(stage1_theta, lai_theta),
            label="capacity_locked",
            extra_json={
                "protected_stage1_phenology.json": stage1_theta,
                "protected_lai_theta.json": lai_theta,
            },
        )
        _export_genotypes(
            capacity_result.cfg,
            _combined_genotype_theta(stage1_theta, {**lai_theta, **capacity_result.best_theta}),
            RESULTS_DIR / capacity_result.cfg["calibrator"]["name"],
            value_source="stage1+lai+capacity_locked",
        )

    if part in {"species_large"}:
        lai_theta = _read_json(Path(os.environ["DSSATCAL_LAI_THETA"]))
        species_result = _run_and_report(
            configure_species_large_stage(stage1_theta, lai_theta),
            label="species_large",
            extra_json={
                "protected_stage1_phenology.json": stage1_theta,
                "lai_theta_used_as_start.json": lai_theta,
            },
        )
        _export_genotypes(
            species_result.cfg,
            _combined_genotype_theta(stage1_theta, {**lai_theta, **species_result.best_theta}),
            RESULTS_DIR / species_result.cfg["calibrator"]["name"],
            value_source="stage1+lai+species_large",
        )

    if part in {"species_lhs"}:
        lai_theta = _read_json(Path(os.environ["DSSATCAL_LAI_THETA"]))
        species_result = _run_and_report(
            configure_species_lhs_stage(stage1_theta, lai_theta),
            label="species_lhs",
            extra_json={
                "protected_stage1_phenology.json": stage1_theta,
                "protected_lai_theta.json": lai_theta,
            },
        )
        _export_genotypes(
            species_result.cfg,
            _combined_genotype_theta(stage1_theta, {**lai_theta, **species_result.best_theta}),
            RESULTS_DIR / species_result.cfg["calibrator"]["name"],
            value_source="stage1+lai+species_lhs",
        )


if __name__ == "__main__":
    main()
