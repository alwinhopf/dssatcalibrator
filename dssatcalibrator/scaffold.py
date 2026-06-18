"""Scaffold a new crop/cultivar from an analog DSSAT module (optional helper).

DSSAT cannot calibrate a species that does not exist — it tunes coefficients inside
an existing crop module. For a **new cultivar** of an existing species, or a **new
species adapted from the most similar analog module** (the carinata-from-canola
pattern), this helper does the manual scaffolding:

* clone the analog ``.CUL``/``.ECO``/``.SPE`` under a new genotype stem + crop code;
* (optionally) duplicate one cultivar row under a new anchor code to calibrate;
* emit a **starter ``parameters:`` config block** — coefficient names with
  ``start`` from the chosen cultivar, ``min``/``max`` from the file's MINIMA/MAXIMA
  calibration rows (or ±``spread`` if absent), informative ``normal`` priors, and a
  rough ``obligatory``/``candidate`` role split (phenology vs growth).

This is a *starting point* a user reviews against literature — it does not invent
physiology. ``.SPE`` editing for new species stays gated (``gating.species``).
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from .config import resolve_template_dir
from .writers import (cultivar_field_map, read_cul_calibration_bounds,
                      read_cultivar_values)

logger = logging.getLogger(__name__)

# Coefficient-name heuristics for the AgMIP obligatory(=phenology)/candidate split.
_PHENOLOGY_HINT = {"CSDL", "PPSEN", "EM-FL", "FL-SD", "SD-PM", "FL-SH", "FL-LF",
                   "PL-EM", "PLEM", "P1", "P2", "P3", "P4", "P5", "P1V", "P1D",
                   "PHINT", "EM-V1", "PHTHRS"}


def _role(name: str) -> str:
    return "obligatory" if name.upper() in _PHENOLOGY_HINT else "candidate"


def scaffold_crop(*, dssat_dir: str | Path, analog_stem: str, new_stem: str,
                  new_code: str, source_anchor: str, new_anchor: str | None = None,
                  out_dir: str | Path | None = None, spread: float = 0.3,
                  copy_spe: bool = True) -> dict:
    """Clone analog genotype files and emit a starter parameter block.

    Returns ``{"files": {...}, "parameters_yaml": str, "coefficients": {...}}``.
    Writes the cloned ``.CUL``/``.ECO``/``.SPE`` into ``out_dir``. Does NOT modify
    the analog files or the DSSAT install.
    """
    dssat_dir = Path(dssat_dir)
    geno = dssat_dir / "Genotype"
    if out_dir is None:
        out_dir = resolve_template_dir(required=True)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    exts = ["CUL", "ECO"] + (["SPE"] if copy_spe else [])
    files = {}
    for e in exts:
        src = geno / f"{analog_stem}.{e}"
        if src.exists():
            dst = out_dir / f"{new_stem}.{e}"
            shutil.copy(src, dst)
            files[e] = str(dst)
        else:
            logger.warning("Analog %s.%s not found under %s", analog_stem, e, geno)

    cul = out_dir / f"{new_stem}.CUL"
    fmap = cultivar_field_map(cul)
    starts = read_cultivar_values(cul, source_anchor)
    bounds = read_cul_calibration_bounds(cul)

    coeffs = {}
    for name in fmap:
        start = starts.get(name)
        if start is None:
            continue
        if name in bounds:
            lo, hi = bounds[name]["min"], bounds[name]["max"]
        else:
            lo, hi = start * (1 - spread), start * (1 + spread)
            if lo > hi:
                lo, hi = hi, lo
        coeffs[name] = {"min": round(lo, 4), "max": round(hi, 4),
                        "start": round(start, 4), "role": _role(name)}

    parameters_yaml = _emit_parameters_yaml(coeffs)
    logger.info("Scaffolded %s (code %s) from analog %s; %d coefficients.",
                new_stem, new_code, analog_stem, len(coeffs))
    return {"files": files, "out_dir": str(out_dir), "parameters_yaml": parameters_yaml,
            "coefficients": coeffs, "new_anchor": new_anchor or source_anchor}


def _emit_parameters_yaml(coeffs: dict) -> str:
    """Render a ready-to-paste ``genetic_cultivar:`` block with normal priors.

    Phenology (obligatory) coefficients are active by default; growth (candidate)
    are declared but inactive — start small and let screening/selection add them.
    """
    lines = ["parameters:", "  genetic_cultivar:"]
    for name, c in coeffs.items():
        active = "true" if c["role"] == "obligatory" else "false"
        sd = round((c["max"] - c["min"]) / 6.0, 4) or 1.0
        lines.append(
            f'    "{name}": {{ active: {active}, role: {c["role"]}, '
            f'min: {c["min"]}, max: {c["max"]}, start: {c["start"]}, '
            f'prior: {{dist: normal, sd: {sd}}} }}'
        )
    return "\n".join(lines) + "\n"
