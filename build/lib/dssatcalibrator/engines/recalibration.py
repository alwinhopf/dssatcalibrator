from __future__ import annotations
import logging
from copy import deepcopy
from datetime import date
import pandas as pd

logger = logging.getLogger(__name__)

class InSeasonRecalibrator:
    """Mid-season parameter re-estimation engine (the coupled in-season path).

    Instead of updating state directly, this engine adjusts parameters (growth rate coefficients,
    stress factors, etc.) so that the model trajectory aligns with the observations received
    up to a specific date. Each call runs the full calibration pipeline on the observations
    available so far, so it IS coupled to DSSAT (unlike the EnKF/forcing prototypes).
    """

    def __init__(self, cfg: dict):
        self.cfg = cfg

    def recalibrate(self, obs_df: pd.DataFrame, current_date: date,
                    warm_start_theta: dict | None = None) -> dict:
        """Filter observations up to ``current_date`` and run a calibration.

        Parameters
        ----------
        warm_start_theta
            If given (and ``recalibration.warm_start`` is on), each active
            parameter's ``start`` value is seeded with the previous checkpoint's
            best fit, so successive checkpoints refine rather than restart.
        """
        ts_limit = pd.Timestamp(current_date)

        # Undated scalars are end-of-season FileA outcomes and are unavailable
        # at an in-season checkpoint.
        filtered_table = obs_df[obs_df["date"].notna() & (obs_df["date"] <= ts_limit)].copy()

        logger.info(f"Recalibrating with {len(filtered_table)} observations up to {current_date}")

        # deepcopy so seeding `start` values never mutates the caller's config.
        cfg_updated = deepcopy(self.cfg)
        # Pass the pre-filtered DataFrame directly in the config
        cfg_updated["source"] = {**cfg_updated.get("source", {}), "table": filtered_table}

        recal_cfg = self.cfg.get("assimilation", {}).get("recalibration", {})
        recal_n = recal_cfg.get("recal_sample_size", 100)

        # Warm start: seed active-parameter start values with the previous best theta.
        if warm_start_theta and recal_cfg.get("warm_start", True):
            for _group, params in (cfg_updated.get("parameters") or {}).items():
                if not isinstance(params, dict):
                    continue
                for name, spec in params.items():
                    if isinstance(spec, dict) and name in warm_start_theta:
                        spec["start"] = float(warm_start_theta[name])

        method = cfg_updated.setdefault("method", {})
        method.setdefault("sample", {})["n"] = recal_n
        # Optionally pin the checkpoint estimator (default keeps the config's engine).
        engine = recal_cfg.get("engine")
        if engine and str(engine).lower() != "none":
            method.setdefault("bayesian", {})["engine"] = engine

        # Run calibration
        from ..orchestrator import calibrate
        result = calibrate(cfg_updated, progress=False)
        return result.best_theta
