from datetime import date
import pandas as pd
import numpy as np
from .base import ObservationSource

def _apply_obs_operator(config: dict, value: float) -> float:
    """Optional satellite-LAI observation operator: effective->true LAI scale/offset.

    Satellite-retrieved LAI is biased vs the model's LAI (clumping, saturation). An
    optional linear operator ``value' = scale*value + offset`` (config
    ``obs_operator: {scale, offset}``) lets the user correct that bias instead of
    letting it leak into the calibrated coefficients. Identity by default.
    """
    op = config.get("obs_operator", {}) or {}
    return float(op.get("scale", 1.0)) * value + float(op.get("offset", 0.0))


class SentinelLAISource(ObservationSource):
    """Sentinel-2 derived LAI via vegetation indices."""
    name = "sentinel2_lai"
    source_type = "satellite"

    def fetch(self, experiment: str, date_range: tuple[date, date], **kwargs) -> pd.DataFrame:
        path = self.config.get("data_path")
        if not path:
            return pd.DataFrame()

        df = pd.read_csv(path)
        df["date"] = pd.to_datetime(df["date"])
        df = df[(df["exp_id"] == experiment) &
                (df["date"].dt.date >= date_range[0]) &
                (df["date"].dt.date <= date_range[1])]

        max_cloud = float(self.config.get("max_cloud_fraction", 1.0))
        out = []
        for _, r in df.iterrows():
            # cloud masking: drop scenes cloudier than the threshold outright (a
            # clouded LAI retrieval is unusable, not merely noisier).
            cloud = float(r.get("cloud_fraction", 0.0))
            if cloud > max_cloud:
                continue
            metadata = {"cloud_fraction": cloud}
            val = _apply_obs_operator(self.config, float(r["value"]))
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": "LAID",
                "kind": "timeseries",
                "date": r["date"],
                "value": val,
                "sigma": self.error_model("LAID", val, metadata),
                "weight": 1.0,
                "source": self.name,
                "quality_flag": int(r.get("quality_flag", 0)),
                "spatial_res_m": 10.0
            })
        return pd.DataFrame(out)

    def quality_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """Drop rows whose quality_flag is non-zero (clouded/bad retrieval)."""
        if df.empty or "quality_flag" not in df.columns:
            return df
        if not self.config.get("drop_bad_quality", True):
            return df
        return df[df["quality_flag"].fillna(0) == 0].reset_index(drop=True)

    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        base_rmse = self.config.get("error_model", {}).get("base_rmse", 0.7)
        sat_lai = self.config.get("error_model", {}).get("saturation_lai", 4.0)
        
        if value > sat_lai:
            base_rmse *= 1.0 + 0.3 * (value - sat_lai)
        cloud_factor = 1.0 + metadata.get("cloud_fraction", 0.0) * 0.5
        return base_rmse * cloud_factor

    def variable_mapping(self) -> dict[str, str]:
        return {"sentinel_lai": "LAID"}


class MODISLAISource(ObservationSource):
    """MODIS LAI (coarser but higher temporal frequency)."""
    name = "modis_lai"
    source_type = "satellite"
    
    def fetch(self, experiment: str, date_range: tuple[date, date], **kwargs) -> pd.DataFrame:
        path = self.config.get("data_path")
        if not path:
            return pd.DataFrame()
        
        df = pd.read_csv(path)
        df["date"] = pd.to_datetime(df["date"])
        df = df[(df["exp_id"] == experiment) & 
                (df["date"].dt.date >= date_range[0]) & 
                (df["date"].dt.date <= date_range[1])]
        
        max_qc = int(self.config.get("max_qc_flag", 99))
        out = []
        for _, r in df.iterrows():
            qc = int(r.get("qc_flag", 0))
            if qc > max_qc:
                continue
            metadata = {"qc_flag": qc}
            val = _apply_obs_operator(self.config, float(r["value"]))
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": "LAID",
                "kind": "timeseries",
                "date": r["date"],
                "value": val,
                "sigma": self.error_model("LAID", val, metadata),
                "weight": 0.5,
                "source": self.name,
                "quality_flag": qc,
                "spatial_res_m": 250.0
            })
        return pd.DataFrame(out)

    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        base_rmse = self.config.get("error_model", {}).get("base_rmse", 0.66)
        qc = metadata.get("qc_flag", 0)
        return base_rmse * (1.0 + 0.5 * qc)

    def variable_mapping(self) -> dict[str, str]:
        return {"modis_lai": "LAID"}
