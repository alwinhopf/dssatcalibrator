from datetime import date
import pandas as pd
import numpy as np
from .base import ObservationSource

class UAVMultispectralSource(ObservationSource):
    """UAV-based multispectral imagery."""
    name = "uav_multispectral"
    source_type = "uav"
    
    def fetch(self, experiment: str, date_range: tuple[date, date], **kwargs) -> pd.DataFrame:
        path = self.config.get("data_path")
        if not path:
            return pd.DataFrame()
            
        df = pd.read_csv(path)
        df["date"] = pd.to_datetime(df["date"])
        df = df[(df["exp_id"] == experiment) & 
                (df["date"].dt.date >= date_range[0]) & 
                (df["date"].dt.date <= date_range[1])]
                
        out = []
        for _, r in df.iterrows():
            var = r["variable"]
            val = r["value"]
            metadata = {"flight_quality": r.get("flight_quality", "good")}
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": var,
                "kind": "timeseries",
                "date": r["date"],
                "value": float(val),
                "sigma": self.error_model(var, float(val), metadata),
                "weight": 1.0,
                "source": self.name,
                "quality_flag": 0,
                "spatial_res_m": 0.05
            })
        return pd.DataFrame(out)
        
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        flight_quality = metadata.get("flight_quality", "good")
        
        cfg_models = self.config.get("error_model", {})
        if variable in cfg_models:
            val = cfg_models[variable].get("value", 0.15)
        else:
            base_errors = {"LAID": 0.4, "canopy_cover": 0.05, "canopy_height": 0.03}
            val = base_errors.get(variable, 0.15 * abs(value))
            
        if flight_quality == "poor":
            val *= 1.5
        return val

    def variable_mapping(self) -> dict[str, str]:
        return {"uav_lai": "LAID", "uav_canopy_cover": "canopy_cover", "uav_canopy_height": "canopy_height"}
