from datetime import date
import pandas as pd
import numpy as np
from .base import ObservationSource

class SoilMoistureSensorSource(ObservationSource):
    """IoT soil moisture sensors."""
    name = "soil_moisture_iot"
    source_type = "iot"
    
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
            metadata = {
                "sensor_type": self.config.get("sensor_type", "capacitance"),
                "calibration_status": self.config.get("calibration_status", "factory")
            }
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": "SW",
                "kind": "timeseries",
                "date": r["date"],
                "value": float(r["value"]),
                "sigma": self.error_model("SW", float(r["value"]), metadata),
                "weight": 1.0,
                "source": self.name,
                "quality_flag": 0,
                "spatial_res_m": np.nan
            })
        return pd.DataFrame(out)
        
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        sensor = metadata.get("sensor_type", "capacitance")
        cal = metadata.get("calibration_status", "factory")
        base = {"capacitance": 0.04, "tdr": 0.02, "tensiometer": 0.03}.get(sensor, 0.04)
        if cal == "field_calibrated":
            base *= 0.6
        return base

    def variable_mapping(self) -> dict[str, str]:
        return {"soil_moisture": "SW"}


class CanopyTemperatureSource(ObservationSource):
    """Thermal canopy temperature sensors."""
    name = "canopy_temperature"
    source_type = "iot"
    
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
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": "TMEAN",
                "kind": "timeseries",
                "date": r["date"],
                "value": float(r["value"]),
                "sigma": self.error_model("TMEAN", float(r["value"]), {}),
                "weight": 1.0,
                "source": self.name,
                "quality_flag": 0,
                "spatial_res_m": np.nan
            })
        return pd.DataFrame(out)
        
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        return self.config.get("error_model", {}).get("value", 1.0)

    def variable_mapping(self) -> dict[str, str]:
        return {"canopy_temp": "TMEAN"}
