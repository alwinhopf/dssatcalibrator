from datetime import date
import pandas as pd
import numpy as np
from .base import ObservationSource

class FarmPhenologySource(ObservationSource):
    """Growth stages from farm management software."""
    name = "farm_phenology"
    source_type = "farm_software"
    
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
            metadata = {"date_precision": self.config.get("date_precision", "exact")}
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": var,
                "kind": "phenology" if var in ("ADAT", "MDAT") else "timeseries",
                "date": r["date"],
                "value": float(val) if var == "GSTD" else val,
                "sigma": self.error_model(var, float(val) if var == "GSTD" else 1.0, metadata),
                "weight": 1.0,
                "source": self.name,
                "quality_flag": 0,
                "spatial_res_m": np.nan
            })
        return pd.DataFrame(out)
        
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        precision = metadata.get("date_precision", "exact")
        if variable in ("GSTD", "growth_stage"):
            return {"exact": 1.0, "weekly": 2.0, "biweekly": 3.0}.get(precision, 1.0)
        return {"exact": 2.0, "weekly": 5.0, "biweekly": 7.0}.get(precision, 2.0)

    def variable_mapping(self) -> dict[str, str]:
        return {"growth_stage": "GSTD", "anthesis_date": "ADAT", "maturity_date": "MDAT"}


class FarmManagementSource(ObservationSource):
    """Management event records (irrigation, fertilization, spray dates)."""
    name = "farm_management"
    source_type = "farm_software"
    
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
            out.append({
                "exp_id": r["exp_id"],
                "treatment": int(r["treatment"]),
                "variable": var,
                "kind": "management_constraint",
                "date": r["date"],
                "value": float(val),
                "sigma": self.error_model(var, float(val), {}),
                "weight": 1.0,
                "source": self.name,
                "quality_flag": 0,
                "spatial_res_m": np.nan
            })
        return pd.DataFrame(out)
        
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        if "date" in variable.lower():
            return 1.0
        return abs(value * 0.05)

    def variable_mapping(self) -> dict[str, str]:
        return {}
