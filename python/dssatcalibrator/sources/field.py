from datetime import date
from pathlib import Path
import pandas as pd
import numpy as np
from .base import ObservationSource
from ..observations import read_filea, read_filet

class FieldMeasurementSource(ObservationSource):
    """Traditional field measurements (FileA / FileT)."""
    name = "field_measurements"
    source_type = "field"
    
    def fetch(self, experiment: str, date_range: tuple[date, date], **kwargs) -> pd.DataFrame:
        hemp_dir = self.config.get("hemp_dir") or self.config.get("data_path")
        if not hemp_dir:
            return pd.DataFrame()
        
        hemp_dir = Path(hemp_dir)
        crop_ext = self.config.get("crop_ext", "HM")
        
        fa = hemp_dir / f"{experiment}.{crop_ext}A"
        ft = hemp_dir / f"{experiment}.{crop_ext}T"
        
        frames = []
        if fa.exists():
            frames.append(read_filea(fa, experiment))
        if ft.exists():
            frames.append(read_filet(ft, experiment))
            
        if not frames:
            return pd.DataFrame()
            
        df = pd.concat(frames, ignore_index=True)
        
        def in_range(row):
            if pd.isna(row["date"]):
                return True
            return date_range[0] <= row["date"].date() <= date_range[1]
            
        df = df[df.apply(in_range, axis=1)].copy()
        
        df["source"] = self.name
        df["quality_flag"] = 0
        df["spatial_res_m"] = np.nan
        
        df["sigma"] = df.apply(lambda r: self.error_model(r["variable"], r["value"], {}), axis=1)
        
        return df
        
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        error_models = {
            "LAID": ("relative", 0.15),
            "CWAD": ("relative", 0.12),
            "HWAM": ("relative", 0.08),
            "GSTD": ("absolute", 1.0),
            "ADAT": ("absolute", 3.0),
            "MDAT": ("absolute", 3.0),
        }
        
        cfg_models = self.config.get("error_model", {})
        if variable in cfg_models:
            kind = cfg_models[variable].get("type", "relative")
            val = cfg_models[variable].get("value", 0.15)
        else:
            kind, val = error_models.get(variable, ("relative", 0.15))
            
        if kind == "relative":
            return max(abs(val * value), 1e-6)
        return val

    def variable_mapping(self) -> dict[str, str]:
        return {}
