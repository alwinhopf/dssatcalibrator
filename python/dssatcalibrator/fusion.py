import logging
from datetime import date
import pandas as pd
import numpy as np
from .sources import ObservationSource

logger = logging.getLogger(__name__)

SCHEMA_EXTENDED = [
    "exp_id", "treatment", "variable", "kind", "date", "value", 
    "sigma", "weight", "source", "quality_flag", "spatial_res_m"
]

class ObservationFuser:
    """Merges observations from multiple sources with conflict resolution."""
    
    def __init__(self, sources: list[ObservationSource], cfg: dict):
        self.sources = {s.name: s for s in sources}
        self.cfg = cfg
    
    def collect(self, experiment: str, date_range: tuple[date, date], **kwargs) -> pd.DataFrame:
        """Gather from all active sources, apply QC, merge."""
        frames = []
        for name, src in self.sources.items():
            try:
                df = src.fetch(experiment, date_range, **kwargs)
                if df.empty:
                    continue
                df = src.quality_filter(df)
                df["source"] = name
                df["source_type"] = src.source_type
                
                # Apply source-specific sigma if not already set or is NaN
                if "sigma" not in df.columns or df["sigma"].isna().any():
                    df["sigma"] = df.apply(
                        lambda r: src.error_model(r["variable"], r["value"], r.get("metadata", {})),
                        axis=1
                    )
                
                # Ensure all SCHEMA_EXTENDED columns exist
                for col in SCHEMA_EXTENDED:
                    if col not in df.columns:
                        df[col] = np.nan
                
                frames.append(df[SCHEMA_EXTENDED])
            except Exception as e:
                logger.warning(f"Source {name} failed to fetch for {experiment}: {e}")
        
        if not frames:
            return pd.DataFrame(columns=SCHEMA_EXTENDED)
        
        merged = pd.concat(frames, ignore_index=True)
        return self.resolve_conflicts(merged)
    
    def resolve_conflicts(self, df: pd.DataFrame) -> pd.DataFrame:
        """Handle overlapping observations from different sources.
        
        Strategy:
        - "keep_all": Keep all observations
        - "inverse_variance": Weight by 1/sigma² (Bayesian optimal)
        - "priority": Use highest-priority source, discard others
        """
        if df.empty:
            return df
            
        strategy = self.cfg.get("fusion", {}).get("conflict_resolution", "keep_all")
        if strategy == "keep_all":
            return df
        elif strategy == "inverse_variance":
            return self._inverse_variance_merge(df)
        elif strategy == "priority":
            return self._priority_merge(df)
        return df
    
    def _inverse_variance_merge(self, df: pd.DataFrame) -> pd.DataFrame:
        """Merge coincident observations using inverse-variance weighting."""
        df_copy = df.copy()
        df_copy["date_key"] = df_copy["date"].dt.date
        
        groups = df_copy.groupby(["exp_id", "treatment", "variable", "date_key"])
        merged = []
        for key, g in groups:
            if len(g) == 1:
                merged.append(g.iloc[0].drop("date_key"))
            else:
                g_clean = g.dropna(subset=["sigma"])
                if g_clean.empty:
                    merged.append(g.iloc[0].drop("date_key"))
                    continue
                
                sigmas = np.where(g_clean["sigma"] == 0, 1e-6, g_clean["sigma"])
                w = 1.0 / (sigmas ** 2)
                val = (g_clean["value"] * w).sum() / w.sum()
                sig = 1.0 / np.sqrt(w.sum())
                
                row = g_clean.iloc[0].copy()
                row["value"] = val
                row["sigma"] = sig
                row["source"] = "+".join(sorted(g_clean["source"].unique()))
                row["quality_flag"] = int(g_clean["quality_flag"].max())
                row["spatial_res_m"] = g_clean["spatial_res_m"].mean()
                merged.append(row.drop("date_key"))
                
        return pd.DataFrame(merged)
        
    def _priority_merge(self, df: pd.DataFrame) -> pd.DataFrame:
        """Keep only the highest priority source for coincident measurements."""
        priority_list = self.cfg.get("fusion", {}).get("source_priority", [])
        if not priority_list:
            return df
            
        p_map = {name: i for i, name in enumerate(priority_list)}
        
        df_copy = df.copy()
        df_copy["date_key"] = df_copy["date"].dt.date
        df_copy["priority_rank"] = df_copy["source"].map(lambda x: p_map.get(x, 9999))
        
        df_copy = df_copy.sort_values("priority_rank")
        
        dedup = df_copy.drop_duplicates(subset=["exp_id", "treatment", "variable", "date_key"], keep="first")
        return dedup.drop(columns=["date_key", "priority_rank"])
