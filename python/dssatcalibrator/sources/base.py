from abc import ABC, abstractmethod
from datetime import date
import pandas as pd

class ObservationSource(ABC):
    """Abstract base class representing a pluggable source of crop observations."""
    
    name: str                    # e.g., "sentinel2_lai", "farm_phenology"
    source_type: str             # "satellite" | "farm_software" | "field" | "iot" | "uav"
    
    def __init__(self, config: dict):
        self.config = config

    @abstractmethod
    def fetch(self, experiment: str, date_range: tuple[date, date], **kwargs) -> pd.DataFrame:
        """Return observations in the standard schema.
        
        Columns: exp_id | treatment | variable | kind | date | value 
                 | sigma | weight | source | quality_flag | spatial_res_m
        """
        pass
    
    @abstractmethod
    def error_model(self, variable: str, value: float, metadata: dict) -> float:
        """Calculate source-specific observation uncertainty (sigma)."""
        pass
    
    @abstractmethod  
    def variable_mapping(self) -> dict[str, str]:
        """Map source variable names to DSSAT model variable names."""
        pass
    
    def quality_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """Optional source-specific quality control filtering."""
        return df
    
    def temporal_aggregate(self, df: pd.DataFrame, window: str = "5D") -> pd.DataFrame:
        """Optional temporal compositing."""
        return df
