from .satellite import SentinelLAISource, MODISLAISource
from .farm_software import FarmPhenologySource, FarmManagementSource
from .field import FieldMeasurementSource
from .iot import SoilMoistureSensorSource, CanopyTemperatureSource
from .uav import UAVMultispectralSource
from .base import ObservationSource

ADAPTER_REGISTRY = {
    "sentinel2_lai":      SentinelLAISource,
    "modis_lai":          MODISLAISource,
    "farm_phenology":     FarmPhenologySource,
    "farm_management":    FarmManagementSource,
    "field_measurements": FieldMeasurementSource,
    "soil_moisture_iot":  SoilMoistureSensorSource,
    "canopy_temperature": CanopyTemperatureSource,
    "uav_multispectral":  UAVMultispectralSource,
}

def build_sources(cfg: dict) -> list[ObservationSource]:
    """Instantiate all active observation sources from config."""
    sources = []
    for name, block in cfg.get("observation_sources", {}).items():
        if not block.get("active", False):
            continue
        adapter_name = block.get("adapter", name)
        if adapter_name not in ADAPTER_REGISTRY:
            raise ValueError(f"Unknown adapter: {adapter_name}. "
                             f"Available: {list(ADAPTER_REGISTRY)}")
        src = ADAPTER_REGISTRY[adapter_name](block)
        sources.append(src)
    return sources
