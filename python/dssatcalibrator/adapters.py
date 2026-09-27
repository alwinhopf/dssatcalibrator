"""Optional model adapters for the shared calibration engine.

Adapters are installed through the ``dssatcalibrator.adapters`` entry-point
group. DSSAT remains built in; selecting an unavailable platform fails closed.
"""
from importlib.metadata import entry_points


def model_platform(cfg):
    raw = str((cfg.get("model") or {}).get("platform", cfg.get("platform", "dssat")))
    value = raw.lower().replace("-", "_")
    return {"dssat_csm": "dssat", "apsim": "apsim_ng", "apsimx": "apsim_ng",
            "apsim_next_generation": "apsim_ng", "apsim_nextgen": "apsim_ng"}.get(value, value)


def model_adapter(cfg):
    platform = model_platform(cfg)
    if platform == "dssat":
        return None
    matches = list(entry_points(group="dssatcalibrator.adapters", name=platform))
    if len(matches) != 1:
        raise ValueError(f"Model platform {platform!r} requires exactly one installed adapter; "
                         "install the corresponding cropmodelcalibrator extension.")
    return matches[0].load()()
