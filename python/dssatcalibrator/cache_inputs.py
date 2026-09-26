"""Content manifests shared by simulation and objective caches."""
from pathlib import Path
import hashlib


def tree_digest(root, suffixes):
    if not root or not Path(root).is_dir():
        return {}
    root = Path(root)
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*"))
            if p.is_file() and p.suffix.lower() in suffixes}


def simulation_inputs(cfg, paths):
    roots = {k: paths.get(k) for k in ("soil", "weather", "genotype")}
    if paths.get("root"):
        roots["standard_data"] = Path(paths["root"]) / "StandardData"
    for key in ("soil", "weather"):
        section = cfg.get(key) or {}
        if section.get("provider", "file") not in {"file", "none", ""}:
            roots[key + "_acquired"] = section.get("cache_dir", key + "_cache")
    suffixes = {".sol", ".wth", ".cul", ".eco", ".spe", ".sda", ".wda", ".cde", ".co2"}
    result = {key: tree_digest(root, suffixes) for key, root in roots.items()}
    result["implementation"] = tree_digest(Path(__file__).parent, {".py"})
    result["external_soils"] = {}
    for key in ("source_sol_file", "external_soil_file"):
        value = (cfg.get("soil") or {}).get(key)
        if value and Path(value).is_file():
            result["external_soils"][key] = hashlib.sha256(Path(value).read_bytes()).hexdigest()
    if (cfg.get("execution") or {}).get("backend") == "dssatengine":
        import dssatengine
        result["engine"] = tree_digest(Path(dssatengine.__file__).parent, {".py"})
    return result
