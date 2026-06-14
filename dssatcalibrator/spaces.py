"""Parameter search space built from the config's active parameters."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import active_parameters


@dataclass
class ParameterSpace:
    names: list[str]
    low: np.ndarray
    high: np.ndarray
    start: np.ndarray
    specs: list[dict]

    @property
    def ndim(self) -> int:
        return len(self.names)

    def to_theta(self, vector) -> dict[str, float]:
        """Map a 1-D parameter vector (in native units) to a named theta dict."""
        return {n: float(v) for n, v in zip(self.names, vector)}

    def clip(self, vector) -> np.ndarray:
        return np.clip(np.asarray(vector, float), self.low, self.high)

    @classmethod
    def from_config(cls, cfg: dict) -> "ParameterSpace":
        specs = active_parameters(cfg)
        if not specs:
            raise ValueError("No active parameters in config (set active: true on some).")
        names = [s["name"] for s in specs]
        low = np.array([float(s["min"]) for s in specs])
        high = np.array([float(s["max"]) for s in specs])
        start = np.array([float(s.get("start", 0.5 * (lo + hi)))
                          for s, lo, hi in zip(specs, low, high)])
        return cls(names=names, low=low, high=high, start=np.clip(start, low, high), specs=specs)
