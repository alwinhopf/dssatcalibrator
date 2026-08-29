import logging

logger = logging.getLogger(__name__)

class ForcingAssimilator:
    """Direct state replacement — UNCOUPLED PROTOTYPE.

    When an observation is received this forces the tracked state value to match
    (or move toward) the observed value. Like the EnKF, it is **not coupled to
    DSSAT**: the ``model_state`` dict it manipulates is never fed back into a
    running simulation, so as wired today it merely records the observed values.
    The orchestrator refuses to run this mode unless
    ``assimilation.allow_uncoupled: true`` is set; use ``mode: recalibration`` for
    the coupled in-season path.
    """
    
    def __init__(self, cfg: dict):
        self.cfg = cfg
        forcing_cfg = cfg.get("assimilation", {}).get("forcing", {})
        self.min_confidence = forcing_cfg.get("min_confidence", 0.8)
        self.smoothing = forcing_cfg.get("smoothing", True)
        
    def apply(self, model_state: dict, observation: dict) -> dict:
        """Replace or adjust model state based on observation."""
        updated = model_state.copy()
        var = observation.get("variable")
        val = observation.get("value")
        confidence = observation.get("confidence", 1.0)
        
        if confidence >= self.min_confidence:
            if self.smoothing:
                if var in updated:
                    updated[var] = confidence * val + (1.0 - confidence) * updated[var]
                else:
                    updated[var] = val
            else:
                updated[var] = val
                
            logger.info(f"Forced state {var} to {updated[var]} (obs value: {val})")
        return updated
