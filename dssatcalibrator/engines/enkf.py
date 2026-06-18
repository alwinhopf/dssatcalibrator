import logging
import numpy as np

logger = logging.getLogger(__name__)

class EnsembleKalmanFilter:
    """Ensemble Kalman Filter (EnKF) — UNCOUPLED PROTOTYPE.

    This implements the stochastic-EnKF *update* math for an ensemble of model
    states when a new observation arrives. The linear algebra below is correct in
    isolation, but it is **not coupled to DSSAT**: there is no forecast step that
    runs DSSAT between updates, and an updated state is never re-injected into a
    running simulation (DSSAT-CSM has no clean state-restart hook). As wired today
    the orchestrator feeds it a *synthetic* ensemble, so the result is illustrative
    only. The orchestrator therefore refuses to run this mode unless
    ``assimilation.allow_uncoupled: true`` is set. For a coupled in-season workflow
    use ``mode: recalibration`` (see :class:`..recalibration.InSeasonRecalibrator`).
    """
    
    def __init__(self, cfg: dict):
        self.cfg = cfg
        enkf_cfg = cfg.get("assimilation", {}).get("enkf", {})
        self.n_ensemble = enkf_cfg.get("n_ensemble", 50)
        self.inflation = enkf_cfg.get("inflation", 1.05)
        self.state_vars = enkf_cfg.get("state_variables", ["LAID", "CWAD"])
        self.rng = np.random.default_rng(cfg.get("calibrator", {}).get("seed", 42))
        
    def assimilate(self, ensemble_states: np.ndarray, obs_var: str, obs_value: float, obs_sigma: float) -> np.ndarray:
        """Run one EnKF update step.
        
        Parameters
        ----------
        ensemble_states : np.ndarray of shape (n_ensemble, n_state_vars)
            The forecast ensemble states before update.
        obs_var : str
            The name of the variable observed.
        obs_value : float
            The observed value.
        obs_sigma : float
            Standard deviation of the observation error.
            
        Returns
        -------
        updated_states : np.ndarray of shape (n_ensemble, n_state_vars)
            The updated ensemble states.
        """
        if obs_var not in self.state_vars:
            logger.warning(f"Observed variable {obs_var} not in EnKF state variables: {self.state_vars}")
            return ensemble_states
            
        obs_var_idx = self.state_vars.index(obs_var)
        
        x_mean = ensemble_states.mean(axis=0)
        X = ensemble_states - x_mean
        
        X *= self.inflation
        
        H_x = ensemble_states[:, obs_var_idx]
        
        HP = X[:, obs_var_idx]
        HPHt = np.var(H_x) + (obs_sigma ** 2)
        if HPHt == 0:
            logger.warning("Innovation covariance HPHt is zero. Skipping EnKF update.")
            return ensemble_states
            
        K = (X.T @ HP) / ((self.n_ensemble - 1) * HPHt)
        
        obs_perturbed = obs_value + obs_sigma * self.rng.standard_normal(self.n_ensemble)
        
        innovation = obs_perturbed - H_x
        updated_states = ensemble_states + np.outer(innovation, K)
        
        logger.info(f"EnKF update for {obs_var}: obs = {obs_value}, mean before = {x_mean[obs_var_idx]:.3f}, mean after = {updated_states.mean(axis=0)[obs_var_idx]:.3f}")
        return updated_states
